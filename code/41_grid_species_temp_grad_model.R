# ============================================================
# Scientific question / 科学问题:
#   v1 province-level hazard model uses temp_grad_z =
#   (temp_anom[province, year] - temp_native_anom[species, year]),
#   a species-specific temperature gradient measuring how much a
#   province's climate anomaly deviates from the species' native-
#   range climate anomaly. The v2 grid models (scripts 28/40)
#   replaced this with time-invariant climate_velocity_z or
#   province-level prov_temp_anom_z, dropping the species-native
#   component entirely. This script restores the correct grid-level
#   analog:
#     grid_temp_grad = grid_temp_anom[grid, year] - temp_native_anom[species, year]
#   and re-fits M0-M4 to test whether the climate x effort
#   interaction survives with the species-specific gradient.
#   v1 省级模型使用物种特异温度梯度 temp_grad_z；v2 网格模型用
#   时不变 climate_velocity_z 替代，丢失了物种原生气候分量。
#   本脚本恢复正确的网格级物种特异温度梯度，重拟合 M0-M4，
#   验证 climate × effort 交互是否稳健。
#
# Objective / 分析目标:
#   1. Produce year-resolved grid temperature anomaly panel.
#      构建网格年际温度异常面板。
#   2. Merge species-native climate into the grid risk set.
#      将物种原生气候合并到网格风险集。
#   3. Compute grid_temp_grad = grid_temp_anom - temp_native_anom.
#      计算网格级物种特异温度梯度。
#   4. Fit M0-M4 with grid_temp_grad_z as climate_z at 100km.
#      以 grid_temp_grad_z 为气候变量拟合 M0-M4。
#   5. Fit M0-M4 with climate_velocity_z and prov_temp_anom_z
#      for comparison.
#      同时拟合 climate_velocity_z 和 prov_temp_anom_z 供对比。
#   6. Emit comparison table of model results.
#      输出模型比较表。
#
# Input data / 输入数据:
#   data/raw/climate_metrics_province_year.csv       (province-year temp_anom)
#   data/derived/community_grid_100km_climate_native.csv  (grid WorldClim bio1)
#   data/derived/risk_set_grid_100km_v2.csv              (existing 9.4M-row risk set)
#   data/raw/grid_100km_base.csv                         (grid-province mapping)
#   bird_new_record_hazard_model/.../species_year_native_climate.csv (species-native)
#   data/raw/species_range_native_anom.csv               (alternative species-native)
#
# Expected output / 预期输出:
#   data/derived/grid_100km_year_resolved_climate.csv
#   data/derived/risk_set_grid_100km_with_species_temp_grad.csv
#   results/tables/table_grid_species_temp_grad_comparison.csv
#   results/tables/table_grid_species_temp_grad_coefs.csv
#   results/diagnostics/table_grid_temp_grad_method_comparison.csv
#   figures/diagnostics/fig_grid_temp_grad_forest_comparison.pdf
#
# Key assumptions / 关键假设:
#   - If CHELSA rasters are absent, grid_temp_anom is approximated
#     via province-year temp_anom + WorldClim bio1 spatial offset
#     (Tier 2 fallback).
#     CHELSA 缺失时，用省级异常 + bio1 空间偏移近似。
#   - temp_native_anom is scale-invariant: it describes the species'
#     historical-range climate anomaly and does not depend on the
#     spatial unit of analysis.
#     物种原生气候异常不依赖空间尺度，与省级模型完全相同。
#   - Baseline period consistent with v1 (implicit in province panel).
#   - Study window 2002-2024.
#
# Main packages / 主要包: data.table, sf, terra, exactextractr,
#   glmmTMB, ggplot2, glue, arrow.
# Output directory / 输出路径: data/derived/, results/tables/,
#   results/diagnostics/, figures/diagnostics/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
  library(ggplot2)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))

# ---- 0. Configuration --------------------------------------------------------
CFG <- list(
  year_min      = 2002L,
  year_max      = 2024L,
  grid_km       = 100,
  subsample_max = 1e6,   # max rows for glmmTMB; keep all events, sample non-events
  # CHELSA path (if available)
  chelsa_dir    = path_spatial("chelsa"),
  # Species-native climate: prefer v1 file, fallback to v2 data/raw
  v1_native_path = file.path(
    "/Users/dingchenchen/Documents/New records/bird-new-distribution-records",
    "tasks/bird_new_record_hazard_model/results/combined_threshold_100_test",
    "derived_inputs/species_year_native_climate.csv")
)

# Z-score helper (same as script 40). Z-score 辅助函数。
zify <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# Logging helper. 日志辅助。
log <- function(...) message(paste0("[41] ", ...))

# ============================================================
# STEP 1: Build year-resolved grid temperature anomaly panel
# 构建网格年际温度异常面板
# ============================================================

log("=== Step 1: Year-resolved grid temperature anomaly ===")

# Three-tier strategy:
# Tier 1: CHELSA v2.1 native (full spatial + temporal variation)
# Tier 2: Province anomaly + WorldClim bio1 spatial offset (default)
# Tier 3: Province anomaly mirror (no within-province variation)

# --- 1a. Check CHELSA availability ---
chelsa_available <- dir.exists(CFG$chelsa_dir) &&
  length(list.files(file.path(CFG$chelsa_dir, "tas"), pattern = "\\.tif$")) > 0

if (chelsa_available) {
  log("CHELSA rasters found -- using Tier 1 (native)")

  suppressPackageStartupMessages({
    library(terra)
    library(sf)
    library(exactextractr)
  })
  source(file.path("code", "utils", "utils_spatial.R"))

  # Load grid SF. 加载网格空间对象。
  grid_sf_path <- path_derived(glue("grid_{CFG$grid_km}km_sf.gpkg"))
  if (!file.exists(grid_sf_path)) {
    stop("[41] Grid SF missing: ", grid_sf_path,
         ". Run script 06/40 first to build grid infrastructure.")
  }
  grid_sf <- sf::st_read(grid_sf_path, quiet = TRUE)
  gid_col <- "grid_id"

  # Extract annual mean temperature per grid cell from CHELSA.
  # 从 CHELSA 逐年月栅格提取网格年均温。
  annual_grid_temp <- rbindlist(lapply(CFG$year_min:CFG$year_max, function(yy) {
    files <- list.files(file.path(CFG$chelsa_dir, "tas"),
                        pattern = glue("CHELSA_tas_.*{yy}.*\\.tif$"),
                        full.names = TRUE)
    if (length(files) == 0) return(NULL)
    rast_yr <- terra::rast(files)
    annual_mean <- terra::app(rast_yr, fun = mean, na.rm = TRUE)
    vals <- exactextractr::exact_extract(annual_mean, grid_sf,
                                         fun = "mean", progress = FALSE)
    data.table(grid_id = grid_sf[[gid_col]], year = yy,
               annual_mean_temp = vals)
  }))

  # Baseline from CHELSA 1981-2010. 计算 1981-2010 基线均值。
  baseline_grid_temp <- rbindlist(lapply(1981:2010, function(yy) {
    files <- list.files(file.path(CFG$chelsa_dir, "tas"),
                        pattern = glue("CHELSA_tas_.*{yy}.*\\.tif$"),
                        full.names = TRUE)
    if (length(files) == 0) return(NULL)
    rast_yr <- terra::rast(files)
    annual_mean <- terra::app(rast_yr, fun = mean, na.rm = TRUE)
    vals <- exactextractr::exact_extract(annual_mean, grid_sf,
                                         fun = "mean", progress = FALSE)
    data.table(grid_id = grid_sf[[gid_col]], year = yy,
               annual_mean_temp = vals)
  }))
  baseline_means <- baseline_grid_temp[,
    .(baseline_temp = mean(annual_mean_temp, na.rm = TRUE)), by = grid_id]

  grid_clim_yr <- merge(annual_grid_temp, baseline_means, by = "grid_id")
  grid_clim_yr[, grid_temp_anom := annual_mean_temp - baseline_temp]
  grid_clim_yr[, grid_temp_anom_source := "chelsa"]

} else {
  log("CHELSA not available -- using Tier 2 (prov_anom + bio1 offset)")

  # --- Tier 2: Province-year anomaly + WorldClim bio1 spatial offset ---
  # 省级年际异常 + WorldClim bio1 空间偏移
  #
  # Formula:
  #   grid_temp_anom = prov_temp_anom + (bio1[grid] - mean(bio1[grids in province]))
  #
  # prov_temp_anom captures the temporal signal (year-to-year variation).
  # bio1_offset captures within-province spatial heterogeneity (elevation/latitude).
  # prov_temp_anom 捕捉年际信号；bio1_offset 捕捉省内空间异质性。

  prov_clim <- fread(path_raw("climate_metrics_province_year.csv"))
  grid_clim <- fread(path_derived("community_grid_100km_climate_native.csv"))
  grid_base <- fread(path_raw("grid_100km_base.csv"))

  # Province-year temperature anomaly. 省级年际温度异常。
  prov_anom <- prov_clim[, .(province, year, temp_anom)]

  # bio1 spatial offset: each grid's bio1 minus province-mean bio1.
  # bio1 空间偏移：每格 bio1 减去省内均值。
  grid_bio1 <- merge(grid_base[, .(grid_id, province)],
                     grid_clim[, .(grid_id, bio1)],
                     by = "grid_id")
  prov_bio1_mean <- grid_bio1[, .(bio1_prov_mean = mean(bio1, na.rm = TRUE)),
                               by = province]
  grid_bio1 <- merge(grid_bio1, prov_bio1_mean, by = "province")
  grid_bio1[, bio1_offset := bio1 - bio1_prov_mean]

  # Merge province anomaly with grid bio1 offset.
  # 合并省级异常与网格 bio1 偏移。
  grid_prov_anom <- merge(grid_base[, .(grid_id, province)],
                          prov_anom, by = "province", allow.cartesian = TRUE)
  grid_clim_yr <- merge(grid_prov_anom[, .(grid_id, year, temp_anom)],
                        grid_bio1[, .(grid_id, bio1_offset)],
                        by = "grid_id")
  grid_clim_yr[, grid_temp_anom := temp_anom + bio1_offset]
  grid_clim_yr[, grid_temp_anom_source := "prov_bio1_offset"]

  # Also compute Tier 3 (province mirror) for comparison.
  # 同时计算 Tier 3（省级复制）供对比。
  grid_clim_yr[, grid_temp_anom_prov_mirror := temp_anom]
}

# Z-score grid_temp_anom across all grid-year rows (global standardization).
# 全局 z-score。
grid_clim_yr[, grid_temp_anom_z := zify(grid_temp_anom)]

# Diagnostic: within-province variation check. 省内变异检查。
prov_var_diag <- grid_clim_yr[, .(
  sd_within = sd(grid_temp_anom, na.rm = TRUE),
  n_grids = uniqueN(grid_id)
), by = .(province, year)]
log("Within-province SD of grid_temp_anom (median across province-years): ",
    sprintf("%.4f", median(prov_var_diag$sd_within, na.rm = TRUE)))

# Save. 保存。
ensure_dir(dirname(path_derived("grid_100km_year_resolved_climate.csv")))
fwrite(grid_clim_yr, path_derived("grid_100km_year_resolved_climate.csv"))
log("Saved grid climate panel: ", nrow(grid_clim_yr), " rows")

# ============================================================
# STEP 2: Load species-native climate
# 加载物种原生气候
# ============================================================

log("=== Step 2: Species-native climate ===")

if (file.exists(CFG$v1_native_path)) {
  species_native <- fread(CFG$v1_native_path)
  # columns: species, year, temp_native_anom, prec_native_anom
  setnames(species_native,
           c("temp_native_anom", "prec_native_anom"),
           c("temp_native_anom_species", "prec_native_anom_species"))
  log("Using v1 species_year_native_climate.csv (", nrow(species_native), " rows)")
} else {
  species_native <- fread(path_raw("species_range_native_anom.csv"))
  # columns: year, temp_native_anom_range, species
  setnames(species_native, "temp_native_anom_range", "temp_native_anom_species")
  log("Using species_range_native_anom.csv (", nrow(species_native), " rows)")
}
species_native <- species_native[year >= CFG$year_min & year <= CFG$year_max]
log("Species with native climate: ", uniqueN(species_native$species),
    " | Year range: ", range(species_native$year))

# ============================================================
# STEP 3: Build enhanced risk set
# 构建增强风险集
# ============================================================

log("=== Step 3: Build enhanced risk set ===")

# Load existing risk set. Only select needed columns to save memory.
# 读取现有风险集，仅选必要列以节省内存。
risk <- fread(
  path_derived("risk_set_grid_100km_v2.csv"),
  select = c("province", "year", "grid_id", "species", "event",
             "climate_velocity_z", "mahalanobis_dist_z",
             "log_n_events_z", "log_n_observers_z", "log_n_dates_z",
             "effort_pc1_z", "prov_temp_anom_z"),
  encoding = "UTF-8"
)
log("Loaded risk set: ", nrow(risk), " rows, ", uniqueN(risk$species), " species")

# Merge year-resolved grid temperature. 合并网格年际温度。
risk <- merge(risk,
              grid_clim_yr[, .(grid_id, year, grid_temp_anom,
                               grid_temp_anom_z, grid_temp_anom_source)],
              by = c("grid_id", "year"), all.x = TRUE)

# Merge species-native climate. 合并物种原生气候。
risk <- merge(risk,
              species_native[, .(species, year, temp_native_anom_species)],
              by = c("species", "year"), all.x = TRUE)

# Compute grid_temp_grad. 计算网格级物种特异温度梯度。
risk[, grid_temp_grad := grid_temp_anom - temp_native_anom_species]
risk[, grid_temp_grad_z := zify(grid_temp_grad)]

# Year centering (for v1-consistent M0 that includes year_c).
# 年份中心化（与 v1 省级模型一致）。
risk[, year_c := year - 2013]

# Diagnostic: grid_temp_grad_z distribution. 分布诊断。
log("grid_temp_grad_z: mean = ", sprintf("%.4f", mean(risk$grid_temp_grad_z, na.rm = TRUE)),
    " | SD = ", sprintf("%.4f", sd(risk$grid_temp_grad_z, na.rm = TRUE)),
    " | NAs = ", sum(is.na(risk$grid_temp_grad_z)))

# Correlation between grid_temp_grad_z and prov_temp_anom_z.
# 两者相关性诊断。
corr_diag <- cor(risk$grid_temp_grad_z, risk$prov_temp_anom_z, use = "complete.obs")
log("Correlation grid_temp_grad_z ~ prov_temp_anom_z: ", sprintf("%.4f", corr_diag))

# Within-province variation of grid_temp_grad. 省内变异。
within_prov_var <- risk[, .(sd_temp_grad = sd(grid_temp_grad_z, na.rm = TRUE)),
                        by = .(province, year)]
log("Within-province SD of grid_temp_grad_z (median): ",
    sprintf("%.4f", median(within_prov_var$sd_temp_grad, na.rm = TRUE)))

# Complete-case filter. 完整案例筛选。
risk_cc <- risk[!is.na(grid_temp_grad_z) & !is.na(log_n_events_z)]
log("Complete cases: ", nrow(risk_cc), " | events: ", sum(risk_cc$event))

# Save enhanced risk set. 保存增强风险集。
fwrite(risk_cc, path_derived("risk_set_grid_100km_with_species_temp_grad.csv"))
log("Saved enhanced risk set")

# ============================================================
# STEP 4: Fit hazard models M0-M4
# 拟合离散时间风险模型
# ============================================================

log("=== Step 4: Fit hazard models ===")

# --- 4a. Stratified subsampling for large risk sets ---
# 大数据集分层子采样：保留全部事件，采样非事件行。
subsample_stratified <- function(dt, max_rows = CFG$subsample_max) {
  if (nrow(dt) <= max_rows) return(dt)
  events    <- dt[event == 1L]
  non_events <- dt[event == 0L]
  n_sample  <- max_rows - nrow(events)
  if (n_sample < nrow(non_events)) {
    non_events <- non_events[sample(.N, n_sample)]
  }
  rbindlist(list(events, non_events))
}

# --- 4b. Define fitting function ---
# 定义模型拟合函数（M0 含 year_c，与 v1 省级一致）。
fit_models <- function(dt, scale_label, re_form) {
  if (!"climate_z" %in% names(dt))
    stop("internal: need climate_z in dt for fit_models()")
  if (!"effort_z" %in% names(dt))
    stop("internal: need effort_z in dt for fit_models()")
  dt <- dt[!is.na(event) & !is.na(climate_z) & !is.na(effort_z)]
  dt <- subsample_stratified(dt)
  if (nrow(dt) < 200L) return(NULL)
  log("  scale = ", scale_label, " | rows = ", nrow(dt),
      " | events = ", sum(dt$event))

  # M0 includes year_c (v1 convention). M0 含 year_c（v1 惯例）。
  forms <- list(
    M0 = sprintf("event ~ year_c + %s", re_form),
    M1 = sprintf("event ~ year_c + effort_z + %s", re_form),
    M2 = sprintf("event ~ year_c + climate_z + %s", re_form),
    M3 = sprintf("event ~ year_c + climate_z + effort_z + %s", re_form),
    M4 = sprintf("event ~ year_c + climate_z * effort_z + %s", re_form))

  fits <- list()
  for (nm in names(forms)) {
    t0 <- Sys.time()
    fit <- tryCatch(
      glmmTMB::glmmTMB(stats::as.formula(forms[[nm]]),
                       data = dt,
                       family = stats::binomial(link = "cloglog"),
                       control = glmmTMBControl(
                         optCtrl = list(iter.max = 1000, eval.max = 1000))),
      error = function(e) {
        log("    ", nm, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(fit)) {
      pdHess <- isTRUE(fit$sdr$pdHess)
      log(sprintf("    %s fitted (%.1fs, AIC=%.1f, pdHess=%s)",
                  nm, as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  AIC(fit), pdHess))
      fits[[nm]] <- fit
    }
  }
  fits
}

# --- 4c. Extract results ---
# 提取模型系数与诊断。
extract_results <- function(fits, scale_label) {
  rows <- list()
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    cf <- tryCatch(glmmTMB::fixef(fit)$cond, error = function(e) NULL)
    se <- tryCatch(sqrt(diag(stats::vcov(fit)$cond)), error = function(e) NULL)
    if (is.null(cf)) next
    pdHess <- isTRUE(fit$sdr$pdHess)
    for (tm in names(cf)) {
      i <- match(tm, names(cf))
      beta <- cf[i]; sei <- if (!is.null(se)) se[i] else NA
      rows[[length(rows) + 1L]] <- data.table::data.table(
        scale = scale_label, model = nm, term = tm,
        beta = beta, se = sei,
        hr = exp(beta),
        hr.low = exp(beta - 1.96 * sei),
        hr.high = exp(beta + 1.96 * sei),
        p.value = 2 * stats::pnorm(-abs(beta / sei)),
        AIC = AIC(fit), n_rows = nobs(fit), pdHess = pdHess)
    }
  }
  data.table::rbindlist(rows, fill = TRUE)
}

# --- 4d. Prepare data for fitting ---
# 准备拟合数据。
risk_fit <- copy(risk_cc)
risk_fit[, species := factor(as.character(species))]
risk_fit[, grid_id := factor(as.character(grid_id))]

re_form <- "(1|species) + (1|grid_id)"

all_results <- list()

# --- Model Set A: grid_temp_grad_z (species-specific gradient) ---
# 模型组 A：物种特异温度梯度（本脚本核心）。
log("--- Model Set A: grid_temp_grad_z ---")
risk_fit[, climate_z := grid_temp_grad_z]
risk_fit[, effort_z  := log_n_events_z]
fits_A <- fit_models(risk_fit, "100km_grid_temp_grad", re_form)
if (!is.null(fits_A)) {
  all_results[["A"]] <- extract_results(fits_A, "100km_grid_temp_grad")
}

# --- Model Set B: climate_velocity_z (current v2 approach) ---
# 模型组 B：气候速率（v2 当前方案）。
log("--- Model Set B: climate_velocity_z ---")
risk_fit[, climate_z := climate_velocity_z]
fits_B <- fit_models(risk_fit, "100km_climate_velocity", re_form)
if (!is.null(fits_B)) {
  all_results[["B"]] <- extract_results(fits_B, "100km_climate_velocity")
}

# --- Model Set C: prov_temp_anom_z (province-level fallback) ---
# 模型组 C：省级温度异常回退。
log("--- Model Set C: prov_temp_anom_z ---")
risk_fit[, climate_z := prov_temp_anom_z]
fits_C <- fit_models(risk_fit, "100km_prov_temp_anom", re_form)
if (!is.null(fits_C)) {
  all_results[["C"]] <- extract_results(fits_C, "100km_prov_temp_anom")
}

# Combine all results. 合并所有结果。
res_all <- rbindlist(all_results, fill = TRUE)

# ============================================================
# STEP 5: Comparison tables and diagnostics
# 结果对比表与诊断
# ============================================================

log("=== Step 5: Comparison tables ===")

ensure_dir(path_tables())
ensure_dir(path_diagnostics())

# --- 5a. Full coefficient table ---
fwrite(res_all, path_tables("table_grid_species_temp_grad_coefs.csv"))
log("Saved coefficient table: ", nrow(res_all), " rows")

# --- 5b. M4 interaction comparison (core scientific finding) ---
# M4 交互项对比（核心科学发现）。
m4_interact <- res_all[model == "M4" & grepl(":", term)]
if (nrow(m4_interact) > 0) {
  m4_interact[, climate_spec := fcase(
    grepl("temp_grad", scale), "grid_temp_grad (species-specific)",
    grepl("climate_velocity", scale), "climate_velocity (time-invariant)",
    grepl("prov_temp_anom", scale), "prov_temp_anom (province-level)",
    default = scale)]
  fwrite(m4_interact, path_tables("table_grid_species_temp_grad_comparison.csv"))
  log("M4 interaction comparison:")
  for (i in seq_len(nrow(m4_interact))) {
    log(sprintf("  %s: HR = %.3f [%.3f, %.3f], p = %.2e",
                m4_interact$climate_spec[i],
                m4_interact$hr[i],
                m4_interact$hr.low[i],
                m4_interact$hr.high[i],
                m4_interact$p.value[i]))
  }
}

# --- 5c. AIC comparison across model sets ---
# AIC 对比。
aic_comp <- res_all[, .(scale, model, AIC, n_rows, pdHess)]
aic_comp[, dAIC := AIC - min(AIC), by = scale]
log("AIC comparison (M4 only):")
aic_m4 <- aic_comp[model == "M4"]
if (nrow(aic_m4) > 0) {
  for (i in seq_len(nrow(aic_m4))) {
    log(sprintf("  %s: AIC = %.1f, dAIC = %.1f",
                aic_m4$scale[i], aic_m4$AIC[i], aic_m4$dAIC[i]))
  }
}

# --- 5d. Method comparison diagnostics ---
# 方法对比诊断。
method_diag <- data.table(
  method = c("grid_temp_grad_z", "climate_velocity_z", "prov_temp_anom_z"),
  description = c("Species-specific gradient (grid anomaly - native anomaly)",
                   "Climate velocity (time-invariant, Loarie proxy)",
                   "Province temperature anomaly (spatial mirror)"),
  species_specific = c(TRUE, FALSE, FALSE),
  year_resolved = c(TRUE, FALSE, TRUE),
  within_province_variation = c(TRUE, FALSE, FALSE),
  correlation_with_grid_temp_grad = c(1.0, NA_real_, corr_diag)
)
# Fill in correlations. 填充相关性。
if ("B" %in% names(all_results)) {
  method_diag[method == "climate_velocity_z",
              correlation_with_grid_temp_grad :=
                cor(risk_cc$grid_temp_grad_z, risk_cc$climate_velocity_z,
                    use = "complete.obs")]
}
fwrite(method_diag, path_diagnostics("table_grid_temp_grad_method_comparison.csv"))

# ============================================================
# STEP 6: Forest plot comparison
# 森林图对比
# ============================================================

log("=== Step 6: Forest plot ===")

if (nrow(m4_interact) > 0) {
  p_forest <- ggplot(m4_interact,
                     aes(x = hr, y = reorder(climate_spec, hr),
                         xmin = hr.low, xmax = hr.high)) +
    geom_point(size = 3, colour = "steelblue") +
    geom_errorbarh(height = 0.2, linewidth = 0.8) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    labs(x = "Hazard Ratio (climate x effort interaction)",
         y = "",
         title = "M4 interaction: grid-scale climate x effort",
         subtitle = "Comparison across climate variable specifications") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(size = 10))

  ensure_dir(path_diag_fig())
  ggsave(path_diag_fig("fig_grid_temp_grad_forest_comparison.pdf"),
         p_forest, width = 18, height = 8, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(path_diag_fig("fig_grid_temp_grad_forest_comparison.png"),
         p_forest, width = 18, height = 8, units = "cm", dpi = 600)
  log("Saved forest plot")
}

# ============================================================
# STEP 7: Session info
# 会话信息
# ============================================================

if (exists("dump_session_info")) {
  dump_session_info(path_logs("41_grid_species_temp_grad_model_sessionInfo.txt"))
}
log("done.")
