# ============================================================
# Scientific question / 科学问题:
#   The hazard model M0-M4 was published at province scale only.
#   Does the climate x effort interaction in new-bird-record
#   detections hold at finer administrative scales (prefecture / 市,
#   county / 县) and at the 100 km regular grid?
#   省级 hazard 模型已发表；本脚本将同一套 M0-M4 模型在市、县、
#   100 km 网格三个尺度上独立重建风险集并重新拟合，检验
#   climate x effort 交互在更细尺度上是否稳健。
#
# Objective / 分析目标:
#   - Build risk sets at prefecture, county, and 100 km grid scales
#   - Fit M0-M4 hazard models (cloglog logistic GLMM) at each scale
#   - Produce coefficient/AIC tables and an M4 forest plot across 4
#     scales (province + 3 new scales), ordered county / prefecture /
#     province / 100km
#   - Write all outputs to outputs_multiscale/ so the canonical
#     province-level results under results/, figures/, data/derived/
#     are NEVER touched
#
# Input data / 输入数据:
#   - data/raw/events_100km_grid_assigned.csv         (lon/lat events)
#   - data/raw/effort_panel_upgraded.csv              (province effort)
#   - data/raw/climate_metrics_province_year.csv      (province climate)
#   - data/raw/hazard_risk_upgraded_complete_case.csv (province SDM risk)
#   - data/derived/risk_set_grid_100km_v2.csv         (grid risk set, prebuilt)
#   - tasks/bird_hazard_model_effort_upgrade/2019中国地图-审图号GS(2019)1822号/
#       市（等积投影）.shp  (371 prefectures)
#       县（等积投影）.shp  (2901 counties)
#       省（等积投影）.shp  (provinces, for centroid lookup)
#
# Main workflow / 主要流程:
#   1. Setup paths, parse CLI args, route every output under
#      outputs_multiscale/
#   2. Load events, province SDM candidates, province climate/effort
#   3. Fit province models (M0-M4) for forest-plot reference
#   4. Spatially join events to prefecture, build prefecture risk set,
#      fit M0-M4
#   5. Same for county
#   6. Load prebuilt 100 km grid risk set, fit M0-M4 with (1|species)
#      + (1|grid_id)
#   7. Save coefficients, AIC table, forest plot, sessionInfo
#
# Expected output / 预期输出:
#   outputs_multiscale/data/derived/
#     risk_set_province.csv
#     risk_set_prefecture.csv
#     risk_set_county.csv
#     risk_set_grid_100km.csv
#   outputs_multiscale/results/tables/
#     table_multiscale_coefficients.csv
#     table_multiscale_aic.csv
#   outputs_multiscale/figures/main/
#     fig_multiscale_M4_forest.pdf / .png
#   outputs_multiscale/logs/
#     51_multiscale.log
#     51_sessionInfo.txt
#
# Key assumptions / 关键假设:
#   - Climate at prefecture/county is INHERITED from province (same
#     convention as scripts 38 and 40). Within-province climate
#     variation requires a separate methodological upgrade.
#   - Effort at prefecture/county is event-share allocation of
#     province totals (matches build_admin_risk in script 40).
#   - SDM candidate units = all prefectures/counties inside the SDM-
#     suitable provinces (no finer SDM threshold available).
#   - Year window: 2002-2024.
#   - Family: binomial(cloglog) -> discrete-time hazard.
#
# Main packages / 主要包: data.table, sf, glmmTMB, ggplot2.
# Output directory / 输出路径: outputs_multiscale/
# ============================================================

# ---- 0. CLI args ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
arg <- list()
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--base-dir")        { arg$base_dir   <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--output-dir") { arg$output_dir <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--v1-dir")     { arg$v1_dir     <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--year-min")   { arg$year_min   <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--year-max")   { arg$year_max   <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--skip-county"){ arg$skip_county <- TRUE; i <- i + 1L }
  else { i <- i + 1L }
}
`%||%` <- function(a, b) if (is.null(a)) b else a

BASE_DIR   <- arg$base_dir   %||% "/Users/dingchenchen/Documents/New records/bird-new-distribution-records/tasks/bird_hazard_model_effort_upgrade_v2"
OUTPUT_DIR <- arg$output_dir %||% "outputs_multiscale"
V1_DIR     <- arg$v1_dir     %||% file.path(dirname(BASE_DIR), "bird_hazard_model_effort_upgrade")
YEAR_MIN   <- arg$year_min   %||% 2002L
YEAR_MAX   <- arg$year_max   %||% 2024L
SKIP_COUNTY <- isTRUE(arg$skip_county)

OUT <- file.path(BASE_DIR, OUTPUT_DIR)
ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
for (sd in c("data/derived", "results/tables", "results/diagnostics",
             "figures/main", "figures/diagnostics", "logs")) {
  ens(file.path(OUT, sd))
}

log_file <- file.path(OUT, "logs", "51_multiscale.log")
log <- function(...) {
  msg <- sprintf("[51 %s] %s\n", format(Sys.time(), "%H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg); cat(msg, file = log_file, append = TRUE)
}
log("=== 51_multiscale_full.R START ===")
log("base_dir   = ", BASE_DIR)
log("output_dir = ", OUT)
log("v1_dir     = ", V1_DIR)
log("year window = ", YEAR_MIN, "-", YEAR_MAX)
log("skip_county = ", SKIP_COUNTY)

# ---- 1. Packages ----------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(glmmTMB)
  library(ggplot2)
})
sf::sf_use_s2(FALSE)   # 平面几何，避免 China shp 的 S2 拓扑报错
options(warn = 1)
set.seed(42)

zify <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# ---- 2. Load core data ----------------------------------------------------
log("=== Step 2: loading events, province SDM candidates, climate, effort ===")

events_file <- file.path(BASE_DIR, "data", "raw", "events_100km_grid_assigned.csv")
ev <- fread(events_file, encoding = "UTF-8")
setnames(ev, tolower(names(ev)))
if (!"year" %in% names(ev) && "pub_year" %in% names(ev))
  setnames(ev, "pub_year", "year")
ev[, year := as.integer(year)]
ev <- ev[!is.na(longitude) & !is.na(latitude) &
         year >= YEAR_MIN & year <= YEAR_MAX]
log("events: ", nrow(ev), " rows, ", uniqueN(ev$species), " species, ",
    uniqueN(ev$province), " provinces")

# Province SDM candidates (canonical source = v1 complete_case risk set)
# Search v1 dir, then v2/data/raw (where rsync places it on server).
prov_risk_candidates <- c(
  file.path(V1_DIR, "data", "hazard_risk_upgraded_complete_case.csv"),
  file.path(BASE_DIR, "data", "raw", "hazard_risk_upgraded_complete_case.csv"),
  file.path(BASE_DIR, "data", "hazard_risk_upgraded_complete_case.csv")
)
prov_risk_file <- prov_risk_candidates[file.exists(prov_risk_candidates)][1]
if (is.na(prov_risk_file))
  stop("Cannot find hazard_risk_upgraded_complete_case.csv. Tried: ",
       paste(prov_risk_candidates, collapse = " | "))
log("province SDM source: ", prov_risk_file)
prov_risk_canon <- fread(prov_risk_file, encoding = "UTF-8")
sdm_candidates <- unique(prov_risk_canon[, .(species, province)])
log("SDM candidate (species x province) = ", nrow(sdm_candidates))

# Province climate panel
prov_clim <- fread(file.path(BASE_DIR, "data", "raw",
                             "climate_metrics_province_year.csv"),
                   encoding = "UTF-8")
setnames(prov_clim, tolower(names(prov_clim)))

# Province effort panel
prov_eff <- fread(file.path(BASE_DIR, "data", "raw",
                            "effort_panel_upgraded.csv"),
                  encoding = "UTF-8")
setnames(prov_eff, tolower(names(prov_eff)))

# ---- 3. Province-level risk set + fit -------------------------------------
log("=== Step 3: fitting PROVINCE models for forest-plot reference ===")
# 使用 v1 已经准备好的省级 risk set；这里只重新跑 M0-M4，不动原始文件。
prov_dt <- copy(prov_risk_canon)
setnames(prov_dt, tolower(names(prov_dt)))

# Choose canonical climate/effort z columns
# v1 risk set uses temp_grad_z, not climate_velocity_z.
# Merge province climate panel to get climate_velocity_z if missing.
if (!"climate_velocity_z" %in% names(prov_dt)) {
  prov_dt <- merge(prov_dt,
    prov_clim[, .(province, year, climate_velocity_z)],
    by = c("province", "year"), all.x = TRUE)
}
if (!"climate_velocity_z" %in% names(prov_dt) && "temp_grad_z" %in% names(prov_dt))
  prov_dt[, climate_velocity_z := temp_grad_z]   # fallback

prov_dt[, effort_z := if ("log_effort_visits_z" %in% names(prov_dt)) log_effort_visits_z
                      else if ("effort_pc1_z" %in% names(prov_dt)) effort_pc1_z
                      else zify(log1p(n_visits))]

prov_dt[, climate_z := climate_velocity_z]

# Save the province risk set we used (not overwriting v1's)
fwrite(prov_dt, file.path(OUT, "data", "derived", "risk_set_province.csv"))
log("province risk set saved: rows=", nrow(prov_dt),
    " events=", sum(prov_dt$event, na.rm = TRUE))

# ---- 4. Spatial setup for prefecture / county -----------------------------
log("=== Step 4: loading prefecture / county shapefiles ===")

shp_root <- NULL
# Per-scale shp search: prefer china_shp (complete .shx), fallback to GS2019
# 每个 scale 独立搜索，优先用 china_shp（有完整 .shx），其次用 GS2019
shp_dirs_ordered <- c(
  file.path(BASE_DIR, "data", "spatial", "china_shp"),
  file.path(BASE_DIR, "data", "spatial", "basemap_GS2019_1822"),
  file.path(V1_DIR, "2019中国地图-审图号GS(2019)1822号"),
  file.path(BASE_DIR, "2019中国地图-审图号GS(2019)1822号"),
  file.path(BASE_DIR, "data", "spatial", "basemap"),
  file.path(BASE_DIR, "data", "spatial")
)

# For each scale, find the best shapefile (must have .shx companion)
find_shp_with_shx <- function(prefix) {
  for (d in shp_dirs_ordered) {
    if (!dir.exists(d)) next
    hits <- list.files(d, pattern = paste0("^", prefix, ".*\\.shp$"),
                       full.names = TRUE)
    hits <- hits[!grepl("境界线", hits)]
    # Prefer 等积投影 version
    eq <- hits[grepl("等积投影", hits)]
    poly <- hits[!grepl("线", hits)]
    candidates <- if (length(eq) > 0L) eq else if (length(poly) > 0L) poly else hits
    for (h in candidates) {
      shx <- sub("\\.shp$", ".shx", h)
      SHX <- sub("\\.shp$", ".SHX", h)
      if (file.exists(shx) || file.exists(SHX)) return(h)
    }
  }
  # Last resort: return first hit even without .shx (GDAL may rebuild)
  for (d in shp_dirs_ordered) {
    if (!dir.exists(d)) next
    hits <- list.files(d, pattern = paste0("^", prefix, ".*\\.shp$"),
                       full.names = TRUE)
    hits <- hits[!grepl("境界线", hits)]
    poly <- hits[!grepl("线", hits)]
    if (length(poly) > 0L) return(poly[1])
    if (length(hits) > 0L) return(hits[1])
  }
  NA_character_
}

pref_shp <- find_shp_with_shx("市")
cnty_shp <- find_shp_with_shx("县")
prov_shp <- find_shp_with_shx("省")
log("prefecture shp: ", pref_shp)
log("county shp:     ", cnty_shp)
log("province shp:   ", prov_shp)

# Province polygon (used to label admin units with their parent province
# from the canonical province name set used elsewhere in v2)
PROV_CN_EN <- c(
  "北京市" = "Beijing", "天津市" = "Tianjin", "河北省" = "Hebei",
  "山西省" = "Shanxi", "内蒙古自治区" = "Inner Mongolia",
  "辽宁省" = "Liaoning", "吉林省" = "Jilin", "黑龙江省" = "Heilongjiang",
  "上海市" = "Shanghai", "江苏省" = "Jiangsu", "浙江省" = "Zhejiang",
  "安徽省" = "Anhui", "福建省" = "Fujian", "江西省" = "Jiangxi",
  "山东省" = "Shandong", "河南省" = "Henan", "湖北省" = "Hubei",
  "湖南省" = "Hunan", "广东省" = "Guangdong",
  "广西壮族自治区" = "Guangxi", "海南省" = "Hainan",
  "重庆市" = "Chongqing", "四川省" = "Sichuan", "贵州省" = "Guizhou",
  "云南省" = "Yunnan", "西藏自治区" = "Tibet",
  "陕西省" = "Shaanxi", "甘肃省" = "Gansu", "青海省" = "Qinghai",
  "宁夏回族自治区" = "Ningxia", "新疆维吾尔自治区" = "Xinjiang",
  "台湾省" = "Taiwan", "香港特别行政区" = "Hong Kong",
  "澳门特别行政区" = "Macau")

prov_sf <- NULL
if (!is.na(prov_shp) && file.exists(prov_shp)) {
  prov_sf <- sf::st_read(prov_shp, quiet = TRUE) |>
    sf::st_make_valid()
  pn_col <- names(prov_sf)[vapply(prov_sf, function(x) {
    is.character(x) && any(grepl("[\u4e00-\u9fff]", x))
  }, logical(1))][1]
  prov_sf$province_en <- unname(PROV_CN_EN[as.character(prov_sf[[pn_col]])])
  prov_sf$province_en[is.na(prov_sf$province_en)] <-
    as.character(prov_sf[[pn_col]])[is.na(prov_sf$province_en)]
}

# Make events sf in shapefile CRS
admin_crs <- {
  if (!is.null(prov_sf)) {
    sf::st_crs(prov_sf)
  } else if (!is.na(pref_shp)) {
    sf::st_crs(sf::st_read(pref_shp, quiet = TRUE))
  } else {
    sf::st_crs(4326)
  }
}
ev_sf <- sf::st_as_sf(ev, coords = c("longitude", "latitude"), crs = 4326)
ev_sf <- sf::st_transform(ev_sf, admin_crs)
# Drop original province column from events to avoid collision with the
# polygon-assigned province after st_join (otherwise both will become
# province.x / province.y, which is fragile across sf versions).
if ("province" %in% names(ev_sf)) ev_sf$province <- NULL

# ---- 5. build_admin_risk: prefecture / county risk-set builder ------------
# Adapted from code/40_execute_five_scale_models.R lines 477-543
# (build_admin_risk()) but with a self-contained province lookup, since
# we cannot assume the original `prov_sf` global is in scope.
build_admin_risk <- function(shp_path, scale_label) {
  if (is.na(shp_path) || !file.exists(shp_path)) {
    log("skip ", scale_label, " (shp missing)"); return(NULL)
  }
  log("loading ", scale_label, " shp: ", basename(shp_path))
  admin <- sf::st_read(shp_path, quiet = TRUE) |>
    sf::st_make_valid() |>
    sf::st_transform(admin_crs)
  admin$unit_id <- paste0(scale_label, "_", seq_len(nrow(admin)))
  log("  ", scale_label, " polygons: ", nrow(admin))

  # Assign parent province via centroid -> province polygon
  if (!is.null(prov_sf)) {
    cent <- sf::st_centroid(admin)
    j <- sf::st_join(cent, prov_sf[, "province_en"], join = sf::st_within)
    admin$province <- j$province_en
  } else {
    # fallback: take Chinese 省 column if present in admin's own attributes
    cn_col <- names(admin)[vapply(admin, function(x) {
      is.character(x) && any(grepl("[\u4e00-\u9fff]", x))
    }, logical(1))]
    cn_col <- intersect(cn_col, c("\u7701", "province"))[1]
    if (!is.na(cn_col)) {
      admin$province <- unname(PROV_CN_EN[as.character(admin[[cn_col]])])
    } else {
      admin$province <- NA_character_
    }
  }
  admin$province[is.na(admin$province)] <- "Unknown"
  log("  province labels assigned: ",
      sum(admin$province != "Unknown"), "/", nrow(admin))

  # Spatial join events -> admin
  ev_admin <- sf::st_join(ev_sf, admin[, c("unit_id", "province")],
                          join = sf::st_within)
  ev_dt <- data.table::as.data.table(sf::st_drop_geometry(ev_admin))
  ev_dt <- ev_dt[!is.na(unit_id)]
  if (nrow(ev_dt) == 0L) {
    log("  no events fell in any ", scale_label, " polygon"); return(NULL)
  }
  log("  events mapped to ", scale_label, ": ", nrow(ev_dt),
      " (", uniqueN(ev_dt$unit_id), " units, ",
      uniqueN(ev_dt$species), " species)")

  # First arrival per (species, unit_id)
  fa <- ev_dt[, .(arrival_year = min(year, na.rm = TRUE)),
              by = .(species, unit_id, province)]

  # Event-share allocation of province effort to each unit-year
  cnt <- ev_dt[, .(n_events_unit = .N),
               by = .(unit_id, year, province)]
  prov_total <- cnt[, .(n_events_prov = sum(n_events_unit)),
                    by = .(province, year)]
  cnt <- merge(cnt, prov_total, by = c("province", "year"))
  cnt[, share := n_events_unit / pmax(n_events_prov, 1)]

  # Pull province effort (n_visits) and allocate
  eff_keep_cols <- intersect(c("province", "year", "n_visits",
                                "n_birding_days", "n_observers",
                                "effort_record"), names(prov_eff))
  cnt <- merge(cnt, prov_eff[, ..eff_keep_cols],
               by = c("province", "year"), all.x = TRUE)
  base_eff <- if ("n_visits" %in% names(cnt)) cnt$n_visits
              else if ("effort_record" %in% names(cnt)) cnt$effort_record
              else 1
  cnt[, n_events_unit_eff := share * pmax(base_eff, 1)]
  cnt[, log_n_events_z := zify(log1p(n_events_unit_eff)), by = year]

  # Build SDM-restricted (species x unit_id) cartesian within each
  # candidate province, then expand to years.
  unit_in_prov <- unique(data.table::as.data.table(sf::st_drop_geometry(admin))[,
    .(unit_id, province)])
  cand <- merge(sdm_candidates, unit_in_prov, by = "province",
                allow.cartesian = TRUE)
  log("  candidate (species x unit): ", nrow(cand))

  risk <- merge(cand,
                data.table::CJ(unit_id = unique(cand$unit_id),
                               year = YEAR_MIN:YEAR_MAX),
                by = "unit_id", allow.cartesian = TRUE)
  # merge back species,province columns lost in CJ
  risk <- merge(risk, cand, by = c("species", "unit_id", "province"))
  risk <- merge(risk, fa[, .(species, unit_id, arrival_year)],
                by = c("species", "unit_id"), all.x = TRUE)
  risk <- risk[is.na(arrival_year) | year <= arrival_year]
  risk[, event := as.integer(year == arrival_year)]
  risk[is.na(event), event := 0L]

  # Attach effort (allocated) and province-year climate
  risk <- merge(risk, cnt[, .(unit_id, year, log_n_events_z,
                              n_events_unit)],
                by = c("unit_id", "year"), all.x = TRUE)
  risk[is.na(log_n_events_z), log_n_events_z := 0]
  risk[is.na(n_events_unit),  n_events_unit  := 0L]
  clim_cols <- intersect(c("province", "year",
                            "climate_velocity_z", "warming_rate_z",
                            "temp_anom_z", "climate_exposure_z",
                            "mahalanobis_dist_z"), names(prov_clim))
  risk <- merge(risk, prov_clim[, ..clim_cols],
                by = c("province", "year"), all.x = TRUE)

  log("  ", scale_label, " risk set: rows=", nrow(risk),
      " events=", sum(risk$event))
  list(risk = risk, n_units = nrow(admin))
}

# ---- 6. Model fitter (M0-M4 cloglog GLMM) ---------------------------------
# Adapted from code/40_execute_five_scale_models.R lines 384-417.
fit_models <- function(dt, scale_label, re_form) {
  stopifnot("climate_z" %in% names(dt), "effort_z" %in% names(dt))
  dt <- dt[!is.na(event) & !is.na(climate_z) & !is.na(effort_z)]
  if (nrow(dt) < 200L) {
    log("  ", scale_label, ": too few rows (", nrow(dt), "), skip"); return(NULL)
  }
  log("  fit_models ", scale_label, ": rows=", nrow(dt),
      " events=", sum(dt$event))
  forms <- list(
    M0 = sprintf("event ~ 1 + %s", re_form),
    M1 = sprintf("event ~ effort_z + %s", re_form),
    M2 = sprintf("event ~ climate_z + %s", re_form),
    M3 = sprintf("event ~ climate_z + effort_z + %s", re_form),
    M4 = sprintf("event ~ climate_z * effort_z + %s", re_form))
  fits <- list()
  for (nm in names(forms)) {
    t0 <- Sys.time()
    fit <- tryCatch(
      glmmTMB::glmmTMB(stats::as.formula(forms[[nm]]),
                       data = dt,
                       family = stats::binomial(link = "cloglog")),
      error = function(e) { log("    ", nm, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(fit)) {
      log(sprintf("    %s fitted (%.1fs, AIC=%.1f)", nm,
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  AIC(fit)))
      fits[[nm]] <- fit
    }
  }
  fits
}

extract_results <- function(fits, scale_label) {
  rows <- list()
  if (is.null(fits)) return(data.table())
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    cf <- tryCatch(glmmTMB::fixef(fit)$cond, error = function(e) NULL)
    se <- tryCatch(sqrt(diag(stats::vcov(fit)$cond)), error = function(e) NULL)
    if (is.null(cf)) next
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
        AIC = AIC(fit), n_rows = nobs(fit))
    }
  }
  data.table::rbindlist(rows, fill = TRUE)
}

# ---- 7. Fit all four scales (incremental save) ----------------------------
# 每个尺度拟合后立即保存结果，避免中途崩溃丢失全部进度
# Save results after each scale so a crash doesn't lose all progress.

all_res <- list()

# --- 7a Province ---
log("=== Step 7a: fitting PROVINCE models ===")
fits_prov <- tryCatch({
  fit_models(
    prov_dt[, .(species, province, year, event, climate_z, effort_z)],
    "province", "(1|species) + (1|province)")
}, error = function(e) {
  log("PROVINCE fitting failed: ", conditionMessage(e)); NULL
})
res_prov <- extract_results(fits_prov, "province")
if (nrow(res_prov) > 0L) {
  fwrite(res_prov, file.path(OUT, "results", "tables",
                            "table_multiscale_province.csv"))
  log("province results saved: ", nrow(res_prov), " rows")
}
all_res[["province"]] <- res_prov

# --- 7b Prefecture ---
log("=== Step 7b: building + fitting PREFECTURE models ===")
fits_pref <- NULL
pref_res_file <- file.path(OUT, "results", "tables",
                           "table_multiscale_prefecture.csv")
pref_risk_file <- file.path(OUT, "data", "derived", "risk_set_prefecture.csv")

if (file.exists(pref_res_file)) {
  log("prefecture results already exist, loading: ", pref_res_file)
  res_pref <- fread(pref_res_file)
} else {
  pref <- tryCatch(build_admin_risk(pref_shp, "prefecture"),
                   error = function(e) {
                     log("prefecture build_admin_risk failed: ",
                         conditionMessage(e)); NULL
                   })
  if (!is.null(pref)) {
    fwrite(pref$risk, pref_risk_file)
    pref$risk[, climate_z := climate_velocity_z]
    pref$risk[, effort_z  := log_n_events_z]
    fits_pref <- tryCatch({
      fit_models(
        pref$risk[, .(species, unit_id, year, event, climate_z, effort_z)],
        "prefecture", "(1|species) + (1|unit_id)")
    }, error = function(e) {
      log("prefecture fit_models failed: ", conditionMessage(e)); NULL
    })
  }
  res_pref <- extract_results(fits_pref, "prefecture")
  if (nrow(res_pref) > 0L) {
    fwrite(res_pref, pref_res_file)
    log("prefecture results saved: ", nrow(res_pref), " rows")
  }
}
all_res[["prefecture"]] <- res_pref

# --- 7c County ---
log("=== Step 7c: building + fitting COUNTY models ===")
fits_cnty <- NULL
cnty_res_file <- file.path(OUT, "results", "tables",
                           "table_multiscale_county.csv")
cnty_risk_file <- file.path(OUT, "data", "derived", "risk_set_county.csv")

if (file.exists(cnty_res_file)) {
  log("county results already exist, loading: ", cnty_res_file)
  res_cnty <- fread(cnty_res_file)
} else if (!SKIP_COUNTY) {
  cnty <- tryCatch(build_admin_risk(cnty_shp, "county"),
                   error = function(e) {
                     log("county build_admin_risk failed: ",
                         conditionMessage(e)); NULL
                   })
  if (!is.null(cnty)) {
    fwrite(cnty$risk, cnty_risk_file)
    cnty$risk[, climate_z := climate_velocity_z]
    cnty$risk[, effort_z  := log_n_events_z]
    fits_cnty <- tryCatch({
      fit_models(
        cnty$risk[, .(species, unit_id, year, event, climate_z, effort_z)],
        "county", "(1|species) + (1|unit_id)")
    }, error = function(e) {
      log("county fit_models failed: ", conditionMessage(e)); NULL
    })
  }
  res_cnty <- extract_results(fits_cnty, "county")
  if (nrow(res_cnty) > 0L) {
    fwrite(res_cnty, cnty_res_file)
    log("county results saved: ", nrow(res_cnty), " rows")
  }
} else {
  log("county fitting skipped via --skip-county flag")
  res_cnty <- data.table()
}
all_res[["county"]] <- res_cnty

# --- 7d 100km Grid ---
log("=== Step 7d: building + fitting 100km GRID models ===")
fits_grid <- NULL
grid_res_file <- file.path(OUT, "results", "tables",
                           "table_multiscale_100km.csv")
grid_risk_out <- file.path(OUT, "data", "derived", "risk_set_grid_100km.csv")

if (file.exists(grid_res_file)) {
  log("100km grid results already exist, loading: ", grid_res_file)
  res_grid <- fread(grid_res_file)
} else {
  grid_file <- file.path(BASE_DIR, "data", "derived",
                         "risk_set_grid_100km_v2.csv")
  if (file.exists(grid_file)) {
    grid_dt <- tryCatch({
      fread(grid_file, encoding = "UTF-8")
    }, error = function(e) {
      log("grid risk set load failed: ", conditionMessage(e)); NULL
    })
    if (!is.null(grid_dt)) {
      setnames(grid_dt, tolower(names(grid_dt)))
      grid_dt <- grid_dt[year >= YEAR_MIN & year <= YEAR_MAX]
      fwrite(grid_dt, grid_risk_out)
      log("100km grid risk set: rows=", nrow(grid_dt),
          " events=", sum(grid_dt$event, na.rm = TRUE))
      grid_dt[, climate_z := climate_velocity_z]
      if (!"log_n_events_z" %in% names(grid_dt)) {
        if ("log_effort_visits_z" %in% names(grid_dt))
          grid_dt[, log_n_events_z := log_effort_visits_z]
        else if ("effort_pc1_z" %in% names(grid_dt))
          grid_dt[, log_n_events_z := effort_pc1_z]
      }
      grid_dt[, effort_z := log_n_events_z]
      fits_grid <- tryCatch({
        fit_models(
          grid_dt[, .(species, grid_id, year, event, climate_z, effort_z)],
          "100km", "(1|species) + (1|grid_id)")
      }, error = function(e) {
        log("100km grid fit_models failed: ", conditionMessage(e)); NULL
      })
    }
  } else {
    log("WARNING: ", grid_file, " not found; skipping 100km grid")
  }
  res_grid <- extract_results(fits_grid, "100km")
  if (nrow(res_grid) > 0L) {
    fwrite(res_grid, grid_res_file)
    log("100km grid results saved: ", nrow(res_grid), " rows")
  }
}
all_res[["100km"]] <- res_grid

# ---- 8. Collect, save coefficient + AIC tables ----------------------------
log("=== Step 8: assembling coefficient + AIC tables ===")
res <- data.table::rbindlist(all_res, fill = TRUE)
fwrite(res, file.path(OUT, "results", "tables",
                      "table_multiscale_coefficients.csv"))
log("coefficients table: ", nrow(res), " rows")

aic_tab <- unique(res[, .(scale, model, AIC, n_rows)])
aic_tab[, dAIC := AIC - min(AIC, na.rm = TRUE), by = scale]
setorder(aic_tab, scale, AIC)
fwrite(aic_tab, file.path(OUT, "results", "tables",
                          "table_multiscale_aic.csv"))
log("AIC table: ", nrow(aic_tab), " rows")

# ---- 9. M4 interaction forest plot ----------------------------------------
log("=== Step 9: M4 interaction forest plot ===")
m4_int <- res[model == "M4" & grepl(":", term)]
if (nrow(m4_int) > 0L) {
  m4_int[, scale := factor(scale,
    levels = c("county", "prefecture", "province", "100km"))]
  p <- ggplot(m4_int, aes(x = hr, y = scale)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                   height = 0.2, linewidth = 0.5) +
    geom_point(size = 2.5, colour = "#B40426") +
    scale_x_continuous(trans = "log",
                       breaks = c(0.25, 0.5, 1, 1.5, 2, 3)) +
    labs(title = "Climate x effort interaction (M4) across 4 spatial scales",
         subtitle = paste0("Hazard model: event ~ climate_z * effort_z + ",
                            "(1|species) + (1|unit). v2 multiscale rebuild."),
         x = "Hazard ratio (log scale)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(OUT, "figures", "main", "fig_multiscale_M4_forest.pdf"),
         p, width = 14, height = 7, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(file.path(OUT, "figures", "main", "fig_multiscale_M4_forest.png"),
         p, width = 14, height = 7, units = "cm", dpi = 600)
  log("forest plot saved")
} else {
  log("WARNING: no M4 interaction rows found, skipping forest plot")
}

# ---- 10. sessionInfo ------------------------------------------------------
sink(file.path(OUT, "logs", "51_sessionInfo.txt"))
print(sessionInfo())
sink()

log("=== 51_multiscale_full.R DONE ===")
log("outputs under: ", OUT)
