# ============================================================
# Scientific question / 科学问题:
#   Does the climate × effort interaction (M4) hold at finer
#   administrative spatial scales (prefecture 市 / county 县) in
#   addition to the headline province scale? v1 only fitted at
#   province + 100 km grid; 100 km grid v1 result already exists in
#   table_grid_model_coefficients.csv. This lite script populates
#   the prefecture + county scales WITHOUT the OOM-prone 100 km
#   grid refit, producing the three new scales for the v2 manuscript
#   five-scale claim.
#   省级是 v1 主结论；100km 网格 v1 已拟合。本脚本仅补市/县两个新
#   尺度，回避 glmmTMB 在 3-6M 行双交叉随机效应上的内存爆炸。
#
# Input data / 输入数据:
#   data/raw/hazard_risk_upgraded_complete_case.csv   (province risk set, v1)
#   data/raw/events_100km_grid_assigned.csv           (coordinate-level events)
#   data/raw/effort_panel_upgraded.csv                (province effort totals)
#   data/raw/climate_metrics_province_year.csv        (province × year climate)
#   data/spatial/basemap_GS2019_1822/{市,县,省}.shp    (admin polygons)
#
# Main workflow / 主要流程:
#   1. Refit province M0-M4 (reproducibility check vs v1's 1.292).
#   2. Spatial-join coordinate events to prefecture & county polygons.
#   3. Allocate effort = share × province_n_visits.
#   4. Build SDM-restricted (species × unit × year) risk sets.
#   5. Fit M0-M4 cloglog hazard at each scale.
#   6. Persist results/tables/table_three_scale_*.csv.
#
# Expected output / 预期输出:
#   results/tables/table_three_scale_coefficients.csv
#   results/tables/table_three_scale_model_comparison.csv
#   results/tables/risk_set_summary_three_scale.csv
#   figures/main/fig_three_scale_interaction_forest.{pdf,png}
#
# Main packages / 主要包: data.table, sf, glmmTMB, ggplot2.
# Output directory / 输出路径: results/tables/, figures/main/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(glmmTMB)
  library(ggplot2)
})
sf::sf_use_s2(FALSE)
options(warn = 1)
set.seed(42)

# ---- 0. Paths (portable) --------------------------------------------------
V2 <- normalizePath(".", mustWork = TRUE)
V1 <- Sys.getenv("V1_ROOT",
                  normalizePath(file.path(V2, "..",
                                          "bird_hazard_model_effort_upgrade"),
                                 mustWork = FALSE))
SHP_DIR <- file.path(V2, "data", "spatial", "basemap_GS2019_1822")
if (!dir.exists(SHP_DIR)) {
  SHP_DIR <- file.path(V1, "2019中国地图-审图号GS(2019)1822号")
}
stopifnot(dir.exists(SHP_DIR))

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ens(file.path(V2, "results", "tables"))
ens(file.path(V2, "figures", "main"))
ens(file.path(V2, "logs"))

log <- function(...) cat(sprintf("[40b %s] ", format(Sys.time(), "%H:%M:%S")),
                          ..., "\n", sep = "")

zify <- function(x) {
  s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# ---- 1. Province name translation (GS2019 中文 → English) -----------------
PROV_CN_EN <- c("北京市"="Beijing","天津市"="Tianjin","河北省"="Hebei",
  "山西省"="Shanxi","内蒙古自治区"="Inner Mongolia","辽宁省"="Liaoning",
  "吉林省"="Jilin","黑龙江省"="Heilongjiang","上海市"="Shanghai",
  "江苏省"="Jiangsu","浙江省"="Zhejiang","安徽省"="Anhui","福建省"="Fujian",
  "江西省"="Jiangxi","山东省"="Shandong","河南省"="Henan","湖北省"="Hubei",
  "湖南省"="Hunan","广东省"="Guangdong","广西壮族自治区"="Guangxi",
  "海南省"="Hainan","重庆市"="Chongqing","四川省"="Sichuan","贵州省"="Guizhou",
  "云南省"="Yunnan","西藏自治区"="Tibet","陕西省"="Shaanxi","甘肃省"="Gansu",
  "青海省"="Qinghai","宁夏回族自治区"="Ningxia","新疆维吾尔自治区"="Xinjiang",
  "台湾省"="Taiwan","香港特别行政区"="Hong Kong","澳门特别行政区"="Macau")

# ---- 2. Load province risk set + sdm candidates --------------------------
log("loading province risk set")
prov_risk <- fread(file.path(V1, "data",
                              "hazard_risk_upgraded_complete_case.csv"),
                    encoding = "UTF-8")
log("province risk: rows=", nrow(prov_risk),
    " species=", uniqueN(prov_risk$species),
    " provinces=", uniqueN(prov_risk$province),
    " events=", sum(prov_risk$event))

sdm_candidates <- unique(prov_risk[, .(species, province)])
log("SDM (species × province) candidates: ", nrow(sdm_candidates))

# Province climate panel
prov_clim <- fread(file.path(V1, "data",
                              "climate_metrics_province_year.csv"),
                    encoding = "UTF-8")

# Province effort panel (totals)
prov_eff <- fread(file.path(V1, "data", "effort_panel_upgraded.csv"),
                   encoding = "UTF-8")

# ---- 3. Coordinate-level events ------------------------------------------
log("loading events_100km_grid_assigned.csv (for coord-level join)")
ev <- fread(file.path(V1, "data", "events_100km_grid_assigned.csv"),
             encoding = "UTF-8")
setnames(ev, tolower(names(ev)))
if (!"year" %in% names(ev) && "pub_year" %in% names(ev))
  setnames(ev, "pub_year", "year")
ev <- ev[!is.na(longitude) & !is.na(latitude) &
          year >= 2002 & year <= 2024]
log("events with valid coords: ", nrow(ev))
ev_sf <- sf::st_as_sf(ev, coords = c("longitude", "latitude"), crs = 4326)

# ---- 4. Province scale — refit M0-M4 for reproducibility -----------------
log("=== fitting province M0-M4 ===")
prov_dt <- prov_risk[, .(species, province, year, event,
                          climate_z = temp_grad_z,
                          effort_z  = log_effort_visits_z)]
prov_dt <- prov_dt[!is.na(climate_z) & !is.na(effort_z) & !is.na(event)]
log("province dt: ", nrow(prov_dt), " rows, ", sum(prov_dt$event), " events")

# ---- 5. Generic fit & extract --------------------------------------------
fit_M_ladder <- function(dt, scale_label, re_term) {
  log("--- ", scale_label, " | rows=", nrow(dt), " events=", sum(dt$event))
  fits <- list()
  forms <- list(
    M0 = sprintf("event ~ 1 + %s", re_term),
    M1 = sprintf("event ~ effort_z + %s", re_term),
    M2 = sprintf("event ~ climate_z + %s", re_term),
    M3 = sprintf("event ~ climate_z + effort_z + %s", re_term),
    M4 = sprintf("event ~ climate_z * effort_z + %s", re_term))
  for (nm in names(forms)) {
    t0 <- Sys.time()
    fit <- tryCatch(
      glmmTMB::glmmTMB(as.formula(forms[[nm]]),
                        data = dt,
                        family = binomial(link = "cloglog")),
      error = function(e) { log("   ", nm, " FAILED: ", conditionMessage(e)); NULL })
    if (!is.null(fit)) {
      log(sprintf("   %s OK (%.1fs, AIC=%.1f)", nm,
                   as.numeric(difftime(Sys.time(), t0, units = "secs")),
                   AIC(fit)))
      fits[[nm]] <- fit
    }
    invisible(gc(verbose = FALSE))
  }
  fits
}

extract_one <- function(fits, scale_label) {
  out <- list()
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    cf <- tryCatch(glmmTMB::fixef(fit)$cond, error = function(e) NULL)
    se <- tryCatch(sqrt(diag(stats::vcov(fit)$cond)), error = function(e) NULL)
    if (is.null(cf)) next
    for (tm in names(cf)) {
      i <- match(tm, names(cf))
      beta <- cf[i]; sei <- if (!is.null(se)) se[i] else NA
      out[[length(out) + 1L]] <- data.table(
        scale = scale_label, model = nm, term = tm,
        beta = beta, se = sei,
        hr = exp(beta),
        hr.low = exp(beta - 1.96 * sei),
        hr.high = exp(beta + 1.96 * sei),
        p.value = 2 * pnorm(-abs(beta / sei)),
        AIC = AIC(fit), n_rows = nobs(fit))
    }
  }
  rbindlist(out, fill = TRUE)
}

fits_prov <- fit_M_ladder(prov_dt, "province",
                           "(1|species) + (1|province)")

# ---- 6. Build admin-scale risk + fit ---------------------------------------
build_admin_risk <- function(scale_cn, scale_label) {
  log("=== ", scale_label, " (", scale_cn, ") ===")
  shp_candidates <- list.files(SHP_DIR, pattern = paste0("^", scale_cn, ".*\\.shp$"),
                                 full.names = TRUE)
  shp_candidates <- shp_candidates[!grepl("境界|线", shp_candidates)]
  if (length(shp_candidates) == 0L) { log("  no shp"); return(NULL) }
  shp_path <- shp_candidates[1]
  log("  shp: ", basename(shp_path))
  admin <- sf::st_read(shp_path, quiet = TRUE) |>
    sf::st_make_valid() |>
    sf::st_transform(4326)
  admin$unit_id <- paste0(scale_label, "_", seq_len(nrow(admin)))
  log("  ", nrow(admin), " ", scale_label, " polygons")

  # Province assignment via centroid join with 省.shp
  prov_shp <- list.files(SHP_DIR, pattern = "^省.*\\.shp$", full.names = TRUE)
  prov_shp <- prov_shp[!grepl("境界|线", prov_shp)][1]
  prov_sf  <- sf::st_read(prov_shp, quiet = TRUE) |>
    sf::st_make_valid() |>
    sf::st_transform(4326)
  prov_name_col <- intersect(c("name","NAME","省","省名","NAME_1"), names(prov_sf))[1]
  if (is.na(prov_name_col)) {
    prov_name_col <- names(prov_sf)[vapply(prov_sf, function(x)
      is.character(x) && any(grepl("[一-龥]", x)), logical(1))][1]
  }
  prov_sf$province <- unname(PROV_CN_EN[as.character(prov_sf[[prov_name_col]])])
  prov_sf$province[is.na(prov_sf$province)] <- as.character(prov_sf[[prov_name_col]])[is.na(prov_sf$province)]
  cen <- sf::st_centroid(admin)
  prov_join <- sf::st_join(cen, prov_sf[, "province"], join = sf::st_within)
  admin$province <- prov_join$province
  admin$province[is.na(admin$province)] <- "Unknown"

  # Map events to admin. ev_sf already carries a `province` column from v1
  # events_100km_grid_assigned.csv; drop it before the join so the admin-
  # derived province is used. 避免列名冲突。
  ev_sf_for_join <- ev_sf
  ev_sf_for_join$province <- NULL
  ev_admin <- sf::st_join(ev_sf_for_join, admin[, c("unit_id", "province")],
                           join = sf::st_within)
  ev_admin_dt <- as.data.table(sf::st_drop_geometry(ev_admin))
  ev_admin_dt <- ev_admin_dt[!is.na(unit_id)]
  if (nrow(ev_admin_dt) == 0L) { log("  no events in ", scale_label); return(NULL) }
  ev_admin_dt[is.na(province), province := "Unknown"]
  log("  events mapped: ", nrow(ev_admin_dt),
      " | provinces: ", uniqueN(ev_admin_dt$province))

  # First arrival per (species, unit_id)
  fa <- ev_admin_dt[, .(arrival_year = min(year, na.rm = TRUE)),
                    by = .(species, unit_id, province)]

  # Effort allocation by record share
  cnt <- ev_admin_dt[, .(n_events_unit = .N),
                     by = .(unit_id, year, province)]
  prov_total <- cnt[, .(n_events_prov = sum(n_events_unit)),
                    by = .(province, year)]
  cnt <- merge(cnt, prov_total, by = c("province", "year"))
  cnt[, share := n_events_unit / pmax(n_events_prov, 1)]
  cnt <- merge(cnt,
                prov_eff[, .(province, year, n_visits)],
                by = c("province", "year"), all.x = TRUE)
  cnt[is.na(n_visits), n_visits := 1L]
  cnt[, allocated_visits := share * n_visits]
  cnt[, log_effort_visits_z := zify(log1p(allocated_visits)), by = year]

  # SDM-candidate-restricted cartesian (species × unit × year)
  unit_in_prov <- unique(as.data.table(sf::st_drop_geometry(admin))[
    , .(unit_id, province)])
  cand <- merge(sdm_candidates, unit_in_prov, by = "province",
                 allow.cartesian = TRUE)
  log("  SDM-restricted (species × ", scale_label, ") pairs: ", nrow(cand))

  # Cartesian with years
  yrs <- 2002:2024
  risk <- cand[rep(seq_len(.N), each = length(yrs))]
  risk[, year := rep(yrs, times = nrow(cand))]
  risk <- merge(risk, fa[, .(species, unit_id, arrival_year)],
                 by = c("species", "unit_id"), all.x = TRUE)
  risk <- risk[is.na(arrival_year) | year <= arrival_year]
  risk[, event := as.integer(year == arrival_year)]
  risk[is.na(event), event := 0L]

  # Attach effort
  risk <- merge(risk, cnt[, .(unit_id, year, log_effort_visits_z,
                               allocated_visits)],
                 by = c("unit_id", "year"), all.x = TRUE)
  risk[is.na(log_effort_visits_z), log_effort_visits_z := 0]

  # Attach province × year temp_grad. v1 climate panel uses
  # `temp_grad_prov_z` (province-level z-score); we rename to temp_grad_z
  # so the model formula matches. 列名兼容。
  pc <- copy(prov_clim[, .(province, year, temp_grad_prov_z)])
  setnames(pc, "temp_grad_prov_z", "temp_grad_z")
  risk <- merge(risk, pc, by = c("province", "year"), all.x = TRUE)

  log("  risk set: rows=", nrow(risk), " events=", sum(risk$event),
      " unique units=", uniqueN(risk$unit_id))
  list(risk = risk, n_units = nrow(admin),
       n_units_with_events = uniqueN(ev_admin_dt$unit_id))
}

# ---- 7. Run prefecture + county -------------------------------------------
pref <- build_admin_risk("市", "prefecture")
fits_pref <- if (!is.null(pref)) {
  pref$risk[, climate_z := temp_grad_z]
  pref$risk[, effort_z  := log_effort_visits_z]
  fit_M_ladder(pref$risk[, .(species, unit_id, year, event,
                              climate_z, effort_z)],
                "prefecture",
                "(1|species) + (1|unit_id)")
} else NULL

invisible(gc(verbose = FALSE))

cnty <- build_admin_risk("县", "county")
fits_cnty <- if (!is.null(cnty)) {
  cnty$risk[, climate_z := temp_grad_z]
  cnty$risk[, effort_z  := log_effort_visits_z]
  fit_M_ladder(cnty$risk[, .(species, unit_id, year, event,
                               climate_z, effort_z)],
                "county",
                "(1|species) + (1|unit_id)")
} else NULL

# ---- 8. Collect + persist -------------------------------------------------
all_results <- rbindlist(list(
  extract_one(fits_prov, "province"),
  extract_one(fits_pref, "prefecture"),
  extract_one(fits_cnty, "county")), fill = TRUE)

fwrite(all_results,
       file.path(V2, "results", "tables",
                  "table_three_scale_coefficients.csv"))
log("wrote table_three_scale_coefficients.csv (", nrow(all_results), " rows)")

# AIC ladder
aic_ladder <- unique(all_results[, .(scale, model, AIC, n_rows)])
aic_ladder[, dAIC := AIC - min(AIC, na.rm = TRUE), by = scale]
setorder(aic_ladder, scale, AIC)
fwrite(aic_ladder,
       file.path(V2, "results", "tables",
                  "table_three_scale_model_comparison.csv"))

# Risk-set summary
risk_summary <- data.table(
  scale = c("province", "prefecture", "county"),
  n_units = c(uniqueN(prov_risk$province),
              if (!is.null(pref)) pref$n_units else NA,
              if (!is.null(cnty)) cnty$n_units else NA),
  n_units_with_events = c(uniqueN(prov_risk[event == 1L, province]),
                           if (!is.null(pref)) pref$n_units_with_events else NA,
                           if (!is.null(cnty)) cnty$n_units_with_events else NA),
  n_rows_in_M4 = c(
    if ("M4" %in% names(fits_prov)) nobs(fits_prov$M4) else NA,
    if (!is.null(fits_pref) && "M4" %in% names(fits_pref)) nobs(fits_pref$M4) else NA,
    if (!is.null(fits_cnty) && "M4" %in% names(fits_cnty)) nobs(fits_cnty$M4) else NA),
  n_events_in_M4 = c(
    sum(prov_risk$event),
    if (!is.null(pref)) sum(pref$risk$event) else NA,
    if (!is.null(cnty)) sum(cnty$risk$event) else NA))
fwrite(risk_summary,
       file.path(V2, "results", "tables",
                  "risk_set_summary_three_scale.csv"))
log("wrote risk_set_summary_three_scale.csv")

# ---- 9. Interaction-term forest plot --------------------------------------
m4_int <- all_results[model == "M4" & grepl(":", term)]
if (nrow(m4_int) > 0L) {
  m4_int[, scale := factor(scale, levels = c("county","prefecture","province"))]
  p_fig <- ggplot(m4_int, aes(x = hr, y = scale)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                    height = 0.2, linewidth = 0.5) +
    geom_point(size = 2.8, colour = "#B40426") +
    geom_text(aes(label = sprintf("HR = %.2f (%.2f, %.2f)", hr, hr.low, hr.high)),
               nudge_y = 0.18, size = 2.8) +
    scale_x_continuous(trans = "log", breaks = c(0.7, 1, 1.3, 1.6, 2)) +
    labs(title = "Climate × effort interaction across 3 administrative scales",
          subtitle = paste0("cloglog M4 with (1|species) + (1|unit). ",
                             "Climate = temp_grad_z; effort = log_effort_visits_z. ",
                             "Province from v1; prefecture/county from coord-level events ",
                             "with record-share effort allocation."),
          x = "Hazard ratio (log scale)", y = NULL) +
    theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(V2, "figures", "main",
                    "fig_three_scale_interaction_forest.pdf"),
         p_fig, width = 16, height = 8, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(file.path(V2, "figures", "main",
                    "fig_three_scale_interaction_forest.png"),
         p_fig, width = 16, height = 8, units = "cm", dpi = 600)
  log("wrote fig_three_scale_interaction_forest.{pdf,png}")
}

log("=== DONE ===")
log("M4 interaction rows:")
print(m4_int[, .(scale, term, beta = round(beta, 3), hr = round(hr, 3),
                 hr.low = round(hr.low, 3), hr.high = round(hr.high, 3),
                 p.value = signif(p.value, 3))])
log("AIC ladder:")
print(aic_ladder)
log("Risk-set summary:")
print(risk_summary)
