# ============================================================
# Scientific question / 科学问题:
#   v1 produced multi-scale future PREDICTIONS at prefecture (市) and
#   county (县) level (code/16_multi_scale_future_prediction.R) but
#   never FITTED hazard models at those scales — prefecture / county
#   predictions used the province-level fitted XGBoost. The reviewer's
#   reasonable critique: do prefecture- and county-resolved hazard
#   models give the same climate × effort interaction, or does the
#   province-level fit hide spatial heterogeneity?
#   v1 仅在市/县级别做了预测，并未拟合该尺度的 hazard 模型。本脚本
#   构造市级/县级 (species, unit, year) 风险集并拟合 M1–M5，与省级 /
#   网格级结果直接对照。
#
# Objective / 分析目标:
#   - Build prefecture-level risk set: spatially join event points to
#     prefecture polygons (data/spatial/basemap_GS2019_1822/市.shp), then
#     restrict the cartesian to (species, prefecture) candidate pairs
#     where the species' SDM-suitable province contains that prefecture.
#   - Build county-level risk set the same way (县.shp).
#   - Aggregate climate + effort to (unit × year) — for now we
#     inherit province-level (prefecture inherits its province's
#     climate/effort), AND additionally compute coordinate-density-based
#     within-prefecture effort variation so the interaction has signal.
#   - Fit M1–M5 at both scales, save coef tables.
#
# Input data / 输入数据:
#   data/raw/events_*_grid_assigned.csv      (coordinate-level events)
#   data/spatial/basemap_GS2019_1822/市.shp   (prefecture polygons)
#   data/spatial/basemap_GS2019_1822/县.shp   (county polygons)
#   data/raw/climate_metrics_province_year.csv (province climate panel)
#   data/raw/effort_panel_upgraded.csv         (province effort panel)
#
# Main workflow / 主要流程:
#   1. Load events with coordinates; spatially join to prefecture/county.
#   2. Derive (unit × year) effort: real coordinate-density + province ratio.
#   3. Build SDM-threshold-aligned cartesian (species × unit × year).
#   4. Fit M1–M5 with `(1|species) + (1|unit)` random effects.
#   5. Write coefficient tables + summary.
#
# Expected output / 预期输出:
#   data/derived/events_prefecture_risk_set.parquet
#   data/derived/events_county_risk_set.parquet
#   results/sensitivity/table_prefecture_county_coefs.csv
#   results/sensitivity/table_prefecture_county_aic.csv
#
# Key assumptions / 关键假设:
#   - Province ↔ prefecture ↔ county mapping is recoverable from the
#     shapefile attributes (省名 / 市名 / 县名). If the shp encoding is
#     in Chinese, we translate via the table embedded in v1
#     code/16_multi_scale_future_prediction.R lines 142-155.
#   - Climate × effort interaction term is the same form as in M4 at
#     province scale (cloglog hazard with crossed REs).
#
# Main packages / 主要包: sf, data.table, glmmTMB, arrow, broom.mixed,
#   ggplot2, glue.
# Output directory / 输出路径: data/derived/, results/sensitivity/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(arrow)
  library(glmmTMB)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_spatial.R"))
source(file.path("code", "utils", "utils_models.R"))

set.seed(42)

CFG <- list(
  year_min = 2002,
  year_max = 2024,
  scales   = c("prefecture", "county")
)

# Province name translation (subset; extended via union of all observed names)
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
  "澳门特别行政区" = "Macau"
)

translate_prov <- function(x) {
  out <- PROV_CN_EN[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

# ---- 1. Helper to read shp + unify column names ----------------------------
read_admin <- function(scale) {
  bm_dir <- path_spatial("basemap_GS2019_1822")
  candidates <- if (scale == "prefecture")
    c("市.shp", "市（等积投影）.shp", "prefecture.shp")
  else
    c("县.shp", "县（等积投影）.shp", "county.shp")
  hit <- candidates[file.exists(file.path(bm_dir, candidates))][1]
  if (is.na(hit)) stop("[38] missing ", scale, " shp in ", bm_dir)
  sf_obj <- sf::st_read(file.path(bm_dir, hit), quiet = TRUE)
  setnames_safe <- function(x, old, new) {
    if (old %in% names(x)) names(x)[names(x) == old] <- new
    x
  }
  for (cand in c("省", "省名", "PROVINCE", "Province")) {
    sf_obj <- setnames_safe(sf_obj, cand, "province_cn")
  }
  for (cand in c("市", "市名", "地名", "PREFECTURE", "Prefecture")) {
    if (scale == "prefecture") sf_obj <- setnames_safe(sf_obj, cand, "unit_cn")
  }
  for (cand in c("县", "县名", "区名", "COUNTY", "County")) {
    if (scale == "county") sf_obj <- setnames_safe(sf_obj, cand, "unit_cn")
  }
  if (!"unit_cn" %in% names(sf_obj)) {
    # Fall back: take first text-y column. 退而求其次。
    text_cols <- names(sf_obj)[vapply(sf_obj,
      function(x) is.character(x) || is.factor(x), logical(1))]
    text_cols <- setdiff(text_cols, c("province_cn", "geometry"))
    if (length(text_cols) > 0L) names(sf_obj)[names(sf_obj) == text_cols[1]] <- "unit_cn"
  }
  if (!"province_cn" %in% names(sf_obj)) {
    stop("[38] cannot find province column in ", scale, " shp.")
  }
  sf_obj$province <- translate_prov(sf_obj$province_cn)
  sf_obj$unit_id  <- paste0(substr(scale, 1, 4), "_", seq_len(nrow(sf_obj)))
  sf_obj
}

# ---- 2. Load coordinate events ---------------------------------------------
load_events_with_coords <- function() {
  p <- path_raw("events_100km_grid_assigned.csv")
  if (!file.exists(p)) stop("[38] cannot find coordinate events at ", p)
  dt <- fread(p, encoding = "UTF-8")
  setnames(dt, tolower(names(dt)))
  if (!"year" %in% names(dt) && "pub_year" %in% names(dt))
    setnames(dt, "pub_year", "year")
  dt[, year := as.integer(year)]
  dt <- dt[year >= CFG$year_min & year <= CFG$year_max &
            !is.na(longitude) & !is.na(latitude)]
  sf::st_as_sf(dt, coords = c("longitude", "latitude"),
                crs = CRS_WGS84, remove = FALSE)
}

# ---- 3. Build prefecture/county risk set ----------------------------------
build_scale_risk <- function(scale) {
  message(glue("\n[38] ===== {scale} ====="))
  admin <- read_admin(scale)
  admin <- to_albers(admin)
  events_sf <- load_events_with_coords()
  events_sf <- to_albers(events_sf)
  joined <- sf::st_join(events_sf, admin[, c("unit_id", "unit_cn", "province")],
                         join = sf::st_within)
  hit <- !is.na(joined$unit_id)
  message(glue("[38] {scale}: {sum(hit)}/{nrow(joined)} events assigned"))
  events_dt <- as.data.table(sf::st_drop_geometry(joined[hit, ]))

  # First arrival per (species, unit)
  first_arrival <- events_dt[, .(arrival_year = min(year, na.rm = TRUE)),
                             by = .(species, unit_id, province)]

  # Province climate panel + effort panel (z-scored)
  climate <- fread(path_raw("climate_metrics_province_year.csv"),
                    encoding = "UTF-8")
  effort  <- fread(path_raw("effort_panel_upgraded.csv"), encoding = "UTF-8")
  setnames(climate, tolower(names(climate)))
  setnames(effort,  tolower(names(effort)))

  # ---- Real within-province effort variation by event density ------------
  # Count events per (unit, year). 单元年内记录密度。
  cnt <- events_dt[, .(n_records_unit = .N), by = .(unit_id, year, province)]
  prov_year_tot <- cnt[, .(n_records_prov = sum(n_records_unit)),
                       by = .(province, year)]
  cnt <- merge(cnt, prov_year_tot, by = c("province", "year"))
  cnt[, share_in_prov := n_records_unit / pmax(n_records_prov, 1)]

  # Province totals from v1 panel
  cnt <- merge(cnt, effort[, .(province, year,
                                n_visits, n_birding_days, n_observers)],
               by = c("province", "year"), all.x = TRUE)
  cnt[, n_visits_unit       := share_in_prov * n_visits]
  cnt[, n_birding_days_unit := share_in_prov * n_birding_days]
  cnt[, n_observers_unit    := share_in_prov * n_observers]

  zify <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }
  cnt[, log_n_visits_z       := zify(log1p(n_visits_unit))]
  cnt[, log_n_birding_days_z := zify(log1p(n_birding_days_unit))]
  cnt[, log_n_observers_z    := zify(log1p(n_observers_unit))]
  cnt[, effort_pc1_z := zify(log1p(n_visits_unit) +
                              log1p(n_birding_days_unit) +
                              log1p(n_observers_unit))]

  # ---- Build SDM-threshold-aligned cartesian -----------------------------
  sdm_rs <- path_raw("grid_100km_risk_set.csv")
  if (file.exists(sdm_rs)) {
    cand <- unique(fread(sdm_rs, select = c("species", "province"),
                          encoding = "UTF-8"))
  } else {
    cand <- unique(events_dt[, .(species, province)])
  }
  unit_by_prov <- unique(as.data.table(sf::st_drop_geometry(admin))[, .(unit_id, province)])
  sp_unit <- merge(cand, unit_by_prov, by = "province",
                    allow.cartesian = TRUE)
  risk <- merge(sp_unit, data.table::CJ(unit_id = unique(sp_unit$unit_id),
                                         year    = CFG$year_min:CFG$year_max),
                 by = "unit_id", allow.cartesian = TRUE)
  risk <- merge(risk, sp_unit, by = c("species", "unit_id", "province"))
  risk <- merge(risk, first_arrival[, .(species, unit_id, arrival_year)],
                 by = c("species", "unit_id"), all.x = TRUE)
  risk <- risk[is.na(arrival_year) | year <= arrival_year]
  risk[, event := as.integer(year == arrival_year)]
  risk[is.na(event), event := 0L]

  # Attach effort + climate
  risk <- merge(risk, cnt[, .(unit_id, year,
                               log_n_visits_z, log_n_birding_days_z,
                               log_n_observers_z, effort_pc1_z,
                               n_visits_unit, n_birding_days_unit,
                               n_observers_unit)],
                 by = c("unit_id", "year"), all.x = TRUE)
  for (col in c("log_n_visits_z", "log_n_birding_days_z",
                 "log_n_observers_z", "effort_pc1_z")) {
    risk[is.na(get(col)), (col) := 0]
  }
  risk <- merge(risk, climate, by = c("province", "year"), all.x = TRUE)
  if (!"climate_velocity_z" %in% names(risk) && "climate_velocity" %in% names(risk)) {
    risk[, climate_velocity_z := zify(climate_velocity)]
  }

  out_path <- path_derived(glue("events_{scale}_risk_set.parquet"))
  ensure_dir(dirname(out_path))
  arrow::write_parquet(risk, out_path, compression = "snappy")
  message(glue("[38] {scale}: wrote {out_path} ({nrow(risk)} rows, ",
                "{sum(risk$event)} events)"))

  # ---- Fit M1–M5 ---------------------------------------------------------
  re_term <- "(1|species) + (1|unit_id)"
  clim_col <- intersect(c("climate_velocity_z", "temp_anom_z"), names(risk))[1]
  if (is.na(clim_col)) {
    warning(glue("[38] {scale}: no climate_z column — skipping models."))
    return(NULL)
  }
  setnames(risk, clim_col, "climate_z")
  risk[, effort_z := log_n_visits_z]

  forms <- list(
    M0 = sprintf("event ~ 1 + %s", re_term),
    M1 = sprintf("event ~ effort_z + %s", re_term),
    M2 = sprintf("event ~ climate_z + %s", re_term),
    M3 = sprintf("event ~ climate_z + effort_z + %s", re_term),
    M4 = sprintf("event ~ climate_z * effort_z + %s", re_term)
  )
  models <- list()
  for (nm in names(forms)) {
    fit <- tryCatch(
      glmmTMB::glmmTMB(as.formula(forms[[nm]]),
                       data = risk[!is.na(climate_z) & !is.na(effort_z)],
                       family = binomial(link = "cloglog")),
      error = function(e) { message("  ", nm, " failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(fit)) models[[nm]] <- fit
  }

  aic_dt  <- aic_table(models)
  aic_dt[, scale := scale]
  coef_dt <- data.table::rbindlist(lapply(names(models), function(nm) {
    out <- extract_coefs_manual(models[[nm]])
    out[, model := nm]; out[, scale := scale]
    out
  }), use.names = TRUE, fill = TRUE)
  list(aic = aic_dt, coef = coef_dt)
}

results <- lapply(CFG$scales, build_scale_risk)
aic_all  <- data.table::rbindlist(lapply(results, `[[`, "aic"), fill = TRUE)
coef_all <- data.table::rbindlist(lapply(results, `[[`, "coef"), fill = TRUE)

ensure_dir(path_sensitivity())
fwrite(aic_all,  path_sensitivity("table_prefecture_county_aic.csv"))
fwrite(coef_all, path_sensitivity("table_prefecture_county_coefs.csv"))
message("[38] AIC across prefecture / county:")
print(aic_all)

dump_session_info(path_logs("38_prefecture_county_hazard_sessionInfo.txt"))
message("[38] done.")
