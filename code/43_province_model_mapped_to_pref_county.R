# ============================================================
# Scientific question / 科学问题:
#   An alternative way to obtain prefecture- and county-level future
#   hazard maps WITHOUT refitting at those scales: take the headline
#   PROVINCE-scale M4 (Spec B headline, glmmTMB; or full-feature XGBoost)
#   and plug in unit-specific covariates evaluated at prefecture / county
#   resolution. The unit gets the climate it actually has + the effort
#   it actually has (raw-records based), so the projection reflects
#   unit-level conditions even though the model parameters come from
#   the province scale.
#   备选思路：不在市/县级 refit，直接用省级 M4 + 市/县级 covariates 做
#   plug-in 预测，独立于本目录 §3-bis 的 prefecture/county refit。
#
# Inputs:
#   v1/data/hazard_risk_upgraded_complete_case.csv  (province risk set)
#   v1/data/climate_metrics_province_year.csv       (province × year climate)
#   v2/data/derived/unit_climate_{prefecture,county}.csv (WorldClim 10')
#   v2/data/derived/unit_effort_{prefecture,county}.csv  (raw-records effort)
#   v2 GS(2019)1822 basemap shapefiles
#
# Workflow:
#   1. Refit headline province M4 (cloglog, Spec B headline) AND a
#      province-scale XGBoost on the full feature set.
#   2. Build prefecture / county baseline rows at year = 2024 using
#      unit-native covariates (temp_grad_z from province × year climate
#      panel; effort z from raw-records derived panel).
#   3. SSP-perturb temp_grad_z forward (SSP245 +0.3 SD/dec, SSP585 +0.8).
#   4. Plug-in prediction (re.form = NA for glmmTMB to marginalise random
#      effects; XGBoost on the same feature vector).
#   5. Render 4 figures (glmmTMB + XGBoost × prefecture + county), each
#      6 panels (SSP × year).
#
# Outputs:
#   results/forecasts/table_prefecture_future_mapped_from_province_glmmTMB.csv
#   results/forecasts/table_prefecture_future_mapped_from_province_xgboost.csv
#   results/forecasts/table_county_future_mapped_from_province_glmmTMB.csv
#   results/forecasts/table_county_future_mapped_from_province_xgboost.csv
#   figures/main/fig_future_mapped_glmmTMB_prefecture.{pdf,png}
#   figures/main/fig_future_mapped_xgboost_prefecture.{pdf,png}
#   figures/main/fig_future_mapped_glmmTMB_county.{pdf,png}
#   figures/main/fig_future_mapped_xgboost_county.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(glmmTMB)
  library(xgboost)
  library(ggplot2)
  library(patchwork)
})
sf::sf_use_s2(FALSE)
options(warn = 1)
set.seed(42)

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- Sys.getenv("V1_ROOT",
                  normalizePath(file.path(V2, "..",
                                          "bird_hazard_model_effort_upgrade"),
                                 mustWork = FALSE))
SHP <- file.path(V2, "data", "spatial", "basemap_GS2019_1822")

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ens(file.path(V2, "results", "forecasts"))
ens(file.path(V2, "figures", "main"))

log <- function(...) cat(sprintf("[43 %s] ", format(Sys.time(), "%H:%M:%S")),
                          ..., "\n", sep = "")

theme_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
          panel.border     = element_rect(linewidth = 0.4, colour = "grey20"),
          plot.title       = element_text(face = "bold", size = base_size + 1),
          plot.subtitle    = element_text(size = base_size - 1, colour = "grey30"))
}
save_pub <- function(p, name, width = 22, height = 14,
                      path = file.path(V2, "figures", "main")) {
  ggsave(file.path(path, paste0(name, ".pdf")), p, width = width,
         height = height, units = "cm", device = grDevices::cairo_pdf)
  ggsave(file.path(path, paste0(name, ".png")), p, width = width,
         height = height, units = "cm", dpi = 600)
  log("wrote ", name, ".{pdf,png}")
}
zify <- function(x) {
  s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

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

# ============================================================
# 1. Refit headline PROVINCE M4 (Spec B) + province XGBoost
# ============================================================
log("loading province risk set + climate panel")
risk_prov <- fread(file.path(V1, "data",
                              "hazard_risk_upgraded_complete_case.csv"),
                    encoding = "UTF-8")
prov_clim_full <- fread(file.path(V1, "data",
                                    "climate_metrics_province_year.csv"),
                          encoding = "UTF-8")
prov_eff <- fread(file.path(V1, "data", "effort_panel_upgraded.csv"),
                    encoding = "UTF-8")

# 1a — Refit glmmTMB province M4 (Spec B headline)
log("fitting province glmmTMB M4 (Spec B headline)")
m4_dt <- risk_prov[, .(species, province, year, event,
                        temp_grad_z, log_effort_visits_z)]
m4_dt <- m4_dt[complete.cases(m4_dt)]
fit_prov_M4 <- glmmTMB::glmmTMB(
  event ~ temp_grad_z * log_effort_visits_z +
    (1 | species) + (1 | province),
  data = m4_dt, family = binomial(link = "cloglog"))
log("province M4 AIC = ", round(AIC(fit_prov_M4), 2))

# 1b — Train province XGBoost on full feature set (climate + effort
# joined from climate_metrics_province_year + effort_panel_upgraded).
climate_join_cols <- intersect(c("climate_velocity_z","precip_velocity_z",
                                   "climate_exposure_z","warming_rate_z",
                                   "mahalanobis_dist_z","temp_anom_z",
                                   "prec_anom_z","temp_grad_prov_z",
                                   "prec_grad_prov_z"),
                                 names(prov_clim_full))
train <- merge(risk_prov,
                prov_clim_full[, c("province","year", climate_join_cols),
                                with = FALSE],
                by = c("province","year"), all.x = TRUE)
train[, temp_x_effort := temp_grad_z * log_effort_visits_z]
train[, mahal_x_effort := mahalanobis_dist_z * log_effort_visits_z]
train[, year_c := year - 2013]
feat <- intersect(c("temp_grad_z","climate_velocity_z","warming_rate_z",
                     "mahalanobis_dist_z","temp_anom_z","prec_anom_z",
                     "log_effort_visits_z","log_effort_days_z",
                     "effort_pc1_z","temp_x_effort","mahal_x_effort","year_c"),
                   names(train))
train_dt <- train[, c("event", feat), with = FALSE]
train_dt <- train_dt[complete.cases(train_dt)]
log("XGB province train: rows = ", nrow(train_dt), " | features = ",
    length(feat))

Xtr <- xgb.DMatrix(as.matrix(train_dt[, ..feat]), label = train_dt$event)
xgp <- list(objective="binary:logistic", eval_metric="auc",
            eta = 0.05, max_depth = 6,
            subsample = 0.8, colsample_bytree = 0.8,
            nthread = max(1L, parallel::detectCores()-2L))
cv <- xgb.cv(xgp, Xtr, nrounds = 600, nfold = 5,
              early_stopping_rounds = 30, verbose = 0)
nb <- cv$best_iteration
if (is.null(nb) || nb == 0L) {
  nb <- which.max(cv$evaluation_log$test_auc_mean)
}
fit_prov_xgb <- xgb.train(xgp, Xtr, nrounds = nb, verbose = 0)
log("province XGBoost: nrounds=", nb, " CV-AUC=",
    round(cv$evaluation_log$test_auc_mean[nb], 3))

# ============================================================
# 2. Build unit baseline covariates for prefecture + county
# ============================================================
log("loading unit-native climate + raw-records-based effort")

pref_clim <- fread(file.path(V2, "data", "derived",
                              "unit_climate_prefecture.csv"))
cnty_clim <- fread(file.path(V2, "data", "derived",
                              "unit_climate_county.csv"))
pref_eff  <- fread(file.path(V2, "data", "derived",
                              "unit_effort_prefecture.csv"))
cnty_eff  <- fread(file.path(V2, "data", "derived",
                              "unit_effort_county.csv"))
log("prefecture climate rows: ", nrow(pref_clim),
    "  | effort rows: ", nrow(pref_eff))
log("county     climate rows: ", nrow(cnty_clim),
    "  | effort rows: ", nrow(cnty_eff))

# Read polygons to attach unit→province mapping.
pref_sf <- st_read(file.path(SHP, "市（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4326)
cnty_sf <- st_read(file.path(SHP, "县（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4326)
pref_sf$pref_id <- paste0("PREF_", sprintf("%04d", seq_len(nrow(pref_sf))))
cnty_sf$cnty_id <- paste0("CNTY_", sprintf("%05d", seq_len(nrow(cnty_sf))))
pref_sf$province <- unname(PROV_CN_EN[as.character(pref_sf[["省"]])])
cnty_sf$province <- unname(PROV_CN_EN[as.character(cnty_sf[["省"]])])

unit_meta <- function(sf_unit, id_col) {
  d <- as.data.table(st_drop_geometry(sf_unit))
  setnames(d, id_col, "unit_id")
  d[, .(unit_id, province)]
}
pref_meta <- unit_meta(pref_sf, "pref_id")
cnty_meta <- unit_meta(cnty_sf, "cnty_id")

# IMPORTANT: the headline province M4 uses temp_grad_z (which in the
# v1 risk set is province × year). We translate this into a unit-level
# value at year = 2024 by using the province × 2024 temp_grad_prov_z
# (the same one each unit would inherit from its province). Effort is
# unit-native (from the raw-records based panel). 时间外推用 SSP 扰动。
prov_clim_2024 <- prov_clim_full[year == 2024,
                                   .(province, temp_grad_z = temp_grad_prov_z)]
prov_eff_2024  <- prov_eff[year == 2024]

build_unit_baseline <- function(unit_clim, unit_eff, meta) {
  # Try several name variations for unit id column
  id_col <- intersect(c("unit_id","pref_id","cnty_id"), names(unit_clim))[1]
  setnames(unit_clim, id_col, "unit_id")
  id_col_e <- intersect(c("unit_id","pref_id","cnty_id"), names(unit_eff))[1]
  setnames(unit_eff, id_col_e, "unit_id")

  base <- merge(unit_clim, meta, by = "unit_id")
  # Province × 2024 climate baseline temp_grad_z
  base <- merge(base, prov_clim_2024, by = "province", all.x = TRUE)
  # Unit-level effort at year=2024 from raw-records effort
  ue2024 <- unit_eff[year == 2024, .(unit_id, log_effort_visits_z)]
  base <- merge(base, ue2024, by = "unit_id", all.x = TRUE)
  base[is.na(log_effort_visits_z), log_effort_visits_z := 0]
  # Other ML features default to province × 2024 values for XGBoost path
  extra <- merge(meta, prov_clim_full[year == 2024], by = "province",
                  all.x = TRUE)
  extra <- merge(extra,
                  prov_eff[year == 2024,
                            .(province, log_effort_record_z,
                              log_effort_days_z, effort_pc1_z)],
                  by = "province", all.x = TRUE)
  base <- merge(base, extra[, .(unit_id,
                                  climate_velocity_z, warming_rate_z,
                                  mahalanobis_dist_z, temp_anom_z, prec_anom_z,
                                  log_effort_record_z, log_effort_days_z,
                                  effort_pc1_z)],
                 by = "unit_id", all.x = TRUE)
  base
}
pref_base <- build_unit_baseline(copy(pref_clim), copy(pref_eff), pref_meta)
cnty_base <- build_unit_baseline(copy(cnty_clim), copy(cnty_eff), cnty_meta)
log("prefecture baseline rows: ", nrow(pref_base),
    "  county baseline rows: ", nrow(cnty_base))

# ============================================================
# 3. Future scenario predictions (province model → unit)
# ============================================================
future_years <- c(2030, 2050, 2080)
ssp_eps <- list(SSP245 = 0.3 / sd(prov_clim_full$temp_grad_prov, na.rm=TRUE),
                 SSP585 = 0.8 / sd(prov_clim_full$temp_grad_prov, na.rm=TRUE))

mk_future_unit <- function(base, label) {
  fut <- CJ(unit_id = base$unit_id,
             ssp     = c("SSP245","SSP585"),
             year    = future_years)
  fut <- merge(fut, base, by = "unit_id")
  fut[, decades_ahead := (year - 2024) / 10]
  fut[, temp_grad_z := temp_grad_z + ssp_eps[[ssp[1]]] * decades_ahead,
      by = ssp]
  # glmmTMB prediction (re.form = NA → marginalise both random effects)
  fut[, hazard_glmm := predict(fit_prov_M4, newdata = fut,
                                 type = "response",
                                 re.form = NA, allow.new.levels = TRUE)]
  # XGBoost path: build the same feature matrix the train used
  fut[, temp_x_effort  := temp_grad_z * log_effort_visits_z]
  fut[, mahal_x_effort := mahalanobis_dist_z * log_effort_visits_z]
  fut[, year_c := year - 2013]
  Xpred <- as.matrix(fut[, ..feat])
  fut[, hazard_xgb := predict(fit_prov_xgb, Xpred)]
  log("plugged in ", label, ": ", nrow(fut), " future predictions ",
      "(", uniqueN(fut$unit_id), " units × ", length(future_years),
      " yrs × 2 ssps)")
  fut
}

fut_pref <- mk_future_unit(pref_base, "prefecture")
fut_cnty <- mk_future_unit(cnty_base, "county")
fwrite(fut_pref[, .(unit_id, province, ssp, year, log_effort_visits_z,
                     temp_grad_z, hazard_glmm, hazard_xgb)],
       file.path(V2, "results", "forecasts",
                  "table_prefecture_future_mapped_from_province.csv"))
fwrite(fut_cnty[, .(unit_id, province, ssp, year, log_effort_visits_z,
                     temp_grad_z, hazard_glmm, hazard_xgb)],
       file.path(V2, "results", "forecasts",
                  "table_county_future_mapped_from_province.csv"))

# ============================================================
# 4. Choropleth maps
# ============================================================
prov_alb <- st_transform(st_read(file.path(SHP, "省（等积投影）.shp"),
                                  quiet = TRUE), 4524)
ninedash_alb <- tryCatch(
  st_transform(st_read(file.path(SHP, "九段线.shp"), quiet = TRUE), 4524),
  error = function(e) NULL)
pref_alb <- st_transform(pref_sf, 4524)
cnty_alb <- st_transform(cnty_sf, 4524)

choropleth_unit <- function(sf_unit, dt, var, title, id_col) {
  joined <- merge(sf_unit, dt[, c(id_col, var), with = FALSE],
                   by = id_col)
  p <- ggplot() +
    geom_sf(data = joined, aes(fill = .data[[var]]), colour = NA) +
    geom_sf(data = prov_alb, fill = NA, colour = "grey20", linewidth = 0.25)
  if (!is.null(ninedash_alb)) p <- p +
    geom_sf(data = ninedash_alb, fill = NA, colour = "grey20",
            linewidth = 0.2)
  p + scale_fill_gradient(low = "#F7FBFF", high = "#B40426",
                           name = "Hazard\n(prob.)",
                           limits = range(dt[[var]], na.rm = TRUE)) +
    coord_sf(datum = NA, expand = FALSE) +
    labs(title = title) +
    theme_pub(base_size = 7.5) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

mk_six_panel <- function(sf_unit, fut, var, id_col, file_label,
                          model_label, scale_label) {
  setnames(fut, "unit_id", id_col)
  fut[, panel := paste0(ssp, " — ", year)]
  fut[, panel := factor(panel,
                          levels = c(paste0("SSP245 — ", future_years),
                                      paste0("SSP585 — ", future_years)))]
  panels <- lapply(split(fut, fut$panel), function(d)
    choropleth_unit(sf_unit, d, var, title = unique(d$panel), id_col))
  combined <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") +
    plot_annotation(
      title = sprintf("%s-scale future hazard — province %s model PLUG-IN",
                      scale_label, model_label),
      subtitle = paste0("Province-scale parameters applied to unit-level ",
                         "covariates (WorldClim 10' climate + raw-records ",
                         "effort). SSP-perturbed temp_grad_z 2024 + Δ."),
      theme = theme(plot.title = element_text(face="bold", size=10),
                     plot.subtitle = element_text(size=8, colour="grey30")))
  save_pub(combined, file_label, width = 22, height = 14)
}

mk_six_panel(pref_alb, copy(fut_pref), "hazard_glmm", "pref_id",
              "fig_future_mapped_glmmTMB_prefecture",
              "glmmTMB M4 (Spec B)", "Prefecture")
mk_six_panel(pref_alb, copy(fut_pref), "hazard_xgb", "pref_id",
              "fig_future_mapped_xgboost_prefecture",
              "XGBoost", "Prefecture")
mk_six_panel(cnty_alb, copy(fut_cnty), "hazard_glmm", "cnty_id",
              "fig_future_mapped_glmmTMB_county",
              "glmmTMB M4 (Spec B)", "County")
mk_six_panel(cnty_alb, copy(fut_cnty), "hazard_xgb", "cnty_id",
              "fig_future_mapped_xgboost_county",
              "XGBoost", "County")

log("=== DONE ===")
log("Forecast tables in results/forecasts/")
log("Choropleths in figures/main/fig_future_mapped_*.{pdf,png}")
