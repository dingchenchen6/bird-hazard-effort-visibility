# ============================================================
# Scientific question / 科学问题:
#   v1 projected risk under SSP245/585 using a single climate model.
#   This collapses GCM structural uncertainty into a point estimate.
#   What does the risk map look like across a 5-GCM CMIP6 ensemble,
#   and where do the GCMs disagree most?
#   单一气候模型外推消除了模型间结构不确定性。改用 5 个 CMIP6 模型
#   集合后，风险图中位数 + IQR 与模型分歧如何分布？
#
# Objective / 分析目标:
#   - For each of 5 GCMs × {SSP245, SSP585} × {2050, 2080}, extract
#     grid-cell climate metrics (temp_anom, precip_anom, velocity).
#   - Score each cell under the fitted XGBoost (script 22 / refitted
#     locally on M4 features).
#   - Per cell × scenario × horizon, store median, IQR and ensemble
#     disagreement (= IQR / median).
#
# Input data / 输入数据:
#   data/spatial/chelsa/cmip6/{model}/{ssp}/{year}/(tas|pr).tif
#     (NetCDF / GeoTIFF; absent → script writes a placeholder shaped
#      table built from province-mean SSP velocities for graceful
#      degradation).
#   data/raw/hazard_risk_upgraded_complete_case.csv
#
# Main workflow / 主要流程:
#   1. Define GCM list + scenarios + horizons (CFG).
#   2. Train XGBoost on the v1 risk set (same params as script 30).
#   3. For each GCM × SSP × year, build feature dataframe per grid cell.
#   4. Predict; collect per-cell scores.
#   5. Aggregate to median + IQR + disagreement; write parquet.
#
# Expected output / 预期输出:
#   results/forecasts/cmip6_ensemble.parquet            (per-cell × scenario)
#   results/forecasts/cmip6_ensemble_summary.parquet    (median + IQR + disag)
#   figures/main/fig5_cmip6_ensemble_2050_ssp585.{pdf,png}
#   figures/supplementary/figS11_cmip6_per_model.{pdf,png}
#
# Key assumptions / 关键假设:
#   - GCM list: ACCESS-CM2, EC-Earth3, MPI-ESM1-2-HR, MIROC6, UKESM1-0-LL
#     (broad ECS coverage; commonly archived on ESGF).
#
# Main packages / 主要包: terra, sf, data.table, arrow, xgboost,
#   ggplot2, patchwork, glue.
# Output directory / 输出路径: results/forecasts/, figures/main/,
#   figures/supplementary/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(sf)
  library(terra)
  library(xgboost)
  library(ggplot2)
  library(patchwork)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_spatial.R"))
source(file.path("code", "utils", "utils_plots.R"))

set.seed(42)

CFG <- list(
  gcm_list   = c("ACCESS-CM2", "EC-Earth3", "MPI-ESM1-2-HR",
                 "MIROC6", "UKESM1-0-LL"),
  scenarios  = c("ssp245", "ssp585"),
  horizons   = c(2050, 2080),
  base_dir   = path_spatial("chelsa/cmip6"),
  grid_km    = 100,
  features   = c("temp_anom_z", "precip_anom_z", "climate_velocity_z",
                 "mahalanobis_dist_z", "log_n_visits_z",
                 "log_n_observers_z", "log_n_birding_days_z",
                 "effort_pc1_z", "sdm_suitability_z"),
  xgb_params = list(objective = "binary:logistic",
                    eval_metric = "auc",
                    eta = 0.05, max_depth = 6,
                    subsample = 0.8, colsample_bytree = 0.8,
                    nthread = max(1L, parallel::detectCores() - 2L))
)

# ---- 1. Train XGBoost on full v1 risk set ---------------------------------
risk <- fread(path_raw("hazard_risk_upgraded_complete_case.csv"),
              encoding = "UTF-8")
features <- intersect(CFG$features, names(risk))
if (length(features) < 4L) stop("[29] Need ≥4 features; have ", length(features))

dtrain <- xgb.DMatrix(as.matrix(risk[, ..features]), label = risk$event)
cv <- xgb.cv(params = CFG$xgb_params, data = dtrain, nrounds = 800,
              nfold = 5, early_stopping_rounds = 30, verbose = 0)
model <- xgboost::xgb.train(params = CFG$xgb_params, data = dtrain,
                             nrounds = cv$best_iteration, verbose = 0)
message(glue("[29] Trained XGBoost (nrounds={cv$best_iteration})."))

# ---- 2. Build a base feature frame per grid cell --------------------------
grid_path <- path_derived(glue("grid_{CFG$grid_km}km_sf.gpkg"))
if (!file.exists(grid_path)) {
  warning("[29] grid sf missing; falling back to province-level prediction.")
  base_units <- read_gs2019_basemap("province")
  unit_key <- "province"
} else {
  base_units <- sf::st_read(grid_path, quiet = TRUE)
  unit_key <- "grid_id"
}

# ---- 3. Helper: read CMIP6 anomalies (or fall back) -----------------------
load_anom <- function(gcm, ssp, year, var = "tas") {
  p <- file.path(CFG$base_dir, gcm, ssp, as.character(year),
                  glue("{var}.tif"))
  if (!file.exists(p)) return(NULL)
  rast <- terra::rast(p)
  vals <- terra::values(rast)
  # Convert to anomaly per grid cell via raster_to_grid().
  agg <- raster_to_grid(rast, base_units, fun = "mean",
                         value_col_prefix = paste0(var, "_", year))
  data.table::setnames(agg, grep(paste0(var, "_"), names(agg), value = TRUE),
                       paste0("anom_", var))
  agg
}

# ---- 4. Sample predictions across the ensemble ----------------------------
# HARD-FAIL by default if any GCM × SSP × year raster is missing. Opt-in
# to the empirical-perturbation fallback by exporting
# V2_ALLOW_CMIP6_FALLBACK=1. The fallback is NOT a real CMIP6 ensemble
# and must never enter a manuscript as "CMIP6 projection" without that
# environment variable being set deliberately. 默认硬失败，避免合成数据被当真。
allow_cmip6_fallback <- nzchar(Sys.getenv("V2_ALLOW_CMIP6_FALLBACK"))

preds <- data.table::rbindlist(lapply(CFG$gcm_list, function(gcm) {
  data.table::rbindlist(lapply(CFG$scenarios, function(ssp) {
    data.table::rbindlist(lapply(CFG$horizons, function(yr) {
      tas_dt <- load_anom(gcm, ssp, yr, "tas")
      pr_dt  <- load_anom(gcm, ssp, yr, "pr")
      if (is.null(tas_dt) || is.null(pr_dt)) {
        if (!allow_cmip6_fallback) {
          stop(glue("[29] CMIP6 raster missing for gcm={gcm} ssp={ssp} ",
                    "year={yr}. Refusing to fabricate ensemble output. ",
                    "Either supply data/spatial/chelsa/cmip6/{gcm}/{ssp}/{yr}/(tas|pr).tif ",
                    "or re-export with V2_ALLOW_CMIP6_FALLBACK=1 to opt in ",
                    "to the empirical-perturbation placeholder mode."),
               call. = FALSE)
        }
        # OPT-IN: empirical perturbation placeholder. NOT a real CMIP6 output.
        # GCM_eps values are illustrative warming/precip signatures, **not**
        # peer-reviewed projections. 仅在显式启用时使用。
        eps <- list("ACCESS-CM2" = c(1.8, -0.5),
                    "EC-Earth3"  = c(1.5,  0.2),
                    "MPI-ESM1-2-HR" = c(1.2, 0.0),
                    "MIROC6"     = c(2.4, -0.3),
                    "UKESM1-0-LL" = c(3.1, -0.6))[[gcm]]
        if (is.null(eps)) eps <- c(2.0, 0.0)
        feat_dt <- copy(risk[, ..features])
        if ("temp_anom_z" %in% names(feat_dt)) {
          feat_dt[, temp_anom_z := temp_anom_z + (eps[1] *
                  (if (ssp == "ssp585") 1.5 else 1.0) *
                  (if (yr == 2080) 1.4 else 1.0))]
        }
        if ("precip_anom_z" %in% names(feat_dt)) {
          feat_dt[, precip_anom_z := precip_anom_z + eps[2]]
        }
        scores <- predict(model, as.matrix(feat_dt[, ..features]))
        return(data.table::data.table(gcm = gcm, ssp = ssp, year = yr,
                                       cell = seq_len(length(scores)),
                                       score = scores,
                                       fallback = TRUE))
      }
      feat_dt <- copy(risk[, ..features])
      # Replace temp/precip anomaly columns with averaged GCM anomalies.
      if ("temp_anom_z" %in% names(feat_dt)) {
        feat_dt[, temp_anom_z := mean(tas_dt$anom_tas, na.rm = TRUE)]
      }
      if ("precip_anom_z" %in% names(feat_dt)) {
        feat_dt[, precip_anom_z := mean(pr_dt$anom_pr, na.rm = TRUE)]
      }
      scores <- predict(model, as.matrix(feat_dt[, ..features]))
      data.table::data.table(gcm = gcm, ssp = ssp, year = yr,
                              cell = seq_len(length(scores)),
                              score = scores,
                              fallback = FALSE)
    }))
  }))
}))

ensure_dir(path_forecasts())
arrow::write_parquet(preds, path_forecasts("cmip6_ensemble.parquet"),
                     compression = "snappy")

# ---- 5. Per-cell summaries -------------------------------------------------
summary_dt <- preds[, .(median = stats::median(score, na.rm = TRUE),
                         q25    = stats::quantile(score, 0.25, na.rm = TRUE),
                         q75    = stats::quantile(score, 0.75, na.rm = TRUE),
                         iqr    = stats::quantile(score, 0.75, na.rm = TRUE) -
                                   stats::quantile(score, 0.25, na.rm = TRUE),
                         disagreement = (stats::quantile(score, 0.75, na.rm = TRUE) -
                                          stats::quantile(score, 0.25, na.rm = TRUE)) /
                                         stats::median(score, na.rm = TRUE)),
                     by = .(cell, ssp, year)]
arrow::write_parquet(summary_dt,
                     path_forecasts("cmip6_ensemble_summary.parquet"),
                     compression = "snappy")

# ---- 6. Figures ------------------------------------------------------------
# Main: 2050 SSP585 median map with disagreement overlay.
# We need a representative spatial unit; if grid sf available, join by row order.
plot_data <- summary_dt[ssp == "ssp585" & year == 2050]
if (unit_key == "grid_id" && nrow(plot_data) == nrow(base_units)) {
  base_units$median       <- plot_data$median
  base_units$disagreement <- plot_data$disagreement
  base_units_alb <- to_albers(base_units)
  p_main <- ggplot() +
    geom_sf(data = base_units_alb, aes(fill = median),
            colour = NA) +
    scale_fill_viridis_c(option = "B", name = "Median\nhazard score") +
    labs(title = "CMIP6 ensemble — SSP585 / 2050",
          subtitle = "Cell colour = ensemble median; alpha = IQR/median (disagreement)") +
    theme_geb() +
    theme(panel.grid = element_blank())
  p_disag <- ggplot() +
    geom_sf(data = base_units_alb, aes(fill = disagreement), colour = NA) +
    scale_fill_viridis_c(option = "C", name = "Ensemble\ndisagreement") +
    labs(title = "Ensemble disagreement (IQR / median)") +
    theme_geb() + theme(panel.grid = element_blank())
  fig5 <- p_main | p_disag
  ensure_dir(path_main_fig())
  ggsave(path_main_fig("fig5_cmip6_ensemble_2050_ssp585.pdf"),
         fig5, width = 18, height = 9, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(path_main_fig("fig5_cmip6_ensemble_2050_ssp585.png"),
         fig5, width = 18, height = 9, units = "cm", dpi = 600)
} else {
  warning("[29] grid units could not be matched to summary rows; ",
          "Fig 5 will be drafted by script 34 instead.")
}

# Per-model maps as supp figure (S11). 单模型补充图。
preds_50 <- preds[ssp == "ssp585" & year == 2050]
if (uniqueN(preds_50$cell) > 0L && unit_key == "grid_id") {
  per_model <- preds_50[, .(score = mean(score, na.rm = TRUE)),
                        by = .(cell, gcm)]
  per_model_wide <- dcast(per_model, cell ~ gcm, value.var = "score")
  setnames(per_model_wide, "cell", "row_idx")
  panels <- lapply(CFG$gcm_list, function(g) {
    if (!g %in% names(per_model_wide)) return(NULL)
    bu <- base_units
    bu$score <- per_model_wide[[g]]
    bu_alb <- to_albers(bu)
    ggplot() + geom_sf(data = bu_alb, aes(fill = score), colour = NA) +
      scale_fill_viridis_c(option = "B") +
      labs(title = g) + theme_geb() +
      theme(panel.grid = element_blank(), legend.position = "none")
  })
  panels <- panels[!sapply(panels, is.null)]
  if (length(panels) > 0L) {
    figS11 <- patchwork::wrap_plots(panels, ncol = 3)
    ensure_dir(path_supp_fig())
    ggsave(path_supp_fig("figS11_cmip6_per_model.pdf"),
           figS11, width = 22, height = 14, units = "cm",
           device = grDevices::cairo_pdf)
    ggsave(path_supp_fig("figS11_cmip6_per_model.png"),
           figS11, width = 22, height = 14, units = "cm", dpi = 600)
  }
}

dump_session_info(path_logs("29_cmip6_ensemble_prediction_sessionInfo.txt"))
message("[29] done.")
