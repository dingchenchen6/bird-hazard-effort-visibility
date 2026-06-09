# ============================================================
# Scientific question / 科学问题:
#   v1 reported XGBoost CV AUC = 0.717 ± 0.027 using row-random folds
#   (code/22:65), which leaks spatial autocorrelation and inflates
#   apparent skill. What does the AUC / PR-AUC / Brier / calibration
#   look like under proper 250 km spatial-block CV?
#   v1 用行随机 5 折 CV，存在空间渗漏。改用 250 km block CV 后表现？
#
# Objective / 分析目标:
#   - Build 250-km spatial blocks across mainland China (k = 5).
#   - Train XGBoost on the same feature set as v1 (script 22).
#   - Report per-fold AUC, PR-AUC, Brier, calibration slope.
#   - Save fold maps and reliability curves.
#
# Input data / 输入数据:
#   data/raw/hazard_risk_upgraded_complete_case.csv  (must contain
#     province + climate + effort columns, latitude/longitude or
#     joinable to data/spatial/basemap_GS2019_1822 for centroid).
#
# Main workflow / 主要流程:
#   1. Load risk set + attach province centroid coordinates.
#   2. Build sf points; generate blockCV folds (250 km, k = 5).
#   3. For each fold, train XGBoost, predict on test, compute metrics.
#   4. Persist results/diagnostics/table_spatial_block_cv.csv and
#      figures/diagnostics/spatial_block_cv.{pdf,png}.
#
# Expected output / 预期输出:
#   results/diagnostics/table_spatial_block_cv.csv
#   results/diagnostics/table_spatial_block_cv_summary.csv
#   figures/diagnostics/spatial_block_cv_folds.pdf/png
#   figures/diagnostics/spatial_block_cv_reliability.pdf/png
#
# Key assumptions / 关键假设:
#   - province centroid is a reasonable proxy for spatial location
#     of each record (used because v1 risk set is at province × year
#     resolution).  For grid-level CV use script 26b (planned).
#
# Main packages / 主要包: data.table, sf, blockCV, xgboost, pROC,
#   PRROC, ggplot2, patchwork, arrow.
# Output directory / 输出路径: results/diagnostics/, figures/diagnostics/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(blockCV)
  library(xgboost)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(patchwork)
  library(arrow)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_spatial.R"))
source(file.path("code", "utils", "utils_plots.R"))

set.seed(42)

CFG <- list(
  block_km = 250,
  k = 5,
  iteration = 50,
  features = c("temp_anom_z", "precip_anom_z", "climate_velocity_z",
               "mahalanobis_dist_z", "log_n_visits_z",
               "log_n_observers_z", "log_n_birding_days_z",
               "effort_pc1_z", "sdm_suitability_z"),
  xgb_params = list(objective = "binary:logistic",
                    eval_metric = "auc",
                    eta = 0.05, max_depth = 6,
                    subsample = 0.8, colsample_bytree = 0.8,
                    nthread = max(1L, parallel::detectCores() - 2L))
)

# ---- 1. Load risk set ------------------------------------------------------
dt <- fread(path_raw("hazard_risk_upgraded_complete_case.csv"),
            encoding = "UTF-8")
features <- intersect(CFG$features, names(dt))
if (length(features) < 4L) stop("[26] Need ≥4 features; have: ", length(features))

# ---- 2. Province centroids -------------------------------------------------
prov_sf <- read_gs2019_basemap("province")
prov_sf$province_norm <- tolower(prov_sf$name %||% prov_sf$NAME)
centroids <- sf::st_centroid(to_albers(prov_sf))
centroid_dt <- data.table::data.table(
  province = prov_sf$province_norm,
  x = sf::st_coordinates(centroids)[, 1],
  y = sf::st_coordinates(centroids)[, 2]
)

dt[, province_norm := tolower(province)]
dt <- merge(dt, centroid_dt, by.x = "province_norm", by.y = "province",
            all.x = TRUE)
if (any(is.na(dt$x))) {
  warning("[26] ", sum(is.na(dt$x)),
          " records without province centroid — dropped.")
  dt <- dt[!is.na(x)]
}

points_sf <- sf::st_as_sf(dt[, .(x, y, event)], coords = c("x", "y"),
                          crs = CRS_ALBERS)

# ---- 3. blockCV folds ------------------------------------------------------
message("[26] Building spatial blocks (250 km, k=5)…")
cv <- blockCV::cv_spatial(
  x = points_sf,
  column    = "event",
  size      = CFG$block_km * 1000,
  k         = CFG$k,
  selection = "random",
  iteration = CFG$iteration,
  plot      = FALSE,
  progress  = FALSE
)

folds_list <- cv$folds_list
metrics <- data.table::rbindlist(lapply(seq_along(folds_list), function(i) {
  tr <- folds_list[[i]]$train
  te <- folds_list[[i]]$test
  train_dt <- dt[tr]
  test_dt  <- dt[te]
  message(glue("  fold {i}: train={nrow(train_dt)} test={nrow(test_dt)} ",
                "pos_train={sum(train_dt$event)} pos_test={sum(test_dt$event)}"))
  dtrain <- xgb.DMatrix(as.matrix(train_dt[, ..features]), label = train_dt$event)
  dtest  <- xgb.DMatrix(as.matrix(test_dt[, ..features]),  label = test_dt$event)
  cv_res <- xgb.cv(params = CFG$xgb_params, data = dtrain,
                    nrounds = 800, nfold = 5,
                    early_stopping_rounds = 30, verbose = 0)
  best <- cv_res$best_iteration
  model <- xgboost::xgb.train(params = CFG$xgb_params, data = dtrain,
                               nrounds = best, verbose = 0)
  pred <- predict(model, dtest)
  obs  <- test_dt$event
  auc  <- tryCatch(as.numeric(pROC::auc(pROC::roc(obs, pred, quiet = TRUE))),
                   error = function(e) NA_real_)
  prauc <- tryCatch(PRROC::pr.curve(scores.class0 = pred[obs == 1],
                                      scores.class1 = pred[obs == 0])$auc.integral,
                    error = function(e) NA_real_)
  brier <- mean((pred - obs)^2)
  cal_slope <- tryCatch(coef(glm(obs ~ qlogis(pmin(0.999, pmax(0.001, pred))),
                                  family = binomial()))[2],
                        error = function(e) NA_real_)
  data.table::data.table(
    fold = i,
    n_train = nrow(train_dt), n_test = nrow(test_dt),
    pos_test = sum(obs),
    auc_roc = auc, auc_pr = prauc, brier = brier,
    calibration_slope = cal_slope
  )
}))

ensure_dir(path_diagnostics())
fwrite(metrics, path_diagnostics("table_spatial_block_cv.csv"))

summary_dt <- data.table::data.table(
  metric = c("auc_roc", "auc_pr", "brier", "calibration_slope"),
  mean = c(mean(metrics$auc_roc, na.rm = TRUE),
            mean(metrics$auc_pr,  na.rm = TRUE),
            mean(metrics$brier,   na.rm = TRUE),
            mean(metrics$calibration_slope, na.rm = TRUE)),
  sd   = c(sd(metrics$auc_roc, na.rm = TRUE),
            sd(metrics$auc_pr,  na.rm = TRUE),
            sd(metrics$brier,   na.rm = TRUE),
            sd(metrics$calibration_slope, na.rm = TRUE))
)
fwrite(summary_dt, path_diagnostics("table_spatial_block_cv_summary.csv"))

# ---- 4. Fold map -----------------------------------------------------------
fold_assign <- data.table::data.table(point_idx = seq_len(nrow(points_sf)),
                                       fold = NA_integer_)
for (i in seq_along(folds_list)) {
  fold_assign[folds_list[[i]]$test, fold := i]
}
points_sf$fold <- factor(fold_assign$fold)

base_provinces <- to_albers(prov_sf)
p_folds <- ggplot() +
  geom_sf(data = base_provinces, fill = "grey98", colour = "grey70",
          linewidth = 0.1) +
  geom_sf(data = points_sf, aes(colour = fold), size = 0.6, alpha = 0.7) +
  scale_colour_manual(values = pal_cat[1:CFG$k], na.value = "grey80") +
  labs(title = glue("Spatial-block CV folds ({CFG$block_km} km, k = {CFG$k})"),
       colour = "Fold") +
  theme_geb() +
  theme(legend.position = "right")

# ---- 5. Reliability curve --------------------------------------------------
# Pool all folds' test predictions for reliability curve. 校准曲线。
# We retrain a single model on all data for a global reliability check;
# fold-level reliability could be added in script 35.
dtrain_full <- xgb.DMatrix(as.matrix(dt[, ..features]), label = dt$event)
cv_full <- xgb.cv(params = CFG$xgb_params, data = dtrain_full, nrounds = 800,
                   nfold = 5, early_stopping_rounds = 30, verbose = 0)
model_full <- xgboost::xgb.train(params = CFG$xgb_params, data = dtrain_full,
                                  nrounds = cv_full$best_iteration, verbose = 0)
pred_full <- predict(model_full, dtrain_full)
cal_dt <- data.table::data.table(
  pred = pred_full,
  obs  = dt$event,
  bin  = cut(pred_full, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
)
cal_summary <- cal_dt[, .(mean_pred = mean(pred), obs_rate = mean(obs),
                          n = .N), by = bin]
p_reliability <- ggplot(cal_summary, aes(x = mean_pred, y = obs_rate)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted",
               colour = "grey50") +
  geom_point(aes(size = n)) + geom_line(linewidth = 0.4) +
  scale_size_continuous(range = c(0.5, 4)) +
  coord_fixed() +
  labs(title = "Reliability (calibration) curve",
        x = "Predicted probability", y = "Observed rate",
        size = "Bin n") +
  theme_geb()

ensure_dir(path_diag_fig())
ggplot2::ggsave(path_diag_fig("spatial_block_cv_folds.pdf"),
                p_folds, width = 16, height = 12, units = "cm",
                device = grDevices::cairo_pdf)
ggplot2::ggsave(path_diag_fig("spatial_block_cv_folds.png"),
                p_folds, width = 16, height = 12, units = "cm", dpi = 600)
ggplot2::ggsave(path_diag_fig("spatial_block_cv_reliability.pdf"),
                p_reliability, width = 9, height = 9, units = "cm",
                device = grDevices::cairo_pdf)
ggplot2::ggsave(path_diag_fig("spatial_block_cv_reliability.png"),
                p_reliability, width = 9, height = 9, units = "cm", dpi = 600)

message("[26] Spatial block CV summary:")
print(summary_dt)

dump_session_info(path_logs("26_spatial_block_cv_sessionInfo.txt"))
message("[26] done.")
