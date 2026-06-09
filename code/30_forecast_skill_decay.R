# ============================================================
# Scientific question / 科学问题:
#   v1 trained XGBoost on 2002–2024 and projected onto SSP2050/2080
#   climates without quantifying covariate-shift risk. How much
#   forecast skill (AUC, calibration) is retained as the prediction
#   horizon lengthens, and which features dominate covariate shift
#   between training and the SSP2050 climate envelope?
#   v1 XGBoost 直接外推到 SSP2050 而不评估预报衰减。本脚本量化预报
#   horizon 上 AUC 衰减与逐特征 covariate shift（PSI）。
#
# Objective / 分析目标:
#   - Train XGBoost on years 2002–2014 (training fold).
#   - Validate on rolling holdouts {2015, 2016, …, 2024}.
#   - Plot AUC vs horizon and Brier vs horizon.
#   - Compute per-feature Population Stability Index (PSI) between
#     training distribution and (a) recent observed, (b) SSP2050
#     ensemble distribution.
#
# Input data / 输入数据:
#   data/raw/hazard_risk_upgraded_complete_case.csv
#   results/forecasts/cmip6_ensemble.parquet (if available; otherwise
#     skip the SSP2050 PSI panel and warn).
#
# Main workflow / 主要流程:
#   1. Load risk set, derive features, time-split.
#   2. Train XGBoost on the training window with early-stopping CV.
#   3. Predict on each holdout year, record AUC/PR-AUC/Brier/calibration.
#   4. Compute PSI per feature: training vs latest holdout year, and
#      training vs SSP2050 ensemble (if available).
#   5. Persist results/forecasts/table_forecast_skill_decay.csv and
#      results/forecasts/table_feature_psi.csv.
#
# Expected output / 预期输出:
#   results/forecasts/table_forecast_skill_decay.csv
#   results/forecasts/table_feature_psi.csv
#   figures/diagnostics/forecast_skill_decay.pdf
#
# Key assumptions / 关键假设:
#   - Training window 2002–2014; holdouts 2015–2024.
#   - Binary event ∈ {0,1}.
#
# Main packages / 主要包: xgboost, data.table, pROC, PRROC, ggplot2,
#   patchwork, arrow.
# Output directory / 输出路径: results/forecasts/, figures/diagnostics/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(xgboost)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(patchwork)
  library(arrow)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_plots.R"))

set.seed(42)

CFG <- list(
  train_years  = 2002:2014,
  test_years   = 2015:2024,
  features     = c("temp_anom_z", "precip_anom_z", "climate_velocity_z",
                   "mahalanobis_dist_z", "log_n_visits_z",
                   "log_n_observers_z", "log_n_birding_days_z",
                   "effort_pc1_z", "sdm_suitability_z"),
  xgb_params   = list(objective = "binary:logistic",
                      eval_metric = "auc",
                      eta = 0.05, max_depth = 6,
                      subsample = 0.8, colsample_bytree = 0.8,
                      lambda = 1, alpha = 0,
                      nthread = max(1L, parallel::detectCores() - 2L))
)

# ---- 1. Load --------------------------------------------------------------
dt <- fread(path_raw("hazard_risk_upgraded_complete_case.csv"),
            encoding = "UTF-8")
feat_present <- intersect(CFG$features, names(dt))
if (length(feat_present) < 4L) {
  stop("[30] Need at least 4 features; found: ",
       paste(feat_present, collapse = ", "))
}
features <- feat_present
message(glue("[30] Using features: {paste(features, collapse=', ')}"))

dt <- dt[!is.na(event)]
train_dt <- dt[year %in% CFG$train_years]
holdouts <- CFG$test_years

# ---- 2. Train --------------------------------------------------------------
mk_dmat <- function(d) xgb.DMatrix(data = as.matrix(d[, ..features]),
                                    label = d$event)
dtrain <- mk_dmat(train_dt)
cv <- xgb.cv(params = CFG$xgb_params, data = dtrain,
              nrounds = 800, nfold = 5, early_stopping_rounds = 30,
              verbose = 0)
best_nrounds <- cv$best_iteration
message(glue("[30] CV best nrounds = {best_nrounds}, ",
              "train AUC = {round(cv$evaluation_log$train_auc_mean[best_nrounds],3)}"))

model <- xgboost::xgb.train(params = CFG$xgb_params,
                             data   = dtrain,
                             nrounds = best_nrounds,
                             verbose = 0)

# ---- 3. Holdout predictions ------------------------------------------------
skill <- data.table::rbindlist(lapply(holdouts, function(yr) {
  test_dt <- dt[year == yr]
  if (nrow(test_dt) < 10L) return(NULL)
  pred <- predict(model, as.matrix(test_dt[, ..features]))
  obs  <- test_dt$event
  auc_roc <- tryCatch(as.numeric(pROC::auc(pROC::roc(obs, pred, quiet = TRUE))),
                       error = function(e) NA_real_)
  pr_auc <- tryCatch(PRROC::pr.curve(scores.class0 = pred[obs == 1],
                                      scores.class1 = pred[obs == 0])$auc.integral,
                     error = function(e) NA_real_)
  brier <- mean((pred - obs)^2)
  cal_intercept <- tryCatch(coef(glm(obs ~ offset(qlogis(pmin(0.999, pmax(0.001, pred)))),
                                      family = binomial()))[1],
                            error = function(e) NA_real_)
  data.table::data.table(
    year          = yr,
    horizon       = yr - max(CFG$train_years),
    n             = nrow(test_dt),
    auc_roc       = auc_roc,
    auc_pr        = pr_auc,
    brier         = brier,
    cal_intercept = cal_intercept
  )
}))

ensure_dir(path_forecasts())
fwrite(skill, path_forecasts("table_forecast_skill_decay.csv"))

# ---- 4. PSI (Population Stability Index) per feature ----------------------
psi <- function(ref, new, bins = 10) {
  brks <- unique(quantile(ref, probs = seq(0, 1, length.out = bins + 1),
                          na.rm = TRUE))
  if (length(brks) < 3L) return(NA_real_)
  p <- table(cut(ref, brks, include.lowest = TRUE)) / length(ref)
  q <- table(cut(new, brks, include.lowest = TRUE)) / length(new)
  eps <- 1e-6
  sum((p - q) * log((p + eps) / (q + eps)))
}

# Training vs latest observed (2024). 训练 vs 最新观测年 PSI。
last_yr <- max(holdouts)
psi_latest <- data.table::data.table(
  feature = features,
  psi_train_vs_2024 = vapply(features, function(f) {
    psi(train_dt[[f]], dt[year == last_yr][[f]])
  }, numeric(1))
)

# Training vs SSP2050 ensemble (if available). 训练 vs SSP2050 PSI。
ssp_path <- path_forecasts("cmip6_ensemble.parquet")
psi_ssp <- if (file.exists(ssp_path)) {
  ssp_dt <- arrow::read_parquet(ssp_path) |> data.table::as.data.table()
  data.table::data.table(
    feature = features,
    psi_train_vs_ssp2050 = vapply(features, function(f) {
      if (!f %in% names(ssp_dt)) return(NA_real_)
      psi(train_dt[[f]], ssp_dt[[f]])
    }, numeric(1))
  )
} else {
  warning("[30] No CMIP6 ensemble at ", ssp_path,
          "; run script 29 to populate SSP2050 PSI.")
  data.table::data.table(feature = features,
                          psi_train_vs_ssp2050 = NA_real_)
}

psi_all <- merge(psi_latest, psi_ssp, by = "feature")
fwrite(psi_all, path_forecasts("table_feature_psi.csv"))

# ---- 5. Diagnostic plots ---------------------------------------------------
p1 <- ggplot(skill, aes(x = horizon, y = auc_roc)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey50") +
  labs(title = "Forecast skill decay (ROC-AUC)",
       x = "Horizon (years beyond training window)",
       y = "AUC") +
  theme_geb()

p2 <- ggplot(skill, aes(x = horizon, y = brier)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  labs(title = "Brier score by horizon", x = "Horizon (years)", y = "Brier") +
  theme_geb()

psi_long <- data.table::melt(psi_all, id.vars = "feature",
                              variable.name = "comparison",
                              value.name = "psi")
p3 <- ggplot(psi_long, aes(x = reorder(feature, psi), y = psi,
                            fill = comparison)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = c(0.1, 0.25), linetype = "dashed",
              colour = c("grey60", "grey30")) +
  coord_flip() +
  scale_fill_manual(values = pal_cat[1:2]) +
  labs(title = "Covariate shift per feature (PSI)",
        subtitle = "Dashed lines: 0.10 = moderate shift, 0.25 = severe shift",
        x = NULL, y = "PSI") +
  theme_geb()

panel <- (p1 | p2) / p3
ensure_dir(path_diag_fig())
ggplot2::ggsave(path_diag_fig("forecast_skill_decay.pdf"),
                panel, width = 18, height = 14, units = "cm",
                device = grDevices::cairo_pdf)
ggplot2::ggsave(path_diag_fig("forecast_skill_decay.png"),
                panel, width = 18, height = 14, units = "cm", dpi = 600)

message("[30] Skill table:")
print(skill)
message("[30] PSI table:")
print(psi_all)

dump_session_info(path_logs("30_forecast_skill_decay_sessionInfo.txt"))
message("[30] done.")
