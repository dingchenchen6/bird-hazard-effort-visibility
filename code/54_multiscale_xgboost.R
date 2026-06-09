# ============================================================
# Scientific question / 科学问题:
#   Can XGBoost models trained at province, prefecture, county,
#   and 100km grid scales replicate the climate × effort
#   interaction found by glmmTMB? Which variables matter most
#   at each scale, and how do future hazards differ?
#   XGBoost 在省/市/县/100km网格尺度上能否复现 climate ×
#   effort 交互？各尺度变量重要性如何？未来风险如何变化？
#
# Objective / 分析目标:
#   - Train XGBoost at each scale with cross-validation
#   - SHAP-based variable importance and interaction effects
#   - Future scenario prediction (SSP245/585 × effort scenarios)
#   - All outputs to outputs_multiscale/ (no overwrite)
#
# Input data / 输入数据:
#   - outputs_multiscale/data/derived/risk_set_{province,prefecture,
#     county,grid_100km}.csv (from script 51)
#
# Expected output / 预期输出:
#   outputs_multiscale/results/tables/
#     table_multiscale_xgboost_performance.csv
#     table_multiscale_shap_importance.csv
#     table_multiscale_xgb_future_{scale}.csv
#   outputs_multiscale/figures/main/
#     fig_multiscale_xgb_shap_importance.{pdf,png}
#     fig_multiscale_xgb_shap_beeswarm.{pdf,png}
#     fig_multiscale_xgb_future_comparison_2050.{pdf,png}
#   outputs_multiscale/logs/54_multiscale_xgboost.log
#
# Key assumptions / 关键假设:
#   - Risk sets from script 51 must exist
#   - xgboost and shapviz packages must be installed
#   - XGBoost uses binary:logistic with scale_pos_weight
#   - 5-fold CV with early stopping for hyperparameter tuning
#   - Future predictions use same scenario logic as script 53
#
# Main packages / 主要包: data.table, xgboost, shapviz, ggplot2.
# Output directory / 输出路径: outputs_multiscale/
# ============================================================

# ---- 0. CLI args ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
arg <- list()
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--base-dir")        { arg$base_dir   <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--output-dir") { arg$output_dir <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--skip-county"){ arg$skip_county <- TRUE; i <- i + 1L }
  else { i <- i + 1L }
}
`%||%` <- function(a, b) if (is.null(a)) b else a

BASE_DIR   <- arg$base_dir   %||% "/Users/dingchenchen/Documents/New records/bird-new-distribution-records/tasks/bird_hazard_model_effort_upgrade_v2"
OUTPUT_DIR <- arg$output_dir %||% "outputs_multiscale"
SKIP_COUNTY <- isTRUE(arg$skip_county)

OUT <- file.path(BASE_DIR, OUTPUT_DIR)
ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
for (sd in c("results/tables", "results/forecasts", "figures/main", "logs"))
  ens(file.path(OUT, sd))

log_file <- file.path(OUT, "logs", "54_multiscale_xgboost.log")
log <- function(...) {
  msg <- sprintf("[54 %s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(msg); cat(msg, file = log_file, append = TRUE)
}
log("=== 54_multiscale_xgboost.R START ===")

# ---- 1. Packages ----------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(xgboost)
  library(ggplot2)
})
options(warn = 1)
set.seed(42)

zify <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# ---- 2. Load risk sets from script 51 ------------------------------------
log("=== Step 2: Loading risk sets ===")

scales <- c("province", "prefecture", "grid_100km")
if (!SKIP_COUNTY) scales <- c(scales, "county")
scales <- intersect(c("province", "prefecture", "county", "grid_100km"), scales)

risk_sets <- list()
for (sc in scales) {
  f <- file.path(OUT, "data", "derived", paste0("risk_set_", sc, ".csv"))
  if (file.exists(f)) {
    risk_sets[[sc]] <- fread(f, encoding = "UTF-8")
    log(sc, ": ", nrow(risk_sets[[sc]]), " rows")
  } else {
    log("WARNING: ", f, " not found")
  }
}
if (length(risk_sets) == 0L) stop("No risk sets. Run 51 first.")

# Province climate/effort for future panels
prov_clim <- fread(file.path(BASE_DIR, "data", "raw",
                             "climate_metrics_province_year.csv"),
                   encoding = "UTF-8")
setnames(prov_clim, tolower(names(prov_clim)))
prov_eff <- fread(file.path(BASE_DIR, "data", "raw",
                            "effort_panel_upgraded.csv"),
                  encoding = "UTF-8")
setnames(prov_eff, tolower(names(prov_eff)))

# ---- 3. Train XGBoost at each scale --------------------------------------
log("=== Step 3: Training XGBoost models ===")

xgb_results <- list()
shap_all <- list()

for (sc in names(risk_sets)) {
  log("  --- XGBoost: ", sc, " ---")
  dt <- copy(risk_sets[[sc]])

  # Determine feature columns available at this scale
  # 候选特征列：尺度特异性的气候+effort列
  climate_features <- intersect(c("climate_velocity_z", "climate_exposure_z",
    "warming_rate_z", "mahalanobis_dist_z", "temp_grad_z", "prec_grad_z",
    "prov_climate_velocity_z", "prov_warming_rate_z", "prov_temp_anom_z",
    "prov_climate_exposure_z", "temp_anom_z", "prec_anom_z"), names(dt))
  effort_features <- intersect(c("effort_z", "log_n_events_z",
    "log_effort_visits_z", "effort_pc1_z", "log_n_visits_z",
    "log_effort_days_z", "log_n_observers_z"), names(dt))
  # year feature
  year_feature <- intersect(c("year_c", "year"), names(dt))[1]
  if (is.na(year_feature)) { dt[, year_c := year - 2013]; year_feature <- "year_c" }

  feature_cols <- c(year_feature, climate_features, effort_features)

  # 添加交互项 (climate × effort)
  # 找第一个 climate 和 effort 特征做交互
  clim_main <- climate_features[1]
  eff_main  <- effort_features[1]
  if (!is.na(clim_main) && !is.na(eff_main)) {
    dt[, clim_x_effort := get(clim_main) * get(eff_main)]
    feature_cols <- c(feature_cols, "clim_x_effort")
  }

  # 完整案例
  xgb_data <- dt[complete.cases(dt[, ..feature_cols])]
  if (nrow(xgb_data) < 500L) {
    log("  too few complete rows: ", nrow(xgb_data), ", skip"); next
  }

  X <- as.matrix(xgb_data[, ..feature_cols])
  y <- xgb_data$event

  pos_weight <- sum(y == 0) / pmax(sum(y == 1), 1L)
  log("  rows=", nrow(X), " events=", sum(y), " pos_weight=", round(pos_weight, 1))

  dtrain <- xgb.DMatrix(data = X, label = y)

  # 5-fold CV 快速调参 (simplified grid for speed)
  log("  tuning via 5-fold CV ...")
  best_score <- Inf
  best_params <- list(max_depth = 6, eta = 0.05, min_child_weight = 5,
                      nrounds = 200L)

  for (md in c(4, 6, 8)) {
    for (eta in c(0.05, 0.1)) {
      params <- list(
        objective = "binary:logistic", eval_metric = "logloss",
        max_depth = md, eta = eta, min_child_weight = 5,
        subsample = 0.8, colsample_bytree = 0.8,
        scale_pos_weight = pos_weight)
      cv_res <- tryCatch(
        xgb.cv(params, dtrain, nrounds = 300L, nfold = 5L,
               early_stopping_rounds = 20L, verbose = 0),
        error = function(e) NULL)
      if (is.null(cv_res)) next
      iter <- cv_res$best_iteration
      if (is.null(iter) || !is.finite(iter)) iter <- 300L
      score <- cv_res$evaluation_log[iter, test_logloss_mean]
      if (!is.na(score) && score < best_score) {
        best_score <- score
        best_params <- list(max_depth = md, eta = eta, min_child_weight = 5,
                            nrounds = iter)
      }
    }
  }
  log("  best CV: md=", best_params$max_depth, " eta=", best_params$eta,
      " nrounds=", best_params$nrounds, " logloss=", round(best_score, 4))

  # 最终模型
  final_params <- list(
    objective = "binary:logistic", eval_metric = "logloss",
    max_depth = best_params$max_depth, eta = best_params$eta,
    min_child_weight = best_params$min_child_weight,
    subsample = 0.8, colsample_bytree = 0.8,
    scale_pos_weight = pos_weight)

  xgb_fit <- xgb.train(final_params, dtrain,
                        nrounds = best_params$nrounds, verbose = 0)

  # AUC on training data
  pred_train <- predict(xgb_fit, X)
  auc <- tryCatch({
    idx <- !is.na(y) & !is.na(pred_train)
    if (sum(idx) > 100L) {
      pROC::auc(y[idx], pred_train)[1]
    } else NA
  }, error = function(e) NA)

  xgb_results[[sc]] <- data.table(
    scale = sc,
    n_rows = nrow(X), n_events = sum(y),
    max_depth = best_params$max_depth, eta = best_params$eta,
    nrounds = best_params$nrounds,
    cv_logloss = best_score, train_auc = auc)
  log("  train AUC=", round(auc, 4))

  # SHAP values (using shapviz if available, otherwise gain importance)
  if (requireNamespace("shapviz", quietly = TRUE)) {
    log("  computing SHAP via shapviz ...")
    sv <- shapviz::shapviz(xgb_fit, X_pred = X)
    sv_imp <- shapviz::importance(sv, kind = "mean_abs")
    sv_dt <- data.table(scale = sc, variable = names(sv_imp),
                         shap_importance = as.numeric(sv_imp))
    sv_dt <- sv_dt[order(-shap_importance)]
    shap_all[[sc]] <- sv_dt
    log("  top 3 SHAP: ", paste(sv_dt$variable[1:3], collapse = ", "))

    # SHAP beeswarm plot (sample 5000 rows for speed)
    n_sample <- min(5000L, nrow(X))
    set.seed(42)
    idx_sample <- sample(seq_len(nrow(X)), n_sample)
    sv_sub <- shapviz::shapviz(xgb_fit, X_pred = X[idx_sample, , drop = FALSE])

    p_bee <- shapviz::sv_importance(sv, kind = "mean_abs", max_display = 12) +
      labs(title = paste0("XGBoost SHAP importance: ", sc)) +
      theme_bw(base_size = 9)

    ggsave(file.path(OUT, "figures", "main",
                     paste0("fig_multiscale_xgb_shap_", sc, ".pdf")),
           p_bee, width = 12, height = 8, units = "cm",
           device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main",
                     paste0("fig_multiscale_xgb_shap_", sc, ".png")),
           p_bee, width = 12, height = 8, units = "cm", dpi = 600)

  } else {
    # fallback: gain-based importance
    log("  shapviz not available, using gain importance")
    imp <- xgb.importance(model = xgb_fit)
    imp_dt <- data.table(scale = sc, variable = imp$Feature,
                          shap_importance = imp$Gain)
    imp_dt <- imp_dt[order(-shap_importance)]
    shap_all[[sc]] <- imp_dt
    log("  top 3 gain: ", paste(imp_dt$variable[1:3], collapse = ", "))
  }

  # Save model
  xgb.save(xgb_fit, file.path(OUT, "results", "forecasts",
                               paste0("xgb_model_", sc, ".model")))

  # ---- Future prediction ---------------------------------------------------
  log("  building future prediction panel ...")
  future_years <- c(2030, 2035, 2040, 2045, 2050)
  climate_scenarios <- c("current", "ssp245", "ssp585")
  effort_scenarios  <- c("baseline", "trend", "doubled")

  # Province mapping for each unit
  unit_col <- if ("unit_id" %in% names(dt)) "unit_id"
              else if ("grid_id" %in% names(dt)) "grid_id"
              else "province"
  # grid_100km may not have a province column
  if (!"province" %in% names(dt)) {
    dt[, province := "ALL"]
    log("  no province column, using overall mean for future scenarios")
  }
  unit_prov <- unique(dt[, c(unit_col, "province"), with = FALSE])

  # Current baseline (2024)
  current_eff <- prov_eff[year == 2024, .(province, log_effort_visits_z, effort_pc1_z)]
  eff_trends <- prov_eff[year >= 2002 & year <= 2024,
    .(effort_trend_z = tryCatch(coef(lm(log_effort_visits_z ~ year))[2],
                                 error = function(e) 0)),
    by = province]
  prov_clim_2024 <- prov_clim[year == 2024,
    .(province, climate_velocity_z, warming_rate_z)]
  clim_vel_sd <- sd(prov_clim$climate_velocity_z, na.rm = TRUE)

  # Build panel
  fut <- as.data.table(expand.grid(
    unit_val = unique(dt[[unit_col]]),
    year = future_years,
    climate_scenario = climate_scenarios,
    effort_scenario = effort_scenarios,
    stringsAsFactors = FALSE))
  setnames(fut, "unit_val", unit_col)
  fut <- merge(fut, unit_prov, by = unit_col, all.x = TRUE)
  fut <- merge(fut, current_eff, by = "province", all.x = TRUE)
  fut <- merge(fut, eff_trends, by = "province", all.x = TRUE)
  fut <- merge(fut, prov_clim_2024, by = "province", all.x = TRUE)

  # Effort scenarios
  if ("log_effort_visits_z" %in% names(fut)) {
    fut[, effort_z_future := fcase(
      effort_scenario == "baseline", log_effort_visits_z,
      effort_scenario == "trend",
        log_effort_visits_z + effort_trend_z * (year - 2024),
      effort_scenario == "doubled", log_effort_visits_z * 2)]
  } else {
    fut[, effort_z_future := 0]
  }
  fut[is.na(effort_z_future), effort_z_future := 0]

  # Climate scenarios
  if ("climate_velocity_z" %in% names(fut)) {
    fut[, climate_velocity_z_future := fcase(
      climate_scenario == "current", climate_velocity_z,
      climate_scenario == "ssp245",
        climate_velocity_z + 0.3 / clim_vel_sd * (year - 2024) / 10,
      climate_scenario == "ssp585",
        climate_velocity_z + 0.6 / clim_vel_sd * (year - 2024) / 10)]
  } else {
    fut[, climate_velocity_z_future := 0]
  }
  fut[is.na(climate_velocity_z_future), climate_velocity_z_future := 0]

  # Build feature matrix matching training columns
  for (fc in feature_cols) {
    if (!fc %in% names(fut)) {
      if (fc == "clim_x_effort") {
        fut[, clim_x_effort := climate_velocity_z_future * effort_z_future]
      } else if (fc == "year_c") {
        fut[, year_c := year - 2013]
      } else if (fc %in% c("climate_velocity_z", "climate_exposure_z",
                           "warming_rate_z", "mahalanobis_dist_z",
                           "temp_grad_z", "prec_grad_z")) {
        # Use future climate for main climate vars
        if (fc == "climate_velocity_z" && "climate_velocity_z_future" %in% names(fut))
          fut[, (fc) := climate_velocity_z_future]
        else
          fut[, (fc) := 0]
      } else if (fc %in% effort_features) {
        if (fc == eff_main)
          fut[, (fc) := effort_z_future]
        else
          fut[, (fc) := 0]
      } else {
        fut[, (fc) := 0]
      }
    }
  }

  X_future <- as.matrix(fut[, ..feature_cols])
  fut[, hazard := predict(xgb_fit, X_future)]
  fut[, scale := sc]

  out_cols <- c(unit_col, "province", "year", "climate_scenario",
                "effort_scenario", "hazard", "scale")
  out_cols <- intersect(out_cols, names(fut))
  fwrite(fut[, ..out_cols],
         file.path(OUT, "results", "forecasts",
                   paste0("table_multiscale_xgb_future_", sc, ".csv")))
  log("  future predictions saved: ", nrow(fut), " rows")
}

# ---- 4. Save results tables -----------------------------------------------
xgb_dt <- rbindlist(xgb_results, fill = TRUE)
fwrite(xgb_dt, file.path(OUT, "results", "tables",
                         "table_multiscale_xgboost_performance.csv"))
log("XGBoost performance table saved")

if (length(shap_all) > 0L) {
  shap_dt <- rbindlist(shap_all, fill = TRUE)
  fwrite(shap_dt, file.path(OUT, "results", "tables",
                            "table_multiscale_shap_importance.csv"))
  log("SHAP importance table saved")

  # Cross-scale SHAP comparison figure
  shap_dt[, scale := factor(scale,
    levels = c("province", "prefecture", "county", "grid_100km"),
    labels = c("Province", "Prefecture", "County", "100km grid"))]
  shap_dt[, imp_std := shap_importance / max(shap_importance, na.rm = TRUE),
          by = scale]

  pal_varimp <- c("Climate" = "#d94801", "Effort" = "#2171b5",
                  "Year" = "#666666", "Other" = "#aaaaaa", "Interaction" = "#6a51a3")
  shap_dt[, category := fcase(
    variable %like% "climate|velocity|exposure|warming|mahalanobis|temp_grad|prec_grad",
    "Climate",
    variable %like% "effort|log_n_|log_effort", "Effort",
    variable == "clim_x_effort", "Interaction",
    variable == "year_c", "Year",
    default = "Other")]

  p_shap <- ggplot(shap_dt, aes(x = reorder(variable, imp_std), y = imp_std,
                                 fill = category)) +
    geom_col(alpha = 0.85) +
    coord_flip() +
    facet_wrap(~ scale, scales = "free_x", ncol = 2) +
    scale_fill_manual(values = pal_varimp, name = "Category") +
    labs(x = NULL, y = "Standardized SHAP importance",
         title = "XGBoost SHAP variable importance across scales") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          strip.text = element_text(face = "bold"))

  ggsave(file.path(OUT, "figures", "main",
                   "fig_multiscale_xgb_shap_importance.pdf"),
         p_shap, width = 22, height = 16, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(file.path(OUT, "figures", "main",
                   "fig_multiscale_xgb_shap_importance.png"),
         p_shap, width = 22, height = 16, units = "cm", dpi = 600)
  log("SHAP importance figure saved")
}

# ---- 5. Future trajectory comparison --------------------------------------
log("=== Step 5: Future trajectory comparison ===")

xgb_fut_files <- list.files(file.path(OUT, "results", "forecasts"),
                             pattern = "xgb_future_.*\\.csv$", full.names = TRUE)
if (length(xgb_fut_files) > 0L) {
  xgb_fut <- rbindlist(lapply(xgb_fut_files, fread), fill = TRUE)
  natl <- xgb_fut[,
    .(hazard_mean = mean(hazard, na.rm = TRUE)),
    by = .(scale, year, climate_scenario, effort_scenario)]
  natl[, scale := factor(scale,
    levels = c("province", "prefecture", "county", "grid_100km"),
    labels = c("Province", "Prefecture", "County", "100km grid"))]

  p_fut <- ggplot(natl[effort_scenario == "baseline"],
                   aes(x = year, y = hazard_mean,
                       colour = climate_scenario, linetype = scale)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    scale_colour_manual(values = c("current" = "#2166AC",
                                    "ssp245" = "#F4A582",
                                    "ssp585" = "#B2182B"),
                        name = "Climate") +
    scale_linetype(name = "Scale") +
    labs(x = "Year", y = "Mean predicted hazard (XGBoost)",
         title = "XGBoost future hazard trajectories across spatial scales") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  ggsave(file.path(OUT, "figures", "main",
                   "fig_multiscale_xgb_future_comparison_2050.pdf"),
         p_fut, width = 18, height = 10, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(file.path(OUT, "figures", "main",
                   "fig_multiscale_xgb_future_comparison_2050.png"),
         p_fut, width = 18, height = 10, units = "cm", dpi = 600)
  log("XGBoost future comparison figure saved")
}

# ---- 6. Session info ------------------------------------------------------
sink(file.path(OUT, "logs", "54_sessionInfo.txt"))
print(sessionInfo())
sink()

log("=== 54_multiscale_xgboost.R DONE ===")
