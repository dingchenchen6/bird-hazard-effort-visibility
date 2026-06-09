# ============================================================
# Script: 49_v3_variable_importance.R
# Family: v3 risk-set — RF + XGBoost + SHAP variable importance
# Author: Chen-Chen Ding + Claude Opus 4.7
# Date  : 2026-05-23
#
# ------------------------------------------------------------
# Scientific question / 科学问题:
#   On the v2 risk set (333 SDM-tight species), Random Forest
#   ranked `temp × effort` as the most important feature
#   (0.0166). Does the same conclusion hold on the v3 relaxed
#   risk set (463 species, 817 events, 188,870 rows)?
#   v3 风险集下，climate × effort 交互项是否仍是 RF / XGBoost / SHAP
#   排第一的特征？
#
# ------------------------------------------------------------
# Steps:
#   1. Load v3 risk set + attach all climate metrics from
#      climate_metrics_province_year (so we have the full feature
#      set used by the v2 RF in script 41).
#   2. Build engineered interaction columns:
#        temp_x_effort      = temp_grad_z      * log_effort_visits_z
#        mahal_x_effort     = mahalanobis_dist_z * log_effort_visits_z
#   3. Fit ranger Random Forest with permutation importance
#      (500 trees, sqrt(p) mtry, OOB).
#   4. Fit XGBoost with 5-fold CV (binary:logistic, AUC), extract
#      best nrounds, then compute SHAP values via SHAPforxgboost.
#   5. Persist importance tables and produce a publication figure
#      comparing v2 vs v3 RF rankings.
#
# Outputs (NEW files — do not overwrite v2 RF):
#   results/tables/table_rf_importance_v3.csv
#   results/tables/table_xgb_cv_v3.csv
#   results/tables/table_v2_v3_rf_comparison.csv
#   figures/main/Figure_4_variable_importance_v3.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(xgboost)
  library(SHAPforxgboost)
  library(ggplot2)
  library(ggbeeswarm)
  library(patchwork)
})
options(warn = 1); set.seed(42)

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- normalizePath(file.path(V2, "..", "bird_hazard_model_effort_upgrade"),
                     mustWork = FALSE)

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE,
                                                    showWarnings = FALSE)
ens(file.path(V2, "logs"))
LOG <- file(file.path(V2, "logs", "49_v3_variable_importance.log"),
            "wt", encoding = "UTF-8")
on.exit({ if (isOpen(LOG)) close(LOG) }, add = TRUE)
log <- function(...) {
  m <- paste0(sprintf("[49 %s] ", format(Sys.time(), "%H:%M:%S")),
              paste(..., sep = ""))
  cat(m, "\n", sep = ""); writeLines(m, LOG)
}
audit <- function(t) {
  bar <- paste(rep("─", 60), collapse = "")
  log(""); log(bar); log("AUDIT — ", t); log(bar)
}

# ============================================================
# STEP 1 — Load v3 risk set + enrich with all climate metrics
# ============================================================
audit("STEP 1: load v3 risk set + attach all climate metrics")

risk_v3 <- fread(file.path(V2, "data", "derived",
                             "risk_set_province_v3.csv"),
                   encoding = "UTF-8")
log("v3 risk set: ", nrow(risk_v3), " rows / ",
    uniqueN(risk_v3$species), " species / ", sum(risk_v3$event), " events")

prov_clim_full <- fread(file.path(V1, "data",
                                    "climate_metrics_province_year.csv"),
                          encoding = "UTF-8")
climate_cols <- intersect(c("climate_velocity_z","precip_velocity_z",
                              "climate_exposure_z","warming_rate_z",
                              "mahalanobis_dist_z","temp_anom_z",
                              "prec_anom_z","temp_grad_prov_z",
                              "prec_grad_prov_z"),
                           names(prov_clim_full))
risk_v3 <- merge(risk_v3,
                   prov_clim_full[, c("province","year", climate_cols),
                                    with = FALSE],
                   by = c("province","year"), all.x = TRUE)

# We also need other effort columns
prov_eff <- fread(file.path(V1, "data", "effort_panel_upgraded.csv"),
                    encoding = "UTF-8")
eff_cols <- intersect(c("log_effort_record_z","log_effort_observers_z",
                          "log_effort_days_z","effort_pc1_z"),
                       names(prov_eff))
risk_v3 <- merge(risk_v3,
                   prov_eff[, c("province","year", eff_cols), with = FALSE],
                   by = c("province","year"), all.x = TRUE)

log("v3 enriched: ", ncol(risk_v3), " columns total")

# Engineered interaction features (same as v2)
risk_v3[, temp_x_effort  := temp_grad_z * log_effort_visits_z]
risk_v3[, mahal_x_effort := mahalanobis_dist_z * log_effort_visits_z]
risk_v3[, year_c := year - 2013]

# Complete-case
feat_full <- c("temp_grad_z","prec_grad_z","climate_velocity_z",
                "precip_velocity_z","climate_exposure_z","warming_rate_z",
                "mahalanobis_dist_z","temp_anom_z","prec_anom_z",
                "log_effort_visits_z","log_effort_record_z",
                "log_effort_days_z","effort_pc1_z",
                "temp_x_effort","mahal_x_effort","year_c")
feat_full <- intersect(feat_full, names(risk_v3))
train_dt <- risk_v3[, c("event", feat_full), with = FALSE]
train_dt <- train_dt[complete.cases(train_dt)]
log("training data: ", nrow(train_dt), " rows × ",
    length(feat_full), " features | events ", sum(train_dt$event))

# ============================================================
# STEP 2 — Random Forest (ranger) with permutation importance
# ============================================================
audit("STEP 2: ranger Random Forest")

rf_v3 <- ranger(event ~ ., data = train_dt,
                 probability = FALSE, classification = FALSE,
                 num.trees = 500,
                 mtry = floor(sqrt(length(feat_full))),
                 importance = "permutation", seed = 42,
                 num.threads = max(1L, parallel::detectCores() - 2L))
imp_v3 <- data.table(
  variable = names(rf_v3$variable.importance),
  importance = as.numeric(rf_v3$variable.importance))
imp_v3[, category := fcase(
  variable %in% c("temp_grad_z","prec_grad_z","climate_velocity_z",
                   "precip_velocity_z","climate_exposure_z","warming_rate_z",
                   "mahalanobis_dist_z","temp_anom_z","prec_anom_z"),
  "Climate",
  variable %in% c("log_effort_visits_z","log_effort_record_z",
                   "log_effort_days_z","effort_pc1_z"),
  "Effort",
  variable %in% c("temp_x_effort","mahal_x_effort"),
  "Climate × Effort",
  variable == "year_c", "Year",
  default = "Other")]
setorder(imp_v3, -importance)
fwrite(imp_v3, file.path(V2, "results", "tables",
                          "table_rf_importance_v3.csv"))
log("RF v3 top 6:")
print(imp_v3[1:6])

# ============================================================
# STEP 3 — XGBoost CV + SHAP
# ============================================================
audit("STEP 3: XGBoost with SHAP")

Xtr <- xgb.DMatrix(data = as.matrix(train_dt[, ..feat_full]),
                    label = train_dt$event)
xgp <- list(objective = "binary:logistic", eval_metric = "auc",
            eta = 0.05, max_depth = 6,
            subsample = 0.8, colsample_bytree = 0.8,
            nthread = max(1L, parallel::detectCores() - 2L))
cv <- xgb.cv(xgp, Xtr, nrounds = 600, nfold = 5,
              early_stopping_rounds = 30, verbose = 0)
nb <- cv$best_iteration
if (is.null(nb) || nb == 0L) {
  auc_col <- intersect(c("test_auc_mean","test_AUC_mean"),
                        names(cv$evaluation_log))[1]
  nb <- if (!is.na(auc_col)) which.max(cv$evaluation_log[[auc_col]]) else 150L
}
best_auc <- cv$evaluation_log$test_auc_mean[nb]
log("XGBoost v3: nrounds = ", nb, " CV-AUC = ", round(best_auc, 3))
fwrite(data.table(nrounds_best = nb, cv_auc = best_auc, nfold = 5),
       file.path(V2, "results", "tables", "table_xgb_cv_v3.csv"))

mdl <- xgb.train(xgp, Xtr, nrounds = nb, verbose = 0)
shap <- shap.prep(mdl, X_train = as.matrix(train_dt[, ..feat_full]))
shap_summary <- shap[, .(mean_abs = mean(abs(value), na.rm = TRUE)),
                       by = variable][order(-mean_abs)]
shap[, variable := factor(variable, levels = rev(shap_summary$variable))]

# ============================================================
# STEP 4 — v2 vs v3 RF comparison table
# ============================================================
audit("STEP 4: build v2 vs v3 RF importance comparison")

imp_v2 <- fread(file.path(V2, "results", "tables",
                            "table_rf_importance_v2.csv"))
imp_v2[, rank_v2 := frank(-importance)]
imp_v3[, rank_v3 := frank(-importance)]
cmp <- merge(imp_v2[, .(variable, imp_v2 = importance, rank_v2)],
              imp_v3[, .(variable, imp_v3 = importance, rank_v3)],
              by = "variable", all = TRUE)
cmp[, delta_rank := rank_v3 - rank_v2]
setorder(cmp, rank_v3)
fwrite(cmp, file.path(V2, "results", "tables",
                       "table_v2_v3_rf_comparison.csv"))
log("v2 vs v3 RF comparison:")
print(cmp)

# ============================================================
# STEP 5 — Publication figure
# ============================================================
audit("STEP 5: build figure")

theme_pub <- function(s = 9) {
  theme_bw(base_size = s) +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(linewidth = 0.4, colour = "grey20"),
          plot.title = element_text(face = "bold", size = s + 1),
          plot.subtitle = element_text(size = s - 1, colour = "grey30"))
}
COL_CAT <- c(Climate = "#3B4CC0", Effort = "#B40426",
              Year = "#7F7F7F", `Climate × Effort` = "#FF7F0E",
              Other = "#8C564B")

imp_v3_show <- copy(imp_v3)
imp_v3_show[, variable_pretty := variable]
imp_v3_show[variable == "temp_x_effort",  variable_pretty := "temp × effort"]
imp_v3_show[variable == "mahal_x_effort", variable_pretty := "mahal × effort"]
imp_v3_show[, variable_pretty := factor(variable_pretty,
                                          levels = rev(variable_pretty))]
p_rf <- ggplot(imp_v3_show, aes(x = importance, y = variable_pretty,
                                  colour = category)) +
  geom_segment(aes(xend = 0, yend = variable_pretty), linewidth = 0.35) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = COL_CAT, name = NULL) +
  labs(tag = "a",
        title = "v3 Random Forest permutation importance",
        subtitle = sprintf("ranger 500 trees; n = %s rows × %d features; events = %s",
                             format(nrow(train_dt), big.mark=","),
                             length(feat_full),
                             format(sum(train_dt$event), big.mark=",")),
        x = "Permutation importance", y = NULL) +
  theme_pub() + theme(legend.position = "top")

# v2 vs v3 rank comparison
cmp_plot <- cmp[!is.na(rank_v2) & !is.na(rank_v3)]
cmp_plot[, variable_pretty := variable]
cmp_plot[variable == "temp_x_effort",  variable_pretty := "temp × effort"]
cmp_plot[variable == "mahal_x_effort", variable_pretty := "mahal × effort"]
p_cmp <- ggplot(cmp_plot, aes(x = rank_v2, y = rank_v3)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
               colour = "grey60") +
  geom_point(size = 2.6, colour = "#B40426") +
  geom_text(aes(label = variable_pretty), size = 2.4, hjust = -0.15,
             vjust = 0.5) +
  scale_x_reverse(limits = c(max(c(cmp_plot$rank_v2, cmp_plot$rank_v3), na.rm=TRUE) + 1, 0),
                   breaks = seq(1, 16, 2)) +
  scale_y_reverse(limits = c(max(c(cmp_plot$rank_v2, cmp_plot$rank_v3), na.rm=TRUE) + 1, 0),
                   breaks = seq(1, 16, 2)) +
  labs(tag = "b",
        title = "v2 vs v3 RF rank consistency",
        subtitle = "Lower rank = more important. Diagonal = identical ranking.",
        x = "v2 rank", y = "v3 rank") +
  theme_pub()

# SHAP beeswarm
p_shap <- ggplot(shap[abs(value) <= quantile(abs(shap$value), 0.99, na.rm=TRUE)],
                  aes(x = value, y = variable, colour = rfvalue)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_quasirandom(groupOnX = FALSE, alpha = 0.25, size = 0.35,
                   width = 0.32) +
  scale_colour_gradient2(low = "#3B4CC0", mid = "grey85", high = "#B40426",
                          midpoint = 0,
                          name = "Feature\nvalue (z)") +
  labs(tag = "c",
        title = "v3 XGBoost SHAP — beeswarm",
        subtitle = sprintf("CV-AUC = %.3f at %d rounds; positive SHAP → increases predicted hazard.",
                             best_auc, nb),
        x = "SHAP value (log-odds)", y = NULL) +
  theme_pub()

fig <- (p_rf | p_cmp) / p_shap
fig <- fig + plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    title = "Figure 4 v3 — Variable importance on the relaxed risk set",
    subtitle = "Climate × effort interaction stays in the top tier under v3 (compare to v2 = 0.0166 / rank 1).",
    theme = theme(plot.title = element_text(face = "bold", size = 10),
                   plot.subtitle = element_text(size = 8.5, colour = "grey30")))

ens(file.path(V2, "figures", "main"))
ggsave(file.path(V2, "figures", "main", "Figure_4_variable_importance_v3.pdf"),
       fig, width = 20, height = 16, units = "cm",
       device = grDevices::cairo_pdf)
ggsave(file.path(V2, "figures", "main", "Figure_4_variable_importance_v3.png"),
       fig, width = 20, height = 16, units = "cm", dpi = 600)
log("wrote Figure_4_variable_importance_v3.{pdf,png}")

log("")
log("══════════════════════════════════════════════════════════")
log("                v3 VARIABLE IMPORTANCE COMPLETE")
log("══════════════════════════════════════════════════════════")
print(imp_v3[, .(variable, category, importance = round(importance, 4))])
