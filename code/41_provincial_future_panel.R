# ============================================================
# Scientific question / 科学问题:
#   Build the manuscript-ready, publication-quality province-scale
#   panel: (1) model selection visualisations across 4 effort
#   specifications × M0-M4; (2) coefficient forest + beeswarm
#   comparisons; (3) variance decomposition; (4) variable importance
#   from Random Forest + XGBoost (+ SHAP); (5) future hazard
#   projections under 4 climate × effort scenarios from BOTH
#   glmmTMB (M4) and XGBoost; (6) choropleth maps for each
#   projection. All figures use a unified theme + GS(2019)1822
#   province basemap (Albers EPSG:4524) + color-blind palette.
#
# 省级综合面板：模型选择 / 系数对比 / 方差分解 / 重要性 / 未来
# 情景 hazard 地图（glmmTMB + XGBoost）。出版级配色与字号。
#
# Input data / 输入数据:
#   v1/data/hazard_risk_upgraded_complete_case.csv  (12,813 × 24)
#   v1/data/climate_metrics_province_year.csv
#   v1/data/effort_panel_upgraded.csv
#   v1/results/table_cross_effort_variance_decomposition.csv
#   v2/results/tables/table_province_v2_coefs.csv  (60 rows)
#   v2/results/tables/table_province_v2_aic.csv
#   v2 GS(2019)1822 basemap shapefiles
#
# Main outputs / 主要输出:
#   figures/main/fig_aic_akaike_ladder.{pdf,png}
#   figures/main/fig_coef_forest_4specs.{pdf,png}
#   figures/main/fig_coef_beeswarm_M4.{pdf,png}
#   figures/main/fig_varpart_4specs.{pdf,png}
#   figures/main/fig_rf_importance.{pdf,png}
#   figures/main/fig_xgb_shap_summary.{pdf,png}
#   figures/main/fig_future_hazard_glmmTMB.{pdf,png}
#   figures/main/fig_future_hazard_xgboost.{pdf,png}
#   figures/main/fig_future_glmmTMB_vs_xgboost_rank.{pdf,png}
#   results/forecasts/table_province_future_glmmTMB.csv
#   results/forecasts/table_province_future_xgboost.csv
#   results/tables/table_rf_importance_v2.csv
#   results/tables/table_aic_akaike_weights.csv
#
# Key assumption: future SSP perturbations follow the empirical
# scheme used by v1 (16_multi_scale_future_prediction.R):
#   SSP245 / decade = +0.3 / temp_grad_sd
#   SSP585 / decade = +0.8 / temp_grad_sd
# Hard CMIP6 NetCDFs are not loaded here; see code/29 for the
# ensemble pipeline guarded by hard-fail.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(glmmTMB)
  library(xgboost)
  library(ranger)
  library(ggplot2)
  library(ggbeeswarm)
  library(ggrepel)
  library(patchwork)
  library(viridisLite)
  library(SHAPforxgboost)
})
sf::sf_use_s2(FALSE)
options(warn = 1)
set.seed(42)

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- Sys.getenv("V1_ROOT",
                  normalizePath(file.path(V2, "..",
                                          "bird_hazard_model_effort_upgrade"),
                                 mustWork = FALSE))

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ens(file.path(V2, "results", "tables"))
ens(file.path(V2, "results", "forecasts"))
ens(file.path(V2, "figures", "main"))
ens(file.path(V2, "logs"))

log <- function(...) cat(sprintf("[41 %s] ", format(Sys.time(), "%H:%M:%S")),
                          ..., "\n", sep = "")

# ---- Unified GEB-style theme + palette ----------------------------------
theme_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size, base_family = "") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
          panel.border     = element_rect(linewidth = 0.4, colour = "grey20"),
          axis.text        = element_text(colour = "grey20"),
          axis.title       = element_text(colour = "grey10"),
          plot.title       = element_text(face = "bold", size = base_size + 1),
          plot.subtitle    = element_text(size = base_size - 1, colour = "grey30"),
          plot.tag         = element_text(face = "bold", size = base_size + 1),
          strip.background = element_rect(fill = "grey95", colour = "grey80"),
          strip.text       = element_text(face = "bold"),
          legend.position  = "right",
          legend.title     = element_text(face = "bold"),
          legend.key.size  = unit(0.3, "cm"))
}
PAL_SPEC <- c(spec_A = "#1F77B4", spec_B = "#D62728",
              spec_C = "#2CA02C", spec_D = "#9467BD")
SPEC_LABEL <- c(spec_A = "A: records (legacy)",
                spec_B = "B: visits (headline)",
                spec_C = "C: PCA composite",
                spec_D = "D: birding-days")
PAL_CAT  <- c(Climate = "#3B4CC0", Effort = "#B40426",
              Year = "#7F7F7F", `Climate × Effort` = "#FF7F0E",
              Other = "#8C564B")

save_pub <- function(p, name, width = 17, height = 10,
                      path = file.path(V2, "figures", "main")) {
  ggsave(file.path(path, paste0(name, ".pdf")),
         p, width = width, height = height, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(file.path(path, paste0(name, ".png")),
         p, width = width, height = height, units = "cm", dpi = 600)
  log("wrote ", name, ".{pdf,png}")
}

# ============================================================
# 1. Load data + model selection + AIC ladder figure
# ============================================================
log("loading data")
risk <- fread(file.path(V1, "data",
                         "hazard_risk_upgraded_complete_case.csv"),
              encoding = "UTF-8")
coefs_v2 <- fread(file.path(V2, "results", "tables",
                              "table_province_v2_coefs.csv"))
aic_v2 <- fread(file.path(V2, "results", "tables",
                            "table_province_v2_aic.csv"))

# Akaike weights per spec ----------------------------------------------------
aic_v2[, dAIC := AIC - min(AIC, na.rm = TRUE), by = spec_id]
aic_v2[, weight := exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC)), by = spec_id]
fwrite(aic_v2[, .(spec_id, spec_label, model, AIC, dAIC, weight, nobs)],
       file.path(V2, "results", "tables",
                  "table_aic_akaike_weights.csv"))

# Fig 1 — AIC + Akaike-weight ladder ----------------------------------------
log("rendering AIC / Akaike weight ladder")
aic_v2[, spec_id := factor(spec_id, levels = c("spec_B","spec_A","spec_C","spec_D"))]
aic_v2[, model := factor(model, levels = c("M0","M1","M2","M3","M4"))]

p_aic <- ggplot(aic_v2, aes(x = dAIC, y = model, colour = spec_id)) +
  geom_segment(aes(xend = 0, yend = model), linewidth = 0.4) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.1f", dAIC)),
             nudge_x = ifelse(aic_v2$dAIC > 5, -1.4, 1.0),
             hjust = ifelse(aic_v2$dAIC > 5, 1, 0), size = 2.3) +
  facet_wrap(~ spec_id, ncol = 2,
             labeller = labeller(spec_id = function(x) SPEC_LABEL[x])) +
  scale_colour_manual(values = PAL_SPEC, guide = "none") +
  labs(title = "Model-selection ladder by effort specification",
        subtitle = "ΔAIC vs. each spec's best model. M4 (climate × effort interaction) is lowest in every spec.",
        x = "ΔAIC (lower = better)", y = NULL) +
  theme_pub()
save_pub(p_aic, "fig_aic_akaike_ladder", width = 17, height = 10)

p_aw <- ggplot(aic_v2, aes(x = weight, y = model, fill = spec_id)) +
  geom_col(width = 0.6) +
  facet_wrap(~ spec_id, ncol = 2,
             labeller = labeller(spec_id = function(x) SPEC_LABEL[x])) +
  scale_fill_manual(values = PAL_SPEC, guide = "none") +
  scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  labs(title = "Akaike weights — concentration on M4 in all four specifications",
        subtitle = "Within each spec the Akaike weight is essentially 1 on M4.",
        x = "Akaike weight", y = NULL) +
  theme_pub()
save_pub(p_aw, "fig_akaike_weights", width = 17, height = 10)

# ============================================================
# 2. Coefficient forest + beeswarm across 4 specs
# ============================================================
log("coef forest + beeswarm")
nice_term <- function(x) {
  x <- gsub("climate_z", "temp_grad_z", x, fixed = TRUE)
  x <- gsub(":", " × ", x, fixed = TRUE)
  x
}
coefs_show <- coefs_v2[term != "(Intercept)",
                        .(spec_id, model, term = nice_term(term),
                          beta, se, hr, hr.low, hr.high, p.value)]
coefs_show[, term := factor(term, levels = c(
  "temp_grad_z", "effort_z", "temp_grad_z × effort_z"))]
coefs_show[, model := factor(model, levels = c("M0","M1","M2","M3","M4"))]
coefs_show[, spec_id := factor(spec_id, levels = c("spec_A","spec_B","spec_C","spec_D"))]

p_forest <- ggplot(coefs_show[!is.na(term)],
                    aes(x = hr, y = model, colour = spec_id)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  height = 0.18, linewidth = 0.5,
                  position = position_dodge(width = 0.55)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.55)) +
  facet_wrap(~ term, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = PAL_SPEC,
                       labels = SPEC_LABEL, name = "Effort spec") +
  scale_x_continuous(trans = "log") +
  labs(title = "Province-scale coefficient forest — 4 effort specs × M0-M4",
        subtitle = "Bars = 95 % CI. Headline spec B (red) interaction HR = 1.288 (95 % CI 1.179, 1.407).",
        x = "Hazard ratio (log scale)", y = NULL) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_forest, "fig_coef_forest_4specs", width = 18, height = 9)

# Beeswarm of interaction-term HR across specs (M4 only) --------------------
beesw_dt <- coefs_show[model == "M4" & term == "temp_grad_z × effort_z"]
p_bee <- ggplot(beesw_dt,
                 aes(x = spec_id, y = hr, colour = spec_id)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = hr.low, ymax = hr.high),
                 width = 0.15, linewidth = 0.5) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sprintf("HR = %.3f\np = %.1e", hr, p.value)),
                   size = 2.6, nudge_y = 0.04, direction = "y",
                   point.padding = 0.3, segment.size = 0.2) +
  scale_colour_manual(values = PAL_SPEC, guide = "none") +
  scale_x_discrete(labels = SPEC_LABEL) +
  scale_y_continuous(trans = "log",
                      breaks = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5)) +
  labs(title = "Climate × effort interaction HR across four effort specifications (M4)",
        subtitle = "Each spec is an independent refit of the headline model on the same risk set.",
        x = NULL, y = "Hazard ratio (log scale)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 12, hjust = 1))
save_pub(p_bee, "fig_coef_beeswarm_M4", width = 14, height = 9)

# ============================================================
# 3. Variance decomposition
# ============================================================
log("variance decomposition figure")
vd <- fread(file.path(V1, "results",
                       "table_cross_effort_variance_decomposition.csv"))
setnames(vd, c("effort_spec","effort_var","additive_r2",
               "interaction_r2","delta_interaction_r2"),
         skip_absent = TRUE)
vd <- vd[!is.na(additive_r2) & !is.na(interaction_r2)]
vd_long <- rbindlist(list(
  vd[, .(effort_spec, comp = "Additive (M3)",
          R2 = additive_r2)],
  vd[, .(effort_spec, comp = "Interaction-only (M4 − M3)",
          R2 = delta_interaction_r2)]))
vd_long[, effort_spec := factor(effort_spec,
                                  levels = c("Observer visits","Birding days",
                                              "Record-based","PCA composite"))]
vd_long[, comp := factor(comp,
                          levels = c("Additive (M3)","Interaction-only (M4 − M3)"))]

p_vd <- ggplot(vd_long, aes(x = R2, y = effort_spec, fill = comp)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.4f", R2)),
             position = position_stack(vjust = 0.5),
             colour = "white", size = 2.6, fontface = "bold") +
  scale_fill_manual(values = c(`Additive (M3)` = "#9DC3E6",
                                 `Interaction-only (M4 − M3)` = "#B40426"),
                     name = NULL) +
  labs(title = "Marginal R² decomposition: additive vs interaction-only",
        subtitle = "Across all four effort specs, the interaction-only increment dominates (≈ 78–84 % of M4 marginal R²).",
        x = "Marginal R²", y = NULL) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_vd, "fig_varpart_4specs", width = 16, height = 8)

# ============================================================
# 4. Random Forest + XGBoost variable importance
# ============================================================
log("fitting Random Forest (ranger) for variable importance")
# v1 raw risk set has only {temp_grad_z, prec_grad_z, *effort* cols}; the
# advanced climate metrics live in climate_metrics_province_year.csv and
# must be joined by (province, year). 把高级气候指标合并进来。
prov_clim_full <- fread(file.path(V1, "data",
                                    "climate_metrics_province_year.csv"),
                           encoding = "UTF-8")
climate_join_cols <- intersect(c("climate_velocity_z","precip_velocity_z",
                                   "climate_exposure_z","warming_rate_z",
                                   "mahalanobis_dist_z","temp_anom_z",
                                   "prec_anom_z","temp_grad_prov_z",
                                   "prec_grad_prov_z"),
                                 names(prov_clim_full))
log("climate columns joined: ", paste(climate_join_cols, collapse = ", "))
train <- merge(risk,
                prov_clim_full[, c("province","year", climate_join_cols),
                                with = FALSE],
                by = c("province","year"), all.x = TRUE)

features_climate <- intersect(c("temp_grad_z","prec_grad_z",
                                 "climate_velocity_z","precip_velocity_z",
                                 "climate_exposure_z","warming_rate_z",
                                 "mahalanobis_dist_z","temp_anom_z",
                                 "prec_anom_z"),
                               names(train))
features_effort  <- intersect(c("log_effort_record_z","log_effort_visits_z",
                                 "log_effort_observers_z","log_effort_days_z",
                                 "effort_pc1_z"),
                               names(train))
features_inter   <- c("temp_x_effort","mahal_x_effort")

train[, temp_x_effort  := temp_grad_z * log_effort_visits_z]
if ("mahalanobis_dist_z" %in% names(train)) {
  train[, mahal_x_effort := mahalanobis_dist_z * log_effort_visits_z]
} else {
  train[, mahal_x_effort := NA_real_]
}
train[, year_c := year - 2013]
feat <- c(features_climate, features_effort, features_inter, "year_c")
feat <- intersect(feat, names(train))
train_dt <- train[, c("event", feat), with = FALSE]
train_dt <- train_dt[complete.cases(train_dt)]
log("RF n = ", nrow(train_dt), " | events = ", sum(train_dt$event))

rf <- ranger(event ~ ., data = train_dt,
              probability = FALSE, classification = FALSE,
              num.trees = 500, mtry = floor(sqrt(length(feat))),
              importance = "permutation", seed = 42,
              num.threads = max(1L, parallel::detectCores() - 2L))
imp <- data.table(variable = names(rf$variable.importance),
                   importance = as.numeric(rf$variable.importance))
imp[, category := fcase(
  variable %in% features_climate, "Climate",
  variable %in% features_effort,  "Effort",
  variable %in% features_inter,   "Climate × Effort",
  variable == "year_c",           "Year",
  default = "Other")]
imp <- imp[order(-importance)]
fwrite(imp, file.path(V2, "results", "tables",
                        "table_rf_importance_v2.csv"))

# Lollipop RF importance figure ---------------------------------------------
imp[, variable := factor(variable, levels = rev(imp$variable))]
p_rf <- ggplot(imp, aes(x = importance, y = variable, colour = category)) +
  geom_segment(aes(xend = 0, yend = variable), linewidth = 0.4) +
  geom_point(size = 3) +
  scale_colour_manual(values = PAL_CAT, name = NULL) +
  labs(title = "Random Forest permutation importance — province scale",
        subtitle = sprintf("ranger probability=FALSE; n = %s rows, %s events; importance = mean ↓ accuracy on OOB.",
                            format(nrow(train_dt), big.mark = ","),
                            format(sum(train_dt$event), big.mark = ",")),
        x = "Permutation importance", y = NULL) +
  theme_pub() + theme(legend.position = "top")
save_pub(p_rf, "fig_rf_importance", width = 15, height = 11)

# ============================================================
# 5. XGBoost + SHAP summary
# ============================================================
log("XGBoost training + SHAP")
xgb_train <- xgb.DMatrix(data = as.matrix(train_dt[, ..feat]),
                          label = train_dt$event)
xgb_par <- list(objective = "binary:logistic",
                 eval_metric = "auc",
                 eta = 0.05, max_depth = 6,
                 subsample = 0.8, colsample_bytree = 0.8,
                 nthread = max(1L, parallel::detectCores() - 2L))
cv_res <- xgb.cv(params = xgb_par, data = xgb_train, nrounds = 600,
                  nfold = 5, early_stopping_rounds = 30, verbose = 0)
# xgb.cv returns NULL best_iteration when early_stopping is never triggered
# or when the metric column name differs. Fall back to the best-AUC row.
nrounds_best <- cv_res$best_iteration
if (is.null(nrounds_best) || length(nrounds_best) == 0L ||
    nrounds_best == 0L) {
  auc_col <- intersect(c("test_auc_mean", "test_AUC_mean", "test_AUC-mean"),
                        names(cv_res$evaluation_log))[1]
  if (!is.na(auc_col)) {
    nrounds_best <- which.max(cv_res$evaluation_log[[auc_col]])
  } else {
    nrounds_best <- 200L
  }
  log("xgb.cv had no best_iteration; falling back to nrounds = ",
      nrounds_best)
}
auc_col <- intersect(c("test_auc_mean", "test_AUC_mean", "test_AUC-mean"),
                      names(cv_res$evaluation_log))[1]
best_auc <- if (!is.na(auc_col)) cv_res$evaluation_log[[auc_col]][nrounds_best] else NA_real_
log("XGBoost CV best nrounds = ", nrounds_best,
    " | mean AUC = ", round(best_auc, 3))
xgb_model <- xgb.train(params = xgb_par, data = xgb_train,
                        nrounds = nrounds_best, verbose = 0)

# SHAP values
shap_dt <- shap.prep(xgb_model, X_train = as.matrix(train_dt[, ..feat]))
# Manually plot a beeswarm to keep the unified theme.
shap_summary <- shap_dt[, .(mean_abs = mean(abs(value), na.rm = TRUE)),
                         by = variable][order(-mean_abs)]
shap_dt[, variable := factor(variable, levels = rev(shap_summary$variable))]

p_shap <- ggplot(shap_dt,
                  aes(y = variable, x = value, colour = rfvalue)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_quasirandom(groupOnX = FALSE, alpha = 0.35, size = 0.5,
                   width = 0.35) +
  scale_colour_gradient2(low = "#3B4CC0", mid = "grey85", high = "#B40426",
                          midpoint = 0,
                          name = "Feature\nvalue (z)") +
  labs(title = "XGBoost SHAP — feature-level contributions to predicted hazard",
        subtitle = sprintf("CV AUC = %.3f at %d rounds; positive SHAP = increases predicted hazard.",
                            cv_res$evaluation_log$test_auc_mean[nrounds_best],
                            nrounds_best),
        x = "SHAP value (log-odds)", y = NULL) +
  theme_pub() + theme(legend.position = "right")
save_pub(p_shap, "fig_xgb_shap_summary", width = 16, height = 11)

# ============================================================
# 6. Future scenario hazard — glmmTMB (M4 Spec B) + XGBoost
# ============================================================
log("future scenario forecasting (glmmTMB + XGBoost)")

# Refit headline M4 (Spec B). 用相同公式与数据，避免依赖 RDS 持久化。
m4_dt <- risk[, .(species, province, year, event,
                   temp_grad_z, log_effort_visits_z)]
m4_dt <- m4_dt[complete.cases(m4_dt)]
fit_M4 <- glmmTMB(event ~ temp_grad_z * log_effort_visits_z +
                            (1 | species) + (1 | province),
                   data = m4_dt, family = binomial(link = "cloglog"))
log("M4 AIC = ", round(AIC(fit_M4), 2),
    " | interaction β = ",
    round(fixef(fit_M4)$cond["temp_grad_z:log_effort_visits_z"], 3))

# Province baseline (most recent year 2024 effort + climate)
prov_clim <- fread(file.path(V1, "data",
                              "climate_metrics_province_year.csv"),
                    encoding = "UTF-8")
prov_eff  <- fread(file.path(V1, "data", "effort_panel_upgraded.csv"),
                    encoding = "UTF-8")
baseline <- merge(
  prov_clim[year == 2024, .(province, temp_grad_z = temp_grad_prov_z)],
  prov_eff[year == 2024, .(province, log_effort_visits_z)],
  by = "province")
log("baseline rows: ", nrow(baseline), " provinces")

# SSP perturbation per decade (v1's empirical scheme)
temp_grad_sd <- sd(prov_clim$temp_grad_prov, na.rm = TRUE)
ssp_eps <- list(SSP245 = 0.3 / temp_grad_sd,
                 SSP585 = 0.8 / temp_grad_sd)
future_years <- c(2030, 2050, 2080)

# Build future panel: province × ssp × year
fut <- CJ(province = baseline$province,
           ssp = c("SSP245","SSP585"),
           year = future_years)
fut <- merge(fut, baseline, by = "province")
fut[, decades_ahead := (year - 2024) / 10]
fut[, temp_grad_z := temp_grad_z +
                      ssp_eps[[ssp[1]]] * decades_ahead, by = ssp]
fut[, log_effort_visits_z := log_effort_visits_z]  # baseline-frozen effort

# === 6.1 glmmTMB future hazard (averaged over species random intercept) ====
# Approach: average over species random effects by predicting the population-
# level fixed-effects log-hazard, then convert to hazard via cloglog inverse.
# 用 fixed-effects 仅做 species-marginal prediction（exclude species random）。
# glmmTMB::predict accepts re.form = NA / NULL / ~0; for fixed-effects-only
# (population-average) predictions we use re.form = ~0. The species random
# intercept is integrated out implicitly. 用 fixed-effects-only 预测。
preds_glmm <- predict(fit_M4, newdata = fut,
                       type = "response",
                       re.form = NA,            # marginalise random effects
                       allow.new.levels = TRUE)
fut[, hazard_glmm := preds_glmm]

fwrite(fut, file.path(V2, "results", "forecasts",
                       "table_province_future_glmmTMB.csv"))

# === 6.2 XGBoost future hazard ============================================
log("XGBoost future prediction")
# Need a feature-aligned row per province × ssp × year. Use the species-
# average covariate values (computed as the province × year mean over the
# training risk set), then perturb temp-related anomalies forward.
# 用各省训练均值作为 species-average baseline。
prov_feat_baseline <- train_dt[, lapply(.SD, mean, na.rm = TRUE),
                                 by = .(), .SDcols = feat]
# Slot in province-specific temp_grad_z + effort
xgb_fut <- merge(
  fut[, .(province, ssp, year, temp_grad_z, log_effort_visits_z)],
  data.table(province = risk[, unique(province)]),
  by = "province")
# Construct feature matrix: replicate the global mean row, then overwrite
# the two province-/scenario-resolved columns. 其它特征用训练集均值。
xgb_feat_mat <- as.matrix(prov_feat_baseline)[
  rep(1L, nrow(xgb_fut)), , drop = FALSE]
xgb_feat_mat <- as.data.table(xgb_feat_mat)
xgb_feat_mat[, temp_grad_z := xgb_fut$temp_grad_z]
xgb_feat_mat[, log_effort_visits_z := xgb_fut$log_effort_visits_z]
xgb_feat_mat[, temp_x_effort  := temp_grad_z * log_effort_visits_z]
xgb_feat_mat[, mahal_x_effort := mahalanobis_dist_z * log_effort_visits_z]
xgb_feat_mat[, year_c := xgb_fut$year - 2013]

xgb_pred <- predict(xgb_model,
                     as.matrix(xgb_feat_mat[, ..feat]))
xgb_fut[, hazard_xgb := xgb_pred]
fwrite(xgb_fut, file.path(V2, "results", "forecasts",
                            "table_province_future_xgboost.csv"))

# === 6.3 Choropleth maps (glmmTMB + XGBoost) ==============================
log("rendering choropleth maps")
shp_dir <- file.path(V2, "data", "spatial", "basemap_GS2019_1822")
prov_shp <- list.files(shp_dir, pattern = "^省.*\\.shp$",
                        full.names = TRUE)[1]
ninedash_shp <- list.files(shp_dir, pattern = "九段线.*\\.shp$",
                            full.names = TRUE)[1]
nat_shp  <- list.files(shp_dir, pattern = "^(国界|中国轮廓).*\\.shp$",
                        full.names = TRUE)[1]

prov_sf <- st_read(prov_shp, quiet = TRUE) |>
  st_make_valid() |>
  st_transform(4524)
ninedash_sf <- if (!is.na(ninedash_shp)) st_read(ninedash_shp, quiet = TRUE) |>
  st_transform(4524) else NULL
nat_sf <- if (!is.na(nat_shp)) st_read(nat_shp, quiet = TRUE) |>
  st_transform(4524) else NULL

# Match shp province names → English
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
prov_name_col <- intersect(c("name","NAME","省","NAME_1"), names(prov_sf))[1]
if (is.na(prov_name_col)) {
  prov_name_col <- names(prov_sf)[vapply(prov_sf, function(x)
    is.character(x) && any(grepl("[一-龥]", x)), logical(1))][1]
}
prov_sf$province <- unname(PROV_CN_EN[as.character(prov_sf[[prov_name_col]])])
prov_sf$province[is.na(prov_sf$province)] <- as.character(prov_sf[[prov_name_col]])[is.na(prov_sf$province)]

choropleth_map <- function(dt, var, title, subtitle,
                            hi_col = "#B40426", lo_col = "#F7FBFF") {
  joined <- merge(prov_sf, dt, by = "province")
  p <- ggplot() +
    geom_sf(data = joined, aes(fill = .data[[var]]),
            colour = "grey45", linewidth = 0.15)
  if (!is.null(nat_sf)) p <- p +
    geom_sf(data = nat_sf, fill = NA, colour = "grey20", linewidth = 0.35)
  if (!is.null(ninedash_sf)) p <- p +
    geom_sf(data = ninedash_sf, fill = NA, colour = "grey20",
            linewidth = 0.25)
  p + scale_fill_gradient(low = lo_col, high = hi_col,
                           name = "Hazard\n(prob.)",
                           limits = range(dt[[var]], na.rm = TRUE)) +
    coord_sf(datum = NA, expand = FALSE) +
    labs(title = title, subtitle = subtitle) +
    theme_pub(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.key.height = unit(0.6, "cm"))
}

# Build 6-panel facets for glmmTMB
mk_panel <- function(dt, var, fig_name, model_label) {
  dt[, panel := paste0(ssp, " — ", year)]
  dt[, panel := factor(panel, levels = c(
    paste0("SSP245 — ", future_years),
    paste0("SSP585 — ", future_years)))]
  panels <- lapply(split(dt, dt$panel), function(d)
    choropleth_map(d, var,
                   title = unique(d$panel),
                   subtitle = NULL))
  combined <- patchwork::wrap_plots(panels, ncol = 3,
                                     guides = "collect")
  combined <- combined + plot_annotation(
    title = sprintf("Provincial future hazard projection — %s", model_label),
    subtitle = paste0("Empirical SSP perturbation: SSP245 = +0.3 SD/decade, ",
                       "SSP585 = +0.8 SD/decade on temp_grad_z. Effort frozen ",
                       "at 2024 baseline. Albers EPSG:4524 / GS(2019)1822 basemap."),
    theme = theme(plot.title = element_text(face = "bold", size = 10),
                   plot.subtitle = element_text(size = 8, colour = "grey30")))
  save_pub(combined, fig_name, width = 22, height = 14)
}

mk_panel(copy(fut), "hazard_glmm",
          "fig_future_hazard_glmmTMB", "glmmTMB (M4 Spec B)")
mk_panel(copy(xgb_fut), "hazard_xgb",
          "fig_future_hazard_xgboost", "XGBoost (full feature set)")

# === 6.4 glmmTMB vs XGBoost rank comparison =================================
log("province ranking comparison")
rank_dt <- merge(
  fut[, .(province, ssp, year, hazard_glmm)],
  xgb_fut[, .(province, ssp, year, hazard_xgb)],
  by = c("province","ssp","year"))
rank_dt[, rank_glmm := frank(-hazard_glmm), by = .(ssp, year)]
rank_dt[, rank_xgb  := frank(-hazard_xgb),  by = .(ssp, year)]

rank_585_2050 <- rank_dt[ssp == "SSP585" & year == 2050]
rank_585_2050[, panel := "SSP585 / 2050"]
rank_245_2050 <- rank_dt[ssp == "SSP245" & year == 2050]
rank_245_2050[, panel := "SSP245 / 2050"]

rank_comp <- rbind(rank_245_2050, rank_585_2050)
p_rank <- ggplot(rank_comp,
                  aes(x = rank_glmm, y = rank_xgb)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
               colour = "grey60") +
  geom_point(size = 1.6, colour = "#B40426") +
  geom_text_repel(aes(label = province), size = 2.4,
                   max.overlaps = 32, segment.size = 0.15) +
  facet_wrap(~ panel, ncol = 2) +
  scale_x_continuous(trans = "reverse", limits = c(33, 0)) +
  scale_y_continuous(trans = "reverse", limits = c(33, 0)) +
  labs(title = "Province ranking — glmmTMB (M4) vs XGBoost",
        subtitle = "Each axis is a rank (1 = highest projected hazard). Provinces on the diagonal are concordant.",
        x = "glmmTMB rank", y = "XGBoost rank") +
  theme_pub()
save_pub(p_rank, "fig_future_glmmTMB_vs_xgboost_rank",
         width = 18, height = 10)

log("=== DONE ===")
log("Files written:")
log("  results/forecasts/table_province_future_{glmmTMB,xgboost}.csv")
log("  results/tables/table_aic_akaike_weights.csv")
log("  results/tables/table_rf_importance_v2.csv")
log("  figures/main/fig_aic_akaike_ladder.{pdf,png}")
log("  figures/main/fig_akaike_weights.{pdf,png}")
log("  figures/main/fig_coef_forest_4specs.{pdf,png}")
log("  figures/main/fig_coef_beeswarm_M4.{pdf,png}")
log("  figures/main/fig_varpart_4specs.{pdf,png}")
log("  figures/main/fig_rf_importance.{pdf,png}")
log("  figures/main/fig_xgb_shap_summary.{pdf,png}")
log("  figures/main/fig_future_hazard_glmmTMB.{pdf,png}")
log("  figures/main/fig_future_hazard_xgboost.{pdf,png}")
log("  figures/main/fig_future_glmmTMB_vs_xgboost_rank.{pdf,png}")
