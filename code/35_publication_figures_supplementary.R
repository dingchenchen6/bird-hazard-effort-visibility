# ============================================================
# Scientific question / 科学问题:
#   Generate the 15 supplementary figures referenced from the v2
#   manuscript and supplementary materials. Each figure is built from
#   a results/diagnostics or results/sensitivity table; when a source
#   is unavailable the script produces a placeholder PDF so reviewers
#   see exactly which dependency is missing.
#   生成 v2 论文引用的 15 张补充图；找不到源数据时输出占位 PDF。
#
# Objective / 分析目标: see table below.
#
# S1  DHARMa global tests (uniformity, dispersion, outliers, ZI)
# S2  DHARMa spatial residual map (placeholder — also produced by 27)
# S3  VIF heatmap across M4 covariates
# S4  Partial-dependence top-8 features (XGBoost)
# S5  SHAP summary (beeswarm) + key dependence plots
# S6  PR-AUC by spatial-block fold
# S7  Reliability per fold
# S8  MAUP elasticity (also produced by 31)
# S9  Offset sensitivity (also produced by 32)
# S10 Grid-event redefinition before/after counts
# S11 CMIP6 per-model maps (also produced by 29)
# S12 Covariate-shift PSI bars
# S13 Effort panel temporal coverage by province
# S14 Species accumulation curves
# S15 Phylogenetic effort heatmap (genus level)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(ggplot2)
  library(patchwork)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_plots.R"))

ensure_dir(path_supp_fig())

placeholder <- function(name, msg) {
  p <- ggplot() + theme_void() +
    annotate("text", x = 1, y = 1, label = msg, size = 3.6) +
    xlim(0, 2) + ylim(0, 2)
  save_pub(p, name, path = path_supp_fig(), width = 14, height = 9)
}

# ---- S1 DHARMa global tests (placeholder if no diag file) -----------------
{
  src <- path_diagnostics("table_dharma_global.csv")
  if (file.exists(src)) {
    g <- fread(src)
    p <- ggplot(g, aes(x = test, y = p.value)) +
      geom_col(width = 0.5, fill = pal_cat[1]) +
      geom_hline(yintercept = 0.05, linetype = "dashed",
                  colour = "grey40") +
      labs(title = "S1 — DHARMa global tests (p-values)") +
      coord_flip() + theme_geb()
    save_pub(p, "figS1_dharma_global", path = path_supp_fig())
  } else {
    placeholder("figS1_dharma_global",
                "S1 placeholder — run DHARMa global tests in script 27.")
  }
}

# ---- S2 DHARMa spatial residual map (script 27 handles primary version) ---
if (!file.exists(path_supp_fig("figS2_dharma_spatial_residuals.pdf"))) {
  placeholder("figS2_dharma_spatial_residuals",
              "S2 placeholder — run script 27 to populate this figure.")
}

# ---- S3 VIF heatmap -------------------------------------------------------
{
  src <- path_diagnostics("table_vif_M4.csv")
  if (file.exists(src)) {
    vif <- fread(src)
    p <- ggplot(vif, aes(x = term, y = "VIF", fill = vif)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      geom_text(aes(label = round(vif, 2)), size = 3) +
      scale_fill_viridis_c(option = "C") +
      coord_fixed() +
      labs(title = "S3 — Variance inflation factors (M4)") +
      theme_geb()
    save_pub(p, "figS3_vif_heatmap", path = path_supp_fig(),
             width = 16, height = 5)
  } else {
    placeholder("figS3_vif_heatmap",
                "S3 placeholder — write table_vif_M4.csv via utils_models::vif_glmmTMB().")
  }
}

# ---- S4 PDP top-8 (XGBoost) ----------------------------------------------
{
  src <- path_diagnostics("table_pdp_top8.csv")
  if (file.exists(src)) {
    pdp <- fread(src)
    p <- ggplot(pdp, aes(x = value, y = yhat)) +
      geom_line(linewidth = 0.4) +
      facet_wrap(~ feature, scales = "free_x", ncol = 4) +
      labs(title = "S4 — Partial dependence (top 8 features)") +
      theme_geb()
    save_pub(p, "figS4_pdp_top8", path = path_supp_fig(),
             width = 18, height = 11)
  } else {
    placeholder("figS4_pdp_top8",
                "S4 placeholder — compute PDPs from xgboost via pdp::partial().")
  }
}

# ---- S5 SHAP summary ------------------------------------------------------
{
  src <- path_diagnostics("table_shap_long.csv")
  if (file.exists(src)) {
    shap <- fread(src)
    p <- ggplot(shap, aes(x = feature, y = shap_value, colour = feature_value)) +
      geom_jitter(width = 0.2, alpha = 0.4, size = 0.4) +
      scale_colour_viridis_c(option = "B") +
      coord_flip() +
      labs(title = "S5 — SHAP summary (beeswarm)") +
      theme_geb()
    save_pub(p, "figS5_shap_summary", path = path_supp_fig(),
             width = 16, height = 11)
  } else {
    placeholder("figS5_shap_summary",
                "S5 placeholder — export SHAP from xgboost::xgb.importance with type='gain'.")
  }
}

# ---- S6 PR-AUC by fold ----------------------------------------------------
{
  src <- path_diagnostics("table_spatial_block_cv.csv")
  if (file.exists(src)) {
    cv <- fread(src)
    p <- ggplot(cv, aes(x = factor(fold), y = auc_pr)) +
      geom_col(width = 0.6, fill = pal_cat[2]) +
      labs(title = "S6 — PR-AUC by spatial-block fold",
            x = "Fold", y = "PR-AUC") +
      theme_geb()
    save_pub(p, "figS6_pr_auc_by_fold", path = path_supp_fig(),
             width = 12, height = 7)
  } else {
    placeholder("figS6_pr_auc_by_fold",
                "S6 placeholder — run script 26 first.")
  }
}

# ---- S7 Reliability per fold ---------------------------------------------
{
  if (file.exists(path_diag_fig("spatial_block_cv_reliability.pdf"))) {
    message("[35] S7 already covered by script 26.")
  } else {
    placeholder("figS7_reliability_per_fold",
                "S7 placeholder — run script 26 to produce reliability curves.")
  }
}

# ---- S8 MAUP elasticity (handled by 31, alias if needed) ------------------
if (!file.exists(path_supp_fig("figS8_maup_elasticity.pdf"))) {
  placeholder("figS8_maup_elasticity",
              "S8 placeholder — run script 31 to produce MAUP elasticity figure.")
}

# ---- S9 Offset sensitivity (handled by 32 diagnostic) ---------------------
if (!file.exists(path_supp_fig("figS9_offset_sensitivity.pdf"))) {
  src_pdf <- path_diag_fig("offset_reformulation_diagnostic.pdf")
  if (file.exists(src_pdf)) {
    file.copy(src_pdf, path_supp_fig("figS9_offset_sensitivity.pdf"),
              overwrite = TRUE)
    src_png <- path_diag_fig("offset_reformulation_diagnostic.png")
    if (file.exists(src_png)) {
      file.copy(src_png, path_supp_fig("figS9_offset_sensitivity.png"),
                overwrite = TRUE)
    }
  } else {
    placeholder("figS9_offset_sensitivity",
                "S9 placeholder — run script 32 to produce offset sensitivity figure.")
  }
}

# ---- S10 Grid event redefinition --------------------------------------
{
  src <- path_diagnostics("table_grid_event_redefinition_summary.csv")
  if (file.exists(src)) {
    dt <- fread(src)
    dt_long <- melt(dt, id.vars = "grid_km",
                     measure.vars = c("n_raw_events", "n_events",
                                       "n_species_grid"),
                     variable.name = "metric", value.name = "value")
    p <- ggplot(dt_long, aes(x = factor(grid_km), y = value, fill = metric)) +
      geom_col(position = "dodge", width = 0.7) +
      scale_fill_manual(values = pal_cat[1:3]) +
      labs(title = "S10 — Grid event redefinition (v1 → v2)",
            x = "Grid resolution (km)", y = "Count") +
      theme_geb()
    save_pub(p, "figS10_grid_event_redefinition",
             path = path_supp_fig(), width = 14, height = 9)
  } else {
    placeholder("figS10_grid_event_redefinition",
                "S10 placeholder — run script 33 to populate.")
  }
}

# ---- S11 CMIP6 per-model maps (handled by 29) ----------------------------
if (!file.exists(path_supp_fig("figS11_cmip6_per_model.pdf"))) {
  placeholder("figS11_cmip6_per_model",
              "S11 placeholder — run script 29 once CMIP6 rasters available.")
}

# ---- S12 Covariate-shift PSI ---------------------------------------------
{
  src <- path_forecasts("table_feature_psi.csv")
  if (file.exists(src)) {
    psi <- fread(src)
    psi_long <- melt(psi, id.vars = "feature",
                      variable.name = "comparison", value.name = "psi")
    p <- ggplot(psi_long, aes(x = reorder(feature, psi), y = psi,
                                fill = comparison)) +
      geom_col(position = "dodge", width = 0.7) +
      geom_hline(yintercept = c(0.1, 0.25), linetype = "dashed",
                  colour = c("grey60", "grey30")) +
      coord_flip() +
      scale_fill_manual(values = pal_cat[1:2]) +
      labs(title = "S12 — Population Stability Index (PSI) per feature",
            x = NULL, y = "PSI") +
      theme_geb()
    save_pub(p, "figS12_covariate_shift_psi",
             path = path_supp_fig(), width = 14, height = 11)
  } else {
    placeholder("figS12_covariate_shift_psi",
                "S12 placeholder — run script 30.")
  }
}

# ---- S13 Effort panel temporal coverage by province -----------------------
{
  src <- path_raw("effort_panel_upgraded.csv")
  if (file.exists(src)) {
    eff <- fread(src, encoding = "UTF-8")
    coverage <- eff[, .(n_years_with_data = sum(!is.na(n_visits))),
                     by = province]
    p <- ggplot(coverage, aes(x = reorder(province, n_years_with_data),
                                y = n_years_with_data)) +
      geom_col(width = 0.7, fill = pal_cat[3]) +
      coord_flip() +
      labs(title = "S13 — Years with effort data by province",
            x = NULL, y = "Years (2002–2024)") +
      theme_geb()
    save_pub(p, "figS13_effort_coverage_by_province",
             path = path_supp_fig(), width = 14, height = 14)
  } else {
    placeholder("figS13_effort_coverage_by_province",
                "S13 placeholder — effort panel missing.")
  }
}

# ---- S14 Species accumulation curve --------------------------------------
{
  src <- path_raw("hazard_risk_upgraded_complete_case.csv")
  if (file.exists(src)) {
    dt <- fread(src, encoding = "UTF-8")
    dt[, year := as.integer(year)]
    acc <- dt[event == 1, .(first = min(year, na.rm = TRUE)), by = species]
    setorder(acc, first)
    acc[, cum_species := seq_len(.N)]
    p <- ggplot(acc, aes(x = first, y = cum_species)) +
      geom_step(linewidth = 0.5) +
      labs(title = "S14 — Cumulative new-record species",
            x = "Year", y = "Cumulative species") +
      theme_geb()
    save_pub(p, "figS14_species_accumulation",
             path = path_supp_fig(), width = 14, height = 9)
  } else {
    placeholder("figS14_species_accumulation",
                "S14 placeholder — risk set missing.")
  }
}

# ---- S15 Phylogenetic effort heatmap (genus) ------------------------------
{
  src <- path_raw("hazard_risk_upgraded_complete_case.csv")
  if (file.exists(src)) {
    dt <- fread(src, encoding = "UTF-8")
    if (!"genus" %in% names(dt) && "species" %in% names(dt)) {
      dt[, genus := tstrsplit(species, " ", fixed = TRUE)[[1]]]
    }
    heat <- dt[, .(events = sum(event, na.rm = TRUE)),
                by = .(genus, province)]
    top_genera <- heat[, .(total = sum(events)), by = genus][order(-total)][1:30, genus]
    heat <- heat[genus %in% top_genera]
    p <- ggplot(heat, aes(x = province, y = genus, fill = events)) +
      geom_tile() + scale_fill_viridis_c(option = "C") +
      labs(title = "S15 — Top-30 genera × province event counts") +
      theme_geb() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6),
            axis.text.y = element_text(size = 6))
    save_pub(p, "figS15_phylogenetic_effort_heatmap",
             path = path_supp_fig(), width = 18, height = 16)
  } else {
    placeholder("figS15_phylogenetic_effort_heatmap",
                "S15 placeholder — risk set missing.")
  }
}

dump_session_info(path_logs("35_publication_figures_supplementary_sessionInfo.txt"))
message("[35] supplementary figures generated under figures/supplementary/.")
