# ============================================================
# Scientific question / 科学问题:
#   Render the 6 main figures of the v2 manuscript at GEB-compliant
#   600 dpi PDF (cairo) + PNG, all with the same theme, fonts and
#   colour-blind safe palette.
#   按 GEB 规范统一生成 6 张主图（600 dpi PDF + PNG）。
#
# Objective / 分析目标:
#   Fig 1  Conceptual / study domain
#   Fig 2  Effort × climate interaction forest (4 specs)
#   Fig 3  Multi-scale variance decomposition (province / 50 / 100 km)
#   Fig 4  Spatial-block CV + Moran's I (composite)
#   Fig 5  CMIP6 2050 SSP585 ensemble (if not already produced by 29)
#   Fig 6  Forecast skill decay + scenario narrative
#
# Input data / 输入数据:
#   results/tables/*.csv, results/diagnostics/*.csv,
#   results/forecasts/*.parquet, figures/diagnostics/*.png
#
# Main workflow / 主要流程: assemble panels, save to figures/main/.
# Expected output / 预期输出: figures/main/figN.{pdf,png}
# Main packages / 主要包: ggplot2, patchwork, ggspatial, scales,
#   data.table, arrow.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(ggplot2)
  library(patchwork)
  library(glue)
  library(sf)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_spatial.R"))
source(file.path("code", "utils", "utils_plots.R"))

ensure_dir(path_main_fig())

# ----------------------------------------------------------------------------
# Fig 1 — Conceptual & study domain
# 概念图 + 研究域。
# ----------------------------------------------------------------------------
make_fig1 <- function() {
  basemap <- china_basemap_gs2019_1822()
  panel_a <- ggplot() +
    geom_sf(data = basemap$province, fill = "grey95", colour = "grey60",
            linewidth = 0.1) +
    geom_sf(data = basemap$national, fill = NA, colour = "grey20",
            linewidth = 0.4) +
    geom_sf(data = basemap$ninedash, fill = NA, colour = "grey20",
            linewidth = 0.3) +
    coord_sf(datum = NA) +
    labs(title = "(a) Study domain: 32 mainland provinces") +
    theme_geb() +
    theme(panel.grid = element_blank())

  # Panel b — conceptual DAG.
  dag <- data.frame(
    x = c(1, 1, 2, 3, 3),
    y = c(3, 1, 2, 3, 1),
    label = c("Climate velocity\n气候速率",
              "Survey effort\n调查努力",
              "Climate × Effort\n交互",
              "Detection prob.\n探测概率",
              "First arrival\n首次到达"),
    type = c("input", "input", "interaction", "mediator", "outcome")
  )
  edges <- data.frame(
    x = c(1, 1, 1, 1, 2, 3),
    y = c(3, 1, 3, 1, 2, 3),
    xend = c(2, 2, 3, 3, 3, 3),
    yend = c(2, 2, 3, 1, 1, 1)
  )
  panel_b <- ggplot() +
    geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
                  arrow = arrow(length = unit(0.15, "cm")), colour = "grey50") +
    geom_label(data = dag, aes(x = x, y = y, label = label, fill = type),
                size = 2.4, label.size = 0.2) +
    scale_fill_manual(values = c(input = "#B8D8F2", interaction = "#F2C8B8",
                                  mediator = "#FFF2CC", outcome = "#C6E0B4"),
                       guide = "none") +
    xlim(0.5, 3.5) + ylim(0.5, 3.5) +
    labs(title = "(b) Causal DAG") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 9))

  # Panel c — sample-size flow.
  panel_c_dt <- data.frame(
    stage = factor(c("Records", "Provinces", "Species", "Years"),
                    levels = c("Records", "Provinces", "Species", "Years")),
    value = c(12813, 32, 333, 23)
  )
  panel_c <- ggplot(panel_c_dt, aes(x = stage, y = value, fill = stage)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = scales::comma(value)), vjust = -0.3, size = 2.6) +
    scale_y_continuous(trans = "log10",
                        labels = scales::comma_format()) +
    scale_fill_manual(values = pal_cat[1:4], guide = "none") +
    labs(title = "(c) Sample size", x = NULL, y = "Count (log scale)") +
    theme_geb()

  fig1 <- (panel_a | panel_b) / panel_c + plot_layout(heights = c(2, 1))
  save_pub(fig1, "fig1_conceptual_and_domain", width = 17.5, height = 16)
}

# ----------------------------------------------------------------------------
# Fig 2 — Effort × climate forest across specs
# ----------------------------------------------------------------------------
make_fig2 <- function() {
  src <- path_tables("table_cross_specification_key_coefficients.csv")
  if (!file.exists(src)) {
    src <- path_raw("../results/table_cross_specification_key_coefficients.csv")
  }
  if (!file.exists(src)) {
    message("[Fig 2] missing source table; skip.")
    return(invisible(NULL))
  }
  coef_tab <- fread(src, encoding = "UTF-8")
  # Tolerant column name handling. 列名兼容。
  std <- function(x, candidates, fb) {
    hit <- intersect(candidates, names(x))
    if (length(hit)) return(setnames(x, hit[1], fb)) else x
  }
  coef_tab <- std(coef_tab, c("hr_estimate", "hazard_ratio", "estimate_hr"), "hr")
  coef_tab <- std(coef_tab, c("hr_low", "lower"), "hr.low")
  coef_tab <- std(coef_tab, c("hr_high", "upper"), "hr.high")
  coef_tab <- std(coef_tab, c("spec", "specification"), "spec")
  coef_tab <- std(coef_tab, c("term", "parameter"), "term")
  interaction_dt <- coef_tab[grepl("interaction|:", term)]
  fig2 <- forest_plot(interaction_dt, term_col = "spec",
                       group_col = if ("model" %in% names(interaction_dt)) "model" else NULL,
                       title = "Climate × effort interaction across 4 specifications") +
            scale_colour_manual(values = pal_cat, na.value = "grey50")
  save_pub(fig2, "fig2_effort_climate_interaction_forest",
           width = 12, height = 8)
}

# ----------------------------------------------------------------------------
# Fig 3 — Multi-scale variance decomposition
# ----------------------------------------------------------------------------
make_fig3 <- function() {
  src <- path_tables("table_variance_decomposition_r2.csv")
  if (!file.exists(src)) {
    src <- path_raw("../results/table_variance_decomposition_r2.csv")
  }
  if (!file.exists(src)) {
    message("[Fig 3] missing source table; skip.")
    return(invisible(NULL))
  }
  vd <- fread(src, encoding = "UTF-8")
  setnames(vd, tolower(names(vd)))
  if (!"scale" %in% names(vd)) vd[, scale := "province"]
  long <- melt(vd, id.vars = c("scale"), variable.name = "component",
               value.name = "share")
  long <- long[grepl("interaction|climate|effort|joint|residual", component)]
  fig3 <- ggplot(long, aes(x = scale, y = share, fill = component)) +
    geom_col(position = "stack", width = 0.7) +
    scale_fill_manual(values = pal_cat[1:uniqueN(long$component)]) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Variance decomposition across spatial scales",
          x = "Spatial scale", y = "Share of conditional R²") +
    theme_geb()
  save_pub(fig3, "fig3_variance_decomposition_multi_scale",
           width = 14, height = 9)
}

# ----------------------------------------------------------------------------
# Fig 4 — Spatial block CV + Moran's I composite
# ----------------------------------------------------------------------------
make_fig4 <- function() {
  cv_path <- path_diagnostics("table_spatial_block_cv.csv")
  mi_path <- path_diagnostics("table_morans_i_residuals.csv")
  if (!file.exists(cv_path) || !file.exists(mi_path)) {
    message("[Fig 4] missing CV / Moran source; skip.")
    return(invisible(NULL))
  }
  cv <- fread(cv_path)
  mi <- fread(mi_path)
  cv_long <- melt(cv, id.vars = c("fold", "n_train", "n_test", "pos_test"),
                  measure.vars = c("auc_roc", "auc_pr"),
                  variable.name = "metric", value.name = "value")
  p_cv <- ggplot(cv_long, aes(x = factor(fold), y = value, fill = metric)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey50") +
    scale_fill_manual(values = pal_cat[1:2]) +
    labs(title = "(a) Spatial-block CV performance",
          x = "Fold", y = "Score") +
    theme_geb()

  p_mi <- ggplot(mi, aes(x = class_km, y = I, colour = model)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
    geom_errorbar(aes(ymin = I - 1.96 * sd, ymax = I + 1.96 * sd),
                   width = 5, linewidth = 0.3) +
    scale_colour_manual(values = pal_cat[seq_len(uniqueN(mi$model))]) +
    labs(title = "(b) Residual Moran's I",
          x = "Distance class (km)", y = "Moran's I") +
    theme_geb()

  fig4 <- p_cv | p_mi
  save_pub(fig4, "fig4_spatial_cv_and_morans",
           width = 17.5, height = 9)
}

# ----------------------------------------------------------------------------
# Fig 5 — defer to 29_cmip6_ensemble_prediction.R unless missing
# ----------------------------------------------------------------------------
make_fig5 <- function() {
  pdf <- path_main_fig("fig5_cmip6_ensemble_2050_ssp585.pdf")
  if (file.exists(pdf)) {
    message("[Fig 5] already produced by script 29.")
    return(invisible(pdf))
  }
  # Fall-back placeholder: just embed forecast skill summary.
  message("[Fig 5] placeholder (run script 29 once CMIP6 rasters available).")
  p <- ggplot() + theme_void() +
    annotate("text", x = 1, y = 1,
             label = "Fig 5: CMIP6 ensemble — run script 29",
             size = 4)
  save_pub(p, "fig5_cmip6_ensemble_2050_ssp585_placeholder",
           width = 17.5, height = 9)
}

# ----------------------------------------------------------------------------
# Fig 6 — Forecast skill decay + scenario narrative
# ----------------------------------------------------------------------------
make_fig6 <- function() {
  skill_path <- path_forecasts("table_forecast_skill_decay.csv")
  psi_path   <- path_forecasts("table_feature_psi.csv")
  if (!file.exists(skill_path)) {
    message("[Fig 6] missing skill table; skip.")
    return(invisible(NULL))
  }
  skill <- fread(skill_path)
  p_skill <- ggplot(skill, aes(x = horizon, y = auc_roc)) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey50") +
    labs(title = "(a) Forecast skill decay",
          x = "Horizon (years)", y = "AUC") +
    theme_geb()
  panels <- list(p_skill)
  if (file.exists(psi_path)) {
    psi <- fread(psi_path)
    psi_long <- melt(psi, id.vars = "feature",
                      variable.name = "comparison",
                      value.name = "psi")
    p_psi <- ggplot(psi_long, aes(x = reorder(feature, psi), y = psi,
                                    fill = comparison)) +
      geom_col(position = "dodge", width = 0.7) +
      geom_hline(yintercept = c(0.1, 0.25), linetype = "dashed",
                  colour = c("grey60", "grey30")) +
      coord_flip() +
      scale_fill_manual(values = pal_cat[1:2]) +
      labs(title = "(b) Covariate shift (PSI)",
            x = NULL, y = "PSI") +
      theme_geb()
    panels <- list(p_skill, p_psi)
  }
  fig6 <- patchwork::wrap_plots(panels, ncol = length(panels))
  save_pub(fig6, "fig6_forecast_skill_and_psi",
           width = 17.5, height = 9)
}

invisible(list(
  fig1 = make_fig1(),
  fig2 = make_fig2(),
  fig3 = make_fig3(),
  fig4 = make_fig4(),
  fig5 = make_fig5(),
  fig6 = make_fig6()
))

dump_session_info(path_logs("34_publication_figures_main_sessionInfo.txt"))
message("[34] main figures generated under figures/main/.")
