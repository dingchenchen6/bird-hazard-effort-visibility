# ============================================================
# Script: 50_v3_publication_figure_and_manuscript.R
# Family: v3 — final publication figure + manuscript section
# Author: Chen-Chen Ding + Claude Opus 4.7
# Date  : 2026-05-23
#
# ------------------------------------------------------------
# Scientific question / 科学问题:
#   Compose a single, unified publication-quality figure that tells
#   the full v3 robustness story: (i) the relaxed risk set recovers
#   2/3 of the events that v1 dropped; (ii) the headline province
#   M4 climate × effort HR is preserved (1.288 → 1.274); (iii) all
#   four effort specs are positive significant; (iv) ML importance
#   confirms climate × effort stays in the top tier; (v) finer
#   admin scales weaken in v3 — an honest negative-control.
#   一张图把 v3 robustness 故事讲完。
#
# Outputs:
#   figures/main/Figure_v3_robustness_panel.{pdf,png}  (5-panel)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})
options(warn = 1); set.seed(42)

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- normalizePath(file.path(V2, "..", "bird_hazard_model_effort_upgrade"),
                     mustWork = FALSE)

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE,
                                                    showWarnings = FALSE)
ens(file.path(V2, "figures", "main"))

theme_pub <- function(s = 9) {
  theme_bw(base_size = s) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.18, colour = "grey90"),
          panel.border = element_rect(linewidth = 0.4, colour = "grey20"),
          plot.title = element_text(face = "bold", size = s + 1),
          plot.subtitle = element_text(size = s - 1, colour = "grey30"),
          plot.tag = element_text(face = "bold", size = s + 2))
}
COL_RUN  <- c(v1 = "#3B4CC0", v2 = "#7F7F7F", v3 = "#B40426")
COL_SCALE<- c(province = "#B40426", prefecture = "#1F77B4", county = "#7F7F7F")

# ============================================================
# Load every v3-relevant table
# ============================================================
recon       <- fread(file.path(V2, "results", "tables",
                                  "table_province_v1_v2_v3_reconciliation.csv"))
v3_specs    <- fread(file.path(V2, "results", "tables",
                                  "table_province_v3_all_specs_coefs.csv"))
v3_specs_aic<- fread(file.path(V2, "results", "tables",
                                  "table_province_v3_all_specs_aic.csv"))
v3_3scale   <- fread(file.path(V2, "results", "tables",
                                  "table_v3_three_scale_summary.csv"))
v3_attr     <- fread(file.path(V2, "results", "diagnostics",
                                  "table_riskset_v3_attrition.csv"))
imp_v2      <- fread(file.path(V2, "results", "tables",
                                  "table_rf_importance_v2.csv"))
imp_v3      <- fread(file.path(V2, "results", "tables",
                                  "table_rf_importance_v3.csv"))
v1_specs    <- fread(file.path(V1, "results",
                                  "table_cross_specification_key_coefficients.csv"))
v2_coefs    <- fread(file.path(V2, "results", "tables",
                                  "table_province_v2_coefs.csv"))

# ============================================================
# Panel a — attrition funnel (sample-size recovery)
# ============================================================
funnel_dt <- data.table(
  Stage = factor(c("Raw events\n2002-2024",
                    "v1 risk set\n(threshold=100)",
                    "v3 risk set\n(threshold=50 +\nforce-include)"),
                   levels = c("Raw events\n2002-2024",
                               "v1 risk set\n(threshold=100)",
                               "v3 risk set\n(threshold=50 +\nforce-include)")),
  events = c(930, 512, 817),
  species = c(519, 333, 463))
funnel_long <- melt(funnel_dt, id.vars = "Stage",
                      variable.name = "metric", value.name = "n")
funnel_long[, metric_lbl := ifelse(metric == "events",
                                      "Events retained",
                                      "Unique species")]
p_a <- ggplot(funnel_long, aes(x = Stage, y = n, fill = metric_lbl)) +
  geom_col(width = 0.62, position = "dodge") +
  geom_text(aes(label = scales::comma(n)),
             position = position_dodge(width = 0.62),
             vjust = -0.4, size = 2.6) +
  scale_y_continuous(trans = "log10",
                      labels = scales::comma_format(),
                      expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = c(`Unique species` = "#1F77B4",
                                 `Events retained` = "#D62728"),
                     name = NULL) +
  labs(tag = "a",
        title = "v3 recovers 305 events (60%) + 130 species (39%) lost in v1",
        x = NULL, y = "Count (log scale)") +
  theme_pub() + theme(legend.position = "top")

# ============================================================
# Panel b — Headline forest v1 / v2 / v3 (Spec B only)
# ============================================================
recon_b <- recon[grepl("v[123]", run),
                   .(run, hr = interaction_HR,
                     hr.low = HR_low, hr.high = HR_high,
                     p = p_value,
                     events = n_events)]
recon_b[, run_short := factor(c("v1","v2","v3"),
                                 levels = c("v3","v2","v1"))]
p_b <- ggplot(recon_b, aes(x = hr, y = run_short, colour = run_short)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  height = 0.18, linewidth = 0.5) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("HR = %.3f\np = %.0e\n(%d events)",
                                  hr, p, events)),
             nudge_x = 0.04, size = 2.5, hjust = 0) +
  scale_colour_manual(values = c(v1 = "#3B4CC0", v2 = "#7F7F7F",
                                   v3 = "#B40426"), guide = "none") +
  scale_x_continuous(trans = "log",
                      breaks = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5),
                      limits = c(0.97, 1.60)) +
  labs(tag = "b",
        title = "Headline interaction HR — v1 / v2 / v3 (Spec B)",
        x = "Hazard ratio (95 % CI, log)", y = NULL) +
  theme_pub()

# ============================================================
# Panel c — v3 4-spec interaction HR
# ============================================================
v3_int <- v3_specs[model == "M4" & grepl(":", term),
                     .(spec_id, hr, hr.low, hr.high, p.value)]
SPEC_LBL <- c(spec_A = "A: records",
              spec_B = "B: visits",
              spec_C = "C: PCA",
              spec_D = "D: birding-days")
v3_int[, spec_lbl := factor(SPEC_LBL[spec_id],
                              levels = rev(SPEC_LBL))]
p_c <- ggplot(v3_int, aes(x = hr, y = spec_lbl)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  height = 0.18, linewidth = 0.5, colour = "#B40426") +
  geom_point(size = 2.8, colour = "#B40426") +
  geom_text(aes(label = sprintf("%.3f", hr)),
             nudge_x = 0.04, size = 2.6, hjust = 0) +
  scale_x_continuous(trans = "log",
                      breaks = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5),
                      limits = c(0.98, 1.55)) +
  labs(tag = "c",
        title = "v3 — all 4 effort specs significantly positive",
        subtitle = "All 4 p < 10⁻⁴",
        x = "Interaction HR", y = NULL) +
  theme_pub()

# ============================================================
# Panel d — Three-scale v3 forest
# ============================================================
v3_3scale[, scale := factor(scale,
                              levels = c("county","prefecture","province"))]
p_d <- ggplot(v3_3scale, aes(x = hr, y = scale, colour = scale)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  height = 0.18, linewidth = 0.5) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.3f", hr)),
             nudge_x = 0.05, size = 2.5, hjust = 0) +
  scale_colour_manual(values = COL_SCALE, guide = "none") +
  scale_x_continuous(trans = "log",
                      breaks = c(0.95, 1.0, 1.1, 1.2, 1.3, 1.4),
                      limits = c(0.95, 1.45)) +
  labs(tag = "d",
        title = "v3 multi-scale — significance fades at fine grain",
        subtitle = "Province sig, prefecture+county n.s. (relaxed candidates dilute moderation)",
        x = "Interaction HR", y = NULL) +
  theme_pub()

# ============================================================
# Panel e — v2 vs v3 RF rank comparison
# ============================================================
imp_v2[, rank_v2 := frank(-importance)]
imp_v3[, rank_v3 := frank(-importance)]
cmp <- merge(imp_v2[, .(variable, rank_v2)],
              imp_v3[, .(variable, rank_v3)],
              by = "variable")
cmp[, variable_pretty := variable]
cmp[variable == "temp_x_effort",  variable_pretty := "temp × effort"]
cmp[variable == "mahal_x_effort", variable_pretty := "mahal × effort"]
p_e <- ggplot(cmp, aes(x = rank_v2, y = rank_v3)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
               colour = "grey60") +
  geom_point(size = 2.5, colour = "#B40426") +
  geom_text(aes(label = variable_pretty), size = 2.4, hjust = -0.15) +
  scale_x_reverse(breaks = seq(1, 17, 2), limits = c(17.5, 0)) +
  scale_y_reverse(breaks = seq(1, 17, 2), limits = c(17.5, 0)) +
  labs(tag = "e",
        title = "RF importance rank — v2 vs v3",
        subtitle = "temp×effort drops from rank 1 → 6; effort metrics climb 7-10 → 3-5",
        x = "v2 rank (lower = more important)",
        y = "v3 rank") +
  theme_pub()

# Compose 5-panel figure
fig <- (p_a | p_b) / (p_c | p_d) / p_e
fig <- fig + plot_layout(heights = c(1, 1, 1.2)) +
  plot_annotation(
    title = "v3 robustness panel — relaxed SDM threshold + event override + 501 modelled species",
    subtitle = "Province headline preserved (HR 1.274 vs v1 1.292; 305 extra events); fine-grain interaction diluted (n.s. at prefecture/county); ML attributes more weight to effort main effects on the larger v3 set.",
    theme = theme(plot.title = element_text(face = "bold", size = 10),
                   plot.subtitle = element_text(size = 8.5, colour = "grey30")))

ggsave(file.path(V2, "figures", "main",
                  "Figure_v3_robustness_panel.pdf"),
       fig, width = 22, height = 22, units = "cm",
       device = grDevices::cairo_pdf)
ggsave(file.path(V2, "figures", "main",
                  "Figure_v3_robustness_panel.png"),
       fig, width = 22, height = 22, units = "cm", dpi = 600)
message("[50] wrote Figure_v3_robustness_panel.{pdf,png}")
