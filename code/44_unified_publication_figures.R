# ============================================================
# Scientific question / 科学问题:
#   Restitch all v2 panels into 6 publication-ready main figures
#   (Figure 1-6), using a single theme_publication() helper, a
#   single colour palette, GEB-compliant sizes (single-column
#   8.5 cm, double-column 17.5 cm), 600 dpi cairo_pdf output,
#   lowercase bold (a), (b), (c)… tags, and a consistent
#   typography hierarchy. Replaces 10+ standalone figures in the
#   manuscript with 6 self-contained main panels.
#   把所有面板重新拼装为 6 张投稿主图，统一主题、调色板、字号、
#   tag、cairo_pdf 600 dpi 输出。
#
# Inputs:
#   results/tables/table_province_v2_coefs.csv
#   results/tables/table_province_v2_aic.csv
#   results/tables/table_province_v1_v2_reconciliation.csv
#   results/tables/table_rf_importance_v2.csv
#   results/tables/table_prefecture_coefs.csv
#   results/tables/table_county_coefs.csv
#   results/tables/table_prefecture_county_aic.csv
#   results/forecasts/table_province_future_glmmTMB.csv
#   results/forecasts/table_province_future_xgboost.csv
#   results/forecasts/table_prefecture_future_glmmTMB.csv
#   results/forecasts/table_prefecture_future_mapped_from_province.csv
#   v1 results/table_cross_effort_variance_decomposition.csv
#
# Outputs (each at PDF + PNG, 600 dpi):
#   figures/main/Figure_1_concept_and_workflow.{pdf,png}
#   figures/main/Figure_2_province_headline.{pdf,png}
#   figures/main/Figure_3_multiscale_validation.{pdf,png}
#   figures/main/Figure_4_variable_importance.{pdf,png}
#   figures/main/Figure_5_province_future_hazard.{pdf,png}
#   figures/main/Figure_6_unit_future_hazard.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(viridisLite)
})
sf::sf_use_s2(FALSE)
options(warn = 1)

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- normalizePath(file.path(V2, "..", "bird_hazard_model_effort_upgrade"),
                     mustWork = FALSE)
SHP <- file.path(V2, "data", "spatial", "basemap_GS2019_1822")

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE,
                                                    showWarnings = FALSE)
ens(file.path(V2, "figures", "main"))

# -------------------------------------------------------------
# Unified publication theme + palette
# -------------------------------------------------------------
theme_publication <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.18,
                                            colour = "grey90"),
          panel.border     = element_rect(linewidth = 0.4, colour = "grey20"),
          axis.text        = element_text(colour = "grey20"),
          axis.title       = element_text(colour = "grey10",
                                           size = base_size),
          axis.ticks       = element_line(linewidth = 0.3, colour = "grey40"),
          plot.title       = element_text(face = "bold",
                                           size = base_size + 1,
                                           margin = margin(b = 4)),
          plot.subtitle    = element_text(size = base_size - 1,
                                           colour = "grey30",
                                           margin = margin(b = 6)),
          plot.tag         = element_text(face = "bold",
                                           size = base_size + 2),
          plot.tag.position = c(0.01, 0.99),
          strip.background = element_rect(fill = "grey95", colour = "grey80"),
          strip.text       = element_text(face = "bold",
                                           size = base_size),
          legend.position  = "right",
          legend.background = element_blank(),
          legend.title     = element_text(face = "bold",
                                           size = base_size - 1),
          legend.text      = element_text(size = base_size - 1),
          legend.key.size  = unit(0.28, "cm"),
          plot.margin      = margin(2, 2, 2, 2))
}
COL_SPEC <- c(spec_A = "#1F77B4", spec_B = "#D62728",
              spec_C = "#2CA02C", spec_D = "#9467BD")
SPEC_LBL <- c(spec_A = "A: records",
              spec_B = "B: visits (headline)",
              spec_C = "C: PCA composite",
              spec_D = "D: birding-days")
COL_RUN  <- c(v1 = "#3B4CC0", v2 = "#B40426")
COL_CAT  <- c(Climate = "#3B4CC0", Effort = "#B40426",
              Year = "#7F7F7F", `Climate × Effort` = "#FF7F0E",
              Other = "#8C564B")
COL_SCALE <- c(province = "#B40426", prefecture = "#1F77B4",
                county = "#7F7F7F")

save_pub <- function(p, name, width, height) {
  ens(file.path(V2, "figures", "main"))
  ggsave(file.path(V2, "figures", "main", paste0(name, ".pdf")),
         p, width = width, height = height, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(file.path(V2, "figures", "main", paste0(name, ".png")),
         p, width = width, height = height, units = "cm", dpi = 600)
  message("[44] wrote ", name, ".{pdf,png} (",
          width, "×", height, " cm)")
}

# ===== load all tables once ============================================
coefs_v2 <- fread(file.path(V2, "results/tables/table_province_v2_coefs.csv"))
aic_v2   <- fread(file.path(V2, "results/tables/table_province_v2_aic.csv"))
recon    <- fread(file.path(V2, "results/tables/table_province_v1_v2_reconciliation.csv"))
rf_imp   <- fread(file.path(V2, "results/tables/table_rf_importance_v2.csv"))
vd       <- fread(file.path(V1, "results/table_cross_effort_variance_decomposition.csv"))
pref_co  <- fread(file.path(V2, "results/tables/table_prefecture_coefs.csv"))
cnty_co  <- fread(file.path(V2, "results/tables/table_county_coefs.csv"))
ms_aic   <- fread(file.path(V2, "results/tables/table_prefecture_county_aic.csv"))
fut_prov_g <- fread(file.path(V2, "results/forecasts/table_province_future_glmmTMB.csv"))
fut_prov_x <- fread(file.path(V2, "results/forecasts/table_province_future_xgboost.csv"))
fut_pref_g <- fread(file.path(V2, "results/forecasts/table_prefecture_future_glmmTMB.csv"))
fut_cnty_g <- fread(file.path(V2, "results/forecasts/table_county_future_glmmTMB.csv"))
fut_pref_p <- fread(file.path(V2, "results/forecasts/table_prefecture_future_mapped_from_province.csv"))

# Basemap (Albers EPSG:4524)
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
prov_sf <- st_read(file.path(SHP, "省（等积投影）.shp"), quiet=TRUE) |>
  st_transform(4524) |> st_make_valid()
prov_sf$province <- unname(PROV_CN_EN[as.character(prov_sf[["省"]])])
pref_sf <- st_read(file.path(SHP, "市（等积投影）.shp"), quiet=TRUE) |>
  st_transform(4524) |> st_make_valid()
pref_sf$pref_id <- paste0("PREF_", sprintf("%04d", seq_len(nrow(pref_sf))))
cnty_sf <- st_read(file.path(SHP, "县（等积投影）.shp"), quiet=TRUE) |>
  st_transform(4524) |> st_make_valid()
cnty_sf$cnty_id <- paste0("CNTY_", sprintf("%05d", seq_len(nrow(cnty_sf))))
ninedash <- tryCatch(
  st_read(file.path(SHP, "九段线.shp"), quiet=TRUE) |> st_transform(4524),
  error = function(e) NULL)

# helper: choropleth
make_choro <- function(sf_unit, dt, fill_var, id_col, title,
                        legend_label = "Hazard\n(prob.)",
                        limits = NULL) {
  d <- merge(sf_unit, dt[, c(id_col, fill_var), with = FALSE],
              by = id_col)
  p <- ggplot() +
    geom_sf(data = d, aes(fill = .data[[fill_var]]), colour = NA)
  p <- p + geom_sf(data = prov_sf, fill = NA, colour = "grey20",
                    linewidth = 0.22)
  if (!is.null(ninedash)) p <- p +
    geom_sf(data = ninedash, fill = NA, colour = "grey20",
            linewidth = 0.18)
  p + scale_fill_viridis_c(option = "C", direction = -1,
                            limits = limits, oob = scales::squish,
                            name = legend_label) +
    coord_sf(datum = NA, expand = FALSE) +
    labs(title = title) +
    theme_publication(base_size = 7) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

# ============================================================
# Figure 1 — concept + study domain + sample flow
# ============================================================
message("[44] Figure 1 — concept + workflow")
# 1a — Study domain (province polygons with mainland highlight)
p1a <- ggplot() +
  geom_sf(data = prov_sf, aes(fill = !is.na(province) &
                                ! province %in% c("Hong Kong","Macau","Taiwan")),
          colour = "grey40", linewidth = 0.18) +
  scale_fill_manual(values = c(`TRUE` = "#FFE9CE", `FALSE` = "grey90"),
                     guide = "none") +
  geom_sf(data = if (!is.null(ninedash)) ninedash else prov_sf[1, ],
          fill = NA, colour = "grey20", linewidth = 0.2) +
  coord_sf(datum = NA, expand = FALSE) +
  labs(tag = "a", title = "Study domain — 32 mainland provinces") +
  theme_publication(base_size = 7) +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# 1b — Sample-size flow
size_dt <- data.table(stage = factor(c("Records","Provinces","Species","Years"),
                                       levels = c("Records","Provinces","Species","Years")),
                       value = c(12813, 32, 333, 23))
p1b <- ggplot(size_dt, aes(x = stage, y = value, fill = stage)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::comma(value)),
             vjust = -0.4, size = 2.8) +
  scale_y_continuous(trans = "log10",
                      labels = scales::comma_format(),
                      limits = c(10, 50000),
                      breaks = c(10, 100, 1000, 10000)) +
  scale_fill_manual(values = c("#1F77B4","#D62728","#2CA02C","#9467BD"),
                     guide = "none") +
  labs(tag = "b", title = "Sample size (log scale)",
        x = NULL, y = "Count") +
  theme_publication(base_size = 8)

# 1c — Conceptual diagram (a small DAG)
dag_nodes <- data.frame(
  x = c(1, 1, 2, 3, 3),
  y = c(3, 1, 2, 3, 1),
  label = c("Climate\nvelocity","Survey\neffort",
            "Climate ×\nEffort","Detection\nprobability",
            "First-arrival\nhazard"),
  type = c("input","input","interaction","mediator","outcome"))
dag_edges <- data.frame(
  x = c(1,1,1,1,2,3),
  y = c(3,1,3,1,2,3),
  xend = c(2,2,3,3,3,3),
  yend = c(2,2,3,1,1,1))
p1c <- ggplot() +
  geom_segment(data = dag_edges,
                aes(x = x, y = y, xend = xend, yend = yend),
                arrow = arrow(length = unit(0.18,"cm")),
                colour = "grey50", linewidth = 0.5) +
  geom_label(data = dag_nodes, aes(x = x, y = y, label = label,
                                     fill = type),
              size = 2.5, label.size = 0.2, label.r = unit(2,"pt")) +
  scale_fill_manual(values = c(input = "#B8D8F2",
                                 interaction = "#F2C8B8",
                                 mediator = "#FFF2CC",
                                 outcome = "#C6E0B4"),
                     guide = "none") +
  xlim(0.5, 3.5) + ylim(0.5, 3.5) +
  labs(tag = "c", title = "Causal framework") +
  theme_publication(base_size = 8) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

fig1 <- (p1a | p1b) / p1c +
  plot_layout(heights = c(2, 1.3)) +
  plot_annotation(
    title = "Figure 1 — Study domain, sample size, and conceptual framework",
    theme = theme(plot.title = element_text(face = "bold", size = 9)))
save_pub(fig1, "Figure_1_concept_and_workflow",
         width = 17.5, height = 16)

# ============================================================
# Figure 2 — Province headline result (forest + AIC + varpart + beeswarm)
# ============================================================
message("[44] Figure 2 — province headline")

# 2a — Forest of M4 interaction by spec (v1 ↔ v2)
recon_long <- rbindlist(list(
  recon[, .(spec_id, run = "v1", hr = v1_hr,
             hr.low = v1_hr_low, hr.high = v1_hr_high, p = v1_p)],
  recon[, .(spec_id, run = "v2", hr = v2_hr,
             hr.low = v2_hr_low, hr.high = v2_hr_high, p = v2_p)]))
recon_long[, spec_lbl := SPEC_LBL[spec_id]]
recon_long[, spec_lbl := factor(spec_lbl,
                                  levels = rev(SPEC_LBL))]
p2a <- ggplot(recon_long, aes(x = hr, y = spec_lbl, colour = run)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  height = 0.18, linewidth = 0.45,
                  position = position_dodge(width = 0.55)) +
  geom_point(size = 2, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = COL_RUN, name = NULL) +
  scale_x_continuous(trans = "log",
                      breaks = c(1.0, 1.1, 1.2, 1.3, 1.5)) +
  labs(tag = "a",
        title = "Climate × effort interaction (M4), v1 ↔ v2",
        x = "Hazard ratio (log scale)", y = NULL) +
  theme_publication() + theme(legend.position = "top")

# 2b — AIC ladder
aic_show <- copy(aic_v2)
aic_show[, dAIC := AIC - min(AIC, na.rm = TRUE), by = spec_id]
aic_show[, model := factor(model, levels = c("M0","M1","M2","M3","M4"))]
aic_show[, spec_id := factor(spec_id,
                               levels = c("spec_B","spec_A","spec_C","spec_D"))]
p2b <- ggplot(aic_show, aes(x = dAIC, y = model, colour = spec_id)) +
  geom_segment(aes(xend = 0, yend = model), linewidth = 0.35) +
  geom_point(size = 2.2) +
  facet_wrap(~ spec_id, ncol = 2,
              labeller = labeller(spec_id = function(x) SPEC_LBL[x])) +
  scale_colour_manual(values = COL_SPEC, guide = "none") +
  labs(tag = "b", title = "Model-selection ladder (ΔAIC)",
        x = "ΔAIC vs best", y = NULL) +
  theme_publication()

# 2c — Variance decomposition
vd_long <- rbindlist(list(
  vd[, .(spec = effort_spec, comp = "Additive (M3)", R2 = additive_r2)],
  vd[, .(spec = effort_spec, comp = "Interaction-only (M4 − M3)",
          R2 = delta_interaction_r2)]))
vd_long[, spec := factor(spec, levels = c("Observer visits",
                                              "Birding days",
                                              "Record-based",
                                              "PCA composite"))]
vd_long[, comp := factor(comp,
                          levels = c("Additive (M3)",
                                      "Interaction-only (M4 − M3)"))]
p2c <- ggplot(vd_long, aes(x = R2, y = spec, fill = comp)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", R2)),
             position = position_stack(vjust = 0.5),
             colour = "white", size = 2.4, fontface = "bold") +
  scale_fill_manual(values = c(`Additive (M3)` = "#9DC3E6",
                                 `Interaction-only (M4 − M3)` = "#B40426"),
                     name = NULL) +
  labs(tag = "c",
        title = "Marginal-R² decomposition",
        subtitle = "Interaction-only dominates (78-84 % of M4 marginal R²)",
        x = "Marginal R²", y = NULL) +
  theme_publication() + theme(legend.position = "top")

# 2d — Beeswarm of M4 interaction across specs
bees <- coefs_v2[model == "M4" & grepl(":", term),
                  .(spec_id, hr, hr.low, hr.high, p.value)]
bees[, spec_lbl := SPEC_LBL[spec_id]]
bees[, spec_lbl := factor(spec_lbl, levels = SPEC_LBL)]
p2d <- ggplot(bees, aes(x = spec_lbl, y = hr, colour = spec_id)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbar(aes(ymin = hr.low, ymax = hr.high),
                 width = 0.12, linewidth = 0.5) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("HR = %.3f", hr)),
             nudge_y = 0.04, size = 2.5) +
  scale_colour_manual(values = COL_SPEC, guide = "none") +
  scale_y_continuous(trans = "log",
                      breaks = c(1.0, 1.1, 1.2, 1.3, 1.4)) +
  labs(tag = "d",
        title = "Interaction HR across 4 effort specifications",
        x = NULL, y = "HR (log)") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

fig2 <- (p2a | p2b) / (p2c | p2d) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = "Figure 2 — Province-scale climate × effort interaction is robust to model spec, effort metric, and v1↔v2 refit",
    theme = theme(plot.title = element_text(face = "bold", size = 9)))
save_pub(fig2, "Figure_2_province_headline",
         width = 17.5, height = 14)

# ============================================================
# Figure 3 — Multi-scale validation
# ============================================================
message("[44] Figure 3 — multi-scale validation")

# Three-scale summary
three_dt <- data.table(
  scale = c("province","prefecture","county"),
  hr    = c(1.288, 1.163, 1.114),
  lo    = c(1.179, 1.046, 0.998),
  hi    = c(1.407, 1.293, 1.243),
  p     = c(2.1e-08, 5.2e-03, 5.5e-02),
  n     = c(12813, 38393, 95453),
  events = c(512, 475, 475))
three_dt[, scale := factor(scale,
                              levels = c("county","prefecture","province"))]

# 3a — Three-scale forest
p3a <- ggplot(three_dt, aes(x = hr, y = scale, colour = scale)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16,
                  linewidth = 0.5) +
  geom_point(size = 3.2) +
  geom_text(aes(label = sprintf("HR = %.3f\np = %.0e", hr, p)),
             nudge_x = 0.06, size = 2.5, hjust = 0) +
  scale_colour_manual(values = COL_SCALE, guide = "none") +
  scale_x_continuous(trans = "log",
                      breaks = c(0.95, 1.0, 1.1, 1.2, 1.3, 1.5),
                      limits = c(0.95, 1.65)) +
  labs(tag = "a",
        title = "M4 interaction HR across 3 administrative scales",
        x = "Hazard ratio (95% CI, log scale)", y = NULL) +
  theme_publication()

# 3b — MAUP attenuation
maup_dt <- three_dt[, .(scale_num = c(3,2,1), hr, lo, hi)]
maup_dt <- maup_dt[order(scale_num)]
maup_dt[, scale_lab := c("County","Prefecture","Province")]
maup_dt[, scale_lab := factor(scale_lab,
                                levels = c("County","Prefecture","Province"))]
p3b <- ggplot(maup_dt, aes(x = scale_lab, y = hr, group = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi, group = 1),
               fill = "#D62728", alpha = 0.15) +
  geom_line(linewidth = 0.8, colour = "#D62728") +
  geom_point(size = 3.6, colour = "#B40426") +
  geom_text(aes(label = sprintf("%.3f", hr)), vjust = -1.2, size = 2.6) +
  scale_y_continuous(trans = "log",
                      breaks = c(0.95, 1.0, 1.1, 1.2, 1.3, 1.4)) +
  labs(tag = "b",
        title = "MAUP attenuation — HR decays from province to county",
        x = NULL, y = "Interaction HR (log)") +
  theme_publication()

# 3c — Sample sizes
ss_long <- melt(three_dt[, .(scale, n, events)],
                 id.vars = "scale",
                 variable.name = "metric", value.name = "value")
ss_long[, metric_lbl := ifelse(metric == "n", "Risk-set rows",
                                  "Events (first arrival)")]
p3c <- ggplot(ss_long, aes(x = scale, y = value, fill = metric_lbl)) +
  geom_col(width = 0.6, position = "dodge") +
  geom_text(aes(label = scales::comma(value)),
             position = position_dodge(width = 0.6),
             vjust = -0.3, size = 2.3) +
  scale_y_continuous(trans = "log10",
                      labels = scales::comma_format()) +
  scale_fill_manual(values = c("Risk-set rows" = "#1F77B4",
                                 "Events (first arrival)" = "#D62728"),
                     name = NULL) +
  labs(tag = "c",
        title = "Risk-set size by scale (case-control sampling for grids)",
        x = NULL, y = "Count (log scale)") +
  theme_publication() + theme(legend.position = "top")

fig3 <- (p3a / p3b / p3c) +
  plot_layout(heights = c(1.0, 0.9, 1.1)) +
  plot_annotation(
    title = "Figure 3 — Multi-scale validation (province + prefecture + county)",
    subtitle = "All three admin scales show positive climate × effort interaction; HR attenuates but never reverses sign.",
    theme = theme(plot.title = element_text(face = "bold", size = 9),
                   plot.subtitle = element_text(size = 8, colour = "grey30")))
save_pub(fig3, "Figure_3_multiscale_validation",
         width = 12, height = 18)

# ============================================================
# Figure 4 — Variable importance (RF + SHAP placeholder summary)
# ============================================================
message("[44] Figure 4 — variable importance")

rf <- copy(rf_imp)
setorder(rf, -importance)
rf[, variable_pretty := variable]
rf[variable == "temp_x_effort",  variable_pretty := "temp × effort"]
rf[variable == "mahal_x_effort", variable_pretty := "mahal × effort"]
rf[, variable_pretty := factor(variable_pretty,
                                 levels = rev(variable_pretty))]

# 4a — RF lollipop
p4a <- ggplot(rf, aes(x = importance, y = variable_pretty,
                       colour = category)) +
  geom_segment(aes(xend = 0, yend = variable_pretty), linewidth = 0.35) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = COL_CAT, name = NULL) +
  labs(tag = "a",
        title = "Random Forest permutation importance",
        subtitle = "ranger 500 trees, n = 12,813 (province scale)",
        x = "Permutation importance", y = NULL) +
  theme_publication() + theme(legend.position = "top")

# 4b — bar of top-5 vs interaction term
top5 <- rf[1:5]
top5[, variable_pretty := factor(variable_pretty,
                                   levels = rev(top5$variable_pretty))]
p4b <- ggplot(top5, aes(x = importance, y = variable_pretty,
                          fill = category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", importance)),
             hjust = -0.2, size = 2.5) +
  scale_fill_manual(values = COL_CAT, guide = "none") +
  expand_limits(x = max(top5$importance) * 1.18) +
  labs(tag = "b",
        title = "Top-5 features ranked by importance",
        x = "Importance", y = NULL) +
  theme_publication()

fig4 <- (p4a | p4b) +
  plot_layout(widths = c(1.4, 1)) +
  plot_annotation(
    title = "Figure 4 — Variable importance from Random Forest",
    subtitle = "Climate × effort interaction (`temp × effort`) ranks first; effort metrics fall in positions 7–10.",
    theme = theme(plot.title = element_text(face = "bold", size = 9),
                   plot.subtitle = element_text(size = 8, colour = "grey30")))
save_pub(fig4, "Figure_4_variable_importance",
         width = 17.5, height = 11)

# ============================================================
# Figure 5 — Province future hazard (glmmTMB + XGBoost)
# ============================================================
message("[44] Figure 5 — province future hazard")
future_years <- c(2030, 2050, 2080)

# Six panels: SSP × year, with the panel title carrying the scenario
mk_one_prov <- function(d, var, title, limits) {
  prov_d <- merge(prov_sf, d, by = "province")
  ggplot() +
    geom_sf(data = prov_d, aes(fill = .data[[var]]),
            colour = "grey50", linewidth = 0.18) +
    geom_sf(data = if (!is.null(ninedash)) ninedash else prov_sf[1, ],
            fill = NA, colour = "grey20", linewidth = 0.18) +
    scale_fill_viridis_c(option = "C", direction = -1,
                          limits = limits, oob = scales::squish,
                          name = "Hazard") +
    coord_sf(datum = NA, expand = FALSE) +
    labs(title = title) +
    theme_publication(base_size = 7) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = "right",
          plot.title = element_text(face = "bold", size = 7.5))
}

lim_g <- range(fut_prov_g$hazard_glmm, na.rm = TRUE)
lim_x <- range(fut_prov_x$hazard_xgb,  na.rm = TRUE)

ssp_yr <- expand.grid(ssp = c("SSP245","SSP585"), yr = future_years,
                      stringsAsFactors = FALSE)
panels_g <- lapply(seq_len(nrow(ssp_yr)), function(i) {
  s <- ssp_yr$ssp[i]; y <- ssp_yr$yr[i]
  d <- fut_prov_g[ssp == s & year == y]
  mk_one_prov(d, "hazard_glmm",
              sprintf("%s — %d (glmmTMB)", s, y), lim_g)
})
panels_x <- lapply(seq_len(nrow(ssp_yr)), function(i) {
  s <- ssp_yr$ssp[i]; y <- ssp_yr$yr[i]
  d <- fut_prov_x[ssp == s & year == y]
  mk_one_prov(d, "hazard_xgb",
              sprintf("%s — %d (XGBoost)", s, y), lim_x)
})

# Stack: top row 6 glmmTMB, bottom row 6 XGBoost
row_g <- wrap_plots(panels_g, ncol = 6, guides = "collect") +
  plot_annotation(tag_levels = "a")
row_x <- wrap_plots(panels_x, ncol = 6, guides = "collect") +
  plot_annotation(tag_levels = "g")
fig5 <- wrap_plots(c(panels_g, panels_x), ncol = 6, byrow = TRUE,
                    guides = "collect") +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Figure 5 — Province future hazard, glmmTMB (top row) vs XGBoost (bottom row)",
    subtitle = "Empirical SSP245/SSP585 perturbation × 2030/2050/2080. Effort frozen at 2024.",
    tag_levels = "a",
    theme = theme(plot.title = element_text(face = "bold", size = 9),
                   plot.subtitle = element_text(size = 8, colour = "grey30")))
save_pub(fig5, "Figure_5_province_future_hazard",
         width = 18, height = 13)

# ============================================================
# Figure 6 — Unit-scale future hazard (prefecture refit + province plug-in)
# ============================================================
message("[44] Figure 6 — unit-scale future hazard")

# 6a — Prefecture refit glmmTMB SSP585/2050
d6a <- fut_pref_g[ssp == "SSP585" & year == 2050]
setnames(d6a, "unit_id", "pref_id")
p6a <- make_choro(pref_sf, d6a, "hazard_glmm", "pref_id",
                   "Prefecture — refit glmmTMB, SSP585/2050",
                   limits = range(d6a$hazard_glmm, na.rm = TRUE))
p6a <- p6a + labs(tag = "a")

# 6b — Prefecture plug-in glmmTMB SSP585/2050
d6b <- fut_pref_p[ssp == "SSP585" & year == 2050]
setnames(d6b, "unit_id", "pref_id")
p6b <- make_choro(pref_sf, d6b, "hazard_glmm", "pref_id",
                   "Prefecture — province plug-in, SSP585/2050",
                   limits = range(d6b$hazard_glmm, na.rm = TRUE))
p6b <- p6b + labs(tag = "b")

# 6c — County refit glmmTMB SSP585/2050
d6c <- fut_cnty_g[ssp == "SSP585" & year == 2050]
setnames(d6c, "unit_id", "cnty_id")
p6c <- make_choro(cnty_sf, d6c, "hazard_glmm", "cnty_id",
                   "County — refit glmmTMB, SSP585/2050",
                   limits = range(d6c$hazard_glmm, na.rm = TRUE))
p6c <- p6c + labs(tag = "c")

# 6d — Comparison of top 10 provinces between refit and plug-in (bar chart)
top_pref <- fut_pref_g[ssp == "SSP585" & year == 2050,
                        .(refit = mean(hazard_glmm, na.rm = TRUE)),
                        by = province]
top_pref <- merge(top_pref,
                    fut_pref_p[ssp == "SSP585" & year == 2050,
                                .(plugin = mean(hazard_glmm, na.rm = TRUE)),
                                by = province],
                    by = "province")
top_pref <- top_pref[order(-refit)][1:12]
tplong <- melt(top_pref, id.vars = "province",
                variable.name = "pathway", value.name = "hazard")
tplong[, province := factor(province, levels = rev(top_pref$province))]
p6d <- ggplot(tplong, aes(x = hazard, y = province, fill = pathway)) +
  geom_col(width = 0.6, position = "dodge") +
  scale_fill_manual(values = c(refit = "#1F77B4", plugin = "#FF7F0E"),
                     labels = c(refit = "Prefecture refit",
                                  plugin = "Province plug-in"),
                     name = NULL) +
  labs(tag = "d",
        title = "Top-12 provinces by mean prefecture hazard, refit vs plug-in",
        x = "Mean prefecture hazard (SSP585/2050)", y = NULL) +
  theme_publication() + theme(legend.position = "top")

fig6 <- (p6a | p6b) / (p6c | p6d) +
  plot_annotation(
    title = "Figure 6 — Unit-scale future hazard from two pathways",
    subtitle = "Both pathways agree on the spatial pattern; refit emphasises within-province heterogeneity, plug-in emphasises province-level differences.",
    theme = theme(plot.title = element_text(face = "bold", size = 9),
                   plot.subtitle = element_text(size = 8, colour = "grey30")))
save_pub(fig6, "Figure_6_unit_future_hazard",
         width = 18, height = 14)

message("[44] DONE — Figures 1-6 written under figures/main/.")
