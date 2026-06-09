# ============================================================
# Script: 53_v4_polished_figures.R
# Family: v4 polished publication figures
# Author: Chen-Chen Ding + Claude Opus 4.7
# Date  : 2026-05-23
#
# ------------------------------------------------------------
# Goals:
#   (a) Coef forest UPGRADED to raincloud + beeswarm layout (Fig 2 v4).
#   (b) AIC ladder + Akaike weight panel with M5 included (Fig 2 v4).
#   (c) Province / prefecture / county future-hazard choropleth maps
#       REBUILT with these requirements:
#         - Each scale uses ITS OWN administrative boundary as the
#           primary geometry (province → 省界, prefecture → 市界,
#           county → 县界); province-level outer frame overlaid.
#         - High-end palette: 'mako' / 'rocket' / 'viridis_c' from
#           viridisLite, with percentile clipping to avoid truncation
#           caused by long-tail outliers.
#         - Unified bottom horizontal legend, no internal axes.
#         - Polished panel layout via patchwork with consistent
#           figure size and tag style.
#   (d) Multi-scale forest (Fig 3 v4) with M0-M5 ladder + 3 scales.
#
# All outputs use the `_v4` suffix and never overwrite earlier figures.
# 输出全部 _v4 后缀，不覆盖。
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
  library(ggdist)
  library(ggbeeswarm)
  library(patchwork)
  library(viridisLite)
  library(scales)
  library(colorspace)
  library(ggrepel)
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

# ------------------------------------------------------------
# Unified publication theme + nature-inspired palettes
# ------------------------------------------------------------
theme_v4 <- function(base_size = 9) {
  theme_bw(base_size = base_size, base_family = "") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.18, colour = "grey92"),
          panel.border     = element_rect(linewidth = 0.4, colour = "grey20"),
          plot.title       = element_text(face = "bold", size = base_size + 1,
                                            margin = margin(b = 4)),
          plot.subtitle    = element_text(size = base_size - 1,
                                            colour = "grey30",
                                            margin = margin(b = 6)),
          plot.tag         = element_text(face = "bold",
                                            size = base_size + 3,
                                            family = ""),
          plot.tag.position = c(0.005, 0.995),
          plot.tag.location = "panel",
          strip.background = element_rect(fill = "grey96", colour = "grey80"),
          strip.text       = element_text(face = "bold", size = base_size),
          axis.text        = element_text(colour = "grey15"),
          axis.title       = element_text(colour = "grey5", size = base_size),
          axis.ticks       = element_line(linewidth = 0.3, colour = "grey40"),
          legend.position  = "right",
          legend.background = element_blank(),
          legend.title     = element_text(face = "bold", size = base_size - 1),
          legend.text      = element_text(size = base_size - 1),
          legend.key.size  = unit(0.32, "cm"))
}
# Spec palette — Nature-Communications-style ordered colour-blind safe
COL_SPEC <- c(spec_A = "#0072B2", spec_B = "#D55E00",
              spec_C = "#009E73", spec_D = "#CC79A7")
SPEC_LBL <- c(spec_A = "A: records",
              spec_B = "B: visits (headline)",
              spec_C = "C: PCA composite",
              spec_D = "D: birding-days")
SPEC_LVL <- c("spec_A","spec_B","spec_C","spec_D")
# Run palette
COL_RUN  <- c(v1 = "#3B4CC0", v2 = "#7F7F7F", v3 = "#B40426")
# Scale palette
COL_SCALE <- c(province = "#B40426", prefecture = "#0072B2",
                county   = "#009E73")
# Model palette (M0-M5)
COL_MODEL <- c(M0 = "#666666", M1 = "#1F77B4", M2 = "#2CA02C",
                M3 = "#9467BD", M4 = "#D62728", M5 = "#FF7F0E")

save_pub <- function(p, name, w, h) {
  ggsave(file.path(V2, "figures", "main", paste0(name, ".pdf")),
         p, width = w, height = h, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(file.path(V2, "figures", "main", paste0(name, ".png")),
         p, width = w, height = h, units = "cm", dpi = 600)
  message("[53] wrote ", name, ".{pdf,png} (", w, "×", h, " cm)")
}

# ============================================================
# PART 1 — Figure 2 v4 (raincloud + beeswarm + AIC w/ M5)
# ============================================================
message("[53] Figure 2 v4 — raincloud + beeswarm + AIC ladder w/ M5")
aic_v2 <- fread(file.path(V2, "results/tables/table_province_v2_with_m5_aic.csv"))
aic_v3 <- fread(file.path(V2, "results/tables/table_province_v3_with_m5_aic.csv"))
coefs_v2 <- fread(file.path(V2, "results/tables/table_province_v2_with_m5_coefs.csv"))
coefs_v3 <- fread(file.path(V2, "results/tables/table_province_v3_with_m5_coefs.csv"))

aic_v2[, run := "v2"]; aic_v3[, run := "v3"]
coefs_v2[, run := "v2"]; coefs_v3[, run := "v3"]
aic_all   <- rbind(aic_v2, aic_v3)
coefs_all <- rbind(coefs_v2, coefs_v3)

# -- Panel A: AIC ladder M0-M5, 4 specs (raincloud-style horizontal lollipops)
aic_show <- aic_all[run == "v3", ]   # use v3 (larger sample)
aic_show[, model := factor(model, levels = c("M0","M1","M2","M3","M4","M5"))]
aic_show[, spec_id := factor(spec_id, levels = SPEC_LVL)]

p_aic <- ggplot(aic_show, aes(x = dAIC, y = model, colour = model)) +
  geom_segment(aes(xend = 0, yend = model), linewidth = 0.5) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.1f", dAIC)),
             hjust = ifelse(aic_show$dAIC > 200, 1.1, -0.2),
             nudge_x = ifelse(aic_show$dAIC > 200, -8, 8),
             size = 2.3) +
  facet_wrap(~ spec_id, ncol = 2,
              labeller = labeller(spec_id = function(x) SPEC_LBL[x])) +
  scale_colour_manual(values = COL_MODEL, guide = "none") +
  labs(tag = "a",
        title = "Model-selection ladder M0-M5 — v3 risk set",
        subtitle = "M4 (climate × effort interaction) consistently wins; M5 (effort-as-offset) underperforms by ≥230 AIC except in PCA spec.",
        x = "ΔAIC vs best model", y = NULL) +
  theme_v4()

# -- Panel B: Interaction-HR forest with all 4 specs, v2 vs v3
int_dt <- coefs_all[model == "M4" & grepl(":", term),
                      .(run, spec_id,
                        beta, se, hr, hr.low, hr.high, p.value)]
int_dt[, spec_lbl := factor(SPEC_LBL[spec_id], levels = rev(SPEC_LBL))]
p_forest <- ggplot(int_dt, aes(x = hr, y = spec_lbl, colour = run)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  height = 0.18, linewidth = 0.55,
                  position = position_dodge(width = 0.5)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = sprintf("HR=%.3f", hr)),
             hjust = 0, nudge_x = 0.02, size = 2.4,
             position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = COL_RUN, name = "Risk set") +
  scale_x_continuous(trans = "log",
                      breaks = c(1.0, 1.1, 1.2, 1.3, 1.5)) +
  labs(tag = "b",
        title = "M4 interaction HR — v2 ↔ v3 across 4 effort specs",
        subtitle = "All 8 intervals are positive; v3 (red) has tighter intervals due to larger event base.",
        x = "Hazard ratio (95 % CI, log scale)", y = NULL) +
  theme_v4() + theme(legend.position = "top")

# -- Panel C: M4 vs M5 AIC dumbbell + raincloud (showing M4 always wins)
m45 <- dcast(aic_all[model %in% c("M4","M5")],
              run + spec_id + spec_label ~ model, value.var = "AIC")
m45[, dAIC := M5 - M4]
m45_long <- melt(m45, measure.vars = c("M4","M5"),
                  variable.name = "model", value.name = "AIC")
m45_long[, spec_lbl := paste0(SPEC_LBL[spec_id], " (", run, ")")]
m45_long[, spec_lbl := factor(spec_lbl,
                                 levels = m45_long[order(run, spec_id), unique(spec_lbl)])]
p_dumbbell <- ggplot(m45_long, aes(x = AIC, y = spec_lbl)) +
  geom_line(aes(group = spec_lbl), linewidth = 0.6, colour = "grey50") +
  geom_point(aes(colour = model, shape = model), size = 3) +
  scale_colour_manual(values = c(M4 = "#D62728", M5 = "#FF7F0E")) +
  scale_shape_manual(values = c(M4 = 16, M5 = 17)) +
  labs(tag = "c",
        title = "M4 (interaction) vs M5 (offset) — M4 wins 7/8",
        subtitle = "Lines connect M4 ↔ M5 within each (run × spec); leftmost point is preferred.",
        x = "AIC (lower = better)", y = NULL,
        colour = "Model", shape = "Model") +
  theme_v4() + theme(legend.position = "top")

# -- Panel D: Raincloud + beeswarm of HR distribution within (run, spec)
# Build a "bootstrap" cloud per (spec, run): n_boot draws from N(hr, se).
set.seed(42); n_boot <- 200
boot <- int_dt[, {
  draws <- rnorm(n_boot, mean = beta, sd = se)
  data.table(hr_draw = exp(draws))
}, by = .(run, spec_id)]
boot[, spec_lbl := factor(SPEC_LBL[spec_id], levels = SPEC_LVL)]
p_rain <- ggplot(boot, aes(x = spec_lbl, y = hr_draw, fill = run, colour = run)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  stat_halfeye(adjust = 0.6, .width = c(0.5, 0.95), justification = -0.15,
                interval_colour = "grey20", point_colour = "grey20",
                alpha = 0.65, position = position_dodge(width = 0.55)) +
  ggbeeswarm::geom_quasirandom(width = 0.05, size = 0.4, alpha = 0.3,
                                  dodge.width = 0.55) +
  scale_colour_manual(values = COL_RUN, name = "Risk set") +
  scale_fill_manual(values = COL_RUN, name = "Risk set") +
  scale_y_continuous(trans = "log",
                      breaks = c(0.95, 1.0, 1.1, 1.2, 1.3, 1.5)) +
  labs(tag = "d",
        title = "M4 interaction HR — bootstrap raincloud (v2 + v3 × 4 specs)",
        subtitle = "Half-violin = posterior density; bar = 50%/95% interval; dots = 200 draws ~ N(β, SE).",
        x = NULL, y = "Hazard ratio (log)") +
  theme_v4() + theme(axis.text.x = element_text(angle = 15, hjust = 1),
                       legend.position = "top")

fig2_v4 <- (p_aic | p_forest) / (p_dumbbell | p_rain) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(
    title = "Figure 2 v4 — Province headline + M5 offset + raincloud panorama",
    theme = theme(plot.title = element_text(face = "bold", size = 10)))
save_pub(fig2_v4, "Figure_2_v4_province_headline_M5_raincloud",
         w = 22, h = 16)

# ============================================================
# PART 2 — Figure 3 v4 (3-scale forest with M0-M5)
# ============================================================
message("[53] Figure 3 v4 — three-scale forest + M5")
prov_m45 <- aic_all[run == "v2" & spec_id == "spec_B",
                      .(scale = "province", model, AIC)]
pref_aic <- fread(file.path(V2, "results/tables/table_prefecture_county_aic.csv"))
pref_m4 <- pref_aic[scale == "prefecture" & model == "M4", .(model, AIC)]
pref_m4[, scale := "prefecture"]
cnty_m4 <- pref_aic[scale == "county" & model == "M4", .(model, AIC)]
cnty_m4[, scale := "county"]
all_3s <- rbind(prov_m45[, .(scale, model, AIC)],
                pref_m4[, .(scale, model, AIC)],
                cnty_m4[, .(scale, model, AIC)])

# Hazard ratios across scales
prov_hr <- coefs_all[run == "v2" & spec_id == "spec_B" & model == "M4" &
                       grepl(":", term),
                       .(scale = "province", hr, hr.low, hr.high, p.value)]
pref_co <- fread(file.path(V2, "results/tables/table_prefecture_coefs.csv"))
cnty_co <- fread(file.path(V2, "results/tables/table_county_coefs.csv"))
pref_hr <- pref_co[model == "M4" & grepl(":", term),
                     .(scale = "prefecture", hr, hr.low, hr.high, p.value)]
cnty_hr <- cnty_co[model == "M4" & grepl(":", term),
                     .(scale = "county", hr, hr.low, hr.high, p.value)]
three <- rbind(prov_hr, pref_hr, cnty_hr)
three[, scale := factor(scale,
                          levels = c("county","prefecture","province"))]

p_three <- ggplot(three, aes(x = hr, y = scale, colour = scale)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  height = 0.16, linewidth = 0.55) +
  geom_point(size = 3.6) +
  geom_text(aes(label = sprintf("HR = %.3f\np = %.0e", hr, p.value)),
             hjust = 0, nudge_x = 0.06, size = 2.7) +
  scale_colour_manual(values = COL_SCALE, guide = "none") +
  scale_x_continuous(trans = "log",
                      breaks = c(0.95, 1.0, 1.1, 1.2, 1.3, 1.5),
                      limits = c(0.95, 1.65)) +
  labs(title = "M4 climate × effort interaction across 3 admin scales",
        subtitle = "Within-province grain dilutes the interaction signal but never reverses it.",
        x = "Hazard ratio (95% CI, log scale)", y = NULL) +
  theme_v4()
save_pub(p_three, "Figure_3_v4_three_scale_forest_M5",
         w = 15, h = 9)

# ============================================================
# PART 3 — Polished choropleth maps (v4): each scale uses its OWN border
# ============================================================
message("[53] Polished choropleth maps — province / prefecture / county")

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

prov_sf <- st_read(file.path(SHP, "省（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4524) |> st_make_valid()
prov_sf$province <- unname(PROV_CN_EN[as.character(prov_sf[["省"]])])
pref_sf <- st_read(file.path(SHP, "市（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4524) |> st_make_valid()
pref_sf$pref_id <- paste0("PREF_", sprintf("%04d", seq_len(nrow(pref_sf))))
cnty_sf <- st_read(file.path(SHP, "县（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4524) |> st_make_valid()
cnty_sf$cnty_id <- paste0("CNTY_", sprintf("%05d", seq_len(nrow(cnty_sf))))
ninedash <- tryCatch(
  st_transform(st_read(file.path(SHP, "九段线.shp"), quiet = TRUE), 4524),
  error = function(e) NULL)

# Common bbox derived from province sf so all maps align.
bb <- st_bbox(prov_sf)
xlim_v <- c(bb["xmin"] - 1e5, bb["xmax"] + 1e5)
ylim_v <- c(bb["ymin"] - 1e5, bb["ymax"] + 1e5)

# Cleaner choropleth helper that:
#  - uses the unit's OWN polygons as fill
#  - overlays the province border as a stronger outer frame
#  - clips top 1 % outliers using percentile rescaling for the fill scale
#  - places the colour-bar at the bottom (horizontal) for clean layout
clean_choro <- function(sf_unit, dt, value_col, id_col,
                          title, palette = "mako", direction = -1,
                          legend_label = "Hazard\n(prob.)") {
  joined <- merge(sf_unit, dt[, c(id_col, value_col), with = FALSE],
                   by = id_col)
  vals <- joined[[value_col]]
  qclip <- stats::quantile(vals, c(0.01, 0.99), na.rm = TRUE)
  joined$.fill <- pmin(pmax(vals, qclip[1]), qclip[2])
  p <- ggplot() +
    geom_sf(data = joined, aes(fill = .fill), colour = "grey60",
            linewidth = 0.08)
  p <- p +
    geom_sf(data = prov_sf, fill = NA, colour = "grey15",
            linewidth = 0.28)
  if (!is.null(ninedash)) p <- p +
    geom_sf(data = ninedash, fill = NA, colour = "grey15",
            linewidth = 0.2)
  p +
    scale_fill_viridis_c(option = palette, direction = direction,
                          name = legend_label,
                          labels = scales::label_number(accuracy = 0.01),
                          guide = guide_colorbar(barwidth = unit(4.2, "cm"),
                                                    barheight = unit(0.32, "cm"),
                                                    title.position = "top",
                                                    title.hjust = 0.5)) +
    coord_sf(xlim = xlim_v, ylim = ylim_v, datum = NA, expand = FALSE) +
    labs(title = title) +
    theme_v4(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          legend.box = "vertical",
          plot.margin = margin(3, 3, 3, 3))
}

# Six-panel composer (SSP × year × model)
mk_six <- function(sf_unit, fut, value_col, id_col, file_label,
                    palette, scale_label) {
  future_years <- c(2030, 2050, 2080)
  fut <- copy(fut)
  fut[, panel := factor(paste0(ssp, " — ", year),
                          levels = c(paste0("SSP245 — ", future_years),
                                      paste0("SSP585 — ", future_years)))]
  panels <- lapply(split(fut, fut$panel), function(d)
    clean_choro(sf_unit, d, value_col, id_col,
                 title = unique(d$panel), palette = palette))
  combined <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") +
    plot_annotation(
      title = sprintf("%s-scale future hazard (v4 polished, %s palette, percentile-clipped 1-99%%)",
                       scale_label, palette),
      subtitle = paste0("Each panel uses its own administrative boundary; ",
                          "province border overlaid for orientation. ",
                          "Empirical SSP perturbation, effort frozen at 2024."),
      theme = theme(plot.title = element_text(face = "bold", size = 10),
                     plot.subtitle = element_text(size = 8, colour = "grey30"))) &
    theme(legend.position = "bottom")
  save_pub(combined, file_label, w = 22, h = 14)
}

# Province glmmTMB + XGBoost
fut_prov_g <- fread(file.path(V2,
  "results/forecasts/table_province_future_glmmTMB.csv"))
fut_prov_x <- fread(file.path(V2,
  "results/forecasts/table_province_future_xgboost.csv"))
mk_six(prov_sf, fut_prov_g, "hazard_glmm", "province",
        "Figure_5_v4_province_future_glmmTMB",
        palette = "mako", scale_label = "Province")
mk_six(prov_sf, fut_prov_x, "hazard_xgb", "province",
        "Figure_5_v4_province_future_xgboost",
        palette = "rocket", scale_label = "Province")

# Prefecture refit
fut_pref_g <- fread(file.path(V2,
  "results/forecasts/table_prefecture_future_glmmTMB.csv"))
setnames(fut_pref_g, intersect(c("unit_id","pref_id"), names(fut_pref_g))[1],
         "pref_id")
fut_pref_x <- fread(file.path(V2,
  "results/forecasts/table_prefecture_future_xgboost.csv"))
setnames(fut_pref_x, intersect(c("unit_id","pref_id"), names(fut_pref_x))[1],
         "pref_id")
mk_six(pref_sf, fut_pref_g, "hazard_glmm", "pref_id",
        "Figure_6_v4_prefecture_future_glmmTMB",
        palette = "mako", scale_label = "Prefecture")
mk_six(pref_sf, fut_pref_x, "hazard_xgb", "pref_id",
        "Figure_6_v4_prefecture_future_xgboost",
        palette = "rocket", scale_label = "Prefecture")

# County refit
fut_cnty_g <- fread(file.path(V2,
  "results/forecasts/table_county_future_glmmTMB.csv"))
setnames(fut_cnty_g, intersect(c("unit_id","cnty_id"), names(fut_cnty_g))[1],
         "cnty_id")
fut_cnty_x <- fread(file.path(V2,
  "results/forecasts/table_county_future_xgboost.csv"))
setnames(fut_cnty_x, intersect(c("unit_id","cnty_id"), names(fut_cnty_x))[1],
         "cnty_id")
mk_six(cnty_sf, fut_cnty_g, "hazard_glmm", "cnty_id",
        "Figure_7_v4_county_future_glmmTMB",
        palette = "mako", scale_label = "County")
mk_six(cnty_sf, fut_cnty_x, "hazard_xgb", "cnty_id",
        "Figure_7_v4_county_future_xgboost",
        palette = "rocket", scale_label = "County")

message("[53] DONE — Figures 2 v4 + 3 v4 + 5 v4 + 6 v4 + 7 v4 saved.")
