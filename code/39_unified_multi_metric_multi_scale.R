# ============================================================
# Scientific question / 科学问题:
#   The reviewer should be able to see the climate × effort interaction
#   sign and magnitude across:
#       6 climate metrics × 4 effort specifications × 5 spatial scales
#       = 120 model rows
#   in ONE comparison table + ONE forest figure. This script wires the
#   already-fitted models from scripts 17/32/31/38 together with grid-
#   level fits from script 33 outputs.
#   把 6 气候 × 4 effort × 5 尺度 = 120 个模型在一张统一表中横向对比。
#
# Objective / 分析目标:
#   - Province scale (M3 + M4) × {temp_anom_z, climate_velocity_z,
#     warming_rate_z, climate_exposure_z, mahalanobis_dist_z,
#     precip_velocity_z} × {visits, days, observers, pc1}.
#   - Prefecture, county scales — same matrix (from script 38).
#   - 50 km, 100 km grid scales — same matrix using grid-native climate
#     (script 28) + grid-native effort (script 28b) + SDM-thresholded
#     risk set (script 33).
#   - For each cell of the matrix report: β, HR, 95 % CI, p, AIC,
#     marginal R², conditional R².
#
# Input data / 输入数据:
#   data/raw/hazard_risk_upgraded_complete_case.csv   (province scale)
#   data/derived/events_{prefecture,county}_risk_set.parquet
#   data/derived/events_{50,100}km_grid_risk_set_v2.parquet
#
# Expected output / 预期输出:
#   results/tables/table_unified_multi_metric_multi_scale.csv
#   results/tables/table_unified_multi_metric_multi_scale.parquet
#   figures/main/fig3b_unified_interaction_matrix.{pdf,png}  (heatmap)
#   figures/supplementary/figS_unified_forest_by_scale.{pdf,png}
#
# Key assumptions / 关键假设:
#   - For very large grid risk sets (>5M rows) we subsample 1M rows per
#     fit to keep glmmTMB within memory budget; the subsample is
#     stratified by event so all positive cases are retained.
#
# Main packages / 主要包: data.table, arrow, glmmTMB, ggplot2,
#   patchwork, broom.mixed.
# Output directory / 输出路径: results/tables/, figures/main/,
#   figures/supplementary/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(glmmTMB)
  library(ggplot2)
  library(patchwork)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_models.R"))
source(file.path("code", "utils", "utils_plots.R"))

set.seed(42)

CFG <- list(
  climate_metrics = c("temp_anom_z", "climate_velocity_z",
                       "warming_rate_z", "climate_exposure_z",
                       "mahalanobis_dist_z", "precip_velocity_z"),
  effort_metrics  = c("log_n_visits_z", "log_n_birding_days_z",
                       "log_n_observers_z", "effort_pc1_z"),
  scales = c("province", "prefecture", "county", "50km", "100km"),
  subsample_max = 1e6
)

# ---- 1. Loader per scale ---------------------------------------------------
load_for_scale <- function(scale) {
  switch(scale,
    "province" = {
      dt <- fread(path_raw("hazard_risk_upgraded_complete_case.csv"),
                  encoding = "UTF-8")
      dt[, unit := province]
      dt
    },
    "prefecture" = {
      p <- path_derived("events_prefecture_risk_set.parquet")
      if (!file.exists(p)) return(NULL)
      dt <- arrow::read_parquet(p) |> data.table::as.data.table()
      dt[, unit := unit_id]; dt
    },
    "county" = {
      p <- path_derived("events_county_risk_set.parquet")
      if (!file.exists(p)) return(NULL)
      dt <- arrow::read_parquet(p) |> data.table::as.data.table()
      dt[, unit := unit_id]; dt
    },
    "50km"  = ,
    "100km" = {
      km <- as.integer(sub("km", "", scale))
      p <- path_derived(glue("events_{km}km_grid_risk_set_v2.parquet"))
      if (!file.exists(p)) return(NULL)
      dt <- arrow::read_parquet(p) |> data.table::as.data.table()
      dt[, unit := grid_id]; dt
    }
  )
}

subsample_stratified <- function(dt, max_n) {
  if (nrow(dt) <= max_n) return(dt)
  pos <- dt[event == 1L]
  neg <- dt[event == 0L]
  n_neg_keep <- max_n - nrow(pos)
  if (n_neg_keep <= 0L) return(pos)
  neg_keep <- neg[sample(.N, n_neg_keep)]
  rbind(pos, neg_keep)
}

# ---- 2. Single (scale × climate × effort) fit ------------------------------
fit_one_cell <- function(scale, clim, eff) {
  dt <- load_for_scale(scale)
  if (is.null(dt) || !(clim %in% names(dt)) || !(eff %in% names(dt))) {
    return(data.table::data.table(scale = scale, climate = clim,
                                    effort = eff, status = "missing"))
  }
  dt <- dt[!is.na(get(clim)) & !is.na(get(eff)) & !is.na(event)]
  if (nrow(dt) < 200L) {
    return(data.table::data.table(scale = scale, climate = clim,
                                    effort = eff, status = "too_small"))
  }
  dt <- subsample_stratified(dt, CFG$subsample_max)
  setnames(dt, c(clim, eff), c("climate_z", "effort_z"))

  re_term <- if (scale == "province") "(1|species) + (1|province)"
             else                       "(1|species) + (1|unit)"
  forms <- list(
    M3 = sprintf("event ~ climate_z + effort_z + %s", re_term),
    M4 = sprintf("event ~ climate_z * effort_z + %s", re_term)
  )
  out <- data.table::data.table()
  for (nm in names(forms)) {
    fit <- tryCatch(
      glmmTMB::glmmTMB(as.formula(forms[[nm]]), data = dt,
                       family = binomial(link = "cloglog")),
      error = function(e) { message("    ", scale, "/", clim, "/", eff, "/", nm,
                                       " failed: ", conditionMessage(e)); NULL }
    )
    if (is.null(fit)) next
    coefs <- extract_coefs_manual(fit)
    interaction_row <- coefs[grepl(":", term)]
    if (nrow(interaction_row) == 0L) interaction_row <-
      data.table::data.table(term = NA, beta = NA, hr = NA,
                              hr.low = NA, hr.high = NA, p.value = NA)
    r2 <- marginal_R2(fit)
    out <- rbind(out, data.table::data.table(
      scale = scale, climate = clim, effort = eff, model = nm,
      term = ifelse(nm == "M4", interaction_row$term[1], "—"),
      beta = ifelse(nm == "M4", interaction_row$beta[1], NA),
      hr   = ifelse(nm == "M4", interaction_row$hr[1],   NA),
      hr.low  = ifelse(nm == "M4", interaction_row$hr.low[1], NA),
      hr.high = ifelse(nm == "M4", interaction_row$hr.high[1], NA),
      p.value = ifelse(nm == "M4", interaction_row$p.value[1], NA),
      AIC    = stats::AIC(fit),
      R2_marginal    = r2$R2_marginal,
      R2_conditional = r2$R2_conditional,
      n_rows = nrow(dt),
      n_events = sum(dt$event),
      status = "ok"
    ))
  }
  out
}

# ---- 3. Grid the matrix ----------------------------------------------------
grid_cells <- as.data.table(expand.grid(
  scale   = CFG$scales,
  climate = CFG$climate_metrics,
  effort  = CFG$effort_metrics,
  stringsAsFactors = FALSE
))
message(glue("[39] {nrow(grid_cells)} (scale × climate × effort) cells to fit"))

results <- data.table::rbindlist(lapply(seq_len(nrow(grid_cells)), function(i) {
  ce <- grid_cells[i]
  message(glue("[39] [{i}/{nrow(grid_cells)}] {ce$scale} | {ce$climate} | {ce$effort}"))
  fit_one_cell(ce$scale, ce$climate, ce$effort)
}), fill = TRUE)

ensure_dir(path_tables())
fwrite(results, path_tables("table_unified_multi_metric_multi_scale.csv"))
arrow::write_parquet(results,
                      path_tables("table_unified_multi_metric_multi_scale.parquet"),
                      compression = "snappy")
message(glue("[39] wrote {nrow(results)} rows to results/tables/"))

# ---- 4. Heatmap of interaction HR -----------------------------------------
heat <- results[model == "M4" & status == "ok",
                .(scale, climate, effort, hr,
                   sig = ifelse(p.value < 0.001, "***",
                          ifelse(p.value < 0.01, "**",
                            ifelse(p.value < 0.05, "*", ""))))]
heat[, scale := factor(scale,
                        levels = c("county", "prefecture", "province",
                                    "50km", "100km"))]

p_heat <- ggplot(heat, aes(x = effort, y = climate, fill = hr)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = paste0(round(hr, 2), sig)), size = 2.4) +
  facet_wrap(~ scale, nrow = 1) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "#F5F5F5", high = "#B40426",
                        midpoint = 1, name = "HR\n(interaction)") +
  labs(title = "Climate × effort interaction HR across metrics and scales",
        subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
        x = "Effort metric", y = "Climate metric") +
  theme_geb() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ensure_dir(path_main_fig())
ggsave(path_main_fig("fig3b_unified_interaction_matrix.pdf"),
       p_heat, width = 22, height = 9, units = "cm",
       device = grDevices::cairo_pdf)
ggsave(path_main_fig("fig3b_unified_interaction_matrix.png"),
       p_heat, width = 22, height = 9, units = "cm", dpi = 600)

# ---- 5. Forest plot per scale (supplementary) -----------------------------
forest_dt <- results[model == "M4" & status == "ok",
                      .(scale, climate, effort, hr, hr.low, hr.high,
                         label = paste(climate, "×", effort))]
forest_dt[, scale := factor(scale,
                              levels = c("county", "prefecture", "province",
                                          "50km", "100km"))]
p_forest <- ggplot(forest_dt,
                    aes(x = hr, y = label, colour = scale)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  position = position_dodge(width = 0.6),
                  height = 0.2, linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.6), size = 1.3) +
  scale_colour_manual(values = pal_cat[1:5]) +
  scale_x_continuous(trans = "log", breaks = c(0.5, 1, 1.5, 2, 3)) +
  labs(title = "Climate × effort interaction across metrics × scales",
        x = "Hazard ratio (log scale)", y = NULL) +
  theme_geb()

ensure_dir(path_supp_fig())
ggsave(path_supp_fig("figS_unified_forest_by_scale.pdf"),
       p_forest, width = 18, height = 22, units = "cm",
       device = grDevices::cairo_pdf)
ggsave(path_supp_fig("figS_unified_forest_by_scale.png"),
       p_forest, width = 18, height = 22, units = "cm", dpi = 600)

dump_session_info(path_logs("39_unified_multi_metric_multi_scale_sessionInfo.txt"))
message("[39] done.")
