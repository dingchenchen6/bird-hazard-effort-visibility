# ============================================================
# Scientific question / 科学问题:
#   Reviewers will reasonably ask whether the climate × effort
#   interaction is an artefact of the chosen spatial unit. We refit
#   M1–M5 at four spatial scales (province, 50 km, 100 km, 200 km)
#   and compare the coefficient of the interaction term — the
#   "MAUP elasticity".
#   把 M1–M5 在 4 个尺度（省/50km/100km/200km）重拟合，比较
#   交互项系数的弹性，回答 MAUP 质疑。
#
# Objective / 分析目标:
#   - For each scale, build a province- or grid-aggregated risk set.
#   - Fit M1–M5 with the canonical formula family.
#   - Tabulate β_interaction, HR_interaction, 95 % CI, AIC.
#   - Report scale-by-scale elasticity in a single CSV.
#
# Input data / 输入数据:
#   data/raw/hazard_risk_upgraded_complete_case.csv  (province risk)
#   data/derived/events_*_grid_risk_set_v2.parquet   (grid risk; script 33)
#
# Expected output / 预期输出:
#   results/sensitivity/table_maup_elasticity.csv
#   results/sensitivity/table_maup_full_coefs.csv
#   figures/supplementary/figS8_maup_elasticity.{pdf,png}
#
# Key assumptions / 关键假设:
#   - 200-km grid is built on-the-fly using utils_spatial::build_china_grid.
#
# Main packages / 主要包: glmmTMB, data.table, arrow, ggplot2.
# Output directory / 输出路径: results/sensitivity/, figures/supplementary/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(glmmTMB)
  library(ggplot2)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_spatial.R"))
source(file.path("code", "utils", "utils_models.R"))
source(file.path("code", "utils", "utils_plots.R"))

set.seed(42)

CFG <- list(
  scales = c("province", "50km", "100km", "200km")
)

load_scale <- function(scale) {
  switch(scale,
    "province" = {
      dt <- fread(path_raw("hazard_risk_upgraded_complete_case.csv"),
                  encoding = "UTF-8")
      dt[, unit := province]
      dt
    },
    "50km" = {
      p <- path_derived("events_50km_grid_risk_set_v2.parquet")
      if (!file.exists(p)) {
        warning("[31] 50km risk set not built; run script 33 first.")
        return(NULL)
      }
      dt <- arrow::read_parquet(p) |> data.table::as.data.table()
      dt[, unit := grid_id]
      dt
    },
    "100km" = {
      p <- path_derived("events_100km_grid_risk_set_v2.parquet")
      if (!file.exists(p)) {
        warning("[31] 100km risk set not built; run script 33 first.")
        return(NULL)
      }
      dt <- arrow::read_parquet(p) |> data.table::as.data.table()
      dt[, unit := grid_id]
      dt
    },
    "200km" = {
      # Build 200km grid on-the-fly. 临时构造 200km 网格。
      boundary <- read_gs2019_basemap("province") |> to_albers()
      g200 <- build_china_grid(boundary, size_km = 200, min_land_frac = 0.4)
      # Re-aggregate 100km risk set into 200km buckets. 把 100km 风险集
      # 重映射到 200km 桶。
      p100 <- path_derived("events_100km_grid_risk_set_v2.parquet")
      if (!file.exists(p100)) {
        warning("[31] 200km requires 100km risk set; run script 33 first.")
        return(NULL)
      }
      r100 <- arrow::read_parquet(p100) |> data.table::as.data.table()
      g100_sf <- sf::st_read(path_derived("grid_100km_sf.gpkg"), quiet = TRUE)
      cen100 <- sf::st_centroid(to_albers(g100_sf))
      map200 <- sf::st_join(cen100, g200, join = sf::st_within)
      map_dt <- data.table::data.table(grid_id = g100_sf$grid_id,
                                        grid200 = map200$grid_id)
      r100 <- merge(r100, map_dt, by = "grid_id")
      r100[, unit := grid200]
      r100
    }
  )
}

extract_for_scale <- function(scale) {
  dt <- load_scale(scale)
  if (is.null(dt)) return(NULL)
  message(glue("[31] Fitting M1–M5 at scale = {scale} (n = {nrow(dt)})"))
  climate_candidates <- c("climate_velocity_z", "temp_anom_z",
                          "temp_anom_native_z", "climate_velocity_native_z")
  climate_col <- intersect(climate_candidates, names(dt))[1]
  effort_candidates <- c("log_n_visits_z", "effort_pc1_z",
                          "log_n_birding_days_z")
  effort_col <- intersect(effort_candidates, names(dt))[1]
  if (is.na(climate_col) || is.na(effort_col)) {
    warning(glue("[31] missing climate/effort at scale {scale}; skip."))
    return(NULL)
  }
  dt[, climate_z := get(climate_col)]
  dt[, effort_z  := get(effort_col)]
  dt <- dt[!is.na(event) & !is.na(climate_z) & !is.na(effort_z)]

  re_term <- if (scale == "province") "(1|species) + (1|province)"
              else                      "(1|species) + (1|unit)"
  forms <- list(
    M1 = sprintf("event ~ effort_z + %s", re_term),
    M2 = sprintf("event ~ climate_z + %s", re_term),
    M3 = sprintf("event ~ climate_z + effort_z + %s", re_term),
    M4 = sprintf("event ~ climate_z * effort_z + %s", re_term),
    M5 = sprintf("event ~ climate_z + %s + offset(log1p(person_hours_raw))", re_term)
  )
  if (!"person_hours_raw" %in% names(dt)) forms$M5 <- NULL

  models <- list()
  for (nm in names(forms)) {
    f <- as.formula(forms[[nm]])
    fit <- tryCatch(
      glmmTMB::glmmTMB(f, data = dt, family = binomial(link = "cloglog")),
      error = function(e) { message("   ", nm, " failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(fit)) models[[nm]] <- fit
  }

  coef_tab <- data.table::rbindlist(lapply(names(models), function(nm) {
    out <- extract_coefs_manual(models[[nm]])
    out[, model := nm]
    out[, scale := scale]
    out
  }), use.names = TRUE, fill = TRUE)
  data.table::setcolorder(coef_tab, c("scale", "model", "term"))
  coef_tab
}

results <- data.table::rbindlist(lapply(CFG$scales, extract_for_scale),
                                  fill = TRUE)
ensure_dir(path_sensitivity())
fwrite(results, path_sensitivity("table_maup_full_coefs.csv"))

# Elasticity table on the interaction term. 弹性表。
interaction_dt <- results[grepl(":", term) & model == "M4",
                          .(scale, beta, hr, hr.low, hr.high, p.value)]
fwrite(interaction_dt, path_sensitivity("table_maup_elasticity.csv"))

# Figure: forest plot across scales. 跨尺度森林图。
if (nrow(interaction_dt) > 0L) {
  fig <- forest_plot(interaction_dt, term_col = "scale",
                      title = "MAUP elasticity — climate × effort interaction (M4)")
  ggsave(path_supp_fig("figS8_maup_elasticity.pdf"),
         fig, width = 12, height = 8, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(path_supp_fig("figS8_maup_elasticity.png"),
         fig, width = 12, height = 8, units = "cm", dpi = 600)
}

message("[31] MAUP elasticity:")
print(interaction_dt)

dump_session_info(path_logs("31_maup_sensitivity_sessionInfo.txt"))
message("[31] done.")
