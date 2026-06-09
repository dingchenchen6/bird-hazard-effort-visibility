# ============================================================
# Scientific question / 科学问题:
#   How will new-bird-record hazard change under future climate
#   (SSP245 / SSP585) and effort scenarios at province, prefecture,
#   county, and 100km grid scales?
#   在省级、市级、县级和100km网格尺度上，未来气候和调查
#   努力情景下新纪录风险如何变化？
#
# Objective / 分析目标:
#   - Refit M4 (climate × effort interaction) at each scale
#   - Build future prediction panels (SSP245/585 × baseline/trend/doubled effort)
#   - Predict hazard and generate choropleth maps at all 4 scales
#   - Output: prediction CSVs + map figures, all under outputs_multiscale/
#
# Input data / 输入数据:
#   - Risk sets from 51_multiscale_full.R (outputs_multiscale/data/derived/)
#   - Province climate: data/raw/climate_metrics_province_year.csv
#   - Province effort:  data/raw/effort_panel_upgraded.csv
#   - Shapefiles: data/spatial/basemap_GS2019_1822/
#
# Expected output / 预期输出:
#   outputs_multiscale/results/tables/
#     table_multiscale_future_{province,prefecture,county,grid_100km}.csv
#   outputs_multiscale/figures/main/
#     fig_multiscale_future_{province,prefecture,county,grid_100km}_2050_ssp585.{pdf,png}
#     fig_multiscale_future_comparison_2050.{pdf,png}
#   outputs_multiscale/logs/53_multiscale_future.log
#
# Key assumptions / 关键假设:
#   - M4 glmmTMB model with cloglog link is the prediction engine
#   - Future climate = current z-score + SSP-adjusted delta per decade
#   - Effort scenarios: baseline (2024 frozen), trend (linear extrapolation),
#     doubled (2× baseline)
#   - Prefecture/county climate = inherited from province (same as v2)
#   - Grid climate = time-invariant bio vars + province climate overlay
#
# Main packages / 主要包: data.table, sf, glmmTMB, ggplot2.
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

log_file <- file.path(OUT, "logs", "53_multiscale_future.log")
log <- function(...) {
  msg <- sprintf("[53 %s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(msg); cat(msg, file = log_file, append = TRUE)
}
log("=== 53_multiscale_future_prediction.R START ===")

# ---- 1. Packages ----------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(glmmTMB)
  library(ggplot2)
})
sf::sf_use_s2(FALSE)
options(warn = 1)
set.seed(42)

zify <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# ---- 2. Load risk sets + province climate/effort --------------------------
log("=== Step 2: Loading data ===")

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

# Province climate and effort for future scenario construction
prov_clim <- fread(file.path(BASE_DIR, "data", "raw",
                             "climate_metrics_province_year.csv"),
                   encoding = "UTF-8")
setnames(prov_clim, tolower(names(prov_clim)))
prov_eff <- fread(file.path(BASE_DIR, "data", "raw",
                            "effort_panel_upgraded.csv"),
                  encoding = "UTF-8")
setnames(prov_eff, tolower(names(prov_eff)))

# Current baseline (2024)
current_eff <- prov_eff[year == 2024, .(province, log_effort_visits_z, effort_pc1_z)]
# Effort trend per province
eff_trends <- prov_eff[year >= 2002 & year <= 2024,
  .(effort_trend_z = tryCatch(coef(lm(log_effort_visits_z ~ year))[2],
                               error = function(e) 0)),
  by = province]

# Climate scenario deltas (per decade, from CMIP6 literature approximations)
# SSP245: ~+0.3°C/decade → ~+0.3/temp_sd z-units/decade
# SSP585: ~+0.6°C/decade → ~+0.6/temp_sd z-units/decade
clim_vel_sd  <- sd(prov_clim$climate_velocity_z, na.rm = TRUE)
warming_sd   <- sd(prov_clim$warming_rate_z, na.rm = TRUE)

# ---- 3. Shapefiles for mapping --------------------------------------------
log("=== Step 3: Loading shapefiles ===")

PROV_CN_EN <- c(
  "北京市" = "Beijing", "天津市" = "Tianjin", "河北省" = "Hebei",
  "山西省" = "Shanxi", "内蒙古自治区" = "Inner Mongolia",
  "辽宁省" = "Liaoning", "吉林省" = "Jilin", "黑龙江省" = "Heilongjiang",
  "上海市" = "Shanghai", "江苏省" = "Jiangsu", "浙江省" = "Zhejiang",
  "安徽省" = "Anhui", "福建省" = "Fujian", "江西省" = "Jiangxi",
  "山东省" = "Shandong", "河南省" = "Henan", "湖北省" = "Hubei",
  "湖南省" = "Hunan", "广东省" = "Guangdong",
  "广西壮族自治区" = "Guangxi", "海南省" = "Hainan",
  "重庆市" = "Chongqing", "四川省" = "Sichuan", "贵州省" = "Guizhou",
  "云南省" = "Yunnan", "西藏自治区" = "Tibet",
  "陕西省" = "Shaanxi", "甘肃省" = "Gansu", "青海省" = "Qinghai",
  "宁夏回族自治区" = "Ningxia", "新疆维吾尔自治区" = "Xinjiang",
  "台湾省" = "Taiwan", "香港特别行政区" = "Hong Kong",
  "澳门特别行政区" = "Macau")

# Prefer china_shp (complete .shx), then GS2019
shp_dirs_ordered <- c(
  file.path(BASE_DIR, "data", "spatial", "china_shp"),
  file.path(BASE_DIR, "data", "spatial", "basemap_GS2019_1822"),
  file.path(dirname(BASE_DIR), "bird_hazard_model_effort_upgrade",
            "2019中国地图-审图号GS(2019)1822号"),
  file.path(BASE_DIR, "2019中国地图-审图号GS(2019)1822号")
)

find_shp <- function(prefix) {
  for (d in shp_dirs_ordered) {
    if (!dir.exists(d)) next
    hits <- list.files(d, pattern = paste0("^", prefix, ".*\\.shp$"),
                       full.names = TRUE)
    hits <- hits[!grepl("境界线", hits)]
    eq <- hits[grepl("等积投影", hits)]
    poly <- hits[!grepl("线", hits)]
    candidates <- if (length(eq) > 0L) eq else if (length(poly) > 0L) poly else hits
    # Check for .shx
    for (h in candidates) {
      shx <- sub("\\.shp$", ".shx", h)
      if (file.exists(shx) || file.exists(sub("\\.shp$", ".SHX", h))) return(h)
    }
    if (length(candidates) > 0L) return(candidates[1])
  }
  NULL
}

# Albers equal-area CRS for China
china_crs <- "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +ellps=GRS80 +units=m +no_defs"

load_shp <- function(prefix) {
  p <- find_shp(prefix)
  if (is.null(p)) return(NULL)
  sf <- st_read(p, quiet = TRUE) |> st_make_valid() |> st_transform(china_crs)
  # Assign province_en via Chinese name column
  cn_col <- names(sf)[vapply(sf, function(x) {
    is.character(x) && any(grepl("[\u4e00-\u9fff]", x))
  }, logical(1))][1]
  if (!is.na(cn_col)) {
    sf$province_en <- unname(PROV_CN_EN[as.character(sf[[cn_col]])])
    sf$province_en[is.na(sf$province_en)] <- as.character(sf[[cn_col]])[is.na(sf$province_en)]
  }
  sf
}

prov_sf <- load_shp("省")
pref_sf <- load_shp("市")
cnty_sf <- load_shp("县")

log("province shp: ", if (!is.null(prov_sf)) nrow(prov_sf) else "missing")
log("prefecture shp: ", if (!is.null(pref_sf)) nrow(pref_sf) else "missing")
log("county shp: ", if (!is.null(cnty_sf)) nrow(cnty_sf) else "missing")

# ---- 4. Fit M4 + predict future at each scale ----------------------------
log("=== Step 4: Fitting M4 + future prediction ===")

future_years <- c(2030, 2035, 2040, 2045, 2050)
climate_scenarios <- c("current", "ssp245", "ssp585")
effort_scenarios  <- c("baseline", "trend", "doubled")

# 存储各尺度预测结果
all_preds <- list()

for (sc in names(risk_sets)) {
  log("  --- ", sc, " ---")
  dt <- copy(risk_sets[[sc]])

  # Determine unit column
  unit_col <- if ("unit_id" %in% names(dt)) "unit_id"
              else if ("grid_id" %in% names(dt)) "grid_id"
              else "province"
  dt[, unit_val := factor(get(unit_col))]
  dt[, species := factor(species)]

  # Ensure climate_z, effort_z exist
  if (!"climate_z" %in% names(dt)) dt[, climate_z := zify(climate_velocity_z)]
  if (!"effort_z"  %in% names(dt)) {
    eff_col <- intersect(c("log_n_events_z", "log_effort_visits_z",
                           "effort_pc1_z"), names(dt))[1]
    if (!is.na(eff_col)) dt[, effort_z := get(eff_col)]
    else next
  }

  re_form <- if (unit_col == "province") "(1|species) + (1|province)"
             else paste0("(1|species) + (1|unit_val)")

  # Fit M4
  fml_m4 <- as.formula(paste("event ~ climate_z * effort_z +", re_form))
  fit_m4 <- tryCatch(
    glmmTMB(fml_m4, data = dt, family = binomial(link = "cloglog")),
    error = function(e) { log("  M4 failed: ", conditionMessage(e)); NULL })
  if (is.null(fit_m4)) { log("  skipping ", sc); next }
  log("  M4 fitted, AIC=", round(AIC(fit_m4), 1))

  # Get province column for effort/climate assignment
  # 对于grid_100km，可能没有province列——用grid_id中的省份编码或跳过
  if (!"province" %in% names(dt) && unit_col != "province") {
    # 对于grid: 尝试从grid_id恢复province信息
    # 如果没有province列，则使用所有省份的均值作为替代
    log("  no province column, using overall mean for future scenarios")
    dt[, province := "ALL"]
  }

  # Build future prediction panel
  units <- unique(dt[, .(get(unit_col), get(if ("province" %in% names(dt)) "province" else unit_col))])
  setnames(units, c(unit_col, "province"))

  fut <- as.data.table(expand.grid(
    unit_val_tmp = unique(dt[[unit_col]]),
    year = future_years,
    climate_scenario = climate_scenarios,
    effort_scenario = effort_scenarios,
    stringsAsFactors = FALSE))
  setnames(fut, "unit_val_tmp", unit_col)

  # Merge province for each unit
  unit_prov <- unique(dt[, c(unit_col, "province"), with = FALSE])
  if (nrow(unit_prov) > 0L && "province" %in% names(dt)) {
    fut <- merge(fut, unit_prov, by = unit_col, all.x = TRUE)
  }

  # Baseline climate and effort per province (2024)
  fut <- merge(fut, current_eff, by = "province", all.x = TRUE)
  fut <- merge(fut, eff_trends, by = "province", all.x = TRUE)

  # Baseline climate z per province (2024)
  prov_clim_2024 <- prov_clim[year == 2024,
    .(province, climate_velocity_z, warming_rate_z)]
  fut <- merge(fut, prov_clim_2024, by = "province", all.x = TRUE)

  # Effort scenarios
  fut[, effort_z := fcase(
    effort_scenario == "baseline", log_effort_visits_z,
    effort_scenario == "trend",
      log_effort_visits_z + effort_trend_z * (year - 2024),
    effort_scenario == "doubled", log_effort_visits_z * 2)]
  fut[is.na(effort_z), effort_z := 0]

  # Climate scenarios (adjust climate_velocity_z)
  fut[, climate_z := fcase(
    climate_scenario == "current", climate_velocity_z,
    climate_scenario == "ssp245",
      climate_velocity_z + 0.3 / clim_vel_sd * (year - 2024) / 10,
    climate_scenario == "ssp585",
      climate_velocity_z + 0.6 / clim_vel_sd * (year - 2024) / 10)]
  fut[is.na(climate_z), climate_z := 0]

  # Set unit_val and species for prediction
  fut[, unit_val := factor(get(unit_col), levels = levels(dt$unit_val))]
  # Use most common species as representative (for population-level prediction)
  top_sp <- names(sort(table(dt$species), decreasing = TRUE))[1]
  fut[, species := factor(top_sp, levels = levels(dt$species))]

  # Predict
  pred <- predict(fit_m4, newdata = fut, type = "response", se.fit = TRUE)
  fut[, hazard := pred$fit]
  fut[, hazard_se := pred$se.fit]
  fut[, scale := sc]

  # Save prediction table
  out_cols <- c(unit_col, "province", "year", "climate_scenario",
                "effort_scenario", "climate_z", "effort_z",
                "hazard", "hazard_se", "scale")
  out_cols <- intersect(out_cols, names(fut))
  fwrite(fut[, ..out_cols],
         file.path(OUT, "results", "forecasts",
                   paste0("table_multiscale_future_", sc, ".csv")))
  log("  predictions saved: ", nrow(fut), " rows")

  all_preds[[sc]] <- fut[, ..out_cols]
}

# ---- 5. Choropleth maps ---------------------------------------------------
log("=== Step 5: Generating choropleth maps ===")

theme_map <- theme_bw(base_size = 9) +
  theme(axis.title = element_blank(), axis.text = element_text(size = 6),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70", linewidth = 0.3),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"), legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8))

make_choropleth <- function(pred_dt, shp, unit_col, title, fill_limits = c(0, 0.1)) {
  if (is.null(shp)) { log("  shp missing, skip map"); return(NULL) }
  # Merge predictions with shapefile
  shp_dt <- st_drop_geometry(shp)
  # Find a join key
  join_key <- intersect(c("province_en", "province"), names(shp_dt))[1]
  if (is.na(join_key)) join_key <- names(shp_dt)[1]

  pred_summary <- pred_dt[year == 2050 & climate_scenario == "ssp585" &
                           effort_scenario == "baseline",
    .(hazard = mean(hazard, na.rm = TRUE)),
    by = setdiff(intersect(names(pred_dt), c(unit_col, "province", "province_en")),
                 c("year", "climate_scenario", "effort_scenario"))]

  if (nrow(pred_summary) == 0L) return(NULL)

  # Try to match: if unit_col is "province", join on province_en
  if (unit_col == "province" && "province_en" %in% names(shp_dt)) {
    setnames(pred_summary, "province", "province_en")
    shp_merged <- merge(shp, pred_summary, by = "province_en", all.x = TRUE)
  } else {
    # For prefecture/county: use row index as fallback
    shp_merged <- cbind(shp, hazard = pred_summary$hazard[seq_len(nrow(shp))])
  }

  ggplot(shp_merged) +
    geom_sf(aes(fill = hazard), colour = NA, linewidth = 0) +
    scale_fill_gradientn(colours = c("#2166AC", "#67A9CF", "#F7F7F7",
                                      "#EF8A62", "#B2182B"),
                         na.value = "grey90",
                         limits = fill_limits,
                         name = "Hazard") +
    labs(title = title) +
    theme_map
}

# Province map
if ("province" %in% names(all_preds) && !is.null(prov_sf)) {
  p_prov <- make_choropleth(all_preds[["province"]], prov_sf,
                            "province", "Province: 2050 SSP585")
  if (!is.null(p_prov)) {
    ggsave(file.path(OUT, "figures", "main",
                     "fig_multiscale_future_province_2050_ssp585.pdf"),
           p_prov, width = 12, height = 8, units = "cm",
           device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main",
                     "fig_multiscale_future_province_2050_ssp585.png"),
           p_prov, width = 12, height = 8, units = "cm", dpi = 600)
    log("province map saved")
  }
}

# Prefecture map
if ("prefecture" %in% names(all_preds) && !is.null(pref_sf)) {
  p_pref <- make_choropleth(all_preds[["prefecture"]], pref_sf,
                            "prefecture", "Prefecture: 2050 SSP585")
  if (!is.null(p_pref)) {
    ggsave(file.path(OUT, "figures", "main",
                     "fig_multiscale_future_prefecture_2050_ssp585.pdf"),
           p_pref, width = 12, height = 8, units = "cm",
           device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main",
                     "fig_multiscale_future_prefecture_2050_ssp585.png"),
           p_pref, width = 12, height = 8, units = "cm", dpi = 600)
    log("prefecture map saved")
  }
}

# County map
if ("county" %in% names(all_preds) && !is.null(cnty_sf)) {
  p_cnty <- make_choropleth(all_preds[["county"]], cnty_sf,
                            "county", "County: 2050 SSP585")
  if (!is.null(p_cnty)) {
    ggsave(file.path(OUT, "figures", "main",
                     "fig_multiscale_future_county_2050_ssp585.pdf"),
           p_cnty, width = 12, height = 8, units = "cm",
           device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main",
                     "fig_multiscale_future_county_2050_ssp585.png"),
           p_cnty, width = 12, height = 8, units = "cm", dpi = 600)
    log("county map saved")
  }
}

# ---- 6. Cross-scale comparison figure -------------------------------------
log("=== Step 6: Cross-scale comparison ===")

all_pred_dt <- rbindlist(all_preds, fill = TRUE)
if (nrow(all_pred_dt) > 0L) {
  # National mean hazard by scale × scenario × year
  natl <- all_pred_dt[,
    .(hazard_mean = mean(hazard, na.rm = TRUE),
      hazard_se = mean(hazard_se, na.rm = TRUE)),
    by = .(scale, year, climate_scenario, effort_scenario)]

  natl[, scale := factor(scale,
    levels = c("province", "prefecture", "county", "grid_100km"),
    labels = c("Province", "Prefecture", "County", "100km grid"))]

  p_comp <- ggplot(natl[effort_scenario == "baseline"],
                    aes(x = year, y = hazard_mean,
                        colour = climate_scenario, linetype = scale)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_colour_manual(values = c("current" = "#2166AC",
                                    "ssp245" = "#F4A582",
                                    "ssp585" = "#B2182B"),
                        name = "Climate") +
    scale_linetype(name = "Scale") +
    labs(x = "Year", y = "Mean predicted hazard",
         title = "Future hazard trajectories across spatial scales",
         subtitle = "Baseline effort scenario; M4 glmmTMB cloglog model") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")

  ggsave(file.path(OUT, "figures", "main",
                   "fig_multiscale_future_comparison_2050.pdf"),
         p_comp, width = 18, height = 10, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(file.path(OUT, "figures", "main",
                   "fig_multiscale_future_comparison_2050.png"),
         p_comp, width = 18, height = 10, units = "cm", dpi = 600)
  log("comparison figure saved")
}

# ---- 7. Session info ------------------------------------------------------
sink(file.path(OUT, "logs", "53_sessionInfo.txt"))
print(sessionInfo())
sink()

log("=== 53_multiscale_future_prediction.R DONE ===")
