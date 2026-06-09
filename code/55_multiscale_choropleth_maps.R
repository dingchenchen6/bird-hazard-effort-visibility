# ============================================================
# Scientific question / 科学问题:
#   How do current and future hazard patterns appear across
#   province, prefecture, county, and 100km grid scales?
#   当前和未来新纪录风险在省/市/县/100km网格尺度上的
#   空间分布格局如何？
#
# Objective / 分析目标:
#   - Read prediction outputs from scripts 51/53/54
#   - Generate publication-quality choropleth maps at all 4 scales
#   - Include: current hazard, 2050 SSP585 (glmmTMB + XGBoost),
#     multi-scale comparison panels
#   - Use 中国shp with 十段线 and 国界 for proper China basemap
#   - All outputs to outputs_multiscale/
#
# Input data / 输入数据:
#   - outputs_multiscale/data/derived/risk_set_*.csv
#   - outputs_multiscale/results/forecasts/table_multiscale_future_*.csv
#   - outputs_multiscale/results/forecasts/table_multiscale_xgb_future_*.csv
#   - data/spatial/china_shp/ (省/市/县.shp + 十段线.shp)
#   - data/spatial/basemap_GS2019_1822/ (等积投影版)
#
# Expected output / 预期输出:
#   outputs_multiscale/figures/main/
#     fig_choropleth_{province,prefecture,county,grid_100km}_current.{pdf,png}
#     fig_choropleth_{province,prefecture,county,grid_100km}_2050_ssp585_glmmTMB.{pdf,png}
#     fig_choropleth_{province,prefecture,county,grid_100km}_2050_ssp585_xgboost.{pdf,png}
#     fig_choropleth_multi_scale_comparison.{pdf,png}
#     fig_choropleth_glmmTMB_vs_xgboost_2050.{pdf,png}
#
# Main packages / 主要包: data.table, sf, ggplot2, patchwork.
# Output directory / 输出路径: outputs_multiscale/
# ============================================================

# ---- 0. CLI args ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
arg <- list()
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--base-dir")        { arg$base_dir   <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--output-dir") { arg$output_dir <- args[i + 1L]; i <- i + 2L }
  else { i <- i + 1L }
}
`%||%` <- function(a, b) if (is.null(a)) b else a

BASE_DIR   <- arg$base_dir   %||% "/Users/dingchenchen/Documents/New records/bird-new-distribution-records/tasks/bird_hazard_model_effort_upgrade_v2"
OUTPUT_DIR <- arg$output_dir %||% "outputs_multiscale"

OUT <- file.path(BASE_DIR, OUTPUT_DIR)
ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ens(file.path(OUT, "figures", "main"))
ens(file.path(OUT, "logs"))

log_file <- file.path(OUT, "logs", "55_multiscale_choropleth.log")
log <- function(...) {
  msg <- sprintf("[55 %s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(msg); cat(msg, file = log_file, append = TRUE)
}
log("=== 55_multiscale_choropleth_maps.R START ===")

# ---- 1. Packages ----------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
})
sf::sf_use_s2(FALSE)
options(warn = 1)

# ---- 2. Load shapefiles ---------------------------------------------------
log("=== Step 2: Loading shapefiles ===")

china_crs <- "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +ellps=GRS80 +units=m +no_defs"

# Search for 中国shp
shp_dirs <- c(
  file.path(BASE_DIR, "data", "spatial", "china_shp"),
  file.path(dirname(BASE_DIR), "bird_hazard_model_effort_upgrade", "中国shp"),
  file.path(BASE_DIR, "中国shp"))
china_shp_dir <- shp_dirs[dir.exists(shp_dirs)][1]

# Also search for GS2019 (等积投影)
gs_dirs <- c(
  file.path(BASE_DIR, "data", "spatial", "basemap_GS2019_1822"),
  file.path(dirname(BASE_DIR), "bird_hazard_model_effort_upgrade",
            "2019中国地图-审图号GS(2019)1822号"))
gs_dir <- gs_dirs[dir.exists(gs_dirs)][1]

log("china_shp: ", if (!is.null(china_shp_dir)) china_shp_dir else "NOT FOUND")
log("GS2019:    ", if (!is.null(gs_dir)) gs_dir else "NOT FOUND")

# Province CN→EN mapping
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

# Load helper: find shp excluding 境界线
find_shp_poly <- function(dir, prefix) {
  if (is.null(dir)) return(NULL)
  hits <- list.files(dir, pattern = paste0("^", prefix, ".*\\.shp$"),
                     full.names = TRUE)
  hits <- hits[!grepl("境界线", hits)]
  if (length(hits) == 0L) return(NULL)
  eq <- hits[grepl("等积投影", hits)]
  if (length(eq) > 0L) return(eq[1])
  poly <- hits[!grepl("线", hits)]
  if (length(poly) > 0L) return(poly[1])
  hits[1]
}

assign_prov_en <- function(sf_obj) {
  cn_col <- names(sf_obj)[vapply(sf_obj, function(x) {
    is.character(x) && any(grepl("[\u4e00-\u9fff]", x))
  }, logical(1))][1]
  if (!is.na(cn_col)) {
    sf_obj$province_en <- unname(PROV_CN_EN[as.character(sf_obj[[cn_col]])])
    sf_obj$province_en[is.na(sf_obj$province_en)] <-
      as.character(sf_obj[[cn_col]])[is.na(sf_obj$province_en)]
  }
  sf_obj
}

# Load shapefiles (prefer china_shp for cleaner names, GS2019 for 等积投影)
prov_sf <- NULL
pref_sf <- NULL
cnty_sf <- NULL
ten_dash_sf <- NULL

# Province: prefer china_shp (has 省名 column for easy mapping)
sp <- find_shp_poly(china_shp_dir, "省")
if (is.null(sp)) sp <- find_shp_poly(gs_dir, "省")
if (!is.null(sp)) {
  prov_sf <- st_read(sp, quiet = TRUE) |> st_make_valid() |> st_transform(china_crs)
  prov_sf <- assign_prov_en(prov_sf)
  log("province: ", nrow(prov_sf), " features")
}

# Prefecture
sp <- find_shp_poly(china_shp_dir, "市")
if (is.null(sp)) sp <- find_shp_poly(gs_dir, "市")
if (!is.null(sp)) {
  pref_sf <- st_read(sp, quiet = TRUE) |> st_make_valid() |> st_transform(china_crs)
  pref_sf <- assign_prov_en(pref_sf)
  log("prefecture: ", nrow(pref_sf), " features")
}

# County
sp <- find_shp_poly(china_shp_dir, "县")
if (is.null(sp)) sp <- find_shp_poly(gs_dir, "县")
if (!is.null(sp)) {
  cnty_sf <- st_read(sp, quiet = TRUE) |> st_make_valid() |> st_transform(china_crs)
  cnty_sf <- assign_prov_en(cnty_sf)
  log("county: ", nrow(cnty_sf), " features")
}

# 十段线
td <- file.path(china_shp_dir, "十段线.shp")
if (!file.exists(td)) td <- find_shp_poly(gs_dir, "九段线")
if (file.exists(td)) {
  ten_dash_sf <- st_read(td, quiet = TRUE) |> st_make_valid() |> st_transform(china_crs)
  log("十段线 loaded: ", nrow(ten_dash_sf), " features")
}

# ---- 3. Load prediction data -----------------------------------------------
log("=== Step 3: Loading prediction data ===")

# Current hazard (from risk sets)
risk_sets <- list()
for (sc in c("province", "prefecture", "county", "grid_100km")) {
  f <- file.path(OUT, "data", "derived", paste0("risk_set_", sc, ".csv"))
  if (file.exists(f)) risk_sets[[sc]] <- fread(f, encoding = "UTF-8")
}

# glmmTMB future
glmm_fut <- list()
for (sc in c("province", "prefecture", "county", "grid_100km")) {
  f <- file.path(OUT, "results", "forecasts",
                 paste0("table_multiscale_future_", sc, ".csv"))
  if (file.exists(f)) glmm_fut[[sc]] <- fread(f, encoding = "UTF-8")
}

# XGBoost future
xgb_fut <- list()
for (sc in c("province", "prefecture", "county", "grid_100km")) {
  f <- file.path(OUT, "results", "forecasts",
                 paste0("table_multiscale_xgb_future_", sc, ".csv"))
  if (file.exists(f)) xgb_fut[[sc]] <- fread(f, encoding = "UTF-8")
}

# ---- 4. Map theme and helper -----------------------------------------------
theme_map <- theme_bw(base_size = 9) +
  theme(axis.title = element_blank(), axis.text = element_text(size = 6),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey70", linewidth = 0.3),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"), legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 7), legend.title = element_text(size = 8))

# Common fill scale: blue-white-red diverging
scale_hazard <- scale_fill_gradientn(
  colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
  na.value = "grey90", name = "Hazard")

# Add 十段线 layer if available
add_ten_dash <- function(p) {
  if (!is.null(ten_dash_sf)) p + geom_sf(data = ten_dash_sf, fill = NA,
    colour = "grey50", linewidth = 0.3) else p
}

# ---- 5. Current hazard maps -----------------------------------------------
log("=== Step 5: Current hazard maps ===")

# Province current hazard
if ("province" %in% names(risk_sets) && !is.null(prov_sf)) {
  dt <- risk_sets[["province"]]
  current_prov <- dt[, .(hazard = mean(event)), by = province]
  setnames(current_prov, "province", "province_en")
  shp_merged <- merge(prov_sf, current_prov, by = "province_en", all.x = TRUE)

  p <- ggplot(shp_merged) +
    geom_sf(aes(fill = hazard), colour = "grey60", linewidth = 0.1) +
    scale_hazard + labs(title = "Province: current hazard rate") + theme_map
  p <- add_ten_dash(p)

  ggsave(file.path(OUT, "figures", "main", "fig_choropleth_province_current.pdf"),
         p, width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf)
  ggsave(file.path(OUT, "figures", "main", "fig_choropleth_province_current.png"),
         p, width = 12, height = 8, units = "cm", dpi = 600)
  log("province current map saved")
}

# Prefecture current hazard
if ("prefecture" %in% names(risk_sets) && !is.null(pref_sf)) {
  dt <- risk_sets[["prefecture"]]
  if ("unit_id" %in% names(dt)) {
    current_pref <- dt[, .(hazard = mean(event)), by = unit_id]
    # Merge by row index (unit_id = prefecture_N)
    pref_sf$unit_id <- paste0("prefecture_", seq_len(nrow(pref_sf)))
    shp_merged <- merge(pref_sf, current_pref, by = "unit_id", all.x = TRUE)

    p <- ggplot(shp_merged) +
      geom_sf(aes(fill = hazard), colour = NA) +
      scale_hazard + labs(title = "Prefecture: current hazard rate") + theme_map
    p <- add_ten_dash(p)

    ggsave(file.path(OUT, "figures", "main", "fig_choropleth_prefecture_current.pdf"),
           p, width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main", "fig_choropleth_prefecture_current.png"),
           p, width = 12, height = 8, units = "cm", dpi = 600)
    log("prefecture current map saved")
  }
}

# County current hazard
if ("county" %in% names(risk_sets) && !is.null(cnty_sf)) {
  dt <- risk_sets[["county"]]
  if ("unit_id" %in% names(dt)) {
    current_cnty <- dt[, .(hazard = mean(event)), by = unit_id]
    cnty_sf$unit_id <- paste0("county_", seq_len(nrow(cnty_sf)))
    shp_merged <- merge(cnty_sf, current_cnty, by = "unit_id", all.x = TRUE)

    p <- ggplot(shp_merged) +
      geom_sf(aes(fill = hazard), colour = NA) +
      scale_hazard + labs(title = "County: current hazard rate") + theme_map
    p <- add_ten_dash(p)

    ggsave(file.path(OUT, "figures", "main", "fig_choropleth_county_current.pdf"),
           p, width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main", "fig_choropleth_county_current.png"),
           p, width = 12, height = 8, units = "cm", dpi = 600)
    log("county current map saved")
  }
}

# ---- 6. Future glmmTMB maps (2050 SSP585) --------------------------------
log("=== Step 6: Future glmmTMB maps ===")

# Helper to make future choropleth at any scale
make_future_choropleth <- function(fut_dt, shp, unit_col, scale_name,
                                    method_label = "glmmTMB") {
  if (is.null(shp) || is.null(fut_dt) || nrow(fut_dt) == 0L) return(NULL)
  pred_summary <- fut_dt[year == 2050 & climate_scenario == "ssp585" &
                          effort_scenario == "baseline",
    .(hazard = mean(hazard, na.rm = TRUE)),
    by = unit_col]
  if (nrow(pred_summary) == 0L) return(NULL)

  # For province scale, join on province_en
  if (unit_col == "province" && "province_en" %in% names(shp)) {
    setnames(pred_summary, "province", "province_en")
    shp_merged <- merge(shp, pred_summary, by = "province_en", all.x = TRUE)
  } else {
    # For prefecture/county: assign unit_id to shp rows by index
    shp_tmp <- shp
    shp_tmp$unit_id <- paste0(scale_name, "_", seq_len(nrow(shp_tmp)))
    shp_merged <- merge(shp_tmp, pred_summary, by = "unit_id", all.x = TRUE)
  }

  fill_name <- paste0("Hazard\n(", method_label, ")")
  p <- ggplot(shp_merged) +
    geom_sf(aes(fill = hazard), colour = if (nrow(shp) > 500L) NA else "grey60",
            linewidth = if (nrow(shp) > 500L) 0 else 0.1) +
    scale_fill_gradientn(colours = c("#2166AC", "#67A9CF", "#F7F7F7",
                                      "#EF8A62", "#B2182B"),
                         na.value = "grey90", name = fill_name) +
    labs(title = paste0(toupper(substr(scale_name, 1, 1)), substr(scale_name, 2, nchar(scale_name)),
                        ": 2050 SSP585 hazard (", method_label, ")")) +
    theme_map
  add_ten_dash(p)
}

# Province
if ("province" %in% names(glmm_fut) && !is.null(prov_sf)) {
  p <- make_future_choropleth(glmm_fut[["province"]], prov_sf, "province", "province")
  if (!is.null(p)) {
    ggsave(file.path(OUT, "figures", "main", "fig_choropleth_province_2050_ssp585_glmmTMB.pdf"),
           p, width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main", "fig_choropleth_province_2050_ssp585_glmmTMB.png"),
           p, width = 12, height = 8, units = "cm", dpi = 600)
    log("province glmmTMB future map saved")
  }
}

# Prefecture
if ("prefecture" %in% names(glmm_fut) && !is.null(pref_sf)) {
  p <- make_future_choropleth(glmm_fut[["prefecture"]], pref_sf, "unit_id", "prefecture")
  if (!is.null(p)) {
    ggsave(file.path(OUT, "figures", "main", "fig_choropleth_prefecture_2050_ssp585_glmmTMB.pdf"),
           p, width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main", "fig_choropleth_prefecture_2050_ssp585_glmmTMB.png"),
           p, width = 12, height = 8, units = "cm", dpi = 600)
    log("prefecture glmmTMB future map saved")
  }
}

# County
if ("county" %in% names(glmm_fut) && !is.null(cnty_sf)) {
  p <- make_future_choropleth(glmm_fut[["county"]], cnty_sf, "unit_id", "county")
  if (!is.null(p)) {
    ggsave(file.path(OUT, "figures", "main", "fig_choropleth_county_2050_ssp585_glmmTMB.pdf"),
           p, width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main", "fig_choropleth_county_2050_ssp585_glmmTMB.png"),
           p, width = 12, height = 8, units = "cm", dpi = 600)
    log("county glmmTMB future map saved")
  }
}

# ---- 7. Future XGBoost maps (2050 SSP585) --------------------------------
log("=== Step 7: Future XGBoost maps ===")

# Province
if ("province" %in% names(xgb_fut) && !is.null(prov_sf)) {
  p <- make_future_choropleth(xgb_fut[["province"]], prov_sf, "province", "province", "XGBoost")
  if (!is.null(p)) {
    ggsave(file.path(OUT, "figures", "main",
                     "fig_choropleth_province_2050_ssp585_xgboost.pdf"),
           p, width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main",
                     "fig_choropleth_province_2050_ssp585_xgboost.png"),
           p, width = 12, height = 8, units = "cm", dpi = 600)
    log("province XGBoost future map saved")
  }
}

# Prefecture
if ("prefecture" %in% names(xgb_fut) && !is.null(pref_sf)) {
  p <- make_future_choropleth(xgb_fut[["prefecture"]], pref_sf, "unit_id", "prefecture", "XGBoost")
  if (!is.null(p)) {
    ggsave(file.path(OUT, "figures", "main",
                     "fig_choropleth_prefecture_2050_ssp585_xgboost.pdf"),
           p, width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main",
                     "fig_choropleth_prefecture_2050_ssp585_xgboost.png"),
           p, width = 12, height = 8, units = "cm", dpi = 600)
    log("prefecture XGBoost future map saved")
  }
}

# County
if ("county" %in% names(xgb_fut) && !is.null(cnty_sf)) {
  p <- make_future_choropleth(xgb_fut[["county"]], cnty_sf, "unit_id", "county", "XGBoost")
  if (!is.null(p)) {
    ggsave(file.path(OUT, "figures", "main",
                     "fig_choropleth_county_2050_ssp585_xgboost.pdf"),
           p, width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf)
    ggsave(file.path(OUT, "figures", "main",
                     "fig_choropleth_county_2050_ssp585_xgboost.png"),
           p, width = 12, height = 8, units = "cm", dpi = 600)
    log("county XGBoost future map saved")
  }
}

# ---- 8. glmmTMB vs XGBoost comparison (province) -------------------------
log("=== Step 8: glmmTMB vs XGBoost comparison ===")

if ("province" %in% names(glmm_fut) && "province" %in% names(xgb_fut) &&
    !is.null(prov_sf)) {
  g_prov <- glmm_fut[["province"]][year == 2050 & climate_scenario == "ssp585" &
    effort_scenario == "baseline", .(province, hazard_glmm = hazard)]
  x_prov <- xgb_fut[["province"]][year == 2050 & climate_scenario == "ssp585" &
    effort_scenario == "baseline", .(province, hazard_xgb = hazard)]

  # Aggregate by province
  g_agg <- g_prov[, .(hazard_glmm = mean(hazard_glmm, na.rm = TRUE)), by = province]
  x_agg <- x_prov[, .(hazard_xgb = mean(hazard_xgb, na.rm = TRUE)), by = province]
  comp <- merge(g_agg, x_agg, by = "province")
  comp[, delta := hazard_xgb - hazard_glmm]
  setnames(comp, "province", "province_en")
  shp_merged <- merge(prov_sf, comp, by = "province_en", all.x = TRUE)

  p <- ggplot(shp_merged) +
    geom_sf(aes(fill = delta), colour = "grey60", linewidth = 0.1) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, na.value = "grey90",
                         name = "XGBoost - glmmTMB") +
    labs(title = "Province: XGBoost vs glmmTMB hazard difference (2050 SSP585)") +
    theme_map
  p <- add_ten_dash(p)

  ggsave(file.path(OUT, "figures", "main",
                   "fig_choropleth_glmmTMB_vs_xgboost_2050.pdf"),
         p, width = 12, height = 8, units = "cm", device = grDevices::cairo_pdf)
  ggsave(file.path(OUT, "figures", "main",
                   "fig_choropleth_glmmTMB_vs_xgboost_2050.png"),
         p, width = 12, height = 8, units = "cm", dpi = 600)
  log("glmmTMB vs XGBoost comparison map saved")
}

# ---- 9. Session info ------------------------------------------------------
sink(file.path(OUT, "logs", "55_sessionInfo.txt"))
print(sessionInfo())
sink()

log("=== 55_multiscale_choropleth_maps.R DONE ===")
