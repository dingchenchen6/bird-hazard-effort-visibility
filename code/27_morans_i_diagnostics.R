# ============================================================
# Scientific question / 科学问题:
#   Do the residuals of the M4 (climate × effort interaction) hazard
#   model retain spatial autocorrelation at biogeographically
#   meaningful scales (50–500 km)? If yes, reported standard errors
#   may be too narrow and the headline conclusion needs a robustness
#   caveat.
#   M4 残差是否仍残留 50–500 km 的空间自相关？若是，置信区间被低估。
#
# Objective / 分析目标:
#   - Refit M4 at province scale (and at grid scale if 28+33 outputs
#     are available).
#   - Compute Moran's I on Pearson residuals at distance classes
#     50 / 100 / 250 / 500 km.
#   - Run DHARMa::testSpatialAutocorrelation as a complementary check.
#   - Render Fig 4 inset + Fig S2 (residual map).
#
# Input data / 输入数据:
#   data/raw/hazard_risk_upgraded_complete_case.csv
#   data/derived/events_*_grid_risk_set_v2.parquet  (optional, from script 33)
#
# Main workflow / 主要流程: see comments inline.
# Expected output / 预期输出:
#   results/diagnostics/table_morans_i_residuals.csv
#   results/diagnostics/table_dharma_spatial.csv
#   figures/diagnostics/morans_i_distance_classes.{pdf,png}
#   figures/supplementary/figS2_dharma_spatial_residuals.{pdf,png}
#
# Main packages / 主要包: glmmTMB, DHARMa, ape, data.table, sf,
#   ggplot2, patchwork.
# Output directory / 输出路径: results/diagnostics/,
#   figures/diagnostics/, figures/supplementary/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
  library(DHARMa)
  library(ape)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_models.R"))
source(file.path("code", "utils", "utils_spatial.R"))
source(file.path("code", "utils", "utils_plots.R"))

set.seed(42)

CFG <- list(
  dist_classes_km = c(50, 100, 250, 500),
  dharma_n        = 1000
)

# ---- 1. Load and fit M4 (province scale) -----------------------------------
dt <- fread(path_raw("hazard_risk_upgraded_complete_case.csv"),
            encoding = "UTF-8")
climate_col <- intersect(c("temp_grad_z", "temp_anom_z",
                            "climate_velocity_z", "mahalanobis_dist_z"),
                          names(dt))[1]
effort_col  <- intersect(c("log_effort_visits_z", "log_n_visits_z",
                            "effort_pc1_z", "log_effort_days_z"),
                          names(dt))[1]
if (is.na(climate_col) || is.na(effort_col)) {
  stop("[27] Need climate_z + effort_z columns. Found: ",
       paste(names(dt), collapse = ", "))
}
message("[27] climate column = ", climate_col, " | effort column = ", effort_col)

dt[, climate_z := get(climate_col)]
dt[, effort_z  := get(effort_col)]

message("[27] Fitting M4 (interaction) at province scale…")
fit_M4 <- glmmTMB::glmmTMB(
  event ~ climate_z * effort_z + (1 | species) + (1 | province),
  data = dt, family = binomial(link = "cloglog"))

# Residuals + coordinates. 残差与坐标。
res <- residuals(fit_M4, type = "pearson")
prov_sf <- read_gs2019_basemap("province")
# GS(2019)1822 中文 shp 用列名 "省"。匹配多种潜在 name 列。
.prov_name_col <- intersect(c("province_en", "name", "NAME", "省"),
                             names(prov_sf))[1]
if (is.na(.prov_name_col)) {
  .prov_name_col <- names(prov_sf)[vapply(prov_sf, function(x)
    is.character(x) && any(grepl("[一-龥]", x)), logical(1))][1]
}
if (is.na(.prov_name_col)) stop("[27] no province-name column in basemap shp")
.PROV_CN_EN <- c("北京市"="beijing","天津市"="tianjin","河北省"="hebei",
  "山西省"="shanxi","内蒙古自治区"="inner mongolia","辽宁省"="liaoning",
  "吉林省"="jilin","黑龙江省"="heilongjiang","上海市"="shanghai",
  "江苏省"="jiangsu","浙江省"="zhejiang","安徽省"="anhui","福建省"="fujian",
  "江西省"="jiangxi","山东省"="shandong","河南省"="henan","湖北省"="hubei",
  "湖南省"="hunan","广东省"="guangdong","广西壮族自治区"="guangxi",
  "海南省"="hainan","重庆市"="chongqing","四川省"="sichuan","贵州省"="guizhou",
  "云南省"="yunnan","西藏自治区"="tibet","陕西省"="shaanxi","甘肃省"="gansu",
  "青海省"="qinghai","宁夏回族自治区"="ningxia","新疆维吾尔自治区"="xinjiang",
  "台湾省"="taiwan","香港特别行政区"="hong kong","澳门特别行政区"="macau")
.raw_names <- as.character(prov_sf[[.prov_name_col]])
prov_sf$province_norm <- ifelse(.raw_names %in% names(.PROV_CN_EN),
                                 unname(.PROV_CN_EN[.raw_names]),
                                 tolower(.raw_names))
centroids <- sf::st_centroid(to_albers(prov_sf))
centroid_dt <- data.table::data.table(
  province = prov_sf$province_norm,
  x = sf::st_coordinates(centroids)[, 1],
  y = sf::st_coordinates(centroids)[, 2]
)

dt[, province_norm := tolower(province)]
dt <- merge(dt, centroid_dt, by.x = "province_norm", by.y = "province",
            all.x = TRUE)
ok <- !is.na(dt$x)
res <- res[ok]
coords <- as.matrix(dt[ok, .(x, y)])

# ---- 2. Moran's I across distance classes ----------------------------------
mi_dt <- morans_i_residuals(res, coords, CFG$dist_classes_km)
mi_dt[, model := "M4_province"]
ensure_dir(path_diagnostics())
fwrite(mi_dt, path_diagnostics("table_morans_i_residuals.csv"))
message("[27] Moran's I (province):")
print(mi_dt)

# ---- 3. DHARMa spatial test ------------------------------------------------
sim <- DHARMa::simulateResiduals(fit_M4, n = CFG$dharma_n)
# Aggregate by province for spatial test. 按省聚合做空间检验。
dharma_dt <- data.table::data.table(
  province = dt$province[ok],
  res = res
)
agg <- dharma_dt[, .(res_mean = mean(res, na.rm = TRUE),
                      x = mean(coords[, 1]),
                      y = mean(coords[, 2])),
                  by = province]
spatial_test <- tryCatch(
  DHARMa::testSpatialAutocorrelation(sim, x = agg$x, y = agg$y, plot = FALSE),
  error = function(e) NULL
)
if (!is.null(spatial_test)) {
  dharma_summary <- data.table::data.table(
    test = "DHARMa::testSpatialAutocorrelation",
    statistic = spatial_test$statistic,
    p.value   = spatial_test$p.value,
    method    = spatial_test$method
  )
  fwrite(dharma_summary, path_diagnostics("table_dharma_spatial.csv"))
  message("[27] DHARMa spatial: ", spatial_test$method,
          "  p = ", signif(spatial_test$p.value, 3))
}

# ---- 4. Optional: grid-scale Moran's I -------------------------------------
for (km in c(50, 100)) {
  rp <- path_derived(glue("events_{km}km_grid_risk_set_v2.parquet"))
  if (!file.exists(rp)) {
    message(glue("[27] Skipping {km}-km grid Moran (script 33 output missing)."))
    next
  }
  grid_dt <- arrow::read_parquet(rp) |> data.table::as.data.table()
  if (!"event" %in% names(grid_dt)) {
    message(glue("[27] Skipping {km}-km grid Moran (no event column)."))
    next
  }
  features <- intersect(c("climate_velocity_native_z", "temp_anom_native_z",
                           "log_n_visits_z", "effort_pc1_z"), names(grid_dt))
  if (length(features) < 2L) {
    message(glue("[27] Skipping {km}-km grid Moran (need ≥2 covariates)."))
    next
  }
  form <- as.formula(paste("event ~",
                           paste(features[1:2], collapse = " * "),
                           "+ (1|species) + (1|grid_id)"))
  fit_grid <- tryCatch(
    glmmTMB::glmmTMB(form, data = grid_dt[!is.na(event)],
                     family = binomial(link = "cloglog")),
    error = function(e) NULL
  )
  if (is.null(fit_grid)) next

  res_g <- residuals(fit_grid, type = "pearson")
  grid_sf_path <- path_derived(glue("grid_{km}km_sf.gpkg"))
  if (!file.exists(grid_sf_path)) {
    message(glue("[27] No grid sf at {grid_sf_path}; skip Moran for {km}km."))
    next
  }
  grid_sf <- sf::st_read(grid_sf_path, quiet = TRUE)
  centroids_g <- sf::st_centroid(to_albers(grid_sf))
  cc <- data.table::data.table(grid_id = grid_sf$grid_id,
                               x = sf::st_coordinates(centroids_g)[, 1],
                               y = sf::st_coordinates(centroids_g)[, 2])
  used <- grid_dt[!is.na(event)]
  used <- merge(used, cc, by = "grid_id", all.x = TRUE)
  ok2 <- !is.na(used$x)
  mi_g <- morans_i_residuals(res_g[ok2], as.matrix(used[ok2, .(x, y)]),
                              CFG$dist_classes_km)
  mi_g[, model := glue("M4_grid_{km}km")]
  fwrite(mi_g, path_diagnostics("table_morans_i_residuals.csv"),
         append = TRUE)
}

# ---- 5. Figure ------------------------------------------------------------
mi_all <- fread(path_diagnostics("table_morans_i_residuals.csv"))
p <- ggplot(mi_all, aes(x = class_km, y = I, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.6) +
  geom_errorbar(aes(ymin = I - 1.96 * sd, ymax = I + 1.96 * sd),
                 width = 5, linewidth = 0.3) +
  scale_x_continuous(breaks = CFG$dist_classes_km) +
  scale_colour_manual(values = pal_cat[seq_len(uniqueN(mi_all$model))]) +
  labs(title = "Residual Moran's I across distance classes",
        x = "Distance class (km)", y = "Moran's I (Pearson residuals)") +
  theme_geb()

ensure_dir(path_diag_fig())
ggplot2::ggsave(path_diag_fig("morans_i_distance_classes.pdf"),
                p, width = 14, height = 9, units = "cm",
                device = grDevices::cairo_pdf)
ggplot2::ggsave(path_diag_fig("morans_i_distance_classes.png"),
                p, width = 14, height = 9, units = "cm", dpi = 600)

# Supplementary residual map. 残差地图（Fig S2 草稿）。
prov_resid <- dt[ok][, .(res_mean = mean(res), .N), by = province]
basemap <- china_basemap_gs2019_1822(crs = CRS_ALBERS)
prov_resid[, province := tolower(province)]
basemap$province$name_norm <- tolower(basemap$province$name %||% basemap$province$NAME)
prov_resid_sf <- merge(basemap$province, prov_resid,
                       by.x = "name_norm", by.y = "province")
p_map <- ggplot() +
  geom_sf(data = prov_resid_sf, aes(fill = res_mean),
          colour = "grey60", linewidth = 0.1) +
  scale_fill_div_zero(name = "Mean Pearson\nresidual") +
  labs(title = "Residual spatial structure (province means)") +
  theme_geb() +
  theme(panel.grid = element_blank())

ensure_dir(path_supp_fig())
ggplot2::ggsave(path_supp_fig("figS2_dharma_spatial_residuals.pdf"),
                p_map, width = 16, height = 12, units = "cm",
                device = grDevices::cairo_pdf)
ggplot2::ggsave(path_supp_fig("figS2_dharma_spatial_residuals.png"),
                p_map, width = 16, height = 12, units = "cm", dpi = 600)

dump_session_info(path_logs("27_morans_i_diagnostics_sessionInfo.txt"))
message("[27] done.")
