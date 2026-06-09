# ============================================================
# Scientific question / 科学问题:
#   In v1 the 100-km grid hazard model inherited province-mean climate
#   velocity, collapsing within-province climatic heterogeneity. Does
#   recomputing climate velocity directly at the native grid scale
#   from CHELSA v2.1 rasters change the magnitude and the sign of the
#   climate × effort interaction?
#   v1 网格模型直接继承省级气候速率，丢失了省内异质性。直接基于
#   CHELSA v2.1 在网格尺度重算气候速率与温降异常，结论是否稳健？
#
# Objective / 分析目标:
#   For each 50- and 100-km grid cell over mainland China, compute
#   four climate metrics at the **native** grid scale:
#     1) climate_velocity_z  (km / decade)
#     2) temp_anom_z         (recent − baseline mean tmean)
#     3) precip_anom_z       (recent − baseline mean prec)
#     4) mahalanobis_dist_z  (multivariate climate displacement)
#
# Input data / 输入数据:
#   - CHELSA v2.1 monthly rasters (tas, pr), 1981-2010 baseline + 2002-2024.
#   - 50-km / 100-km grids from utils_spatial::build_china_grid().
#
# Main workflow / 主要流程:
#   1. Source utils; read GS(2019)1822 mainland boundary.
#   2. Build (or reload) Albers grids at 50 + 100 km.
#   3. Aggregate CHELSA monthly rasters → annual mean per grid cell.
#   4. Compute climate velocity (spatial gradient ÷ temporal trend).
#   5. Compute temp / precip anomalies (recent vs baseline).
#   6. Compute Mahalanobis distance in temp × precip space.
#   7. z-score and write parquet.
#
# Expected output / 预期输出:
#   data/derived/grid_{50,100}km_climate_native.parquet
#   results/diagnostics/grid_climate_native_summary.csv
#
# Key assumptions / 关键假设:
#   - CHELSA v2.1 rasters stored under data/spatial/chelsa/ (NetCDF or GeoTIFF).
#     If absent, the script falls back to the province-level proxy and
#     warns loudly, so v1 results remain reproducible.
#   - Baseline = 1981–2010; recent = 2015–2024.
#
# Main packages / 主要包: terra, sf, exactextractr, data.table, arrow,
#   future.apply, glue.
# Output directory / 输出路径: data/derived/, results/diagnostics/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(terra)
  library(arrow)
  library(glue)
  library(future.apply)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_spatial.R"))

# ---- 0. Config -------------------------------------------------------------
CFG <- list(
  grid_sizes_km   = c(50, 100),
  baseline_years  = 1981:2010,
  recent_years    = 2015:2024,
  chelsa_dir      = path_spatial("chelsa"),        # /tas/, /pr/
  fallback_csv    = path_raw("grid_100km_climate.csv"),
  workers         = max(1L, parallel::detectCores() - 2L)
)

if (Sys.getenv("V2_WORKERS") != "") {
  CFG$workers <- as.integer(Sys.getenv("V2_WORKERS"))
}

message(glue("[28] Workers = {CFG$workers}"))
plan_ok <- requireNamespace("future", quietly = TRUE)
if (plan_ok) future::plan(future::multisession, workers = CFG$workers)

# ---- 1. Province boundary --------------------------------------------------
boundary <- tryCatch(read_gs2019_basemap("national"),
                     error = function(e) read_gs2019_basemap("province"))
boundary <- to_albers(boundary)

# ---- 2. Grids --------------------------------------------------------------
build_or_load_grid <- function(km) {
  out <- path_derived(glue("grid_{km}km_sf.gpkg"))
  if (file.exists(out)) {
    g <- sf::st_read(out, quiet = TRUE)
  } else {
    g <- build_china_grid(boundary, size_km = km, min_land_frac = 0.4)
    ensure_dir(dirname(out))
    sf::st_write(g, out, delete_dsn = TRUE, quiet = TRUE)
  }
  g
}

# ---- 3. Aggregate CHELSA to grid cell × year -------------------------------
aggregate_chelsa_year <- function(year, var, grid_sf) {
  pattern <- glue("CHELSA_{var}_.*{year}.*\\.tif$")
  files <- list.files(file.path(CFG$chelsa_dir, var), pattern = pattern,
                       full.names = TRUE, recursive = TRUE)
  if (length(files) == 0L) return(NULL)
  rast <- terra::rast(files)
  yr_mean <- terra::app(rast, fun = mean, na.rm = TRUE)
  out <- raster_to_grid(yr_mean, grid_sf, fun = "mean",
                        value_col_prefix = paste0(var, "_", year))
  out[, year := year]
  out[, var  := var]
  out
}

# ---- 4. Climate velocity = spatial gradient / temporal trend ---------------
# Implementation of Loarie et al. 2009. 时间趋势 / 空间梯度。
compute_velocity_for_grid <- function(annual_dt, grid_sf, var) {
  setDT(annual_dt)
  # temporal trend per cell. 每格的时间趋势。
  trend <- annual_dt[, .(slope = tryCatch(coef(lm(value ~ year))[2],
                                            error = function(e) NA_real_)),
                    by = grid_id]
  # spatial gradient via neighbour terra::terrain on baseline mean raster.
  baseline_dt <- annual_dt[year %in% CFG$baseline_years,
                           .(mean_val = mean(value, na.rm = TRUE)),
                           by = grid_id]
  # Build raster from cell centroids. 从网格中心构建栅格。
  centroids <- sf::st_centroid(grid_sf)
  coord_dt <- data.table::data.table(
    grid_id = grid_sf$grid_id,
    x = sf::st_coordinates(centroids)[, 1],
    y = sf::st_coordinates(centroids)[, 2]
  )
  baseline_dt <- merge(baseline_dt, coord_dt, by = "grid_id")
  template <- terra::rast(terra::ext(baseline_dt$x, baseline_dt$y),
                          resolution = unique(diff(sort(unique(baseline_dt$x))))[1])
  vals <- terra::rasterize(as.matrix(baseline_dt[, .(x, y)]),
                            template, values = baseline_dt$mean_val)
  slope <- terra::terrain(vals, v = "slope", unit = "tangent")
  slope_dt <- terra::extract(slope, sf::st_centroid(grid_sf))
  setDT(slope_dt)
  slope_dt[, grid_id := grid_sf$grid_id]
  vel <- merge(trend, slope_dt[, .(grid_id, gradient = slope)], by = "grid_id")
  vel[, velocity_km_per_decade := (slope * 10) / pmax(gradient, 1e-6) / 1000]
  vel[, var := var]
  vel
}

# ---- 5. Fallback: use province-level proxy if CHELSA absent ----------------
# HARD FAIL by default. Set V2_ALLOW_CLIMATE_FALLBACK=1 to opt in to the
# province-mirror proxy. The opt-in path writes fallback=TRUE in the
# output so downstream consumers can detect it. P0-1 is only resolved
# when this function is NOT called. 默认硬失败，避免静默回退到省级气候。
fallback_warning <- function(km) {
  allow <- nzchar(Sys.getenv("V2_ALLOW_CLIMATE_FALLBACK"))
  if (!allow) {
    stop(glue("[28] CHELSA dir missing ({CFG$chelsa_dir}). ",
              "Refusing to fall back silently to province-mean climate ",
              "because that secretly re-introduces v1's P0-1 MAUP artefact. ",
              "Either: (a) place CHELSA v2.1 monthly rasters under {CFG$chelsa_dir}/{{tas,pr}}/, ",
              "or (b) re-export this script with V2_ALLOW_CLIMATE_FALLBACK=1 ",
              "to opt in to the proxy. Aborting."), call. = FALSE)
  }
  warning(glue("[28] CHELSA dir missing ({CFG$chelsa_dir}); FALLBACK MODE ON ",
                "via V2_ALLOW_CLIMATE_FALLBACK=1. Output rows will be marked ",
                "fallback=TRUE. P0-1 NOT fixed in this run."))
  prov <- fread(path_raw("climate_metrics_province_year.csv"), encoding = "UTF-8")
  base <- fread(path_raw(glue("grid_{km}km_base.csv")), encoding = "UTF-8")
  base[, province := tolower(province)]
  prov[, province := tolower(province)]
  recent <- prov[year %in% CFG$recent_years,
                 lapply(.SD, mean, na.rm = TRUE),
                 by = province,
                 .SDcols = patterns("velocity|anom|mahalanobis")]
  out <- merge(base, recent, by = "province", all.x = TRUE)
  zify <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  for (col in grep("velocity|anom|mahalanobis", names(out), value = TRUE)) {
    out[, paste0(col, "_z") := zify(get(col))]
  }
  out[, fallback := TRUE]
  out
}

# ---- 6. Main driver --------------------------------------------------------
run_for_grid <- function(km) {
  message(glue("\n========== Grid {km} km =========="))
  grid_sf <- build_or_load_grid(km)
  out_parq <- path_derived(glue("grid_{km}km_climate_native.parquet"))

  if (!dir.exists(CFG$chelsa_dir)) {
    out <- fallback_warning(km)
    ensure_dir(dirname(out_parq))
    arrow::write_parquet(out, out_parq, compression = "snappy")
    return(invisible(out_parq))
  }

  # ---- Aggregate CHELSA per year × variable -------------------------------
  years_all <- c(CFG$baseline_years, CFG$recent_years)
  vars_ <- c("tas", "pr")
  annual <- data.table::rbindlist(future.apply::future_lapply(years_all, function(yy) {
    parts <- lapply(vars_, function(vv) {
      one <- aggregate_chelsa_year(yy, vv, grid_sf)
      if (is.null(one)) return(NULL)
      one_long <- data.table::melt(one, id.vars = c("grid_id", "year", "var"),
                                   variable.name = "month", value.name = "value")
      one_long
    })
    parts <- parts[!sapply(parts, is.null)]
    if (length(parts) == 0L) return(NULL)
    data.table::rbindlist(parts, fill = TRUE)
  }), fill = TRUE)

  if (nrow(annual) == 0L) {
    out <- fallback_warning(km)
  } else {
    # ---- Velocity for each variable ----------------------------------------
    velocity <- data.table::rbindlist(lapply(vars_, function(vv) {
      sub <- annual[var == vv]
      compute_velocity_for_grid(sub, grid_sf, vv)
    }), fill = TRUE)
    velocity_wide <- data.table::dcast(velocity, grid_id ~ var,
                                       value.var = "velocity_km_per_decade")
    setnames(velocity_wide, c("tas", "pr"),
             c("velocity_tas_km_decade", "velocity_pr_km_decade"))

    # ---- Temp / precip anomalies ------------------------------------------
    baseline <- annual[year %in% CFG$baseline_years,
                       .(baseline = mean(value, na.rm = TRUE)),
                       by = .(grid_id, var)]
    recent <- annual[year %in% CFG$recent_years,
                     .(recent = mean(value, na.rm = TRUE)),
                     by = .(grid_id, var)]
    anom <- merge(baseline, recent, by = c("grid_id", "var"))
    anom[, anom := recent - baseline]
    anom_wide <- data.table::dcast(anom, grid_id ~ var, value.var = "anom")
    setnames(anom_wide, c("tas", "pr"),
             c("temp_anom_native", "precip_anom_native"))

    # ---- Mahalanobis distance --------------------------------------------
    M <- as.matrix(anom_wide[, .(temp_anom_native, precip_anom_native)])
    mu <- colMeans(M, na.rm = TRUE)
    S  <- cov(M, use = "pairwise.complete.obs")
    mahala <- as.numeric(mahalanobis(M, mu, S))
    anom_wide[, mahalanobis_dist_native := mahala]

    # ---- Merge + z-score -------------------------------------------------
    out <- merge(velocity_wide, anom_wide, by = "grid_id")
    zify <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    for (col in setdiff(names(out), "grid_id")) {
      out[, paste0(col, "_z") := zify(get(col))]
    }
    out[, fallback := FALSE]
  }

  ensure_dir(dirname(out_parq))
  arrow::write_parquet(out, out_parq, compression = "snappy")

  diag <- data.table::data.table(
    grid_km        = km,
    n_cells        = nrow(out),
    fallback       = out$fallback[1],
    velocity_mean  = mean(out$velocity_tas_km_decade, na.rm = TRUE),
    velocity_sd    = sd(out$velocity_tas_km_decade,   na.rm = TRUE),
    temp_anom_mean = mean(out$temp_anom_native, na.rm = TRUE),
    temp_anom_sd   = sd(out$temp_anom_native,   na.rm = TRUE)
  )
  ensure_dir(path_diagnostics())
  diag_path <- path_diagnostics("grid_climate_native_summary.csv")
  existing <- if (file.exists(diag_path)) fread(diag_path, encoding = "UTF-8") else NULL
  fwrite(data.table::rbindlist(list(existing, diag), fill = TRUE,
                               use.names = TRUE), diag_path)
  message(glue("[28] Wrote {out_parq} ({nrow(out)} cells, fallback={out$fallback[1]})"))
  invisible(out_parq)
}

invisible(lapply(CFG$grid_sizes_km, run_for_grid))

# Tidy down. 收尾。
if (plan_ok) future::plan(future::sequential)
dump_session_info(path_logs("28_grid_native_climate_sessionInfo.txt"))
message("[28] done.")
