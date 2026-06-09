# ============================================================
# Scientific question / 科学问题:
#   Provide validated spatial primitives (Albers reprojection,
#   blockCV spatial folds, Moran's I residual diagnostics, raster→grid
#   aggregation) so spatial reasoning is consistent across v2 scripts.
#   提供统一的空间原语，让所有 v2 脚本共享同一空间逻辑。
#
# Objective / 分析目标:
#   - Build mainland-China grids in EPSG:4524 (Albers equal-area).
#   - Generate spatial-block CV folds via blockCV::cv_spatial.
#   - Compute Moran's I on Pearson residuals across distance classes.
#   - Aggregate rasters (CHELSA monthly tmean / prec) to grids.
#
# Input data / 输入数据:
#   data/spatial/basemap_GS2019_1822/省（等积投影）.shp etc.
#   Any sf / SpatVector / SpatRaster passed by callers.
#
# Main workflow / 主要流程:
#   1. CRS helpers: WGS84 (4326) ↔ Albers (4524).
#   2. Build grid_sf(size_km).
#   3. blockCV folds.
#   4. morans_i_residuals(): Moran's I per distance class.
#   5. raster_to_grid(): area-weighted aggregation.
#
# Expected output / 预期输出: NA — library functions.
# Key assumptions / 关键假设:
#   - sf >= 1.0, terra >= 1.7, blockCV >= 3.1, exactextractr >= 0.10.
# Main packages / 主要包: sf, terra, blockCV, exactextractr, ape,
#   spdep (fallback), data.table.
# Output directory / 输出路径: NA.
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(data.table)
})

# CRS constants. 坐标系常量。
CRS_WGS84  <- "EPSG:4326"
CRS_ALBERS <- "EPSG:4524"   # CGCS2000 / Albers Equal-Area, China

# ---- 1. CRS helpers ---------------------------------------------------------
to_albers <- function(x) sf::st_transform(x, CRS_ALBERS)
to_wgs84  <- function(x) sf::st_transform(x, CRS_WGS84)

# ---- 2. Build a regular China grid ----------------------------------------
# size_km — nominal cell size; only cells overlapping mainland China are kept.
# 仅保留与大陆相交的网格，省级与海洋格子被剔除。
build_china_grid <- function(boundary_sf, size_km = 100, min_land_frac = 0.5) {
  stopifnot(inherits(boundary_sf, "sf"))
  b <- to_albers(boundary_sf) |> sf::st_union()
  cell_m <- size_km * 1000
  bbox <- sf::st_bbox(b)
  grid <- sf::st_make_grid(b, cellsize = cell_m, square = TRUE,
                           offset = c(floor(bbox$xmin / cell_m) * cell_m,
                                      floor(bbox$ymin / cell_m) * cell_m))
  grid_sf <- sf::st_sf(grid_id = paste0("g", size_km, "_", seq_along(grid)),
                       geometry = grid)
  # area-weighted land fraction. 估计陆地占比。
  inter <- sf::st_intersection(grid_sf, b)
  area_inter <- as.numeric(sf::st_area(inter))
  inter_df <- data.table::as.data.table(sf::st_drop_geometry(inter))
  inter_df[, area_inter := area_inter]
  inter_df[, area_cell := cell_m * cell_m]
  inter_df[, land_frac := pmin(1, area_inter / area_cell)]
  keep <- inter_df[land_frac >= min_land_frac, grid_id]
  grid_sf <- grid_sf[grid_sf$grid_id %in% keep, ]
  grid_sf$cell_km <- size_km
  grid_sf$land_frac <- inter_df[match(grid_sf$grid_id, grid_id), land_frac]
  grid_sf
}

# ---- 3. Spatial block CV (blockCV::cv_spatial) -----------------------------
make_spatial_blocks <- function(points_sf, block_km = 250, k = 5,
                                response_col = NULL, seed = 42,
                                iteration = 50) {
  stopifnot(inherits(points_sf, "sf"))
  points_alb <- to_albers(points_sf)
  set.seed(seed)
  if (!requireNamespace("blockCV", quietly = TRUE)) {
    stop("Install blockCV (>= 3.1) before calling make_spatial_blocks().")
  }
  cv <- blockCV::cv_spatial(
    x          = points_alb,
    column     = response_col,
    size       = block_km * 1000,
    k          = k,
    selection  = "random",
    iteration  = iteration,
    plot       = FALSE,
    progress   = FALSE
  )
  list(
    folds  = cv$folds_list,
    blocks = cv$blocks,
    summary = data.table::data.table(
      fold  = seq_len(k),
      train = vapply(cv$folds_list, function(f) length(f$train), integer(1)),
      test  = vapply(cv$folds_list, function(f) length(f$test),  integer(1))
    )
  )
}

# ---- 4. Moran's I on Pearson residuals ------------------------------------
# Uses ape::Moran.I on inverse-distance weights binned by distance class.
# 距离分级的 Moran's I。
morans_i_residuals <- function(residuals, coords,
                                dist_classes_km = c(50, 100, 250, 500)) {
  stopifnot(length(residuals) == nrow(coords))
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Install ape before calling morans_i_residuals().")
  }
  coords <- as.matrix(coords)
  d <- as.matrix(stats::dist(coords))
  out <- data.table::rbindlist(lapply(dist_classes_km, function(dk) {
    thresh <- dk * 1000
    w <- 1 / d
    w[d > thresh | d == 0] <- 0
    if (all(w == 0)) {
      return(data.table::data.table(class_km = dk, I = NA_real_,
                                    expected = NA_real_, sd = NA_real_,
                                    p.value = NA_real_, n_pairs = 0L))
    }
    rs <- rowSums(w)
    w[rs > 0, ] <- w[rs > 0, ] / rs[rs > 0]
    mi <- ape::Moran.I(residuals, w)
    data.table::data.table(class_km = dk,
                           I = mi$observed,
                           expected = mi$expected,
                           sd = mi$sd,
                           p.value = mi$p.value,
                           n_pairs = sum(w > 0))
  }))
  out[]
}

# ---- 5. Raster → grid aggregation ------------------------------------------
# Area-weighted mean of `rast` values within each grid_sf cell.
# 面积加权汇总栅格到网格。
raster_to_grid <- function(rast, grid_sf, fun = "mean",
                            value_col_prefix = "v") {
  if (!requireNamespace("exactextractr", quietly = TRUE)) {
    stop("Install exactextractr for raster_to_grid().")
  }
  if (terra::nlyr(rast) == 0L) stop("Empty raster passed to raster_to_grid().")
  # Reproject grid to raster CRS for extract. 网格投影到栅格 CRS。
  rast_crs <- terra::crs(rast, describe = TRUE)$code
  if (!is.na(rast_crs) && nzchar(rast_crs)) {
    target <- paste0("EPSG:", rast_crs)
    grid_for_extract <- sf::st_transform(grid_sf, target)
  } else {
    grid_for_extract <- grid_sf
  }
  ex <- exactextractr::exact_extract(rast, grid_for_extract, fun = fun,
                                     progress = FALSE)
  if (is.matrix(ex) || is.data.frame(ex)) {
    ex <- as.data.frame(ex, stringsAsFactors = FALSE)
    colnames(ex) <- paste0(value_col_prefix, "_", seq_len(ncol(ex)))
  } else {
    ex <- data.frame(v_1 = ex)
  }
  out <- cbind(sf::st_drop_geometry(grid_sf), ex)
  data.table::setDT(out)
  out
}

# ---- 6. Read GS(2019)1822 basemap ------------------------------------------
read_gs2019_basemap <- function(layer = c("province", "prefecture", "county",
                                          "national", "ninedash")) {
  layer <- match.arg(layer)
  dir <- file.path("data", "spatial", "basemap_GS2019_1822")
  # Files in v1 follow Chinese names; map them. 文件名映射。
  candidates <- list(
    province   = c("省（等积投影）.shp", "province.shp"),
    prefecture = c("市（等积投影）.shp", "prefecture.shp"),
    county     = c("县（等积投影）.shp", "county.shp"),
    national   = c("国界（等积投影）.shp", "national.shp"),
    ninedash   = c("九段线.shp", "ninedash.shp")
  )
  for (nm in candidates[[layer]]) {
    path <- file.path(dir, nm)
    if (file.exists(path)) return(sf::st_read(path, quiet = TRUE))
  }
  stop("Basemap layer '", layer, "' not found under ", dir,
       "; check symlink data/spatial/basemap_GS2019_1822/")
}
