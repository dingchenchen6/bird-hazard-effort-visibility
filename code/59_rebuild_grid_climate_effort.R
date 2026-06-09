# ============================================================
# Script: 59_rebuild_grid_climate_effort.R
# 目的/Goal: 修复网格协变量两处问题, 为真正的网格原生分析铺路。
#   (A) 努力源修正: 旧 28b 误用"新记录坐标"(events_*_grid_assigned, 2024 仅 64 格)
#       作努力。这里改用【真正合并调查努力】table_effort_by_grid_year_source_100km.csv
#       的 Combined 源(eBird_GBIF + China_Birdwatch), 2024 有 ~1120 格 → 密集且正确。
#   (B) 时变气候: 用 CRU TS 4.09 月值(1901-2024)逐年提取每个社区网格的年均温,
#       计算温度异常(vs 2002-2010 基线), 得到 grid×year 时变气候(覆盖 2024)。
#   两者都对齐到同一套社区 100km 网格(grid_id 54-3665, community_grid_100km_climate_native)。
# 输入: New project/bird_dynamic_occupancy_analysis/results_v2/table_effort_by_grid_year_source_100km.csv
#       New project/bird_dynamic_occupancy_analysis/data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc
#       data/derived/community_grid_100km_climate_native.csv  (grid_id, lon, lat, bio1, mahalanobis_dist)
# 输出: data/derived/grid_100km_effort_combined.parquet   (grid_id, year, effort + z)
#       data/derived/grid_100km_climate_cru.parquet        (grid_id, year, tmean, temp_anom, ...)
#       results/diagnostics/table_grid_effort_density_fix.csv  (旧vs新 2024 覆盖对比)
# 运行: Rscript --no-init-file code/59_rebuild_grid_climate_effort.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(terra); library(arrow)
})
options(warn = 1)
V2 <- normalizePath(".", mustWork = TRUE)
NP <- "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis"
log <- function(...) cat(sprintf("[59 %s] ", format(Sys.time(),"%H:%M:%S")),...,"\n")

# grid geometry / coordinates (community 100km grid)
gc <- fread(file.path(V2,"data/derived/community_grid_100km_climate_native.csv"))
grid_xy <- gc[, .(grid_id, lon, lat, bio1, mahalanobis_dist)]
log("community grid cells:", nrow(grid_xy))

# ---------- (A) corrected effort from real merged Combined source ----------
eff <- fread(file.path(NP,"results_v2/table_effort_by_grid_year_source_100km.csv"))
ec <- eff[source_short == "Combined" & year >= 2002 & year <= 2024,
          .(grid_id = grid_cell, year, n_events, n_observers, n_dates)]
# fill zero cells (grid x year full panel over community grid)
full <- CJ(grid_id = grid_xy$grid_id, year = 2002:2024)
ec <- merge(full, ec, by = c("grid_id","year"), all.x = TRUE)
for (c in c("n_events","n_observers","n_dates")) ec[is.na(get(c)), (c) := 0]
# log + within-year grid-native z (effort = checklists/events ~ "visits")
ec[, `:=`(log_events = log1p(n_events), log_obs = log1p(n_observers),
          log_dates = log1p(n_dates))]
zf <- function(x){ s<-sd(x); if(!is.finite(s)||s==0) rep(0,length(x)) else (x-mean(x))/s }
ec[, `:=`(log_effort_visits_z = zf(log_events),
          log_effort_obs_z    = zf(log_obs),
          log_effort_dates_z  = zf(log_dates)), by = year]
write_parquet(ec, file.path(V2,"data/derived/grid_100km_effort_combined.parquet"))
log("wrote grid_100km_effort_combined.parquet:", nrow(ec),"rows")
log(sprintf("2024 grids with events>0 = %d (was 64 in old wrong panel)",
            ec[year==2024 & n_events>0, .N]))

# density-fix comparison table
oldp <- tryCatch(as.data.table(read_parquet(file.path(V2,
          "data/derived/grid_100km_effort_native.parquet"))), error=function(e) NULL)
cmp <- data.table(
  panel = c("OLD (28b, from new-record coords)","NEW (59, merged Combined effort)"),
  source = c("events_100km_grid_assigned.csv","table_effort_by_grid_year_source_100km.csv (Combined)"),
  cells_2024_nonzero = c(if(!is.null(oldp)) oldp[year==2024 & n_visits_grid>0,.N] else NA_integer_,
                          ec[year==2024 & n_events>0,.N]),
  total_cells = c(if(!is.null(oldp)) uniqueN(oldp$grid_id) else NA_integer_, uniqueN(ec$grid_id)))
fwrite(cmp, file.path(V2,"results/diagnostics/table_grid_effort_density_fix.csv"))
print(cmp)

# ---------- (B) CRU TS time-varying climate per grid x year ----------
cru_f <- file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc")
r <- tryCatch(rast(cru_f, subds="tmp"), error=function(e) rast(cru_f))
log("CRU layers:", nlyr(r), "| taking 2002-2024 monthly (tmp)")
# CRU monthly from 1901-01; layer index for year-month:
base_year <- 1901
idx <- function(y,m) (y-base_year)*12 + m
lyr_idx <- unlist(lapply(2002:2024, function(y) idx(y,1):idx(y,12)))
lyr_idx <- lyr_idx[lyr_idx <= nlyr(r)]
rr <- r[[lyr_idx]]
# annual mean per year
pts <- vect(grid_xy[, .(lon, lat)], geom=c("lon","lat"), crs="EPSG:4326")
ann <- list()
for (y in 2002:2024) {
  li <- idx(y,1):idx(y,12); li <- li[li <= nlyr(r)]
  if (length(li) < 12) next
  ym <- mean(r[[li]], na.rm=TRUE)
  v <- terra::extract(ym, pts)[,2]
  ann[[as.character(y)]] <- data.table(grid_id = grid_xy$grid_id, year = y, tmean = v)
}
clim <- rbindlist(ann)
# temperature anomaly vs 2002-2010 baseline per grid; spatial mean-centred too
base <- clim[year %in% 2002:2010, .(tbase = mean(tmean, na.rm=TRUE)), by = grid_id]
clim <- merge(clim, base, by = "grid_id", all.x = TRUE)
clim[, temp_anom := tmean - tbase]
clim[, temp_anom_z := { s<-sd(temp_anom,na.rm=TRUE); if(!is.finite(s)||s==0) 0 else (temp_anom-mean(temp_anom,na.rm=TRUE))/s }, by = year]
write_parquet(clim, file.path(V2,"data/derived/grid_100km_climate_cru.parquet"))
log("wrote grid_100km_climate_cru.parquet:", nrow(clim),"rows (grid x year 2002-2024)")
log(sprintf("CRU sanity: 2024 mean tmean=%.2f C | mean temp_anom(2024)=%.3f C",
            clim[year==2024, mean(tmean,na.rm=TRUE)], clim[year==2024, mean(temp_anom,na.rm=TRUE)]))
log("DONE")
