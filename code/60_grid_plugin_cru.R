# ============================================================
# Script: 60_grid_plugin_cru.R
# 设计/Design: 真正的网格原生 plug-in (终版)。省级拟合 temp_anom_z × effort
#   (HR=1.273, p=9e-11, 与头条 temp_grad 1.288 等强); 把【CRU 时变温度异常】+
#   【合并 Combined 调查努力(eBird+China_Birdwatch)】两者都在 100km 社区网格
#   原生、密集、时变, 以网格原生 z 表达, 套省级斜率做 plug-in 相对 hazard 面。
#   注: 网格事件仅 48 个 → 网格"重拟合"不可行, 故用省级拟合 + 网格 plug-in。
#   未来: 温度异常按 SSP 前推(+0.3/+0.8 SD/decade, 沿用既有方法; 无本地 CMIP6)。
# 输入: data/raw/hazard_risk_upgraded_complete_case.csv
#       data/raw/climate_metrics_province_year.csv (province temp_anom_z)
#       data/derived/grid_100km_climate_cru.parquet  (grid×year temp_anom)  [59]
#       data/derived/grid_100km_effort_combined.parquet (grid×year Combined) [59]
#       data/derived/community_grid_100km_climate_native.csv (grid lon/lat)
# 输出: results/forecasts/table_grid_100km_plugin_cru.csv
#       figures/main/Figure_8_grid100_native_plugin_hazard.{pdf,png}
# 运行: Rscript --no-init-file code/60_grid_plugin_cru.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(glmmTMB); library(arrow)
  library(sf); library(ggplot2); library(viridisLite)
})
sf::sf_use_s2(FALSE); options(warn=1)
V2 <- normalizePath(".", mustWork=TRUE)
SHP <- file.path(V2,"data","spatial","basemap_GS2019_1822")
log <- function(...) cat(sprintf("[60 %s] ",format(Sys.time(),"%H:%M:%S")),...,"\n")
inv_cloglog <- function(eta) 1 - exp(-exp(eta))
zc <- function(x){s<-sd(x,na.rm=TRUE); if(!is.finite(s)||s==0) rep(0,length(x)) else (x-mean(x,na.rm=TRUE))/s}

# ---- 1. province fit: temp_anom_z * effort ----
r <- fread(file.path(V2,"data/raw/hazard_risk_upgraded_complete_case.csv"))
cm <- fread(file.path(V2,"data/raw/climate_metrics_province_year.csv"))
d <- merge(r, cm[,.(province,year,temp_anom_z)], by=c("province","year"), all.x=TRUE)
d <- d[is.finite(temp_anom_z)]
m4 <- glmmTMB(event~temp_anom_z*log_effort_visits_z+(1|species)+(1|province),
              data=d, family=binomial("cloglog"))
cf <- fixef(m4)$cond
b0<-cf["(Intercept)"]; bC<-cf["temp_anom_z"]; bE<-cf["log_effort_visits_z"]
bI<-cf["temp_anom_z:log_effort_visits_z"]
log(sprintf("province M4(temp_anom): b0=%.3f bC=%.3f bE=%.3f bI=%.3f (int HR=%.3f)",
            b0,bC,bE,bI,exp(bI)))

# ---- 2. grid-native covariates (baseline = 2020-2024 mean), grid-native z ----
clim <- as.data.table(read_parquet(file.path(V2,"data/derived/grid_100km_climate_cru.parquet")))
eff  <- as.data.table(read_parquet(file.path(V2,"data/derived/grid_100km_effort_combined.parquet")))
xy   <- fread(file.path(V2,"data/derived/community_grid_100km_climate_native.csv"))[,.(grid_id,lon,lat)]
cb <- clim[year %in% 2020:2024, .(temp_anom = mean(temp_anom,na.rm=TRUE)), by=grid_id]
eb <- eff[year %in% 2020:2024, .(n_events = mean(n_events,na.rm=TRUE)), by=grid_id]
g <- Reduce(function(a,b) merge(a,b,by="grid_id",all=FALSE), list(xy, cb, eb))
g <- g[is.finite(temp_anom) & is.finite(n_events) & is.finite(lon) & is.finite(lat)]
g[, climate_z := zc(temp_anom)]                 # grid-native z (relative across grids)
g[, effort_z  := zc(log1p(n_events))]
log(sprintf("grid cells modelled=%d | effort>0 cells=%d | climate_z[%.2f,%.2f] effort_z[%.2f,%.2f]",
            nrow(g), sum(g$n_events>0), min(g$climate_z),max(g$climate_z),
            min(g$effort_z),max(g$effort_z)))

# ---- 3. predict current + SSP futures (perturb climate_z forward) ----
incr <- c(SSP245=0.3, SSP585=0.8); yrs<-c(2030,2050,2080); dec<-(yrs-2024)/10
pr <- function(cz,ez) inv_cloglog(b0 + bC*cz + bE*ez + bI*cz*ez)
out <- list(g[, .(grid_id,lon,lat,ssp="Current",year=2024L,climate_z,effort_z,
                  hazard=pr(climate_z,effort_z))])
for(sp in names(incr)) for(k in seq_along(yrs))
  out[[paste0(sp,yrs[k])]] <- g[, .(grid_id,lon,lat,ssp=sp,year=yrs[k],
     climate_z=climate_z+incr[[sp]]*dec[k], effort_z,
     hazard=pr(climate_z+incr[[sp]]*dec[k],effort_z))]
res <- rbindlist(out)
fwrite(res, file.path(V2,"results/forecasts/table_grid_100km_plugin_cru.csv"))
log(sprintf("wrote table_grid_100km_plugin_cru.csv (%d rows) | Current mean=%.4f max=%.4f | SSP585-2050 mean=%.4f",
            nrow(res), mean(out[[1]]$hazard), max(out[[1]]$hazard),
            mean(out[["SSP5852050"]]$hazard)))

# ---- 4. map ----
prov <- tryCatch(st_transform(st_make_valid(st_read(file.path(SHP,"省（等积投影）.shp"),quiet=TRUE)),4326),
                 error=function(e) NULL)
res[, panel := factor(ifelse(ssp=="Current","Current (2024)",paste0(ssp," — ",year)),
     levels=c("Current (2024)",paste0("SSP245 — ",yrs),paste0("SSP585 — ",yrs)))]
qc <- quantile(res$hazard,0.99,na.rm=TRUE); res[, hz:=pmin(hazard,qc)]
p <- ggplot()
if(!is.null(prov)) p <- p + geom_sf(data=prov,fill="grey97",colour="grey80",linewidth=0.1)
p <- p + geom_tile(data=res, aes(lon,lat,fill=hz), width=1.0, height=1.0) +
  scale_fill_viridis_c(option="mako",direction=-1,name="Relative hazard") +
  facet_wrap(~panel,ncol=4) +
  coord_sf(xlim=c(73,135),ylim=c(17,54),expand=FALSE,datum=NA) +
  labs(title="100 km grid-native plug-in hazard — CRU temperature anomaly × merged survey effort",
       subtitle=paste0("Province M4 (temp_anom_z × visits, HR=1.27) applied to GRID-NATIVE covariates: ",
       "CRU TS time-varying temperature anomaly + merged eBird/GBIF + China-Birdwatch effort (Combined). ",
       "Relative surface (grid-native z). Future = climate +0.3/+0.8 SD/decade; effort held. 99th-pct clip."),
       x=NULL,y=NULL) +
  theme_bw(base_size=8) +
  theme(legend.position="bottom",panel.grid=element_blank(),
        plot.title=element_text(face="bold",size=9.5),
        plot.subtitle=element_text(size=6.6,colour="grey35"),
        strip.text=element_text(face="bold",size=7))
ggsave(file.path(V2,"figures/main/Figure_8_grid100_native_plugin_hazard.pdf"),p,
       width=24,height=14,units="cm",device=grDevices::cairo_pdf)
ggsave(file.path(V2,"figures/main/Figure_8_grid100_native_plugin_hazard.png"),p,
       width=24,height=14,units="cm",dpi=300)
log("wrote Figure_8_grid100_native_plugin_hazard.{pdf,png}")
log("DONE")
