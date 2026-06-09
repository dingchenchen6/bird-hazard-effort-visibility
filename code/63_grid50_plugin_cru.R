# ============================================================
# Script: 63_grid50_plugin_cru.R  —— 任务1: 50km 网格 plug-in
# 设计: 与 100km(code 60)同口径的纯 plug-in。省级拟合 temp_anom_z×visits,
#   套到 50km 网格原生气候(CRU 时变温度异常)×原生努力(合并 Combined),
#   网格原生 z + 省级斜率 → 相对 hazard 面(现状 + SSP 未来)。
# 50km 努力来源: 把已合并的 10km 努力表(Combined)按 50km 分箱聚合(轻量,
#   复用 occupancy 项目已去重的合并努力, 而非重聚合 1450 万原始事件)。
# 输入: 10km 努力表 + 10km 网格质心(occupancy 项目); CRU TS tmp .nc; 省界 shp;
#       data/raw/{hazard_risk_upgraded_complete_case,climate_metrics_province_year}.csv
# 输出: results/forecasts/table_grid_50km_plugin_cru.csv
#       figures/main/Figure_8b_grid50_native_plugin_hazard.{pdf,png}
# 运行: Rscript --no-init-file code/63_grid50_plugin_cru.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(glmmTMB); library(terra); library(sf); library(ggplot2); library(viridisLite)
})
sf::sf_use_s2(FALSE); options(warn=1)
V2 <- normalizePath(".", mustWork=TRUE)
NP <- "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis"
SHP <- file.path(V2,"data","spatial","basemap_GS2019_1822")
log <- function(...) cat(sprintf("[63 %s] ",format(Sys.time(),"%H:%M:%S")),...,"\n")
inv_cloglog <- function(eta) 1 - exp(-exp(eta))
zc <- function(x){s<-sd(x,na.rm=TRUE); if(!is.finite(s)||s==0) rep(0,length(x)) else (x-mean(x,na.rm=TRUE))/s}
PROV<-c("北京市"="Beijing","天津市"="Tianjin","河北省"="Hebei","山西省"="Shanxi","内蒙古自治区"="Inner Mongolia","辽宁省"="Liaoning","吉林省"="Jilin","黑龙江省"="Heilongjiang","上海市"="Shanghai","江苏省"="Jiangsu","浙江省"="Zhejiang","安徽省"="Anhui","福建省"="Fujian","江西省"="Jiangxi","山东省"="Shandong","河南省"="Henan","湖北省"="Hubei","湖南省"="Hunan","广东省"="Guangdong","广西壮族自治区"="Guangxi","海南省"="Hainan","重庆市"="Chongqing","四川省"="Sichuan","贵州省"="Guizhou","云南省"="Yunnan","西藏自治区"="Tibet","陕西省"="Shaanxi","甘肃省"="Gansu","青海省"="Qinghai","宁夏回族自治区"="Ningxia","新疆维吾尔自治区"="Xinjiang","台湾省"="Taiwan","香港特别行政区"="Hong Kong","澳门特别行政区"="Macau")

# ---- 1. province headline plug-in model: temp_anom_z * visits ----
r <- fread(file.path(V2,"data/raw/hazard_risk_upgraded_complete_case.csv"))
cm <- fread(file.path(V2,"data/raw/climate_metrics_province_year.csv"))
d <- merge(r, cm[,.(province,year,temp_anom_z)], by=c("province","year"), all.x=TRUE)[is.finite(temp_anom_z)]
m4 <- glmmTMB(event~temp_anom_z*log_effort_visits_z+(1|species)+(1|province),data=d,family=binomial("cloglog"))
cf<-fixef(m4)$cond; b0<-cf[1]; bC<-cf["temp_anom_z"]; bE<-cf["log_effort_visits_z"]; bI<-cf["temp_anom_z:log_effort_visits_z"]
log(sprintf("province M4: b0=%.3f bC=%.3f bE=%.3f bI=%.3f (int HR=%.3f)",b0,bC,bE,bI,exp(bI)))

# ---- 2. 10km Combined effort -> 50km bins ----
e10 <- fread(file.path(NP,"results_v2/table_effort_by_grid_year_source_10km.csv"))[source_short=="Combined"]
g10 <- as.data.table(sf::st_drop_geometry(readRDS(file.path(NP,"data/derived_v2/china_grid_10km_v2.rds"))))[,.(grid_cell,centroid_lon,centroid_lat)]
e10 <- merge(e10, g10, by="grid_cell")
# project centroids to Albers, bin to 50km
xy <- sf::sf_project("EPSG:4326","EPSG:4524", as.matrix(e10[,.(centroid_lon,centroid_lat)]))
e10[, `:=`(ax=xy[,1], ay=xy[,2])]
cell_m <- 50000
x0 <- floor(min(e10$ax)/cell_m)*cell_m; y0 <- floor(min(e10$ay)/cell_m)*cell_m
e10[, `:=`(ix=floor((ax-x0)/cell_m), iy=floor((ay-y0)/cell_m))]
e10[, cell50 := paste(ix,iy,sep="_")]
eff50 <- e10[, .(n_events=sum(n_events,na.rm=TRUE)), by=.(cell50,ix,iy,year)]
# 50km cell centroids (Albers bin centre -> WGS84)
cells <- unique(eff50[,.(cell50,ix,iy)])
cc <- sf::sf_project("EPSG:4524","EPSG:4326", as.matrix(cells[,.(x0+(ix+0.5)*cell_m, y0+(iy+0.5)*cell_m)]))
cells[, `:=`(lon=cc[,1], lat=cc[,2])]
log("50km cells with effort:", nrow(cells), "| effort rows:", nrow(eff50))

# ---- 3. province per 50km cell (point-in-polygon) ----
prov <- st_make_valid(st_transform(st_read(file.path(SHP,"省（等积投影）.shp"),quiet=TRUE),4326))
prov$province <- unname(PROV[as.character(prov[["省"]])])
ji <- as.data.table(st_drop_geometry(st_join(st_as_sf(cells,coords=c("lon","lat"),crs=4326,remove=FALSE),prov["province"],join=st_intersects)))
cells <- merge(cells, ji[!is.na(province),.(cell50,province)], by="cell50")
cells <- cells[!duplicated(cell50)]
log("50km cells assigned to province:", nrow(cells))

# ---- 4. CRU temp_anom per 50km centroid x year ----
r_cru <- tryCatch(rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc"),subds="tmp"),
                  error=function(e) rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc")))
idx<-function(y,m)(y-1901)*12+m
pts <- vect(cells[,.(lon,lat)], geom=c("lon","lat"), crs="EPSG:4326")
ann <- list()
for(y in 2002:2024){ li<-idx(y,1):idx(y,12); li<-li[li<=nlyr(r_cru)]; if(length(li)<12) next
  ym<-mean(r_cru[[li]],na.rm=TRUE); ann[[as.character(y)]]<-data.table(cell50=cells$cell50,year=y,tmean=terra::extract(ym,pts)[,2]) }
clim <- rbindlist(ann)
base <- clim[year%in%2002:2010,.(tbase=mean(tmean,na.rm=TRUE)),by=cell50]
clim <- merge(clim,base,by="cell50"); clim[,temp_anom:=tmean-tbase]

# ---- 5. baseline (2020-2024 mean) + grid-native z ----
cb <- clim[year%in%2020:2024,.(temp_anom=mean(temp_anom,na.rm=TRUE)),by=cell50]
eb <- eff50[year%in%2020:2024,.(n_events=mean(n_events,na.rm=TRUE)),by=cell50]
g <- Reduce(function(a,b) merge(a,b,by="cell50",all=FALSE), list(cells,cb,eb))
g <- g[is.finite(temp_anom)&is.finite(n_events)&is.finite(lon)&is.finite(lat)]
g[,climate_z:=zc(temp_anom)]; g[,effort_z:=zc(log1p(n_events))]
log(sprintf("50km cells modelled=%d | climate_z[%.2f,%.2f] effort_z[%.2f,%.2f]",nrow(g),min(g$climate_z),max(g$climate_z),min(g$effort_z),max(g$effort_z)))

# ---- 6. predict current + SSP ----
incr<-c(SSP245=0.3,SSP585=0.8); yrs<-c(2030,2050,2080); dec<-(yrs-2024)/10
pr<-function(cz,ez) inv_cloglog(b0+bC*cz+bE*ez+bI*cz*ez)
out<-list(g[,.(cell50,lon,lat,ssp="Current",year=2024L,hazard=pr(climate_z,effort_z))])
for(sp in names(incr)) for(k in seq_along(yrs)) out[[paste0(sp,yrs[k])]]<-g[,.(cell50,lon,lat,ssp=sp,year=yrs[k],hazard=pr(climate_z+incr[[sp]]*dec[k],effort_z))]
res<-rbindlist(out)
fwrite(res,file.path(V2,"results/forecasts/table_grid_50km_plugin_cru.csv"))
log(sprintf("wrote table_grid_50km_plugin_cru.csv (%d rows) | Current mean=%.4f max=%.4f | SSP585-2050 mean=%.4f",
    nrow(res),mean(out[[1]]$hazard),max(out[[1]]$hazard),mean(out[["SSP5852050"]]$hazard)))

# ---- 7. map ----
res[,panel:=factor(ifelse(ssp=="Current","Current (2024)",paste0(ssp," — ",year)),levels=c("Current (2024)",paste0("SSP245 — ",yrs),paste0("SSP585 — ",yrs)))]
qc<-quantile(res$hazard,0.99,na.rm=TRUE); res[,hz:=pmin(hazard,qc)]
p<-ggplot()+geom_sf(data=prov,fill="grey97",colour="grey80",linewidth=0.1)+
  geom_tile(data=res,aes(lon,lat,fill=hz),width=0.55,height=0.55)+
  scale_fill_viridis_c(option="mako",direction=-1,name="Relative hazard")+
  facet_wrap(~panel,ncol=4)+coord_sf(xlim=c(73,135),ylim=c(17,54),expand=FALSE,datum=NA)+
  labs(title="50 km grid plug-in hazard — CRU temperature anomaly × merged survey effort",
       subtitle=paste0("Province M4 (temp_anom_z × visits, HR=",sprintf("%.2f",exp(bI)),") applied to 50 km grid-native CRU climate + merged ",
       "eBird/GBIF+China-Birdwatch effort (aggregated from 10 km Combined). Relative surface; future = climate +0.3/+0.8 SD/decade; effort held."),
       x=NULL,y=NULL)+theme_bw(base_size=8)+
  theme(legend.position="bottom",panel.grid=element_blank(),plot.title=element_text(face="bold",size=9.5),
        plot.subtitle=element_text(size=6.6,colour="grey35"),strip.text=element_text(face="bold",size=7))
ggsave(file.path(V2,"figures/main/Figure_8b_grid50_native_plugin_hazard.pdf"),p,width=24,height=14,units="cm",device=grDevices::cairo_pdf)
ggsave(file.path(V2,"figures/main/Figure_8b_grid50_native_plugin_hazard.png"),p,width=24,height=14,units="cm",dpi=300)
log("wrote Figure_8b_grid50_native_plugin_hazard.{pdf,png}"); log("DONE")
