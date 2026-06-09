# ============================================================
# Script: 64_cmip6_grid_futures.R —— 任务2: CMIP6 多模型集成替代 SSP 扰动
# 设计: 用 WorldClim CMIP6 未来 bio1(年均温) 多 GCM × SSP × 时段, 相对 WorldClim
#   历史 bio1 算【升温增量 delta = 未来 − 历史】(同产品→无跨产品偏差), 5 GCM 取
#   集成中位数; future 温度异常 = CRU 现状异常 + delta; 换算到模型 z 单位 plug-in。
#   替代 code 60/63 的 ±0.3/0.8 SD/decade 扰动。同时更新 100km + 50km 未来面。
# 输入: WorldClim 历史 bio1 (本地); data/external/cmip6_worldclim/*.tif (5 GCM×SSP×3 时段);
#       CRU 现状异常 + Combined 努力(100km: 派生 parquet; 50km: 由 10km 聚合);
#       省级 M4 (temp_anom_z × visits)。
# 输出: results/forecasts/table_grid_{100,50}km_plugin_cmip6.csv
#       figures/main/Figure_8_grid100_native_plugin_hazard.{pdf,png}  (CMIP6 版, 覆盖)
#       figures/main/Figure_8b_grid50_native_plugin_hazard.{pdf,png}  (CMIP6 版, 覆盖)
#       results/diagnostics/table_cmip6_ensemble_delta_summary.csv
# 运行: Rscript --no-init-file code/64_cmip6_grid_futures.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(glmmTMB); library(terra); library(sf); library(ggplot2); library(viridisLite); library(arrow)
})
sf::sf_use_s2(FALSE); options(warn=1)
V2 <- normalizePath(".", mustWork=TRUE)
NP <- "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis"
WC_HIST <- "/Users/dingchenchen/Documents/New project/bird_grid_community_analysis/data/external/worldclim/climate/wc2.1_10m/wc2.1_10m_bio_1.tif"
CMIP_DIR <- file.path(V2,"data/external/cmip6_worldclim")
SHP <- file.path(V2,"data","spatial","basemap_GS2019_1822")
log <- function(...) cat(sprintf("[64 %s] ",format(Sys.time(),"%H:%M:%S")),...,"\n")
inv_cloglog <- function(eta) 1 - exp(-exp(eta))
GCMS<-c("ACCESS-CM2","MPI-ESM1-2-HR","MIROC6","UKESM1-0-LL")  # EC-Earth3 unavailable at WorldClim 10m
SSPS<-c("ssp245","ssp585"); PERIODS<-c("2021-2040","2041-2060","2061-2080"); PYEAR<-c(2030,2050,2080)
PROV<-c("北京市"="Beijing","天津市"="Tianjin","河北省"="Hebei","山西省"="Shanxi","内蒙古自治区"="Inner Mongolia","辽宁省"="Liaoning","吉林省"="Jilin","黑龙江省"="Heilongjiang","上海市"="Shanghai","江苏省"="Jiangsu","浙江省"="Zhejiang","安徽省"="Anhui","福建省"="Fujian","江西省"="Jiangxi","山东省"="Shandong","河南省"="Henan","湖北省"="Hubei","湖南省"="Hunan","广东省"="Guangdong","广西壮族自治区"="Guangxi","海南省"="Hainan","重庆市"="Chongqing","四川省"="Sichuan","贵州省"="Guizhou","云南省"="Yunnan","西藏自治区"="Tibet","陕西省"="Shaanxi","甘肃省"="Gansu","青海省"="Qinghai","宁夏回族自治区"="Ningxia","新疆维吾尔自治区"="Xinjiang","台湾省"="Taiwan","香港特别行政区"="Hong Kong","澳门特别行政区"="Macau")

# ---- province plug-in model ----
r<-fread(file.path(V2,"data/raw/hazard_risk_upgraded_complete_case.csv"))
cm<-fread(file.path(V2,"data/raw/climate_metrics_province_year.csv"))
d<-merge(r,cm[,.(province,year,temp_anom_z)],by=c("province","year"),all.x=TRUE)[is.finite(temp_anom_z)]
m4<-glmmTMB(event~temp_anom_z*log_effort_visits_z+(1|species)+(1|province),data=d,family=binomial("cloglog"))
cf<-fixef(m4)$cond; b0<-cf[1];bC<-cf["temp_anom_z"];bE<-cf["log_effort_visits_z"];bI<-cf["temp_anom_z:log_effort_visits_z"]
log(sprintf("province M4 int HR=%.3f",exp(bI)))

bio_hist <- rast(WC_HIST)
prov<-st_make_valid(st_transform(st_read(file.path(SHP,"省（等积投影）.shp"),quiet=TRUE),4326)); prov$province<-unname(PROV[as.character(prov[["省"]])])

# ensemble CMIP6 delta (future bio1 - hist bio1) at given points -> data.table(ssp,year,delta)
ensemble_delta <- function(pts_dt){
  P<-vect(as.data.frame(pts_dt[,.(lon,lat)]),geom=c("lon","lat"),crs="EPSG:4326")
  hist0<-terra::extract(bio_hist,P)[,2]
  out<-list()
  for(si in seq_along(SSPS)) for(pi in seq_along(PERIODS)){
    dl<-matrix(NA_real_,nrow=nrow(pts_dt),ncol=length(GCMS))
    for(gi in seq_along(GCMS)){
      f<-file.path(CMIP_DIR,sprintf("wc2.1_10m_bioc_%s_%s_%s.tif",GCMS[gi],SSPS[si],PERIODS[pi]))
      if(!file.exists(f)) next
      fut<-rast(f,lyrs=1)   # bio1 = band 1
      dl[,gi]<-terra::extract(fut,P)[,2]-hist0
    }
    delta<-apply(dl,1,median,na.rm=TRUE)
    out[[paste0(SSPS[si],pi)]]<-data.table(cell=pts_dt$cell,ssp=toupper(SSPS[si]),year=PYEAR[pi],delta=delta)
  }
  rbindlist(out)
}

# generic grid plug-in with CMIP6 futures
run_grid <- function(g, label, fig, tbl, tile_w){
  # g: data.table(cell, lon, lat, temp_anom (present baseline), n_events)
  g<-g[is.finite(temp_anom)&is.finite(n_events)&is.finite(lon)&is.finite(lat)]
  mC<-mean(g$temp_anom); sC<-sd(g$temp_anom); mE<-mean(log1p(g$n_events)); sE<-sd(log1p(g$n_events))
  g[,climate_z:=(temp_anom-mC)/sC]; g[,effort_z:=(log1p(n_events)-mE)/sE]
  pr<-function(cz,ez) inv_cloglog(b0+bC*cz+bE*ez+bI*cz*ez)
  ed<-ensemble_delta(g[,.(cell,lon,lat)])
  out<-list(g[,.(cell,lon,lat,ssp="Current",year=2024L,hazard=pr(climate_z,effort_z))])
  for(si in c("SSP245","SSP585")) for(yy in PYEAR){
    dd<-ed[ssp==si&year==yy]; gg<-merge(g,dd[,.(cell,delta)],by="cell")
    gg[,cz_fut:=climate_z + delta/sC]   # CMIP6 warming in z units
    out[[paste0(si,yy)]]<-gg[,.(cell,lon,lat,ssp=si,year=yy,hazard=pr(cz_fut,effort_z))]
  }
  res<-rbindlist(out); fwrite(res,file.path(V2,tbl))
  log(sprintf("%s: cells=%d Current mean=%.4f SSP585-2050 mean=%.4f SSP585-2080 mean=%.4f",
      label,nrow(g),mean(out[[1]]$hazard),mean(out[["SSP5852050"]]$hazard),mean(out[["SSP5852080"]]$hazard)))
  res[,panel:=factor(ifelse(ssp=="Current","Current (2024)",paste0(ssp," — ",year)),levels=c("Current (2024)",paste0("SSP245 — ",PYEAR),paste0("SSP585 — ",PYEAR)))]
  qc<-quantile(res$hazard,0.99,na.rm=TRUE); res[,hz:=pmin(hazard,qc)]
  p<-ggplot()+geom_sf(data=prov,fill="grey97",colour="grey80",linewidth=0.1)+
    geom_tile(data=res,aes(lon,lat,fill=hz),width=tile_w,height=tile_w)+
    scale_fill_viridis_c(option="mako",direction=-1,name="Relative hazard")+
    facet_wrap(~panel,ncol=4)+coord_sf(xlim=c(73,135),ylim=c(17,54),expand=FALSE,datum=NA)+
    labs(title=sprintf("%s grid plug-in hazard — CRU climate × merged effort, CMIP6 ensemble futures",label),
         subtitle=paste0("Province M4 (temp_anom_z × visits, HR=",sprintf("%.2f",exp(bI)),"). Futures = median warming of 4 CMIP6 GCMs ",
         "(ACCESS-CM2/MPI-ESM1-2-HR/MIROC6/UKESM1-0-LL) under SSP245/585 vs WorldClim baseline, added to CRU anomaly. Effort held."),
         x=NULL,y=NULL)+theme_bw(base_size=8)+
    theme(legend.position="bottom",panel.grid=element_blank(),plot.title=element_text(face="bold",size=9.5),
          plot.subtitle=element_text(size=6.4,colour="grey35"),strip.text=element_text(face="bold",size=7))
  ggsave(file.path(V2,paste0(fig,".pdf")),p,width=24,height=14,units="cm",device=grDevices::cairo_pdf)
  ggsave(file.path(V2,paste0(fig,".png")),p,width=24,height=14,units="cm",dpi=300)
  log("wrote",fig)
  invisible(ed)
}

# ---- 100km grid panel (community grid) ----
xy100<-fread(file.path(V2,"data/derived/community_grid_100km_climate_native.csv"))[,.(cell=grid_id,lon,lat)]
clim100<-as.data.table(read_parquet(file.path(V2,"data/derived/grid_100km_climate_cru.parquet")))[year%in%2020:2024,.(temp_anom=mean(temp_anom,na.rm=TRUE)),by=.(cell=grid_id)]
eff100<-as.data.table(read_parquet(file.path(V2,"data/derived/grid_100km_effort_combined.parquet")))[year%in%2020:2024,.(n_events=mean(n_events,na.rm=TRUE)),by=.(cell=grid_id)]
g100<-Reduce(function(a,b)merge(a,b,by="cell"),list(xy100,clim100,eff100))
ed100<-run_grid(g100,"100 km","figures/main/Figure_8_grid100_native_plugin_hazard","results/forecasts/table_grid_100km_plugin_cmip6.csv",1.0)

# ---- 50km grid panel (rebuild from 10km Combined, same as code 63) ----
e10<-fread(file.path(NP,"results_v2/table_effort_by_grid_year_source_10km.csv"))[source_short=="Combined"]
g10<-as.data.table(sf::st_drop_geometry(readRDS(file.path(NP,"data/derived_v2/china_grid_10km_v2.rds"))))[,.(grid_cell,centroid_lon,centroid_lat)]
e10<-merge(e10,g10,by="grid_cell"); xyA<-sf::sf_project("EPSG:4326","EPSG:4524",as.matrix(e10[,.(centroid_lon,centroid_lat)]))
e10[,`:=`(ax=xyA[,1],ay=xyA[,2])]; cm_<-50000; x0<-floor(min(e10$ax)/cm_)*cm_; y0<-floor(min(e10$ay)/cm_)*cm_
e10[,`:=`(ix=floor((ax-x0)/cm_),iy=floor((ay-y0)/cm_))]; e10[,cell:=paste(ix,iy,sep="_")]
eff50<-e10[year%in%2020:2024,.(n_events=mean(n_events,na.rm=TRUE)),by=.(cell,ix,iy)]
cc<-sf::sf_project("EPSG:4524","EPSG:4326",as.matrix(eff50[,.(x0+(ix+0.5)*cm_,y0+(iy+0.5)*cm_)])); eff50[,`:=`(lon=cc[,1],lat=cc[,2])]
P50<-vect(eff50[,.(lon,lat)],geom=c("lon","lat"),crs="EPSG:4326")
r_cru<-tryCatch(rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc"),subds="tmp"),error=function(e)rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc")))
idx<-function(y,m)(y-1901)*12+m
annf<-function(y){li<-idx(y,1):idx(y,12);li<-li[li<=nlyr(r_cru)];if(length(li)<12)return(NULL);terra::extract(mean(r_cru[[li]],na.rm=TRUE),P50)[,2]}
t20<-rowMeans(sapply(2020:2024,annf),na.rm=TRUE); tb<-rowMeans(sapply(2002:2010,annf),na.rm=TRUE)
eff50[,temp_anom:=t20-tb]
g50<-eff50[,.(cell,lon,lat,temp_anom,n_events)]
ed50<-run_grid(g50,"50 km","figures/main/Figure_8b_grid50_native_plugin_hazard","results/forecasts/table_grid_50km_plugin_cmip6.csv",0.55)

# ensemble delta summary
ds<-rbind(ed100[,.(grid="100km",ssp,year,delta_C=delta)],ed50[,.(grid="50km",ssp,year,delta_C=delta)])
summ<-ds[,.(mean_delta_C=round(mean(delta_C,na.rm=TRUE),2),median_delta_C=round(median(delta_C,na.rm=TRUE),2)),by=.(grid,ssp,year)]
fwrite(summ,file.path(V2,"results/diagnostics/table_cmip6_ensemble_delta_summary.csv"))
print(summ); log("DONE")
