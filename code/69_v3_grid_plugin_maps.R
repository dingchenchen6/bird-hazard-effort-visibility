# ============================================================
# Script: 69_v3_grid_plugin_maps.R —— v3 宽松集的网格 plug-in 未来图(与 v2 Fig 8/8b 并列)
# 设计: 与 code 66 同口径, 但省级 plug-in 模型在【v3 宽松集】上拟合
#   (temp_anom_z × visits, HR=1.274), 套到相同的 100/50km 网格原生 CRU 气候 ×
#   合并 Combined 努力; 未来 = 4-GCM CMIP6 集成升温 + SSP 差异化努力增长。
#   九段线官方底图 + 冷→暖配色 + 四周留边距。
# 输出: figures/main/Figure_8_v3_grid100_plugin_hazard.{pdf,png}
#       figures/main/Figure_8b_v3_grid50_plugin_hazard.{pdf,png}
#       results/forecasts/table_grid_100km_plugin_v3.csv
# 运行: Rscript --no-init-file code/69_v3_grid_plugin_maps.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(glmmTMB); library(terra); library(sf); library(ggplot2); library(arrow)
})
sf::sf_use_s2(FALSE); options(warn=1)
V2<-normalizePath(".",mustWork=TRUE); NP<-"/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis"
WC_HIST<-"/Users/dingchenchen/Documents/New project/bird_grid_community_analysis/data/external/worldclim/climate/wc2.1_10m/wc2.1_10m_bio_1.tif"
CMIP_DIR<-file.path(V2,"data/external/cmip6_worldclim"); SHP<-file.path(V2,"data","spatial","basemap_GS2019_1822")
log<-function(...) cat(sprintf("[69 %s] ",format(Sys.time(),"%H:%M:%S")),...,"\n")
inv_cloglog<-function(eta) 1-exp(-exp(eta))
GCMS<-c("ACCESS-CM2","MPI-ESM1-2-HR","MIROC6","UKESM1-0-LL"); SSPS<-c("ssp245","ssp585")
PERIODS<-c("2021-2040","2041-2060","2061-2080"); PYEAR<-c(2030,2050,2080)
EFF_GROW<-c(SSP245=0.30,SSP585=0.60)
COOL_WARM<-c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF","#FEE090","#FDAE61","#F46D43","#D73027","#A50026")
XLIM<-c(71,137); YLIM<-c(2.5,55)
PROV<-c("北京市"="Beijing","天津市"="Tianjin","河北省"="Hebei","山西省"="Shanxi","内蒙古自治区"="Inner Mongolia","辽宁省"="Liaoning","吉林省"="Jilin","黑龙江省"="Heilongjiang","上海市"="Shanghai","江苏省"="Jiangsu","浙江省"="Zhejiang","安徽省"="Anhui","福建省"="Fujian","江西省"="Jiangxi","山东省"="Shandong","河南省"="Henan","湖北省"="Hubei","湖南省"="Hunan","广东省"="Guangdong","广西壮族自治区"="Guangxi","海南省"="Hainan","重庆市"="Chongqing","四川省"="Sichuan","贵州省"="Guizhou","云南省"="Yunnan","西藏自治区"="Tibet","陕西省"="Shaanxi","甘肃省"="Gansu","青海省"="Qinghai","宁夏回族自治区"="Ningxia","新疆维吾尔自治区"="Xinjiang","台湾省"="Taiwan","香港特别行政区"="Hong Kong","澳门特别行政区"="Macau")

# ---- v3 province plug-in model ----
v3<-fread(file.path(V2,"data/derived/risk_set_province_v3.csv"))
cm<-fread(file.path(V2,"data/raw/climate_metrics_province_year.csv"))[,.(province,year,temp_anom_z)]
d<-merge(v3,cm,by=c("province","year"),all.x=TRUE)[is.finite(temp_anom_z)&is.finite(log_effort_visits_z)]
m4<-glmmTMB(event~temp_anom_z*log_effort_visits_z+(1|species)+(1|province),d,family=binomial("cloglog"))
cf<-fixef(m4)$cond; b0<-cf[1];bC<-cf["temp_anom_z"];bE<-cf["log_effort_visits_z"];bI<-cf["temp_anom_z:log_effort_visits_z"]
log(sprintf("v3 province M4 int HR=%.3f (b0=%.2f bE=%.2f)",exp(bI),b0,bE))
bio_hist<-rast(WC_HIST)
prov<-st_make_valid(st_transform(st_read(file.path(SHP,"省（等积投影）.shp"),quiet=TRUE),4326))
nation<-st_make_valid(st_transform(st_read(file.path(SHP,"中国轮廓线.shp"),quiet=TRUE),4326))
dash<-st_transform(st_read(file.path(SHP,"九段线.shp"),quiet=TRUE),4326)
theme_map<-function() theme_bw(base_size=8)+theme(legend.position="bottom",panel.grid=element_blank(),panel.spacing=unit(0.45,"lines"),plot.margin=margin(8,10,8,10),plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=6.6,colour="grey35",margin=margin(b=4)),strip.background=element_rect(fill="grey96",colour="grey80"),strip.text=element_text(face="bold",size=7.5),legend.key.width=unit(1.1,"cm"),legend.key.height=unit(0.3,"cm"),legend.title=element_text(size=8,face="bold"))

ensemble_delta<-function(pts){P<-vect(as.data.frame(pts[,.(lon,lat)]),geom=c("lon","lat"),crs="EPSG:4326");h0<-terra::extract(bio_hist,P)[,2];out<-list()
  for(si in seq_along(SSPS))for(pi in seq_along(PERIODS)){dl<-matrix(NA_real_,nrow(pts),length(GCMS))
    for(gi in seq_along(GCMS)){f<-file.path(CMIP_DIR,sprintf("wc2.1_10m_bioc_%s_%s_%s.tif",GCMS[gi],SSPS[si],PERIODS[pi]));if(file.exists(f))dl[,gi]<-terra::extract(rast(f,lyrs=1),P)[,2]-h0}
    out[[paste0(si,pi)]]<-data.table(unit=pts$unit,ssp=toupper(SSPS[si]),year=PYEAR[pi],delta=apply(dl,1,median,na.rm=TRUE))};rbindlist(out)}
predict_long<-function(g){g<-g[is.finite(temp_anom)&is.finite(n_events)&is.finite(lon)&is.finite(lat)]
  mC<-mean(g$temp_anom);sC<-sd(g$temp_anom);mE<-mean(log1p(g$n_events));sE<-sd(log1p(g$n_events))
  g[,climate_z:=(temp_anom-mC)/sC];g[,effort_z:=(log1p(n_events)-mE)/sE];pr<-function(cz,ez)inv_cloglog(b0+bC*cz+bE*ez+bI*cz*ez)
  ed<-ensemble_delta(g[,.(unit,lon,lat)]);out<-list(g[,.(unit,lon,lat,ssp="Current",year=2024L,hazard=pr(climate_z,effort_z))])
  for(si in c("SSP245","SSP585"))for(yy in PYEAR){dec<-(yy-2024)/10;dd<-ed[ssp==si&year==yy];gg<-merge(g,dd[,.(unit,delta)],by="unit")
    gg[,cz:=climate_z+delta/sC];gg[,ez:=effort_z+EFF_GROW[[si]]*dec];out[[paste0(si,yy)]]<-gg[,.(unit,lon,lat,ssp=si,year=yy,hazard=pr(cz,ez))]}
  res<-rbindlist(out);res[,panel:=factor(ifelse(ssp=="Current","Current (2024)",paste0(ssp," — ",year)),levels=c("Current (2024)",paste0("SSP245 — ",PYEAR),paste0("SSP585 — ",PYEAR)))]
  qc<-quantile(res$hazard,0.98,na.rm=TRUE);res[,hz:=pmin(hazard,qc)];res}
subtitle<-sprintf("v3 relaxed-set province M4 (temp_anom_z × visits, HR=%.2f) plug-in to grid-native CRU climate + merged effort. Futures: 4-GCM CMIP6 warming + SSP-differentiated effort growth.",exp(bI))
map_tile<-function(g,label,fig,tile){res<-predict_long(g)
  p<-ggplot()+geom_sf(data=prov,fill="grey98",colour="grey80",linewidth=0.12)+
    geom_tile(data=res,aes(lon,lat,fill=hz),width=tile,height=tile)+
    geom_sf(data=prov,fill=NA,colour="grey35",linewidth=0.15)+geom_sf(data=nation,fill=NA,colour="grey20",linewidth=0.22)+geom_sf(data=dash,colour="grey15",linewidth=0.4)+
    scale_fill_gradientn(colours=COOL_WARM,name="New-record risk (relative)",labels=scales::label_number(accuracy=0.01))+
    facet_wrap(~panel,ncol=4)+coord_sf(xlim=XLIM,ylim=YLIM,expand=TRUE,datum=NA)+
    labs(title=sprintf("%s plug-in new-record risk (v3 relaxed set) — cool (low) to warm (high)",label),subtitle=subtitle,x=NULL,y=NULL)+theme_map()
  ggsave(file.path(V2,paste0(fig,".pdf")),p,width=25,height=15,units="cm",device=grDevices::cairo_pdf)
  ggsave(file.path(V2,paste0(fig,".png")),p,width=25,height=15,units="cm",dpi=300);log("wrote",fig)}

# ---- 100km grid panel ----
xy100<-fread(file.path(V2,"data/derived/community_grid_100km_climate_native.csv"))[,.(unit=grid_id,lon,lat)]
clim100<-as.data.table(read_parquet(file.path(V2,"data/derived/grid_100km_climate_cru.parquet")))[year%in%2020:2024,.(temp_anom=mean(temp_anom,na.rm=TRUE)),by=.(unit=grid_id)]
eff100<-as.data.table(read_parquet(file.path(V2,"data/derived/grid_100km_effort_combined.parquet")))[year%in%2020:2024,.(n_events=mean(n_events,na.rm=TRUE)),by=.(unit=grid_id)]
g100<-Reduce(function(a,b)merge(a,b,by="unit"),list(xy100,clim100,eff100))
fwrite(predict_long(g100),file.path(V2,"results/forecasts/table_grid_100km_plugin_v3.csv"))
map_tile(g100,"100 km grid","figures/main/Figure_8_v3_grid100_plugin_hazard",1.0)

# ---- 50km grid panel (from 10km Combined) ----
e10<-fread(file.path(NP,"results_v2/table_effort_by_grid_year_source_10km.csv"))[source_short=="Combined"]
g10<-as.data.table(sf::st_drop_geometry(readRDS(file.path(NP,"data/derived_v2/china_grid_10km_v2.rds"))))[,.(grid_cell,centroid_lon,centroid_lat)]
e10<-merge(e10,g10,by="grid_cell");xyA<-sf::sf_project("EPSG:4326","EPSG:4524",as.matrix(e10[,.(centroid_lon,centroid_lat)]))
e10[,`:=`(ax=xyA[,1],ay=xyA[,2])];cmm<-50000;x0<-floor(min(e10$ax)/cmm)*cmm;y0<-floor(min(e10$ay)/cmm)*cmm
e10[,`:=`(ix=floor((ax-x0)/cmm),iy=floor((ay-y0)/cmm))];e10[,unit:=paste(ix,iy,sep="_")]
eff50<-e10[year%in%2020:2024,.(n_events=mean(n_events,na.rm=TRUE)),by=.(unit,ix,iy)]
ccc<-sf::sf_project("EPSG:4524","EPSG:4326",as.matrix(eff50[,.(x0+(ix+0.5)*cmm,y0+(iy+0.5)*cmm)]));eff50[,`:=`(lon=ccc[,1],lat=ccc[,2])]
P50<-vect(as.data.frame(eff50[,.(lon,lat)]),geom=c("lon","lat"),crs="EPSG:4326")
rc<-tryCatch(rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc"),subds="tmp"),error=function(e)rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc")))
idx<-function(y,m)(y-1901)*12+m;af<-function(y){li<-idx(y,1):idx(y,12);li<-li[li<=nlyr(rc)];if(length(li)<12)return(rep(NA,nrow(eff50)));terra::extract(mean(rc[[li]],na.rm=TRUE),P50)[,2]}
eff50[,temp_anom:=rowMeans(sapply(2020:2024,af),na.rm=TRUE)-rowMeans(sapply(2002:2010,af),na.rm=TRUE)]
map_tile(eff50[,.(unit,lon,lat,temp_anom,n_events)],"50 km grid","figures/main/Figure_8b_v3_grid50_plugin_hazard",0.55)
log("DONE")
