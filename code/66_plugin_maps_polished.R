# ============================================================
# Script: 66_plugin_maps_polished.R —— 预测图专业精修(统一重做 Fig 6/7/8/8b)
# 改进:
#  (1) 未来加【调查努力增长情景】: effort_z 随年代上升(+EFF_GROW SD/decade),
#      不再冻结 2024; 与 CMIP6 气候同时演进。
#  (2) 市/县用【真实矢量多边形 + 边界】(geom_sf), 不再用质心方块。
#  (3) 配色【冷→暖】(低风险冷蓝 → 高风险暖红, RdYlBu 反向)。
#  (4) 地图四周【留边距】(coord_sf 缓冲 + expand + panel/ plot margin), 不贴边框。
#  (5) 统一专业排版主题。
# 模型: 省级 temp_anom_z × visits (HR=1.273) plug-in; 网格/单元原生 CRU 气候 +
#       合并努力; 未来 = 4-GCM CMIP6 集成升温 + 努力增长情景。
# 输出(覆盖): figures/main/Figure_6_prefecture_plugin_unified.{pdf,png}
#             figures/main/Figure_7_county_plugin_unified.{pdf,png}
#             figures/main/Figure_8_grid100_native_plugin_hazard.{pdf,png}
#             figures/main/Figure_8b_grid50_native_plugin_hazard.{pdf,png}
# 运行: Rscript --no-init-file code/66_plugin_maps_polished.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(glmmTMB); library(terra); library(sf); library(ggplot2); library(arrow)
})
sf::sf_use_s2(FALSE); options(warn=1)
V2<-normalizePath(".",mustWork=TRUE); NP<-"/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis"
WC_HIST<-"/Users/dingchenchen/Documents/New project/bird_grid_community_analysis/data/external/worldclim/climate/wc2.1_10m/wc2.1_10m_bio_1.tif"
CMIP_DIR<-file.path(V2,"data/external/cmip6_worldclim"); SHP<-file.path(V2,"data","spatial","basemap_GS2019_1822")
log<-function(...) cat(sprintf("[66 %s] ",format(Sys.time(),"%H:%M:%S")),...,"\n")
inv_cloglog<-function(eta) 1-exp(-exp(eta))
GCMS<-c("ACCESS-CM2","MPI-ESM1-2-HR","MIROC6","UKESM1-0-LL"); SSPS<-c("ssp245","ssp585")
PERIODS<-c("2021-2040","2041-2060","2061-2080"); PYEAR<-c(2030,2050,2080)
# SSP-differentiated survey-effort growth scenarios (SD/decade): effort tracks
# socioeconomic development — SSP2 (middle-of-road) moderate, SSP5 (fossil-fuelled,
# high human-capital/connectivity) faster citizen-science growth.
EFF_GROW<-c(SSP245=0.30, SSP585=0.60)
PROV<-c("北京市"="Beijing","天津市"="Tianjin","河北省"="Hebei","山西省"="Shanxi","内蒙古自治区"="Inner Mongolia","辽宁省"="Liaoning","吉林省"="Jilin","黑龙江省"="Heilongjiang","上海市"="Shanghai","江苏省"="Jiangsu","浙江省"="Zhejiang","安徽省"="Anhui","福建省"="Fujian","江西省"="Jiangxi","山东省"="Shandong","河南省"="Henan","湖北省"="Hubei","湖南省"="Hunan","广东省"="Guangdong","广西壮族自治区"="Guangxi","海南省"="Hainan","重庆市"="Chongqing","四川省"="Sichuan","贵州省"="Guizhou","云南省"="Yunnan","西藏自治区"="Tibet","陕西省"="Shaanxi","甘肃省"="Gansu","青海省"="Qinghai","宁夏回族自治区"="Ningxia","新疆维吾尔自治区"="Xinjiang","台湾省"="Taiwan","香港特别行政区"="Hong Kong","澳门特别行政区"="Macau")

# ----- cool->warm palette + polished theme + framed coord -----
COOL_WARM <- c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF",
               "#FEE090","#FDAE61","#F46D43","#D73027","#A50026")  # RdYlBu reversed (cool->warm)
XLIM<-c(71,137); YLIM<-c(2.5,55)  # extend south to show the nine-dash line (南海段线); margin to frame
theme_map <- function() theme_bw(base_size=8) +
  theme(legend.position="bottom", panel.grid=element_blank(),
        panel.spacing=unit(0.45,"lines"),
        plot.margin=margin(8,10,8,10),
        plot.title=element_text(face="bold",size=10),
        plot.subtitle=element_text(size=7,colour="grey35",margin=margin(b=4)),
        strip.background=element_rect(fill="grey96",colour="grey80"),
        strip.text=element_text(face="bold",size=7.5),
        legend.key.width=unit(1.1,"cm"), legend.key.height=unit(0.3,"cm"),
        legend.title=element_text(size=8,face="bold"))

# ----- province plug-in model -----
r<-fread(file.path(V2,"data/raw/hazard_risk_upgraded_complete_case.csv"))
cm<-fread(file.path(V2,"data/raw/climate_metrics_province_year.csv"))
d<-merge(r,cm[,.(province,year,temp_anom_z)],by=c("province","year"),all.x=TRUE)[is.finite(temp_anom_z)]
m4<-glmmTMB(event~temp_anom_z*log_effort_visits_z+(1|species)+(1|province),data=d,family=binomial("cloglog"))
cf<-fixef(m4)$cond; b0<-cf[1];bC<-cf["temp_anom_z"];bE<-cf["log_effort_visits_z"];bI<-cf["temp_anom_z:log_effort_visits_z"]
log(sprintf("province M4 int HR=%.3f",exp(bI)))
bio_hist<-rast(WC_HIST)
prov<-st_make_valid(st_transform(st_read(file.path(SHP,"省（等积投影）.shp"),quiet=TRUE),4326))
# official GS(2019)1822 basemap elements incl. nine-dash line (南海段线) for compliance
nation<-st_make_valid(st_transform(st_read(file.path(SHP,"中国轮廓线.shp"),quiet=TRUE),4326))
dash<-st_transform(st_read(file.path(SHP,"九段线.shp"),quiet=TRUE),4326)
basemap_layers<-function() list(
  geom_sf(data=nation,fill=NA,colour="grey20",linewidth=0.22),
  geom_sf(data=dash,colour="grey15",linewidth=0.4))

ensemble_delta<-function(pts){ # pts: data.table(unit,lon,lat)
  P<-vect(as.data.frame(pts[,.(lon,lat)]),geom=c("lon","lat"),crs="EPSG:4326"); h0<-terra::extract(bio_hist,P)[,2]
  out<-list()
  for(si in seq_along(SSPS)) for(pi in seq_along(PERIODS)){
    dl<-matrix(NA_real_,nrow(pts),length(GCMS))
    for(gi in seq_along(GCMS)){f<-file.path(CMIP_DIR,sprintf("wc2.1_10m_bioc_%s_%s_%s.tif",GCMS[gi],SSPS[si],PERIODS[pi])); if(file.exists(f)) dl[,gi]<-terra::extract(rast(f,lyrs=1),P)[,2]-h0}
    out[[paste0(si,pi)]]<-data.table(unit=pts$unit,ssp=toupper(SSPS[si]),year=PYEAR[pi],delta=apply(dl,1,median,na.rm=TRUE))}
  rbindlist(out)}

# build long prediction (unit, lon, lat, panel, hazard) given baseline g(unit,lon,lat,temp_anom,n_events)
predict_long<-function(g){
  g<-g[is.finite(temp_anom)&is.finite(n_events)&is.finite(lon)&is.finite(lat)]
  mC<-mean(g$temp_anom);sC<-sd(g$temp_anom);mE<-mean(log1p(g$n_events));sE<-sd(log1p(g$n_events))
  g[,climate_z:=(temp_anom-mC)/sC]; g[,effort_z:=(log1p(n_events)-mE)/sE]
  pr<-function(cz,ez) inv_cloglog(b0+bC*cz+bE*ez+bI*cz*ez)
  ed<-ensemble_delta(g[,.(unit,lon,lat)])
  out<-list(g[,.(unit,lon,lat,ssp="Current",year=2024L,hazard=pr(climate_z,effort_z))])
  for(si in c("SSP245","SSP585")) for(yy in PYEAR){
    dec<-(yy-2024)/10; dd<-ed[ssp==si&year==yy]; gg<-merge(g,dd[,.(unit,delta)],by="unit")
    gg[,cz:=climate_z+delta/sC]; gg[,ez:=effort_z+EFF_GROW[[si]]*dec]  # 气候(CMIP6)+SSP差异化努力增长
    out[[paste0(si,yy)]]<-gg[,.(unit,lon,lat,ssp=si,year=yy,hazard=pr(cz,ez))]}
  res<-rbindlist(out)
  res[,panel:=factor(ifelse(ssp=="Current","Current (2024)",paste0(ssp," — ",year)),
       levels=c("Current (2024)",paste0("SSP245 — ",PYEAR),paste0("SSP585 — ",PYEAR)))]
  qc<-quantile(res$hazard,0.98,na.rm=TRUE); res[,hz:=pmin(hazard,qc)]; res}

subtitle<-paste0("Province M4 (temp_anom_z × visits, HR=",sprintf("%.2f",exp(bI)),
  ") plug-in to native CRU climate + merged eBird/GBIF+China-Birdwatch effort. ",
  "Futures: 4-GCM CMIP6 median warming + SSP-differentiated effort growth (SSP245 +0.3, SSP585 +0.6 SD/decade).")

# ---------- polygon map (prefecture / county) ----------
map_polygon<-function(shp_file,idpref,pad,eff_file,idcol,label,fig,lw){
  s<-st_make_valid(st_transform(st_read(shp_file,quiet=TRUE),4326))
  s$unit<-paste0(idpref,sprintf(paste0("%0",pad,"d"),seq_len(nrow(s))))
  ct<-st_coordinates(st_centroid(s)); cen<-data.table(unit=s$unit,lon=ct[,1],lat=ct[,2])
  eff<-fread(eff_file); setnames(eff,idcol,"unit",skip_absent=TRUE)
  eb<-eff[year%in%2020:2024,.(n_events=mean(comm_n_events,na.rm=TRUE)),by=unit]
  P<-vect(as.data.frame(cen[,.(lon,lat)]),geom=c("lon","lat"),crs="EPSG:4326")
  rc<-tryCatch(rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc"),subds="tmp"),error=function(e)rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc")))
  idx<-function(y,m)(y-1901)*12+m; af<-function(y){li<-idx(y,1):idx(y,12);li<-li[li<=nlyr(rc)];if(length(li)<12)return(rep(NA,nrow(cen)));terra::extract(mean(rc[[li]],na.rm=TRUE),P)[,2]}
  cen[,temp_anom:=rowMeans(sapply(2020:2024,af),na.rm=TRUE)-rowMeans(sapply(2002:2010,af),na.rm=TRUE)]
  g<-merge(cen,eb,by="unit")
  res<-predict_long(g)
  ps<-merge(s[,c("unit")],res,by="unit")   # sf replicated per panel
  p<-ggplot()+
    geom_sf(data=ps,aes(fill=hz),colour="grey75",linewidth=lw)+
    geom_sf(data=prov,fill=NA,colour="grey25",linewidth=0.18)+
    basemap_layers()+
    scale_fill_gradientn(colours=COOL_WARM,name="New-record risk (relative)",
                         labels=scales::label_number(accuracy=0.01))+
    facet_wrap(~panel,ncol=4)+
    coord_sf(xlim=XLIM,ylim=YLIM,expand=TRUE,datum=NA)+
    labs(title=sprintf("%s plug-in new-record risk — cool (low) to warm (high)",label),subtitle=subtitle,x=NULL,y=NULL)+
    theme_map()
  ggsave(file.path(V2,paste0(fig,".pdf")),p,width=25,height=15,units="cm",device=grDevices::cairo_pdf)
  ggsave(file.path(V2,paste0(fig,".png")),p,width=25,height=15,units="cm",dpi=300)
  log("wrote",fig,"(",nrow(s),"polygons )")}

# ---------- tile map (grid) ----------
map_tile<-function(g,label,fig,tile){
  res<-predict_long(g)
  p<-ggplot()+
    geom_sf(data=prov,fill="grey98",colour="grey80",linewidth=0.12)+
    geom_tile(data=res,aes(lon,lat,fill=hz),width=tile,height=tile)+
    geom_sf(data=prov,fill=NA,colour="grey35",linewidth=0.15)+
    basemap_layers()+
    scale_fill_gradientn(colours=COOL_WARM,name="New-record risk (relative)",
                         labels=scales::label_number(accuracy=0.01))+
    facet_wrap(~panel,ncol=4)+
    coord_sf(xlim=XLIM,ylim=YLIM,expand=TRUE,datum=NA)+
    labs(title=sprintf("%s plug-in new-record risk — cool (low) to warm (high)",label),subtitle=subtitle,x=NULL,y=NULL)+
    theme_map()
  ggsave(file.path(V2,paste0(fig,".pdf")),p,width=25,height=15,units="cm",device=grDevices::cairo_pdf)
  ggsave(file.path(V2,paste0(fig,".png")),p,width=25,height=15,units="cm",dpi=300)
  log("wrote",fig)}

# ===== prefecture / county (polygons) =====
map_polygon(file.path(SHP,"市（等积投影）.shp"),"PREF_",4,file.path(V2,"data/derived/unit_effort_prefecture.csv"),"pref_id","Prefecture (市)","figures/main/Figure_6_prefecture_plugin_unified",0.06)
map_polygon(file.path(SHP,"县（等积投影）.shp"),"CNTY_",5,file.path(V2,"data/derived/unit_effort_county.csv"),"cnty_id","County (县)","figures/main/Figure_7_county_plugin_unified",0.03)

# ===== 100km / 50km grids (tiles) =====
xy100<-fread(file.path(V2,"data/derived/community_grid_100km_climate_native.csv"))[,.(unit=grid_id,lon,lat)]
clim100<-as.data.table(read_parquet(file.path(V2,"data/derived/grid_100km_climate_cru.parquet")))[year%in%2020:2024,.(temp_anom=mean(temp_anom,na.rm=TRUE)),by=.(unit=grid_id)]
eff100<-as.data.table(read_parquet(file.path(V2,"data/derived/grid_100km_effort_combined.parquet")))[year%in%2020:2024,.(n_events=mean(n_events,na.rm=TRUE)),by=.(unit=grid_id)]
g100<-Reduce(function(a,b)merge(a,b,by="unit"),list(xy100,clim100,eff100))
map_tile(g100,"100 km grid","figures/main/Figure_8_grid100_native_plugin_hazard",1.0)

# 50km rebuild from 10km Combined
e10<-fread(file.path(NP,"results_v2/table_effort_by_grid_year_source_10km.csv"))[source_short=="Combined"]
g10<-as.data.table(sf::st_drop_geometry(readRDS(file.path(NP,"data/derived_v2/china_grid_10km_v2.rds"))))[,.(grid_cell,centroid_lon,centroid_lat)]
e10<-merge(e10,g10,by="grid_cell"); xyA<-sf::sf_project("EPSG:4326","EPSG:4524",as.matrix(e10[,.(centroid_lon,centroid_lat)]))
e10[,`:=`(ax=xyA[,1],ay=xyA[,2])]; cmm<-50000; x0<-floor(min(e10$ax)/cmm)*cmm; y0<-floor(min(e10$ay)/cmm)*cmm
e10[,`:=`(ix=floor((ax-x0)/cmm),iy=floor((ay-y0)/cmm))]; e10[,unit:=paste(ix,iy,sep="_")]
eff50<-e10[year%in%2020:2024,.(n_events=mean(n_events,na.rm=TRUE)),by=.(unit,ix,iy)]
ccc<-sf::sf_project("EPSG:4524","EPSG:4326",as.matrix(eff50[,.(x0+(ix+0.5)*cmm,y0+(iy+0.5)*cmm)])); eff50[,`:=`(lon=ccc[,1],lat=ccc[,2])]
P50<-vect(as.data.frame(eff50[,.(lon,lat)]),geom=c("lon","lat"),crs="EPSG:4326")
rc<-tryCatch(rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc"),subds="tmp"),error=function(e)rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc")))
idx<-function(y,m)(y-1901)*12+m; af<-function(y){li<-idx(y,1):idx(y,12);li<-li[li<=nlyr(rc)];if(length(li)<12)return(rep(NA,nrow(eff50)));terra::extract(mean(rc[[li]],na.rm=TRUE),P50)[,2]}
eff50[,temp_anom:=rowMeans(sapply(2020:2024,af),na.rm=TRUE)-rowMeans(sapply(2002:2010,af),na.rm=TRUE)]
g50<-eff50[,.(unit,lon,lat,temp_anom,n_events)]
map_tile(g50,"50 km grid","figures/main/Figure_8b_grid50_native_plugin_hazard",0.55)
log("DONE")
