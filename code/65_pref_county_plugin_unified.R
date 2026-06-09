# ============================================================
# Script: 65_pref_county_plugin_unified.R —— 任务3: 市/县统一为 plug-in
# 设计: 与网格(code 64)完全同口径。只省级拟合 temp_anom_z×visits; 市/县用各自
#   单元的 CRU 时变温度异常(质心提取) × 合并社区努力(comm_n_events) 做 plug-in
#   预测; 未来用 4-GCM CMIP6 集成升温。替代现稿市/县【重拟合】(Table 4)。
# 输入: 市/县 shp(GS2019); data/derived/unit_effort_{prefecture,county}.csv;
#       CRU TS tmp .nc; WorldClim 历史/CMIP6 bio1; 省级风险集+气候。
# 输出: results/forecasts/table_{prefecture,county}_plugin_unified.csv
#       figures/main/Figure_6_prefecture_plugin_unified.{pdf,png}
#       figures/main/Figure_7_county_plugin_unified.{pdf,png}
#       results/tables/table_multiscale_plugin_summary.csv  (省→市→县→网格 均值hazard)
# 运行: Rscript --no-init-file code/65_pref_county_plugin_unified.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(glmmTMB); library(terra); library(sf); library(ggplot2); library(viridisLite)
})
sf::sf_use_s2(FALSE); options(warn=1)
V2<-normalizePath(".",mustWork=TRUE); NP<-"/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis"
WC_HIST<-"/Users/dingchenchen/Documents/New project/bird_grid_community_analysis/data/external/worldclim/climate/wc2.1_10m/wc2.1_10m_bio_1.tif"
CMIP_DIR<-file.path(V2,"data/external/cmip6_worldclim"); SHP<-file.path(V2,"data","spatial","basemap_GS2019_1822")
log<-function(...) cat(sprintf("[65 %s] ",format(Sys.time(),"%H:%M:%S")),...,"\n")
inv_cloglog<-function(eta) 1-exp(-exp(eta))
GCMS<-c("ACCESS-CM2","MPI-ESM1-2-HR","MIROC6","UKESM1-0-LL"); SSPS<-c("ssp245","ssp585"); PERIODS<-c("2021-2040","2041-2060","2061-2080"); PYEAR<-c(2030,2050,2080)

# province model
r<-fread(file.path(V2,"data/raw/hazard_risk_upgraded_complete_case.csv"))
cm<-fread(file.path(V2,"data/raw/climate_metrics_province_year.csv"))
d<-merge(r,cm[,.(province,year,temp_anom_z)],by=c("province","year"),all.x=TRUE)[is.finite(temp_anom_z)]
m4<-glmmTMB(event~temp_anom_z*log_effort_visits_z+(1|species)+(1|province),data=d,family=binomial("cloglog"))
cf<-fixef(m4)$cond; b0<-cf[1];bC<-cf["temp_anom_z"];bE<-cf["log_effort_visits_z"];bI<-cf["temp_anom_z:log_effort_visits_z"]
log(sprintf("province M4 int HR=%.3f",exp(bI)))
bio_hist<-rast(WC_HIST)
prov<-st_make_valid(st_transform(st_read(file.path(SHP,"省（等积投影）.shp"),quiet=TRUE),4326))

ensemble_delta<-function(pts){ # pts: data.table(unit,lon,lat)
  P<-vect(as.data.frame(pts[,.(lon,lat)]),geom=c("lon","lat"),crs="EPSG:4326"); h0<-terra::extract(bio_hist,P)[,2]
  out<-list()
  for(si in seq_along(SSPS)) for(pi in seq_along(PERIODS)){
    dl<-matrix(NA_real_,nrow(pts),length(GCMS))
    for(gi in seq_along(GCMS)){f<-file.path(CMIP_DIR,sprintf("wc2.1_10m_bioc_%s_%s_%s.tif",GCMS[gi],SSPS[si],PERIODS[pi])); if(file.exists(f)) dl[,gi]<-terra::extract(rast(f,lyrs=1),P)[,2]-h0}
    out[[paste0(si,pi)]]<-data.table(unit=pts$unit,ssp=toupper(SSPS[si]),year=PYEAR[pi],delta=apply(dl,1,median,na.rm=TRUE))}
  rbindlist(out)}

run_unit<-function(shp_file,idpref,eff_file,idcol,label,fig,tbl,tile,pad=4){
  s<-st_make_valid(st_transform(st_read(shp_file,quiet=TRUE),4326))
  s$unit<-paste0(idpref,sprintf(paste0("%0",pad,"d"),seq_len(nrow(s))))
  ct<-st_coordinates(st_centroid(s)); cen<-data.table(unit=s$unit,lon=ct[,1],lat=ct[,2])
  # effort: comm_n_events 2020-2024 mean per unit
  eff<-fread(eff_file); setnames(eff,idcol,"unit",skip_absent=TRUE)
  eb<-eff[year%in%2020:2024,.(n_events=mean(comm_n_events,na.rm=TRUE)),by=unit]
  # CRU temp_anom at centroids
  P<-vect(as.data.frame(cen[,.(lon,lat)]),geom=c("lon","lat"),crs="EPSG:4326")
  rc<-tryCatch(rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc"),subds="tmp"),error=function(e)rast(file.path(NP,"data/external/cru_ts/cru_ts4.09.1901.2024.tmp.dat.nc")))
  idx<-function(y,m)(y-1901)*12+m; annf<-function(y){li<-idx(y,1):idx(y,12);li<-li[li<=nlyr(rc)];if(length(li)<12)return(rep(NA,nrow(cen)));terra::extract(mean(rc[[li]],na.rm=TRUE),P)[,2]}
  t20<-rowMeans(sapply(2020:2024,annf),na.rm=TRUE); tb<-rowMeans(sapply(2002:2010,annf),na.rm=TRUE)
  cen[,temp_anom:=t20-tb]
  g<-merge(cen,eb,by="unit"); g<-g[is.finite(temp_anom)&is.finite(n_events)&is.finite(lon)&is.finite(lat)]
  mC<-mean(g$temp_anom);sC<-sd(g$temp_anom);mE<-mean(log1p(g$n_events));sE<-sd(log1p(g$n_events))
  g[,climate_z:=(temp_anom-mC)/sC];g[,effort_z:=(log1p(n_events)-mE)/sE]
  pr<-function(cz,ez) inv_cloglog(b0+bC*cz+bE*ez+bI*cz*ez)
  ed<-ensemble_delta(g[,.(unit,lon,lat)])
  out<-list(g[,.(unit,lon,lat,ssp="Current",year=2024L,hazard=pr(climate_z,effort_z))])
  for(si in c("SSP245","SSP585")) for(yy in PYEAR){dd<-ed[ssp==si&year==yy];gg<-merge(g,dd[,.(unit,delta)],by="unit");out[[paste0(si,yy)]]<-gg[,.(unit,lon,lat,ssp=si,year=yy,hazard=pr(climate_z+delta/sC,effort_z))]}
  res<-rbindlist(out); fwrite(res,file.path(V2,tbl))
  cur<-mean(out[[1]]$hazard,na.rm=TRUE); f50<-mean(out[["SSP5852050"]]$hazard,na.rm=TRUE)
  log(sprintf("%s: units=%d Current mean=%.4f SSP585-2050 mean=%.4f",label,nrow(g),cur,f50))
  res[,panel:=factor(ifelse(ssp=="Current","Current (2024)",paste0(ssp," — ",year)),levels=c("Current (2024)",paste0("SSP245 — ",PYEAR),paste0("SSP585 — ",PYEAR)))]
  qc<-quantile(res$hazard,0.99,na.rm=TRUE);res[,hz:=pmin(hazard,qc)]
  p<-ggplot()+geom_sf(data=prov,fill="grey97",colour="grey80",linewidth=0.1)+geom_tile(data=res,aes(lon,lat,fill=hz),width=tile,height=tile)+
    scale_fill_viridis_c(option="mako",direction=-1,name="Relative hazard")+facet_wrap(~panel,ncol=4)+coord_sf(xlim=c(73,135),ylim=c(17,54),expand=FALSE,datum=NA)+
    labs(title=sprintf("%s plug-in hazard — CRU climate × merged effort, CMIP6 ensemble futures",label),
         subtitle="Province M4 (temp_anom_z × visits) applied to unit-native CRU anomaly + merged eBird/GBIF+China-Birdwatch effort; futures = 4-GCM CMIP6 median warming. Plug-in projection.",
         x=NULL,y=NULL)+theme_bw(base_size=8)+theme(legend.position="bottom",panel.grid=element_blank(),plot.title=element_text(face="bold",size=9.5),plot.subtitle=element_text(size=6.4,colour="grey35"),strip.text=element_text(face="bold",size=7))
  ggsave(file.path(V2,paste0(fig,".pdf")),p,width=24,height=14,units="cm",device=grDevices::cairo_pdf)
  ggsave(file.path(V2,paste0(fig,".png")),p,width=24,height=14,units="cm",dpi=300); log("wrote",fig)
  data.table(scale=label,n_units=nrow(g),cur_mean=round(cur,4),ssp585_2050_mean=round(f50,4))
}
s1<-run_unit(file.path(SHP,"市（等积投影）.shp"),"PREF_",file.path(V2,"data/derived/unit_effort_prefecture.csv"),"pref_id","Prefecture (市)","figures/main/Figure_6_prefecture_plugin_unified","results/forecasts/table_prefecture_plugin_unified.csv",0.9)
s2<-run_unit(file.path(SHP,"县（等积投影）.shp"),"CNTY_",file.path(V2,"data/derived/unit_effort_county.csv"),"cnty_id","County (县)","figures/main/Figure_7_county_plugin_unified","results/forecasts/table_county_plugin_unified.csv",0.5,pad=5)
fwrite(rbind(s1,s2),file.path(V2,"results/tables/table_multiscale_plugin_summary.csv")); print(rbind(s1,s2)); log("DONE")
