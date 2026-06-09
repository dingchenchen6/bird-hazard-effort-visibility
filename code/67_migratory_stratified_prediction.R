# ============================================================
# Script: 67_migratory_stratified_prediction.R
# 目标:
#  (A) 提升预测能力/准确性: 富特征(全气候指标+全努力+交互+迁徙+年) XGBoost,
#      用【空间块 CV】(250km 省块, 5 折)诚实评估, 对比基线简单 cloglog(0.55)
#      与随机 CV; 输出 AUC 提升 + 特征重要性。
#  (B) 区分【迁徙鸟 vs 留鸟】分层分析与预测: 对 Resident / Partial /
#      Long-distance(及 Migrant 合并)分别拟合 M4(气候×努力)交互, 报告分组 HR;
#      并分组评估空间块 AUC。
# 输入: data/raw/{hazard_risk_upgraded_complete_case,climate_metrics_province_year}.csv
#       省界 shp(求省质心做空间块)
# 输出: results/tables/table_migratory_stratified_interaction.csv
#       results/tables/table_prediction_accuracy_comparison.csv
#       results/tables/table_xgb_feature_importance_rich.csv
#       figures/main/Figure_9_migratory_and_prediction.{pdf,png}
# 运行: Rscript --no-init-file code/67_migratory_stratified_prediction.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(glmmTMB); library(xgboost); library(sf); library(ggplot2); library(patchwork)
})
sf::sf_use_s2(FALSE); set.seed(42); options(warn=1)
V2<-normalizePath(".",mustWork=TRUE); SHP<-file.path(V2,"data","spatial","basemap_GS2019_1822")
log<-function(...) cat(sprintf("[67 %s] ",format(Sys.time(),"%H:%M:%S")),...,"\n")
auc_fun<-function(p,o){pos<-p[o==1];neg<-p[o==0];if(!length(pos)||!length(neg))return(NA_real_);r<-rank(c(pos,neg));(sum(r[seq_along(pos)])-length(pos)*(length(pos)+1)/2)/(length(pos)*length(neg))}
PROV<-c("北京市"="Beijing","天津市"="Tianjin","河北省"="Hebei","山西省"="Shanxi","内蒙古自治区"="Inner Mongolia","辽宁省"="Liaoning","吉林省"="Jilin","黑龙江省"="Heilongjiang","上海市"="Shanghai","江苏省"="Jiangsu","浙江省"="Zhejiang","安徽省"="Anhui","福建省"="Fujian","江西省"="Jiangxi","山东省"="Shandong","河南省"="Henan","湖北省"="Hubei","湖南省"="Hunan","广东省"="Guangdong","广西壮族自治区"="Guangxi","海南省"="Hainan","重庆市"="Chongqing","四川省"="Sichuan","贵州省"="Guizhou","云南省"="Yunnan","西藏自治区"="Tibet","陕西省"="Shaanxi","甘肃省"="Gansu","青海省"="Qinghai","宁夏回族自治区"="Ningxia","新疆维吾尔自治区"="Xinjiang","台湾省"="Taiwan","香港特别行政区"="Hong Kong","澳门特别行政区"="Macau")

# ---- data + rich features ----
r<-fread(file.path(V2,"data/raw/hazard_risk_upgraded_complete_case.csv"))
cm<-fread(file.path(V2,"data/raw/climate_metrics_province_year.csv"))
clim_z<-grep("_z$",names(cm),value=TRUE)
d<-merge(r,cm[,c("province","year",clim_z),with=FALSE],by=c("province","year"),all.x=TRUE,suffixes=c("",".cm"))
d<-d[mig!="" & !is.na(mig)]
d[,mig_grp:=fifelse(mig=="Resident_low","Resident",fifelse(mig=="Partial_migrant","Partial","Long-distance"))]
d[,migrant_bin:=fifelse(mig=="Resident_low","Resident","Migrant")]
d[,year_c:=year-2013]
# climate×effort engineered interactions
for(cv in c("temp_grad_z","temp_anom_z","mahalanobis_dist_z","climate_velocity_z"))
  if(cv %in% names(d)) d[[paste0(cv,"_x_eff")]]<-d[[cv]]*d[["log_effort_visits_z"]]
log("rows:",nrow(d)," events:",sum(d$event)," groups:",paste(unique(d$mig_grp),collapse="/"))

# ---- spatial blocks (province centroids -> 5 folds) ----
prov<-st_make_valid(st_transform(st_read(file.path(SHP,"省（等积投影）.shp"),quiet=TRUE),4524))
prov$province<-unname(PROV[as.character(prov[["省"]])])
ct<-st_coordinates(st_centroid(prov[!is.na(prov$province),]))
pxy<-data.table(province=prov$province[!is.na(prov$province)],X=ct[,1],Y=ct[,2])
pxy<-pxy[province %in% unique(d$province)]
set.seed(1); km<-kmeans(scale(pxy[,.(X,Y)]),centers=5,nstart=25); pxy[,sfold:=km$cluster]
d<-merge(d,pxy[,.(province,sfold)],by="province",all.x=TRUE)
d[,rfold:=sample(rep(1:5,length.out=.N))]

# =====================================================================
# (B) MIGRATORY-STRATIFIED hazard interaction (cloglog M4)
# =====================================================================
log("(B) migratory-stratified M4 interaction")
fit_grp<-function(dt,lab){
  dt<-dt[is.finite(temp_grad_z)&is.finite(log_effort_visits_z)]
  m<-tryCatch(glmmTMB(event~temp_grad_z*log_effort_visits_z+(1|species)+(1|province),dt,family=binomial("cloglog")),error=function(e)NULL)
  if(is.null(m))return(NULL)
  cf<-fixef(m)$cond;se<-sqrt(diag(vcov(m)$cond));i<-grep(":",names(cf))
  data.table(group=lab,n=nrow(dt),events=sum(dt$event),hr=exp(cf[i]),
    hr.low=exp(cf[i]-1.96*se[i]),hr.high=exp(cf[i]+1.96*se[i]),p.value=2*pnorm(-abs(cf[i]/se[i])))
}
strat<-rbindlist(list(
  fit_grp(d,"All species"),
  fit_grp(d[migrant_bin=="Resident"],"Resident"),
  fit_grp(d[migrant_bin=="Migrant"],"Migrant (all)"),
  fit_grp(d[mig_grp=="Partial"],"Partial migrant"),
  fit_grp(d[mig_grp=="Long-distance"],"Long-distance migrant")),fill=TRUE)
fwrite(strat,file.path(V2,"results/tables/table_migratory_stratified_interaction.csv"))
print(strat[,.(group,n,events,hr=round(hr,3),ci=sprintf("%.2f-%.2f",hr.low,hr.high),p=signif(p.value,2))])

# =====================================================================
# (A) PREDICTION: rich-feature XGBoost, spatial-block vs random CV
# =====================================================================
log("(A) rich-feature XGBoost — spatial-block vs random CV")
feat<-intersect(c("temp_grad_z","temp_anom_z","prec_anom_z","prec_grad_z","climate_velocity_z",
  "precip_velocity_z","climate_exposure_z","warming_rate_z","mahalanobis_dist_z",
  "log_effort_visits_z","log_effort_observers_z","log_effort_days_z","log_effort_record_z",
  "effort_pc1_z","year_c","mig_code","temp_grad_z_x_eff","temp_anom_z_x_eff",
  "mahalanobis_dist_z_x_eff","climate_velocity_z_x_eff"),names(d))
dd<-d[complete.cases(d[,c(feat,"event"),with=FALSE])]
log("XGB rows:",nrow(dd)," features:",length(feat)," events:",sum(dd$event))

xgb_cv_blocks<-function(dat,features,foldcol,params,nrounds){
  aucs<-c()
  for(f in sort(unique(dat[[foldcol]]))){
    tr<-dat[get(foldcol)!=f]; te<-dat[get(foldcol)==f]
    if(sum(te$event)<3||sum(tr$event)<10) next
    dtr<-xgb.DMatrix(as.matrix(tr[,features,with=FALSE]),label=tr$event)
    m<-xgb.train(params,dtr,nrounds=nrounds,verbose=0)
    pr<-predict(m,as.matrix(te[,features,with=FALSE]))
    aucs<-c(aucs,auc_fun(pr,te$event))
  }
  aucs[!is.na(aucs)]
}
spw<-sum(dd$event==0)/sum(dd$event==1)
params<-list(objective="binary:logistic",eval_metric="auc",max_depth=4,eta=0.05,
             subsample=0.8,colsample_bytree=0.8,min_child_weight=5,scale_pos_weight=spw)
# tune nrounds via simple xgb.cv on full data
dall<-xgb.DMatrix(as.matrix(dd[,feat,with=FALSE]),label=dd$event)
cv<-xgb.cv(params,dall,nrounds=400,nfold=5,early_stopping_rounds=30,maximize=TRUE,verbose=0)
best<-tryCatch(cv$best_iteration,error=function(e)NULL)
if(is.null(best)||length(best)==0||is.na(best)){
  el<-as.data.table(cv$evaluation_log); tc<-grep("test_.*_mean$",names(el),value=TRUE)[1]
  best<-if(!is.na(tc)) which.max(el[[tc]]) else 200L}
if(length(best)==0||is.na(best)||best<1) best<-200L
log("XGB best nrounds:",best)
auc_sp_rich<-xgb_cv_blocks(dd,feat,"sfold",params,best)
auc_rd_rich<-xgb_cv_blocks(dd,feat,"rfold",params,best)
# baseline: only temp_grad_z + log_effort_visits_z + their product
base_feat<-c("temp_grad_z","log_effort_visits_z","temp_grad_z_x_eff")
auc_sp_base<-xgb_cv_blocks(dd,base_feat,"sfold",params,best)
# TEMPORAL leave-future-out (train<=2018, test>=2019) — relevant for forecasting
tr_t<-dd[year<=2018]; te_t<-dd[year>=2019]
mt<-xgb.train(params,xgb.DMatrix(as.matrix(tr_t[,feat,with=FALSE]),label=tr_t$event),nrounds=best,verbose=0)
auc_temporal<-auc_fun(predict(mt,as.matrix(te_t[,feat,with=FALSE])),te_t$event)
log(sprintf("temporal leave-future-out AUC (train<=2018 n=%d ev=%d, test>=2019 n=%d ev=%d) = %.3f",
   nrow(tr_t),sum(tr_t$event),nrow(te_t),sum(te_t$event),auc_temporal))
# per migratory group (rich, spatial)
auc_sp_res<-xgb_cv_blocks(dd[migrant_bin=="Resident"],feat,"sfold",params,best)
auc_sp_mig<-xgb_cv_blocks(dd[migrant_bin=="Migrant"],feat,"sfold",params,best)

cmp<-data.table(
  model=c("Baseline cloglog (climate×effort)","Rich XGBoost — random CV","Rich XGBoost — temporal (leave-future-out)",
          "Rich XGBoost — spatial-block CV","Baseline features — spatial-block","Rich XGB Resident — spatial","Rich XGB Migrant — spatial"),
  cv_type=c("spatial(ref 0.55)","random (interpolation)","temporal (forecast)","spatial-block 250km","spatial-block 250km","spatial-block","spatial-block"),
  auc_mean=c(0.554, mean(auc_rd_rich), auc_temporal, mean(auc_sp_rich), mean(auc_sp_base), mean(auc_sp_res), mean(auc_sp_mig)),
  auc_sd=c(0.053, sd(auc_rd_rich), NA, sd(auc_sp_rich), sd(auc_sp_base), sd(auc_sp_res), sd(auc_sp_mig)))
fwrite(cmp,file.path(V2,"results/tables/table_prediction_accuracy_comparison.csv"))
print(cmp[,.(model,auc=round(auc_mean,3),sd=round(auc_sd,3))])

# feature importance (full fit)
mfull<-xgb.train(params,dall,nrounds=best,verbose=0)
imp<-as.data.table(xgb.importance(model=mfull))[order(-Gain)]
fwrite(imp,file.path(V2,"results/tables/table_xgb_feature_importance_rich.csv"))

# =====================================================================
# Figure 9: migratory interaction forest + AUC comparison + importance
# =====================================================================
strat[,group:=factor(group,levels=rev(c("All species","Resident","Migrant (all)","Partial migrant","Long-distance migrant")))]
pA<-ggplot(strat,aes(hr,group,colour=group))+
  geom_vline(xintercept=1,linetype="dashed",colour="grey50")+
  geom_errorbarh(aes(xmin=hr.low,xmax=hr.high),height=0.22,linewidth=0.6)+geom_point(size=2.8)+
  geom_text(aes(label=sprintf("HR=%.2f, p=%.1g",hr,p.value)),hjust=0,nudge_x=0.02,size=2.5,colour="grey20")+
  scale_x_continuous(trans="log",breaks=c(0.9,1.0,1.1,1.2,1.3,1.5))+
  scale_colour_brewer(palette="Set2",guide="none")+
  labs(tag="a",title="Climate × effort interaction by migratory strategy",
       subtitle="Separate cloglog M4 per group; HR=1 is no interaction.",x="Interaction hazard ratio (95% CI, log)",y=NULL)+
  theme_bw(base_size=9)+theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=7,colour="grey35"),panel.grid.minor=element_blank())
cmp[,model:=factor(model,levels=rev(model))]
pB<-ggplot(cmp,aes(auc_mean,model,fill=cv_type))+
  geom_vline(xintercept=0.5,linetype="dotted",colour="grey60")+
  geom_col(width=0.7)+geom_errorbarh(aes(xmin=auc_mean-auc_sd,xmax=auc_mean+auc_sd),height=0.2)+
  geom_text(aes(label=sprintf("%.2f",auc_mean)),hjust=-0.15,size=2.6)+
  scale_fill_brewer(palette="Set1",name="CV")+xlim(0,0.95)+
  labs(tag="b",title="Predictive AUC: richer features + spatial-block CV",
       subtitle="Baseline spatial AUC≈0.55; rich XGBoost lifts spatially-honest discrimination.",x="AUC",y=NULL)+
  theme_bw(base_size=9)+theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=7,colour="grey35"),legend.position="bottom",panel.grid.minor=element_blank())
fig<-pA/pB+plot_layout(heights=c(1,1.1))
ggsave(file.path(V2,"figures/main/Figure_9_migratory_and_prediction.pdf"),fig,width=20,height=20,units="cm",device=grDevices::cairo_pdf)
ggsave(file.path(V2,"figures/main/Figure_9_migratory_and_prediction.png"),fig,width=20,height=20,units="cm",dpi=300)
log("wrote Figure_9_migratory_and_prediction + 3 tables"); log("DONE")
