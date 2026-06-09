# ============================================================
# Script: 68_v3_full_parallel_analysis.R
# 目标: 在【v3 宽松集 188,870 行 / 463 种 / 817 事件】上做与 v2 平行的【完整分析】:
#   富集 v3 风险集(join 物种→mig + 省×年气候/努力指标) 后运行:
#   (1) 内生性: 努力滞后1年(省×年) 重拟合 M4 交互
#   (2) 迁徙分层: Resident/Migrant/Partial/Long-distance 交互 HR (mig 覆盖子集)
#   (3) 预测三轴: 富特征 XGBoost 的 随机/时间/空间块 AUC + 特征重要性
#   (头条 4-spec、M5 offset、多尺度 已在既有 v3 表中)
# 输入: data/derived/risk_set_province_v3.csv (province,year,species,event,temp_grad_z,log_effort_visits_z)
#       data/raw/{hazard_risk_upgraded_complete_case (species->mig + 省×年其他努力),
#                 climate_metrics_province_year} ; 省界 shp
# 输出: results/tables/table_v3_endogeneity_lag.csv
#       results/tables/table_v3_migratory_stratified.csv
#       results/tables/table_v3_prediction_accuracy.csv
#       results/tables/table_v3_xgb_importance.csv
#       figures/main/Figure_10_v3_full_analysis.{pdf,png}
# 运行: Rscript --no-init-file code/68_v3_full_parallel_analysis.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(glmmTMB); library(xgboost); library(sf); library(ggplot2); library(patchwork)
})
sf::sf_use_s2(FALSE); set.seed(42); options(warn=1)
V2<-normalizePath(".",mustWork=TRUE); SHP<-file.path(V2,"data","spatial","basemap_GS2019_1822")
log<-function(...) cat(sprintf("[68 %s] ",format(Sys.time(),"%H:%M:%S")),...,"\n")
auc_fun<-function(p,o){pos<-p[o==1];neg<-p[o==0];if(!length(pos)||!length(neg))return(NA_real_);r<-rank(c(pos,neg));(sum(r[seq_along(pos)])-length(pos)*(length(pos)+1)/2)/(length(pos)*length(neg))}
PROV<-c("北京市"="Beijing","天津市"="Tianjin","河北省"="Hebei","山西省"="Shanxi","内蒙古自治区"="Inner Mongolia","辽宁省"="Liaoning","吉林省"="Jilin","黑龙江省"="Heilongjiang","上海市"="Shanghai","江苏省"="Jiangsu","浙江省"="Zhejiang","安徽省"="Anhui","福建省"="Fujian","江西省"="Jiangxi","山东省"="Shandong","河南省"="Henan","湖北省"="Hubei","湖南省"="Hunan","广东省"="Guangdong","广西壮族自治区"="Guangxi","海南省"="Hainan","重庆市"="Chongqing","四川省"="Sichuan","贵州省"="Guizhou","云南省"="Yunnan","西藏自治区"="Tibet","陕西省"="Shaanxi","甘肃省"="Gansu","青海省"="Qinghai","宁夏回族自治区"="Ningxia","新疆维吾尔自治区"="Xinjiang","台湾省"="Taiwan","香港特别行政区"="Hong Kong","澳门特别行政区"="Macau")

# ---- load + enrich v3 ----
v3<-fread(file.path(V2,"data/derived/risk_set_province_v3.csv"))
v2<-fread(file.path(V2,"data/raw/hazard_risk_upgraded_complete_case.csv"))
cm<-fread(file.path(V2,"data/raw/climate_metrics_province_year.csv"))
sm<-unique(v2[mig!=""&!is.na(mig),.(species,mig)])
eff_py<-unique(v2[,.(province,year,log_effort_observers_z,log_effort_days_z,log_effort_record_z,effort_pc1_z)])
clim_z<-grep("_z$",names(cm),value=TRUE)
v3<-merge(v3,sm,by="species",all.x=TRUE)
v3<-merge(v3,eff_py,by=c("province","year"),all.x=TRUE)
v3<-merge(v3,cm[,c("province","year",setdiff(clim_z,"temp_grad_prov_z")),with=FALSE],by=c("province","year"),all.x=TRUE)
v3[,migrant_bin:=fifelse(mig=="Resident_low","Resident",fifelse(mig%in%c("Partial_migrant","Long_distance_migrant"),"Migrant",NA_character_))]
v3[,mig_grp:=fifelse(mig=="Resident_low","Resident",fifelse(mig=="Partial_migrant","Partial",fifelse(mig=="Long_distance_migrant","Long-distance",NA_character_)))]
v3[,mig_code:=fifelse(mig=="Resident_low",0L,fifelse(mig=="Partial_migrant",1L,fifelse(mig=="Long_distance_migrant",2L,NA_integer_)))]
v3[,year_c:=year-2013]
for(cv in c("temp_grad_z","temp_anom_z","mahalanobis_dist_z","climate_velocity_z"))
  if(cv%in%names(v3)) v3[[paste0(cv,"_x_eff")]]<-v3[[cv]]*v3[["log_effort_visits_z"]]
log("v3 enriched:",nrow(v3),"rows,",sum(v3$event),"events,",uniqueN(v3$species),"species; mig-covered events:",v3[!is.na(mig)&event==1,.N])

# ---- spatial blocks ----
prov<-st_make_valid(st_transform(st_read(file.path(SHP,"省（等积投影）.shp"),quiet=TRUE),4524)); prov$province<-unname(PROV[as.character(prov[["省"]])])
ct<-st_coordinates(st_centroid(prov[!is.na(prov$province),])); pxy<-data.table(province=prov$province[!is.na(prov$province)],X=ct[,1],Y=ct[,2])[province%in%unique(v3$province)]
set.seed(1); km<-kmeans(scale(pxy[,.(X,Y)]),5,nstart=25); pxy[,sfold:=km$cluster]
v3<-merge(v3,pxy[,.(province,sfold)],by="province",all.x=TRUE); v3[,rfold:=sample(rep(1:5,length.out=.N))]

# =========== (1) endogeneity: province-year effort lag ===========
log("(1) v3 endogeneity — effort lag (t-1)")
ep<-unique(v3[,.(province,year,eff=log_effort_visits_z)]); ep[,year:=year+1L]; setnames(ep,"eff","effort_lag")
dl<-merge(v3,ep,by=c("province","year"),all.x=TRUE)[is.finite(effort_lag)&is.finite(temp_grad_z)]
m4l<-glmmTMB(event~temp_grad_z*effort_lag+(1|species)+(1|province),dl,family=binomial("cloglog"))
m3l<-glmmTMB(event~temp_grad_z+effort_lag+(1|species)+(1|province),dl,family=binomial("cloglog"))
cf<-fixef(m4l)$cond;se<-sqrt(diag(vcov(m4l)$cond));i<-grep(":",names(cf))
lag_tab<-data.table(run="v3",term="temp_grad_z:effort_lag(t-1)",n=nrow(dl),events=sum(dl$event),
  hr=exp(cf[i]),hr.low=exp(cf[i]-1.96*se[i]),hr.high=exp(cf[i]+1.96*se[i]),p.value=2*pnorm(-abs(cf[i]/se[i])),dAIC_M3_M4=AIC(m3l)-AIC(m4l))
fwrite(lag_tab,file.path(V2,"results/tables/table_v3_endogeneity_lag.csv"))
log(sprintf("  v3 lag interaction HR=%.3f (%.3f-%.3f) p=%.2e dAIC=%.1f",lag_tab$hr,lag_tab$hr.low,lag_tab$hr.high,lag_tab$p.value,lag_tab$dAIC_M3_M4))

# =========== (2) migratory stratification ===========
log("(2) v3 migratory-stratified interaction")
fit_grp<-function(dt,lab){dt<-dt[is.finite(temp_grad_z)&is.finite(log_effort_visits_z)];if(sum(dt$event)<25)return(NULL)
  m<-tryCatch(glmmTMB(event~temp_grad_z*log_effort_visits_z+(1|species)+(1|province),dt,family=binomial("cloglog")),error=function(e)NULL);if(is.null(m))return(NULL)
  cf<-fixef(m)$cond;se<-sqrt(diag(vcov(m)$cond));i<-grep(":",names(cf))
  data.table(group=lab,n=nrow(dt),events=sum(dt$event),hr=exp(cf[i]),hr.low=exp(cf[i]-1.96*se[i]),hr.high=exp(cf[i]+1.96*se[i]),p.value=2*pnorm(-abs(cf[i]/se[i])))}
strat<-rbindlist(list(fit_grp(v3,"All (v3)"),fit_grp(v3[migrant_bin=="Resident"],"Resident"),
  fit_grp(v3[migrant_bin=="Migrant"],"Migrant (all)"),fit_grp(v3[mig_grp=="Partial"],"Partial migrant"),
  fit_grp(v3[mig_grp=="Long-distance"],"Long-distance migrant")),fill=TRUE)
fwrite(strat,file.path(V2,"results/tables/table_v3_migratory_stratified.csv"))
print(strat[,.(group,n,events,hr=round(hr,3),ci=sprintf("%.2f-%.2f",hr.low,hr.high),p=signif(p.value,2))])

# =========== (3) prediction: rich XGBoost, 3 CV axes ===========
log("(3) v3 prediction — rich XGBoost (random/temporal/spatial)")
feat<-intersect(c("temp_grad_z","temp_anom_z","prec_anom_z","prec_grad_prov_z","climate_velocity_z","precip_velocity_z",
  "climate_exposure_z","warming_rate_z","mahalanobis_dist_z","log_effort_visits_z","log_effort_observers_z",
  "log_effort_days_z","log_effort_record_z","effort_pc1_z","year_c","mig_code",
  "temp_grad_z_x_eff","temp_anom_z_x_eff","mahalanobis_dist_z_x_eff","climate_velocity_z_x_eff"),names(v3))
dd<-v3[complete.cases(v3[,c(feat,"event"),with=FALSE])]
log("  XGB rows:",nrow(dd),"events:",sum(dd$event),"features:",length(feat))
spw<-sum(dd$event==0)/sum(dd$event==1)
params<-list(objective="binary:logistic",eval_metric="auc",max_depth=4,eta=0.05,subsample=0.8,colsample_bytree=0.8,min_child_weight=5,scale_pos_weight=spw)
dall<-xgb.DMatrix(as.matrix(dd[,feat,with=FALSE]),label=dd$event)
cv<-xgb.cv(params,dall,nrounds=400,nfold=5,early_stopping_rounds=30,maximize=TRUE,verbose=0)
best<-tryCatch(cv$best_iteration,error=function(e)NULL); if(is.null(best)||length(best)==0||is.na(best)){el<-as.data.table(cv$evaluation_log);tc<-grep("test_.*_mean$",names(el),value=TRUE)[1];best<-if(!is.na(tc))which.max(el[[tc]]) else 200L}; if(length(best)==0||is.na(best)||best<1)best<-200L
log("  best nrounds:",best)
cvb<-function(dat,fc){a<-c();for(f in sort(unique(dat[[fc]]))){tr<-dat[get(fc)!=f];te<-dat[get(fc)==f];if(sum(te$event)<3||sum(tr$event)<10)next
  m<-xgb.train(params,xgb.DMatrix(as.matrix(tr[,feat,with=FALSE]),label=tr$event),nrounds=best,verbose=0)
  a<-c(a,auc_fun(predict(m,as.matrix(te[,feat,with=FALSE])),te$event))};a[!is.na(a)]}
auc_rd<-cvb(dd,"rfold"); auc_sp<-cvb(dd,"sfold")
trt<-dd[year<=2018];tet<-dd[year>=2019]
mt<-xgb.train(params,xgb.DMatrix(as.matrix(trt[,feat,with=FALSE]),label=trt$event),nrounds=best,verbose=0)
auc_tmp<-auc_fun(predict(mt,as.matrix(tet[,feat,with=FALSE])),tet$event)
cmp<-data.table(model=c("Rich XGBoost — random CV (interpolation)","Rich XGBoost — temporal (leave-future-out)","Rich XGBoost — spatial-block CV"),
  auc_mean=c(mean(auc_rd),auc_tmp,mean(auc_sp)),auc_sd=c(sd(auc_rd),NA,sd(auc_sp)),run="v3")
fwrite(cmp,file.path(V2,"results/tables/table_v3_prediction_accuracy.csv"))
print(cmp[,.(model,auc=round(auc_mean,3))])
mfull<-xgb.train(params,dall,nrounds=best,verbose=0); fwrite(as.data.table(xgb.importance(model=mfull))[order(-Gain)],file.path(V2,"results/tables/table_v3_xgb_importance.csv"))

# =========== Figure 10: v3 panorama ===========
strat[,group:=factor(group,levels=rev(c("All (v3)","Resident","Migrant (all)","Partial migrant","Long-distance migrant")))]
pA<-ggplot(strat,aes(hr,group,colour=group))+geom_vline(xintercept=1,linetype="dashed",colour="grey50")+
  geom_errorbarh(aes(xmin=hr.low,xmax=hr.high),height=0.22,linewidth=0.6)+geom_point(size=2.8)+
  geom_text(aes(label=sprintf("HR=%.2f, p=%.1g",hr,p.value)),hjust=0,nudge_x=0.015,size=2.5,colour="grey20")+
  scale_x_continuous(trans="log",breaks=c(0.9,1.0,1.1,1.2,1.3))+scale_colour_brewer(palette="Set2",guide="none")+
  labs(tag="a",title="v3 relaxed set: climate × effort interaction by migratory strategy",subtitle="188,870 rows / 463 species / 817 events; separate cloglog M4 per group (mig-covered subset).",x="Interaction hazard ratio (95% CI, log)",y=NULL)+
  theme_bw(base_size=9)+theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=7,colour="grey35"),panel.grid.minor=element_blank())
cmp2<-rbind(data.table(model="Baseline headline cloglog (v3, ref)",auc_mean=0.561,auc_sd=NA,run="v3"),cmp)
cmp2[,model:=factor(model,levels=rev(model))]
pB<-ggplot(cmp2,aes(auc_mean,model))+geom_vline(xintercept=0.5,linetype="dotted",colour="grey60")+
  geom_col(width=0.65,fill="#4575B4")+geom_errorbarh(aes(xmin=auc_mean-auc_sd,xmax=auc_mean+auc_sd),height=0.2)+
  geom_text(aes(label=sprintf("%.2f",auc_mean)),hjust=-0.15,size=2.7)+xlim(0,0.9)+
  labs(tag="b",title="v3 predictive AUC across validation regimes",subtitle="Interpolation > temporal forecast > spatial extrapolation — same scope-dependence as the conservative set.",x="AUC",y=NULL)+
  theme_bw(base_size=9)+theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=7,colour="grey35"),panel.grid.minor=element_blank())
ggsave(file.path(V2,"figures/main/Figure_10_v3_full_analysis.pdf"),pA/pB+plot_layout(heights=c(1,0.95)),width=20,height=19,units="cm",device=grDevices::cairo_pdf)
ggsave(file.path(V2,"figures/main/Figure_10_v3_full_analysis.png"),pA/pB+plot_layout(heights=c(1,0.95)),width=20,height=19,units="cm",dpi=300)
log("wrote Figure_10_v3_full_analysis + 4 tables"); log("DONE")
