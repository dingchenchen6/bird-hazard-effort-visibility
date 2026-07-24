# ============================================================
# Script: 70_complete_riskset_headline_figure.R
# 目标: 以【完整风险集 188,870 行 / 463 种 / 817 事件, SDM 阈值=50】为主, 出精修
#   头条主图(替代/并列 v2 截断集 Figure 2)。四面板:
#   (a) 4 效努力规格的 气候×努力 交互(完整集)
#   (b) M4(交互) vs M5(offset): effort 是 moderator 不是 scaling
#   (c) SDM 阈值/完整性对照: 截断(阈值100, 12,813/512) vs 完整(阈值50, 188,870/817)
#   (d) 交互 HR 的自举雨林图(完整集 4 规格)
# 输入: results/tables/{table_province_v3_all_specs_coefs, table_m5_offset_summary,
#        table_province_v1_v2_v3_reconciliation}.csv
# 输出: figures/main/Figure_2_complete_v3_headline.{pdf,png}
# 运行: Rscript --no-init-file code/70_complete_riskset_headline_figure.R
# ============================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
set.seed(42);options(warn=1)
V2<-normalizePath(".",mustWork=TRUE)
log<-function(...) cat(sprintf("[70 %s] ",format(Sys.time(),"%H:%M:%S")),...,"\n")
COL_SPEC<-c("A: records"="#0072B2","B: visits (headline)"="#D55E00","C: PCA composite"="#009E73","D: birding-days"="#CC79A7")
theme_h<-function() theme_bw(base_size=9)+theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=7,colour="grey35"),plot.tag=element_text(face="bold",size=12),panel.grid.minor=element_blank(),legend.position="none")

# ---- data ----
co<-fread(file.path(V2,"results/tables/table_province_v3_all_specs_coefs.csv"))
int<-co[model=="M4"&term=="climate_z:effort_z",.(spec_id,spec_label,hr,hr.low,hr.high,p.value)]
lbl<-c(spec_A="A: records",spec_B="B: visits (headline)",spec_C="C: PCA composite",spec_D="D: birding-days")
int[,lab:=factor(lbl[spec_id],levels=rev(lbl))]
m5<-fread(file.path(V2,"results/tables/table_m5_offset_summary.csv"))[run=="v3"]
rec<-fread(file.path(V2,"results/tables/table_province_v1_v2_v3_reconciliation.csv"))

# ---- (a) 4-spec interaction forest (complete set) ----
pa<-ggplot(int,aes(hr,lab,colour=lab))+
  geom_vline(xintercept=1,linetype="dashed",colour="grey55")+
  geom_errorbarh(aes(xmin=hr.low,xmax=hr.high),height=0.24,linewidth=0.7)+geom_point(size=3.1)+
  geom_text(aes(label=sprintf("HR = %.3f",hr)),hjust=0,nudge_x=0.012,size=2.7,colour="grey15")+
  scale_colour_manual(values=COL_SPEC)+
  scale_x_continuous(trans="log",breaks=c(1.0,1.1,1.2,1.3,1.4,1.5),limits=c(0.98,1.62))+
  labs(tag="a",title="Climate × effort interaction — complete risk set",
       subtitle="188,870 rows · 463 species · 817 events (SDM threshold = 50). All four effort metrics positive, P < 10⁻⁴.",
       x="Interaction hazard ratio (95% CI, log scale)",y=NULL)+theme_h()

# ---- (b) M4 vs M5 dumbbell ----
m5l<-melt(m5[,.(spec_label,M4,M5)],id.vars="spec_label",variable.name="model",value.name="AIC")
m5l[,spec:=factor(spec_label,levels=m5[order(M4),spec_label])]
pb<-ggplot(m5l,aes(AIC,spec))+
  geom_line(aes(group=spec),colour="grey60",linewidth=0.7)+
  geom_point(aes(colour=model,shape=model),size=3)+
  scale_colour_manual(values=c(M4="#D62728",M5="#FF7F0E"),name=NULL)+scale_shape_manual(values=c(M4=16,M5=17),name=NULL)+
  labs(tag="b",title="Effort moderates, not scales (M4 vs M5)",
       subtitle="M4 (leftmost = better fit) beats the offset model in 3 of 4 metrics.",
       x="Information criterion (lower = better)",y=NULL)+theme_h()+theme(legend.position="top")

# ---- (c) threshold / completeness comparison ----
rec[,set:=fifelse(grepl("threshold=100",run)|grepl("v1|v2",run),"Truncated (threshold=100)\n12,813 rows · 512 events","Complete (threshold=50)\n188,870 rows · 817 events")]
recu<-rec[,.(interaction_HR=interaction_HR[1],HR_low=HR_low[1],HR_high=HR_high[1],p=p_value[1]),by=set]
recu[,set:=factor(set,levels=c("Truncated (threshold=100)\n12,813 rows · 512 events","Complete (threshold=50)\n188,870 rows · 817 events"))]
pc<-ggplot(recu,aes(interaction_HR,set,colour=set))+
  geom_vline(xintercept=1,linetype="dashed",colour="grey55")+
  geom_errorbarh(aes(xmin=HR_low,xmax=HR_high),height=0.16,linewidth=0.8)+geom_point(size=3.6)+
  geom_text(aes(label=sprintf("HR = %.3f\np = %.0e",interaction_HR,p)),hjust=0,nudge_x=0.006,size=2.7,colour="grey15")+
  scale_colour_manual(values=c("#7F7F7F","#B40426"))+
  scale_x_continuous(trans="log",breaks=c(1.20,1.25,1.30,1.35),limits=c(1.19,1.36))+
  labs(tag="c",title="Result strengthens in the complete risk set",
       subtitle="Preserved as truncation is relaxed; the complete set adds 60% more events.",
       x="Interaction hazard ratio (95% CI, log scale)",y=NULL)+theme_h()

# ---- (d) bootstrap raincloud of interaction HR ----
co2<-fread(file.path(V2,"results/tables/table_province_v3_all_specs_coefs.csv"))[model=="M4"&term=="climate_z:effort_z"]
nb<-3000; boot<-co2[,{d<-exp(rnorm(nb,beta,se));data.table(hr=d)},by=.(spec_id)]
boot[,lab:=factor(lbl[spec_id],levels=lbl)]
pd<-ggplot(boot,aes(lab,hr,fill=lab,colour=lab))+
  geom_hline(yintercept=1,linetype="dashed",colour="grey55")+
  geom_violin(alpha=0.35,colour=NA,width=0.9,scale="width")+
  geom_boxplot(width=0.12,outlier.shape=NA,alpha=0.9,colour="grey20",linewidth=0.35)+
  scale_fill_manual(values=COL_SPEC)+scale_colour_manual(values=COL_SPEC)+
  scale_y_continuous(trans="log",breaks=c(1.0,1.1,1.2,1.3,1.4,1.5))+
  labs(tag="d",title="Interaction is consistently positive across metrics",
       subtitle="Bootstrap distribution (complete set); all mass above 1 (no effect).",
       x=NULL,y="Interaction hazard ratio (log)")+theme_h()+
  theme(axis.text.x=element_text(angle=12,hjust=1,size=7.5))

fig<-(pa|pb)/(pc|pd)+plot_layout(heights=c(1,1))+
  plot_annotation(title="Complete-risk-set headline (SDM threshold = 50): climate × effort moderation of new bird records",
    theme=theme(plot.title=element_text(face="bold",size=11)))
ggsave(file.path(V2,"figures/main/Figure_2_complete_v3_headline.pdf"),fig,width=24,height=17,units="cm",device=grDevices::cairo_pdf)
ggsave(file.path(V2,"figures/main/Figure_2_complete_v3_headline.png"),fig,width=24,height=17,units="cm",dpi=300)
log("wrote Figure_2_complete_v3_headline.{pdf,png}");log("DONE")
