# ============================================================
# Script: 58_grid100_native_plugin.R  (v2, corrected join)
# 设计/Design: 纯 plug-in 多尺度预测(用户确认)。仅省级拟合头条 M4
#   (temp_grad_z × log_effort_visits_z); 把【努力】真正降尺度到 100km 社区
#   网格(941 格系统, grid_100km_base 提供质心坐标), 按【省级训练口径】重标准化;
#   气候用省级 temp_grad_z 广播(复用已验证的省级未来轨迹 SSP×年)。
#   Province-fit headline M4; grid-native EFFORT + province-broadcast climate;
#   plug-in prediction over SSP245/585 × 2030/2050/2080.
# 诚实说明/Honesty: 真正的"网格原生时变气候"需原始 CHELSA/WorldClim 栅格,
#   本机不可得; 故气候停留在省级分辨率(广播), 网格尺度的真实变化来自努力。
#   100km 社区努力网格稀疏(2024 仅 64/941 格有访问)。50km 缺坐标与匹配指标,
#   本机无法等严谨完成(另行说明)。
# 输入: data/raw/hazard_risk_upgraded_complete_case.csv
#       data/raw/grid_100km_base.csv  (grid_id 1-941, province, centroid_lon/lat)
#       data/derived/grid_100km_effort_native.parquet
#       results/forecasts/table_province_future_glmmTMB.csv  (province temp_grad_z SSP×yr)
# 输出: results/forecasts/table_grid_100km_effort_plugin.csv
#       figures/main/Figure_8_grid100_effort_plugin_hazard.{pdf,png}
# 运行: Rscript --no-init-file code/58_grid100_native_plugin.R
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(glmmTMB); library(arrow)
  library(sf); library(ggplot2); library(viridisLite)
})
sf::sf_use_s2(FALSE); options(warn = 1)
V2 <- normalizePath(".", mustWork = TRUE)
SHP <- file.path(V2, "data", "spatial", "basemap_GS2019_1822")
log <- function(...) cat(sprintf("[58 %s] ", format(Sys.time(),"%H:%M:%S")),...,"\n")
inv_cloglog <- function(eta) 1 - exp(-exp(eta))

# ---- 1. province headline M4 fit ----
risk <- fread(file.path(V2,"data/raw/hazard_risk_upgraded_complete_case.csv"))
m4 <- glmmTMB(event ~ temp_grad_z * log_effort_visits_z + (1|species) + (1|province),
              data = risk, family = binomial("cloglog"))
cf <- fixef(m4)$cond
b0 <- cf["(Intercept)"]; bC <- cf["temp_grad_z"]
bE <- cf["log_effort_visits_z"]; bI <- cf["temp_grad_z:log_effort_visits_z"]
log(sprintf("province headline M4: b0=%.3f bC=%.3f bE=%.3f bI=%.3f (int HR=%.3f)",
            b0,bC,bE,bI,exp(bI)))

# ---- 2. recover effort transform (log1p(n_visits) -> province z) ----
ev <- unique(risk[is.finite(log_effort_visits_z), .(n_visits, log_effort_visits_z)])
ev[, lv := log1p(n_visits)]
fitE <- lm(lv ~ log_effort_visits_z, ev); m_eff <- coef(fitE)[1]; s_eff <- coef(fitE)[2]
log(sprintf("effort transform: log1p mean=%.3f sd=%.3f", m_eff, s_eff))

# ---- 3. grid-native effort (941 system) + coordinates ----
base <- fread(file.path(V2,"data/raw/grid_100km_base.csv"))   # grid_id, province, centroid_lon/lat
ge <- as.data.table(read_parquet(file.path(V2,"data/derived/grid_100km_effort_native.parquet")))
eff24 <- ge[year == 2024, .(grid_id, n_visits_grid)]
g <- merge(base, eff24, by = "grid_id", all.x = TRUE)
g[is.na(n_visits_grid), n_visits_grid := 0]
g[, eff_z := (log1p(n_visits_grid) - m_eff)/s_eff]   # province-scale standardisation
log(sprintf("grid cells=%d | nonzero-effort=%d | eff_z range %.2f..%.2f",
            nrow(g), sum(g$n_visits_grid>0), min(g$eff_z), max(g$eff_z)))

# ---- 4. province temp_grad_z trajectory (reuse validated province future) ----
pf <- fread(file.path(V2,"results/forecasts/table_province_future_glmmTMB.csv"))
prov_clim <- pf[, .(province, ssp, year, temp_grad_z)]
# add a "Current (2024)" baseline = de-perturbed (SSP245 2030 minus 0.3*0.6 SD)
cur <- pf[ssp=="SSP245" & year==2030, .(province, temp_grad_z = temp_grad_z - 0.3*0.6)]
cur[, `:=`(ssp="Current", year=2024L)]
prov_clim <- rbind(prov_clim, cur[, .(province, ssp, year, temp_grad_z)])

# ---- 5. plug-in prediction: grid effort × province climate ----
pred <- merge(g[, .(grid_id, province, centroid_lon, centroid_lat, eff_z)],
              prov_clim, by = "province", allow.cartesian = TRUE)
pred[, hazard := inv_cloglog(b0 + bC*temp_grad_z + bE*eff_z + bI*temp_grad_z*eff_z)]
pred[, `:=`(lon = centroid_lon, lat = centroid_lat)]
fwrite(pred[, .(grid_id, province, lon=centroid_lon, lat=centroid_lat,
                ssp, year, temp_grad_z, eff_z, hazard)],
       file.path(V2,"results/forecasts/table_grid_100km_effort_plugin.csv"))
log("wrote table_grid_100km_effort_plugin.csv (", nrow(pred)," rows )")
log(sprintf("Current mean=%.4f | SSP585-2050 mean=%.4f max=%.4f",
    mean(pred[ssp=="Current"]$hazard),
    mean(pred[ssp=="SSP585" & year==2050]$hazard),
    max(pred[ssp=="SSP585" & year==2050]$hazard)))

# ---- 6. map ----
prov <- tryCatch(st_transform(st_make_valid(st_read(file.path(SHP,"省（等积投影）.shp"),
          quiet=TRUE)), 4326), error=function(e) NULL)
yrs <- c(2030,2050,2080)
pred[, panel := factor(ifelse(ssp=="Current","Current (2024)", paste0(ssp," — ",year)),
     levels=c("Current (2024)", paste0("SSP245 — ",yrs), paste0("SSP585 — ",yrs)))]
qcap <- quantile(pred$hazard, 0.99, na.rm=TRUE); pred[, hz := pmin(hazard, qcap)]
p <- ggplot()
if (!is.null(prov)) p <- p + geom_sf(data=prov, fill="grey96", colour="grey80", linewidth=0.12)
p <- p +
  geom_point(data=pred, aes(lon, lat, colour=hz), size=0.9, shape=15) +
  scale_colour_viridis_c(option="mako", direction=-1, name="Hazard (prob.)") +
  facet_wrap(~panel, ncol=4) +
  coord_sf(xlim=c(73,135), ylim=c(17,54), expand=FALSE, datum=NA) +
  labs(title="100 km grid plug-in hazard — effort downscaled to grid, climate at province resolution",
       subtitle=paste0("Headline province M4 (temp_grad_z × visits) applied to GRID-NATIVE effort ",
       "(941-cell community grid, re-standardised to province scale) × province temp_grad_z ",
       "(SSP +0.3/+0.8 SD/decade). Grid-native CLIMATE needs rasters (unavailable). Fill clipped 99th pct."),
       x=NULL, y=NULL) +
  theme_bw(base_size=8) +
  theme(legend.position="bottom", panel.grid=element_blank(),
        plot.title=element_text(face="bold", size=9.5),
        plot.subtitle=element_text(size=6.8, colour="grey35"),
        strip.text=element_text(face="bold", size=7))
ggsave(file.path(V2,"figures/main/Figure_8_grid100_effort_plugin_hazard.pdf"), p,
       width=24, height=14, units="cm", device=grDevices::cairo_pdf)
ggsave(file.path(V2,"figures/main/Figure_8_grid100_effort_plugin_hazard.png"), p,
       width=24, height=14, units="cm", dpi=300)
log("wrote Figure_8_grid100_effort_plugin_hazard.{pdf,png}")
log("DONE")
