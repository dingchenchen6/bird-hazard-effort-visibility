# ============================================================
# Script: 54_endogeneity_and_spatial_cv.R
# 科学问题 / Scientific question:
#   (A) 内生性 / Endogeneity: 用滞后 1 年的 effort 预测当年新记录，
#       若 climate × effort 交互仍显著为正，则"反向因果(effort 追逐记录)"
#       难以解释结果 —— 支持 effort 是 moderator。
#       Does the climate × effort interaction survive when effort is
#       LAGGED by one year (effort that predates the event)?
#   (B) 空间诚实泛化 / Spatially honest CV: blockCV 250km 省级空间块
#       5 折交叉验证 vs 随机 5 折，量化空间自相关下的真实判别力。
#
# 输入 / Input : data/raw/hazard_risk_upgraded_complete_case.csv
#                data/spatial/basemap_GS2019_1822/省（等积投影）.shp
# 输出 / Output: results/tables/table_effort_lag_refit.csv
#                results/tables/table_spatial_block_cv.csv
# 模型 / Model : event ~ climate_z * effort_z + (1|species)+(1|province)
#                cloglog, climate_z = temp_grad_z, effort_z = Spec
# 运行 / Run   : Rscript --no-init-file code/54_endogeneity_and_spatial_cv.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table); library(glmmTMB); library(sf)
})
set.seed(42)
options(warn = 1)
V2 <- normalizePath(".", mustWork = TRUE)
SHP <- file.path(V2, "data", "spatial", "basemap_GS2019_1822")
ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ens(file.path(V2, "results", "tables"))
log <- function(...) cat(sprintf("[54 %s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

risk <- fread(file.path(V2, "data/raw/hazard_risk_upgraded_complete_case.csv"))
log("risk set:", nrow(risk), "rows,", uniqueN(risk$species), "species,",
    sum(risk$event), "events")

specs <- list(
  spec_A = list(label = "Record-based (legacy)", col = "log_effort_record_z"),
  spec_B = list(label = "Observer visits (headline)", col = "log_effort_visits_z"),
  spec_C = list(label = "PCA composite", col = "effort_pc1_z"),
  spec_D = list(label = "Birding days", col = "log_effort_days_z"))

# AUC via Mann-Whitney (rank-based), robust & dependency-free
auc_fun <- function(pred, obs) {
  pos <- pred[obs == 1]; neg <- pred[obs == 0]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

# ============================================================
# PART A — Effort-lag refit (endogeneity defense)
# ============================================================
log("PART A — effort-lag (t-1) refit")
fit_lag <- function(dt, spec_id, lab, eff_col) {
  d <- dt[, .(species, province, year, event,
              climate_z = temp_grad_z, effort_now = get(eff_col))]
  setorder(d, species, province, year)
  # lag effort by 1 yr within species×province, only when year is consecutive
  d[, effort_lag := shift(effort_now, 1L), by = .(species, province)]
  d[, year_prev  := shift(year, 1L),       by = .(species, province)]
  d[, consec := !is.na(year_prev) & (year - year_prev == 1L)]
  dl <- d[consec == TRUE & !is.na(effort_lag) & !is.na(climate_z) & !is.na(event)]
  out <- list()
  for (nm in c("M3", "M4")) {
    f <- if (nm == "M3")
      event ~ climate_z + effort_lag + (1 | species) + (1 | province)
    else
      event ~ climate_z * effort_lag + (1 | species) + (1 | province)
    fit <- tryCatch(glmmTMB(f, data = dl, family = binomial("cloglog")),
                    error = function(e) {log("  FAIL", spec_id, nm, conditionMessage(e)); NULL})
    if (is.null(fit)) next
    cf <- fixef(fit)$cond; se <- sqrt(diag(vcov(fit)$cond))
    tm <- if (nm == "M4") "climate_z:effort_lag" else NA
    aic <- AIC(fit)
    if (nm == "M4") {
      i <- match(tm, names(cf)); b <- cf[i]; s <- se[i]
      out$row <- data.table(spec_id, spec_label = lab, term = "climate_z:effort_lag",
        n_rows = nrow(dl), n_events = sum(dl$event),
        beta = b, se = s, hr = exp(b),
        hr.low = exp(b - 1.96 * s), hr.high = exp(b + 1.96 * s),
        p.value = 2 * pnorm(-abs(b / s)), AIC_M4 = aic)
    } else out$aic_m3 <- aic
  }
  if (!is.null(out$row) && !is.null(out$aic_m3))
    out$row[, dAIC_M3_minus_M4 := out$aic_m3 - AIC_M4]
  out$row
}
lagA <- rbindlist(lapply(names(specs), function(s)
  fit_lag(risk, s, specs[[s]]$label, specs[[s]]$col)), fill = TRUE)
fwrite(lagA, file.path(V2, "results/tables/table_effort_lag_refit.csv"))
log("  wrote table_effort_lag_refit.csv")
print(lagA[, .(spec_id, n_rows, n_events, hr = round(hr,3),
               ci = sprintf("%.3f-%.3f", hr.low, hr.high),
               p = signif(p.value,3), dAIC = round(dAIC_M3_minus_M4,1))])

# ============================================================
# PART B — Spatial-block CV (province centroids, blockCV 250 km)
# ============================================================
log("PART B — spatial-block CV vs random CV (headline Spec B)")
PROV_CN_EN <- c("北京市"="Beijing","天津市"="Tianjin","河北省"="Hebei",
  "山西省"="Shanxi","内蒙古自治区"="Inner Mongolia","辽宁省"="Liaoning",
  "吉林省"="Jilin","黑龙江省"="Heilongjiang","上海市"="Shanghai",
  "江苏省"="Jiangsu","浙江省"="Zhejiang","安徽省"="Anhui","福建省"="Fujian",
  "江西省"="Jiangxi","山东省"="Shandong","河南省"="Henan","湖北省"="Hubei",
  "湖南省"="Hunan","广东省"="Guangdong","广西壮族自治区"="Guangxi",
  "海南省"="Hainan","重庆市"="Chongqing","四川省"="Sichuan","贵州省"="Guizhou",
  "云南省"="Yunnan","西藏自治区"="Tibet","陕西省"="Shaanxi","甘肃省"="Gansu",
  "青海省"="Qinghai","宁夏回族自治区"="Ningxia","新疆维吾尔自治区"="Xinjiang",
  "台湾省"="Taiwan","香港特别行政区"="Hong Kong","澳门特别行政区"="Macau")

prov_sf <- st_read(file.path(SHP, "省（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4524) |> st_make_valid()
prov_sf$province <- unname(PROV_CN_EN[as.character(prov_sf[["省"]])])
cent <- st_centroid(prov_sf[!is.na(prov_sf$province), ])
cc <- st_coordinates(cent)
prov_xy <- data.table(province = cent$province, X = cc[,1], Y = cc[,2])
prov_xy <- prov_xy[province %in% unique(risk$province)]

# spatial folds: blockCV 250 km if available, else k-means(5) on centroids
make_spatial_folds <- function(xy, k = 5) {
  ok <- requireNamespace("blockCV", quietly = TRUE)
  if (ok) {
    sf_pts <- st_as_sf(xy, coords = c("X","Y"), crs = 4524)
    res <- tryCatch(blockCV::cv_spatial(x = sf_pts, size = 250000, k = k,
             selection = "random", iteration = 50, progress = FALSE,
             plot = FALSE, report = FALSE), error = function(e) NULL)
    if (!is.null(res)) { xy$fold <- res$folds_ids; return(xy) }
  }
  km <- kmeans(scale(xy[, .(X, Y)]), centers = k, nstart = 20)
  xy$fold <- km$cluster; xy
}
prov_xy <- make_spatial_folds(prov_xy, k = 5)
risk2 <- merge(risk, prov_xy[, .(province, fold)], by = "province", all.x = TRUE)

cv_run <- function(dt, eff_col, fold_col) {
  d <- dt[, .(species, province, event, climate_z = temp_grad_z,
              effort_z = get(eff_col), fold = get(fold_col))]
  d <- d[!is.na(fold)]
  aucs <- c()
  for (f in sort(unique(d$fold))) {
    tr <- d[fold != f]; te <- d[fold == f]
    if (sum(te$event) < 1 || sum(tr$event) < 5) next
    fit <- tryCatch(glmmTMB(event ~ climate_z * effort_z +
              (1|species) + (1|province), data = tr,
              family = binomial("cloglog")), error = function(e) NULL)
    if (is.null(fit)) next
    pr <- tryCatch(predict(fit, newdata = te, type = "response",
              allow.new.levels = TRUE), error = function(e) NULL)
    if (is.null(pr)) next
    aucs <- c(aucs, auc_fun(pr, te$event))
  }
  aucs[!is.na(aucs)]
}

# random 5-fold for comparison
risk2[, rfold := sample(rep(1:5, length.out = .N))]
auc_spatial <- cv_run(risk2, "log_effort_visits_z", "fold")
auc_random  <- cv_run(risk2, "log_effort_visits_z", "rfold")

cvtab <- data.table(
  cv_type = c("spatial_block_250km", "random_5fold"),
  n_folds = c(length(auc_spatial), length(auc_random)),
  auc_mean = c(mean(auc_spatial), mean(auc_random)),
  auc_sd   = c(sd(auc_spatial),   sd(auc_random)),
  auc_min  = c(min(auc_spatial),  min(auc_random)),
  auc_max  = c(max(auc_spatial),  max(auc_random)),
  spec = "B_observer_visits", model = "M4_cloglog")
fwrite(cvtab, file.path(V2, "results/tables/table_spatial_block_cv.csv"))
log("  wrote table_spatial_block_cv.csv")
print(cvtab)
log("DONE")
