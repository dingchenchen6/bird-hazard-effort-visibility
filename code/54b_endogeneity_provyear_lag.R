# ============================================================
# Script: 54b_endogeneity_provyear_lag.R
# 修正 / Fix: effort 是"省×年"层级(同省同年所有物种共享)，故在
#   省×年层级把 effort 滞后 1 年(t-1)，再回连到每条风险集行，
#   几乎保留全部行(仅丢最早年)。气候(temp_grad_z, 物种×省×年)保持当年。
#   Lag effort at the province-year level (it is constant across species
#   within a province-year), preserving nearly all risk-set rows.
# 检验 / Test: M4 climate_z * effort_lag —— 若交互仍显著为正，则
#   "effort 追逐已发生记录"的反向因果难以解释结果。
# 输出 / Output: results/tables/table_effort_lag_refit.csv (overwrite, correct)
# ============================================================
suppressPackageStartupMessages({ library(data.table); library(glmmTMB) })
set.seed(42); options(warn = 1)
V2 <- normalizePath(".", mustWork = TRUE)
log <- function(...) cat(sprintf("[54b %s] ", format(Sys.time(),"%H:%M:%S")),...,"\n")
risk <- fread(file.path(V2,"data/raw/hazard_risk_upgraded_complete_case.csv"))
log("risk:", nrow(risk),"rows,",sum(risk$event),"events")

specs <- list(
  spec_A = list(label="Record-based (legacy)",      col="log_effort_record_z"),
  spec_B = list(label="Observer visits (headline)", col="log_effort_visits_z"),
  spec_C = list(label="PCA composite",              col="effort_pc1_z"),
  spec_D = list(label="Birding days",               col="log_effort_days_z"))

fit_lag <- function(dt, spec_id, lab, eff_col) {
  d <- dt[, .(species, province, year, event,
              climate_z = temp_grad_z, effort_now = get(eff_col))]
  # province-year effort panel, shifted +1 yr so it joins as PREVIOUS year
  ep <- unique(d[, .(province, year, effort_lag = effort_now)])
  ep[, year := year + 1L]
  d <- merge(d, ep, by = c("province","year"), all.x = TRUE)
  dl <- d[!is.na(effort_lag) & !is.na(climate_z) & !is.na(event)]
  rows <- list()
  m3 <- tryCatch(glmmTMB(event ~ climate_z + effort_lag + (1|species)+(1|province),
            data=dl, family=binomial("cloglog")), error=function(e) NULL)
  m4 <- tryCatch(glmmTMB(event ~ climate_z * effort_lag + (1|species)+(1|province),
            data=dl, family=binomial("cloglog")), error=function(e) NULL)
  if (is.null(m4)) { log("  M4 FAIL",spec_id); return(NULL) }
  cf <- fixef(m4)$cond; se <- sqrt(diag(vcov(m4)$cond))
  i <- match("climate_z:effort_lag", names(cf)); b <- cf[i]; s <- se[i]
  data.table(spec_id, spec_label=lab, term="climate_z:effort_lag(t-1)",
    n_rows=nrow(dl), n_events=sum(dl$event),
    beta=b, se=s, hr=exp(b), hr.low=exp(b-1.96*s), hr.high=exp(b+1.96*s),
    p.value=2*pnorm(-abs(b/s)),
    dAIC_M3_minus_M4 = if(!is.null(m3)) AIC(m3)-AIC(m4) else NA_real_)
}
lag_tab <- rbindlist(lapply(names(specs), function(s)
  fit_lag(risk, s, specs[[s]]$label, specs[[s]]$col)), fill=TRUE)
fwrite(lag_tab, file.path(V2,"results/tables/table_effort_lag_refit.csv"))
log("wrote table_effort_lag_refit.csv (province-year lag, corrected)")
print(lag_tab[, .(spec_id, n_rows, n_events, hr=round(hr,3),
  ci=sprintf("%.3f-%.3f",hr.low,hr.high), p=signif(p.value,3),
  dAIC=round(dAIC_M3_minus_M4,1))])
log("DONE")
