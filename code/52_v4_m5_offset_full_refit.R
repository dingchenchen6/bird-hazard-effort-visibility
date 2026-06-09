# ============================================================
# Script: 52_v4_m5_offset_full_refit.R
# Family: v4 — add M5 (effort-as-offset) to the M0-M4 ladder
# Author: Chen-Chen Ding + Claude Opus 4.7
# Date  : 2026-05-23
#
# ------------------------------------------------------------
# Scientific question / 科学问题:
#   The previous M0-M4 ladder lacks an explicit "effort-as-offset"
#   model. M5 treats log(effort) as an offset on the cloglog hazard
#   linear predictor — i.e., the hazard is hypothesised to scale
#   one-for-one with log(effort). If M5 fits better than M4, effort
#   is a SCALING FACTOR; if M4 wins, effort is a MODERATOR (the
#   manuscript's preferred interpretation).
#   把 effort 作为 offset 的 M5 加入 ladder，与 M4 交互比较。
#
#   Formula:
#     M5 :  event ~ climate_z
#                   + offset(log(effort_raw + 1))
#                   + (1|species) + (1|unit_id)
#
# ------------------------------------------------------------
# Files written:
#   results/tables/table_province_v2_with_m5_coefs.csv
#   results/tables/table_province_v2_with_m5_aic.csv
#   results/tables/table_province_v3_with_m5_coefs.csv
#   results/tables/table_province_v3_with_m5_aic.csv
#   results/tables/table_m5_offset_summary.csv  (M4 vs M5 ΔAIC across runs/specs)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
})
options(warn = 1)
set.seed(42)

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- normalizePath(file.path(V2, "..", "bird_hazard_model_effort_upgrade"),
                     mustWork = FALSE)
ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE,
                                                    showWarnings = FALSE)
ens(file.path(V2, "results", "tables"))
ens(file.path(V2, "logs"))

LOG_PATH <- file.path(V2, "logs", "52_m5_offset_audit.log")
LOG_CON <- file(LOG_PATH, open = "wt", encoding = "UTF-8")
on.exit({ if (isOpen(LOG_CON)) close(LOG_CON) }, add = TRUE)
.log <- function(...) {
  msg <- paste0(sprintf("[52 %s] ", format(Sys.time(), "%H:%M:%S")),
                 paste0(...))   # use paste0 to concatenate args correctly
  cat(msg, "\n", sep = "")
  writeLines(msg, LOG_CON)
}

# ============================================================
# STEP A — Build risk-set inputs (v2 + v3 in parallel)
# ============================================================
.log("loading inputs")
risk_v2 <- fread(file.path(V1, "data",
                            "hazard_risk_upgraded_complete_case.csv"),
                  encoding = "UTF-8")
risk_v3 <- fread(file.path(V2, "data", "derived",
                            "risk_set_province_v3.csv"),
                  encoding = "UTF-8")
# Need the raw effort_raw counts (n_visits, n_records, n_observers, n_days)
# for offset construction; the v1 risk set carries the *_z versions, so we
# rejoin from effort_panel_upgraded.csv.
prov_eff <- fread(file.path(V1, "data", "effort_panel_upgraded.csv"),
                   encoding = "UTF-8")
eff_raw <- prov_eff[, .(province, year,
                          n_records = effort_record,
                          n_visits, n_birding_days, effort_pc1)]
# Drop pre-existing effort raw columns to avoid .x/.y duplication on merge.
# v1 risk_v2 carries n_visits/n_birding_days/effort_pc1 already; remove them.
drop_existing <- intersect(c("n_visits","n_birding_days","effort_pc1",
                              "n_records","effort_record"),
                            names(risk_v2))
if (length(drop_existing) > 0) {
  cat("[52] dropping pre-existing effort-raw columns from v2: ",
      paste(drop_existing, collapse = ", "), "\n", sep = "")
  risk_v2[, (drop_existing) := NULL]
}
drop_existing3 <- intersect(c("n_visits","n_birding_days","effort_pc1",
                                "n_records","effort_record"),
                              names(risk_v3))
if (length(drop_existing3) > 0) risk_v3[, (drop_existing3) := NULL]
risk_v2 <- merge(risk_v2, eff_raw, by = c("province", "year"), all.x = TRUE)
risk_v3 <- merge(risk_v3, eff_raw, by = c("province", "year"), all.x = TRUE)
.log("risk_v2 rows = ", nrow(risk_v2), " | events = ", sum(risk_v2$event))
.log("risk_v3 rows = ", nrow(risk_v3), " | events = ", sum(risk_v3$event))

specs <- list(
  spec_A = list(label = "Record-based (legacy)",
                effort_z = "log_effort_record_z",
                effort_raw = "n_records"),
  spec_B = list(label = "Observer visits",
                effort_z = "log_effort_visits_z",
                effort_raw = "n_visits"),
  spec_C = list(label = "PCA composite",
                effort_z = "effort_pc1_z",
                effort_raw = "effort_pc1"),
  spec_D = list(label = "Birding days",
                effort_z = "log_effort_days_z",
                effort_raw = "n_birding_days"))

# ============================================================
# STEP B — Generic fit-extract for one (risk-set, spec)
# ============================================================
fit_one <- function(d, spec_id, spec_label, eff_z_col, eff_raw_col,
                     run_label) {
  d <- copy(d)
  d[, climate_z := temp_grad_z]
  if (!eff_z_col %in% names(d)) { .log("  ", spec_id, " missing effort_z col: ", eff_z_col); return(NULL) }
  if (!eff_raw_col %in% names(d)) { .log("  ", spec_id, " missing effort_raw col: ", eff_raw_col); return(NULL) }
  d[, effort_z := get(eff_z_col)]
  d[, effort_raw_lp := log(pmax(get(eff_raw_col), 0) + 1)]  # offset value
  d <- d[!is.na(climate_z) & !is.na(effort_z) & !is.na(effort_raw_lp) &
          !is.na(event)]
  .log("  ", run_label, " ", spec_id, " (", spec_label, ") n = ", nrow(d),
      " events = ", sum(d$event))
  forms <- list(
    M0 = "event ~ 1                          + (1|species) + (1|province)",
    M1 = "event ~ effort_z                    + (1|species) + (1|province)",
    M2 = "event ~ climate_z                   + (1|species) + (1|province)",
    M3 = "event ~ climate_z + effort_z        + (1|species) + (1|province)",
    M4 = "event ~ climate_z * effort_z        + (1|species) + (1|province)",
    M5 = "event ~ climate_z + offset(effort_raw_lp) + (1|species) + (1|province)")
  rows_c <- list(); rows_a <- list()
  for (nm in names(forms)) {
    t0 <- Sys.time()
    fit <- tryCatch(
      glmmTMB(as.formula(forms[[nm]]), data = d,
              family = binomial(link = "cloglog")),
      error = function(e) { .log("    ", nm, " FAILED: ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    .log(sprintf("    %s OK (%.1fs, AIC = %.2f, n = %d)", nm,
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  AIC(fit), nobs(fit)))
    cf <- fixef(fit)$cond
    se <- sqrt(diag(stats::vcov(fit)$cond))
    for (tm in names(cf)) {
      i <- match(tm, names(cf))
      rows_c[[length(rows_c) + 1L]] <- data.table(
        run = run_label, spec_id = spec_id, spec_label = spec_label,
        model = nm, term = tm,
        beta = cf[i], se = se[i],
        hr = exp(cf[i]),
        hr.low  = exp(cf[i] - 1.96 * se[i]),
        hr.high = exp(cf[i] + 1.96 * se[i]),
        p.value = 2 * pnorm(-abs(cf[i] / se[i])))
    }
    rows_a[[length(rows_a) + 1L]] <- data.table(
      run = run_label, spec_id = spec_id, spec_label = spec_label,
      model = nm,
      AIC = AIC(fit), BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)),
      df = attr(logLik(fit), "df"),
      nobs = nobs(fit))
    invisible(gc(verbose = FALSE))
  }
  list(coefs = rbindlist(rows_c, fill = TRUE),
       aic   = rbindlist(rows_a, fill = TRUE))
}

# ============================================================
# STEP C — Refit v2 (4 spec × 6 models = 24 fits)
# ============================================================
.log("")
.log("====== v2 (SDM-tight, 12,813 rows) — 4 spec × M0-M5 ======")
v2_all <- lapply(names(specs), function(sid) {
  fit_one(risk_v2, sid, specs[[sid]]$label,
          specs[[sid]]$effort_z, specs[[sid]]$effort_raw, "v2")
})
v2_coefs <- rbindlist(lapply(v2_all, function(x) x$coefs), fill = TRUE)
v2_aic   <- rbindlist(lapply(v2_all, function(x) x$aic),   fill = TRUE)
v2_aic[, dAIC := AIC - min(AIC, na.rm = TRUE), by = spec_id]
v2_aic[, akaike_weight := exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC)), by = spec_id]
fwrite(v2_coefs, file.path(V2, "results", "tables",
                             "table_province_v2_with_m5_coefs.csv"))
fwrite(v2_aic,   file.path(V2, "results", "tables",
                             "table_province_v2_with_m5_aic.csv"))

# ============================================================
# STEP D — Refit v3 (4 spec × 6 models = 24 fits)
# ============================================================
.log("")
.log("====== v3 (relaxed, 188,870 rows) — 4 spec × M0-M5 ======")
# v3 risk set was built with temp_grad_z + log_effort_visits_z only — we
# need to enrich it with the other effort z columns too, taken from the
# v1 effort panel. 把另外 3 个 effort z 列也加进 v3。
extras <- prov_eff[, .(province, year, log_effort_record_z,
                         log_effort_observers_z, log_effort_days_z, effort_pc1_z)]
risk_v3 <- merge(risk_v3, extras, by = c("province", "year"), all.x = TRUE)
v3_all <- lapply(names(specs), function(sid) {
  fit_one(risk_v3, sid, specs[[sid]]$label,
          specs[[sid]]$effort_z, specs[[sid]]$effort_raw, "v3")
})
v3_coefs <- rbindlist(lapply(v3_all, function(x) x$coefs), fill = TRUE)
v3_aic   <- rbindlist(lapply(v3_all, function(x) x$aic),   fill = TRUE)
v3_aic[, dAIC := AIC - min(AIC, na.rm = TRUE), by = spec_id]
v3_aic[, akaike_weight := exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC)), by = spec_id]
fwrite(v3_coefs, file.path(V2, "results", "tables",
                             "table_province_v3_with_m5_coefs.csv"))
fwrite(v3_aic,   file.path(V2, "results", "tables",
                             "table_province_v3_with_m5_aic.csv"))

# ============================================================
# STEP E — M4 vs M5 summary
# ============================================================
.log("")
.log("====== M4 vs M5 ΔAIC summary ======")
summ_v2 <- dcast(v2_aic[model %in% c("M4","M5")],
                  spec_id + spec_label ~ model, value.var = "AIC")
summ_v2[, dAIC_M4_minus_M5 := M4 - M5]
summ_v2[, run := "v2"]
summ_v3 <- dcast(v3_aic[model %in% c("M4","M5")],
                  spec_id + spec_label ~ model, value.var = "AIC")
summ_v3[, dAIC_M4_minus_M5 := M4 - M5]
summ_v3[, run := "v3"]
m5_summary <- rbind(summ_v2, summ_v3)
m5_summary[, interpretation := ifelse(dAIC_M4_minus_M5 < 0,
                                         "M4 wins (interaction; effort is moderator)",
                                         "M5 wins (offset; effort is scaling factor)")]
fwrite(m5_summary, file.path(V2, "results", "tables",
                                "table_m5_offset_summary.csv"))
print(m5_summary)

.log("")
.log("DONE — v2 + v3 M0-M5 fits persisted.")
