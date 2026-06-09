# ============================================================
# Scientific question / 科学问题:
#   The province-level hazard model is the manuscript's headline.
#   This script refits M0-M4 across all four effort specifications
#   (record, visits, observers/PCA composite, birding-days) at
#   province scale, persists every coefficient + AIC, and emits an
#   automatic reconciliation table against v1's published numbers
#   (so manuscript_v2 can quote v2 figures with full traceability).
#   省级是论文主结论；本脚本只在省级跑 M0-M4 × 4 个 effort 指标，
#   写表保存，并自动跟 v1 已发布表对照。
#
# Input data / 输入数据:
#   data/raw/hazard_risk_upgraded_complete_case.csv
#   v1 reference tables in V1_ROOT/results/
#
# Outputs / 输出:
#   results/tables/table_province_v2_coefs.csv
#   results/tables/table_province_v2_aic.csv
#   results/tables/table_province_v1_v2_reconciliation.csv
#   figures/main/fig_province_interaction_forest_v1_vs_v2.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
  library(ggplot2)
})
set.seed(42)
options(warn = 1)

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- Sys.getenv("V1_ROOT",
                  normalizePath(file.path(V2, "..",
                                          "bird_hazard_model_effort_upgrade"),
                                 mustWork = FALSE))

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ens(file.path(V2, "results", "tables"))
ens(file.path(V2, "figures", "main"))
ens(file.path(V2, "logs"))

log <- function(...) cat(sprintf("[40c %s] ", format(Sys.time(), "%H:%M:%S")),
                          ..., "\n", sep = "")

# ---- 1. Load province risk set --------------------------------------------
risk <- fread(file.path(V1, "data",
                         "hazard_risk_upgraded_complete_case.csv"),
              encoding = "UTF-8")
log("province risk set: rows=", nrow(risk),
    " species=", uniqueN(risk$species),
    " provinces=", uniqueN(risk$province),
    " events=", sum(risk$event))

# 4 effort specs (climate kept as temp_grad_z — v1's headline climate proxy).
specs <- list(
  spec_A = list(label = "Record-based (legacy)",
                effort_col = "log_effort_record_z"),
  spec_B = list(label = "Observer visits",
                effort_col = "log_effort_visits_z"),
  spec_C = list(label = "PCA composite",
                effort_col = "effort_pc1_z"),
  spec_D = list(label = "Birding days",
                effort_col = "log_effort_days_z"))

# ---- 2. Fit & extract -----------------------------------------------------
fit_extract <- function(dt, spec_id, spec_label, effort_col) {
  if (!effort_col %in% names(dt)) {
    log("MISS spec=", spec_id, " (column ", effort_col, " absent)")
    return(NULL)
  }
  d <- dt[, .(species, province, year, event,
               climate_z = temp_grad_z,
               effort_z  = get(effort_col))]
  d <- d[!is.na(event) & !is.na(climate_z) & !is.na(effort_z)]
  log("=== ", spec_id, " (", spec_label, ") | n=", nrow(d),
      " events=", sum(d$event))
  forms <- list(
    M0 = "event ~ 1                          + (1|species) + (1|province)",
    M1 = "event ~ effort_z                    + (1|species) + (1|province)",
    M2 = "event ~ climate_z                   + (1|species) + (1|province)",
    M3 = "event ~ climate_z + effort_z        + (1|species) + (1|province)",
    M4 = "event ~ climate_z * effort_z        + (1|species) + (1|province)")
  coef_rows <- list()
  aic_rows  <- list()
  for (nm in names(forms)) {
    t0 <- Sys.time()
    fit <- tryCatch(
      glmmTMB::glmmTMB(as.formula(forms[[nm]]),
                        data = d,
                        family = binomial(link = "cloglog")),
      error = function(e) { log("  ", nm, " FAILED: ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    dt_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    log(sprintf("  %s OK (%.1fs, AIC=%.2f, nobs=%d)", nm, dt_sec,
                 AIC(fit), nobs(fit)))
    cf <- glmmTMB::fixef(fit)$cond
    se <- sqrt(diag(stats::vcov(fit)$cond))
    for (tm in names(cf)) {
      i <- match(tm, names(cf))
      coef_rows[[length(coef_rows) + 1L]] <- data.table(
        spec_id = spec_id, spec_label = spec_label,
        effort_var = effort_col, model = nm, term = tm,
        beta = cf[i], se = se[i],
        hr = exp(cf[i]),
        hr.low  = exp(cf[i] - 1.96 * se[i]),
        hr.high = exp(cf[i] + 1.96 * se[i]),
        p.value = 2 * pnorm(-abs(cf[i] / se[i])))
    }
    aic_rows[[length(aic_rows) + 1L]] <- data.table(
      spec_id = spec_id, spec_label = spec_label,
      effort_var = effort_col, model = nm,
      AIC = AIC(fit), BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)),
      df = attr(logLik(fit), "df"),
      nobs = nobs(fit),
      converged = fit$sdr$pdHess %||% NA)
    invisible(gc(verbose = FALSE))
  }
  list(coefs = rbindlist(coef_rows, fill = TRUE),
       aic   = rbindlist(aic_rows,  fill = TRUE))
}
"%||%" <- function(a, b) if (is.null(a)) b else a

all_coefs <- list()
all_aic   <- list()
for (sid in names(specs)) {
  res <- fit_extract(risk, sid, specs[[sid]]$label, specs[[sid]]$effort_col)
  if (!is.null(res)) {
    all_coefs[[sid]] <- res$coefs
    all_aic[[sid]]   <- res$aic
  }
}
coefs_v2 <- rbindlist(all_coefs, fill = TRUE)
aic_v2   <- rbindlist(all_aic,   fill = TRUE)

# Compute dAIC vs best within each spec.
aic_v2[, dAIC := AIC - min(AIC, na.rm = TRUE), by = spec_id]
setorder(aic_v2, spec_id, AIC)

fwrite(coefs_v2, file.path(V2, "results", "tables",
                            "table_province_v2_coefs.csv"))
fwrite(aic_v2,   file.path(V2, "results", "tables",
                            "table_province_v2_aic.csv"))
log("wrote table_province_v2_{coefs,aic}.csv")

# ---- 3. Reconcile against v1 ----------------------------------------------
v1_coefs_path <- file.path(V1, "results",
                           "table_cross_specification_key_coefficients.csv")
v1_aic_path   <- file.path(V1, "results",
                           "table_unified_model_comparison.csv")
recon <- NULL
if (file.exists(v1_coefs_path)) {
  v1_c <- fread(v1_coefs_path)
  # Extract M4 interaction rows from both
  v1_int <- v1_c[grepl(":", term) & model == "M4",
                  .(spec_id = spec, term, v1_beta = estimate, v1_se = std.error,
                    v1_hr = hr, v1_hr_low = hr_lower, v1_hr_high = hr_upper,
                    v1_p = p.value)]
  v2_int <- coefs_v2[grepl(":", term) & model == "M4",
                      .(spec_id, term, v2_beta = beta, v2_se = se,
                        v2_hr = hr, v2_hr_low = hr.low, v2_hr_high = hr.high,
                        v2_p = p.value)]
  recon <- merge(v1_int, v2_int, by = "spec_id", suffixes = c("_v1","_v2"),
                  all = TRUE)
  recon[, delta_beta := v2_beta - v1_beta]
  recon[, delta_hr   := v2_hr - v1_hr]
  recon[, abs_delta_hr := abs(delta_hr)]
  setcolorder(recon, c("spec_id", "term.x", "v1_beta", "v2_beta", "delta_beta",
                        "v1_hr", "v2_hr", "delta_hr", "abs_delta_hr",
                        "v1_p", "v2_p"))
  fwrite(recon, file.path(V2, "results", "tables",
                           "table_province_v1_v2_reconciliation.csv"))
  log("wrote v1-vs-v2 reconciliation. Max |Δ HR| = ",
      round(max(recon$abs_delta_hr, na.rm = TRUE), 4))
}

# ---- 4. Forest figure: v1 vs v2 interaction HR ----------------------------
if (!is.null(recon) && nrow(recon) > 0L) {
  long <- rbindlist(list(
    recon[, .(spec_id, label = paste0(spec_id, " (v1)"),
              hr = v1_hr, hr.low = v1_hr_low, hr.high = v1_hr_high,
              source = "v1")],
    recon[, .(spec_id, label = paste0(spec_id, " (v2)"),
              hr = v2_hr, hr.low = v2_hr_low, hr.high = v2_hr_high,
              source = "v2")]))
  long[, label := factor(label, levels = unique(label))]
  p_fig <- ggplot(long, aes(x = hr, y = label, colour = source)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                    height = 0.18, linewidth = 0.5) +
    geom_point(size = 2.4) +
    geom_text(aes(label = sprintf("HR = %.3f", hr)),
               nudge_x = 0.04, size = 2.6, hjust = 0) +
    scale_colour_manual(values = c(v1 = "#3B4CC0", v2 = "#B40426"),
                         name = "Run") +
    scale_x_continuous(trans = "log",
                        breaks = c(0.9, 1.0, 1.1, 1.2, 1.3, 1.5)) +
    labs(title = "Province-scale climate × effort interaction (M4) — v1 vs v2 refit",
          subtitle = "All 4 effort specifications. Bars = 95 % CI; dashed line = HR 1.0 (no effect).",
          x = "Hazard ratio (log scale)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          plot.margin = margin(8, 22, 8, 8))
  ggsave(file.path(V2, "figures", "main",
                    "fig_province_interaction_forest_v1_vs_v2.pdf"),
         p_fig, width = 16, height = 9, units = "cm",
         device = grDevices::cairo_pdf)
  ggsave(file.path(V2, "figures", "main",
                    "fig_province_interaction_forest_v1_vs_v2.png"),
         p_fig, width = 16, height = 9, units = "cm", dpi = 600)
  log("wrote fig_province_interaction_forest_v1_vs_v2.{pdf,png}")
}

log("=== DONE ===")
log("M4 interaction summary (v1 vs v2):")
if (!is.null(recon)) print(recon[, .(spec_id, v1_hr = round(v1_hr,3),
                                      v2_hr = round(v2_hr,3),
                                      delta_hr = round(delta_hr, 4),
                                      v1_p = signif(v1_p,3),
                                      v2_p = signif(v2_p,3))])
log("AIC ladder (v2):")
print(aic_v2[, .(spec_id, model, AIC = round(AIC,2), dAIC = round(dAIC,2))])
