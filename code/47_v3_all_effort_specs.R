# ============================================================
# Script: 47_v3_all_effort_specs.R
# Family: v3 risk-set — 4 effort specifications × M0-M4
# Author: Chen-Chen Ding + Claude Opus 4.7
# Date  : 2026-05-23
#
# ------------------------------------------------------------
# Scientific question / 科学问题:
#   The v3 relaxed risk set (script 46) has been refit only with
#   Spec B (observer-visits) as the effort variable. Reviewers will
#   want the same robustness check that v1 + v2 did: do all four
#   effort specifications still yield a positive, significant
#   climate × effort interaction on the larger v3 risk set?
#
#   v3 风险集 (188,870 行 / 463 物种 / 817 事件) 在 4 个 effort 指标
#   下，climate × effort 交互是否依旧稳健？
#
# Reuses sdm_province_v3_relaxed.csv (8,932 candidate pairs) as
# the inclusion source. Effort variables come from v1's
# effort_panel_upgraded.csv (province × year).
#
# Outputs (NEW files, do not overwrite v1/v2/v3-spec-B):
#   results/tables/table_province_v3_all_specs_coefs.csv
#   results/tables/table_province_v3_all_specs_aic.csv
#   figures/main/Figure_2_province_headline_v3_all_specs.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
options(warn = 1); set.seed(42)

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- normalizePath(file.path(V2, "..", "bird_hazard_model_effort_upgrade"),
                     mustWork = FALSE)

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE,
                                                    showWarnings = FALSE)
ens(file.path(V2, "logs"))
LOG <- file(file.path(V2, "logs", "47_v3_all_specs_audit.log"), "wt",
            encoding = "UTF-8")
on.exit({ if (isOpen(LOG)) close(LOG) }, add = TRUE)
log <- function(...) {
  m <- paste0(sprintf("[47 %s] ", format(Sys.time(), "%H:%M:%S")),
              paste(..., sep = ""))
  cat(m, "\n", sep = ""); writeLines(m, LOG)
}
audit <- function(t) { bar <- paste(rep("─", 60), collapse = "")
  log(""); log(bar); log("AUDIT — ", t); log(bar) }

# ============================================================
# STEP 1 — Load v3 risk set + effort panel + species filter
# ============================================================
audit("STEP 1: load v3 risk-set scaffold + effort + climate")

# 1.1 v3 candidate set (from script 46)
v3_cand <- fread(file.path(V2, "data", "derived",
                             "sdm_province_v3_relaxed.csv"),
                   encoding = "UTF-8")
log("1.1 v3 candidate set: ", nrow(v3_cand),
    " (sp × prov) pairs / ", uniqueN(v3_cand$species), " species")

# 1.2 events 2002-2024 for first-arrival calculation
ev <- fread(file.path(V1, "data", "events_100km_grid_assigned.csv"),
            encoding = "UTF-8")
setnames(ev, tolower(names(ev)))
if (!"year" %in% names(ev) && "pub_year" %in% names(ev))
  setnames(ev, "pub_year", "year")
ev <- ev[year >= 2002 & year <= 2024]
first_arrival <- ev[, .(arrival_year = min(year, na.rm = TRUE)),
                    by = .(species, province)]
log("1.2 first-arrival pairs: ", nrow(first_arrival))

# 1.3 effort panel — all 4 spec columns
prov_eff <- fread(file.path(V1, "data", "effort_panel_upgraded.csv"),
                    encoding = "UTF-8")
eff_cols <- intersect(c("log_effort_record_z","log_effort_visits_z",
                          "effort_pc1_z","log_effort_days_z"),
                       names(prov_eff))
log("1.3 effort spec columns available: ", paste(eff_cols, collapse=", "))

# 1.4 climate panel — temp_grad_prov_z (renamed to temp_grad_z)
prov_clim <- fread(file.path(V1, "data", "climate_metrics_province_year.csv"),
                    encoding = "UTF-8")
clim_use <- prov_clim[, .(province, year,
                            temp_grad_z = temp_grad_prov_z)]

# ============================================================
# STEP 2 — Build full v3 (species, province, year) risk frame
# ============================================================
audit("STEP 2: expand candidate × year, mark events, drop post-arrival")

yrs <- 2002:2024
risk_full <- v3_cand[, .(species, province)][
  rep(seq_len(.N), each = length(yrs))]
risk_full[, year := rep(yrs, times = nrow(v3_cand))]
risk_full <- merge(risk_full, first_arrival,
                    by = c("species","province"), all.x = TRUE)
risk_full <- risk_full[is.na(arrival_year) | year <= arrival_year]
risk_full[, event := as.integer(year == arrival_year)]
risk_full[is.na(event), event := 0L]

# Attach climate
risk_full <- merge(risk_full, clim_use,
                    by = c("province","year"), all.x = TRUE)
# Attach all four effort metrics
risk_full <- merge(risk_full,
                    prov_eff[, c("province","year", eff_cols), with = FALSE],
                    by = c("province","year"), all.x = TRUE)

log("2.1 risk_full: ", nrow(risk_full), " rows | events ",
    sum(risk_full$event))

# Complete-case across temp_grad_z + every effort column we need.
keep_cols <- c("species","province","year","event","temp_grad_z", eff_cols)
risk_full <- risk_full[complete.cases(risk_full[, ..keep_cols])]
log("2.2 after complete.case: ", nrow(risk_full), " rows | events ",
    sum(risk_full$event), " | species ", uniqueN(risk_full$species))

# ============================================================
# STEP 3 — Fit M0-M4 across 4 effort specifications
# ============================================================
audit("STEP 3: fit M0-M4 for each of the 4 effort specs (20 fits)")

specs <- list(
  spec_A = list(label = "Record-based (legacy)",
                col   = "log_effort_record_z"),
  spec_B = list(label = "Observer visits (headline)",
                col   = "log_effort_visits_z"),
  spec_C = list(label = "PCA composite",
                col   = "effort_pc1_z"),
  spec_D = list(label = "Birding days",
                col   = "log_effort_days_z"))

fit_one_spec <- function(spec_id, info) {
  d <- risk_full[, .(species, province, year, event,
                       climate_z = temp_grad_z,
                       effort_z  = get(info$col))]
  d <- d[complete.cases(d)]
  log("=== ", spec_id, " (", info$label, ") | n = ", nrow(d),
      " | events ", sum(d$event), " ===")
  forms <- list(
    M0 = "event ~ 1                          + (1|species) + (1|province)",
    M1 = "event ~ effort_z                    + (1|species) + (1|province)",
    M2 = "event ~ climate_z                   + (1|species) + (1|province)",
    M3 = "event ~ climate_z + effort_z        + (1|species) + (1|province)",
    M4 = "event ~ climate_z * effort_z        + (1|species) + (1|province)")
  cf_rows <- list(); ai_rows <- list()
  for (nm in names(forms)) {
    t0 <- Sys.time()
    fit <- tryCatch(
      glmmTMB(as.formula(forms[[nm]]), data = d,
              family = binomial(link = "cloglog")),
      error = function(e) { log("  ", nm, " FAILED: ",
                                  conditionMessage(e)); NULL })
    if (is.null(fit)) next
    s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    log(sprintf("  %s OK (%.1fs, AIC = %.2f, nobs = %d)",
                nm, s, AIC(fit), nobs(fit)))
    cfx <- fixef(fit)$cond; se <- sqrt(diag(stats::vcov(fit)$cond))
    for (tm in names(cfx)) {
      i <- match(tm, names(cfx))
      cf_rows[[length(cf_rows)+1L]] <- data.table(
        spec_id = spec_id, spec_label = info$label,
        effort_var = info$col, model = nm, term = tm,
        beta = cfx[i], se = se[i],
        hr = exp(cfx[i]),
        hr.low = exp(cfx[i] - 1.96 * se[i]),
        hr.high = exp(cfx[i] + 1.96 * se[i]),
        p.value = 2 * pnorm(-abs(cfx[i] / se[i])))
    }
    ai_rows[[length(ai_rows)+1L]] <- data.table(
      spec_id = spec_id, spec_label = info$label,
      effort_var = info$col, model = nm,
      AIC = AIC(fit), BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)), nobs = nobs(fit))
    invisible(gc(verbose = FALSE))
  }
  list(coefs = rbindlist(cf_rows, fill = TRUE),
       aic   = rbindlist(ai_rows, fill = TRUE))
}

all_co <- list(); all_ai <- list()
for (sid in names(specs)) {
  out <- fit_one_spec(sid, specs[[sid]])
  if (!is.null(out)) { all_co[[sid]] <- out$coefs;  all_ai[[sid]] <- out$aic }
}
coefs_all <- rbindlist(all_co, fill = TRUE)
aic_all   <- rbindlist(all_ai, fill = TRUE)
aic_all[, dAIC := AIC - min(AIC, na.rm = TRUE), by = spec_id]
aic_all[, akaike_weight := exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC)),
        by = spec_id]
setorder(aic_all, spec_id, AIC)

fwrite(coefs_all,
       file.path(V2, "results", "tables",
                  "table_province_v3_all_specs_coefs.csv"))
fwrite(aic_all,
       file.path(V2, "results", "tables",
                  "table_province_v3_all_specs_aic.csv"))
log("wrote table_province_v3_all_specs_{coefs,aic}.csv")

# ============================================================
# STEP 4 — Forest comparing v1 / v2 / v3 across 4 specs
# ============================================================
audit("STEP 4: build 4-spec × 3-run forest plot")

v1_coefs <- fread(file.path(V1, "results",
                              "table_cross_specification_key_coefficients.csv"))
v2_coefs <- fread(file.path(V2, "results", "tables",
                              "table_province_v2_coefs.csv"))
v3_int <- coefs_all[model == "M4" & grepl(":", term),
                      .(spec_id, run = "v3", hr, hr.low, hr.high, p.value)]
v2_int <- v2_coefs[model == "M4" & grepl(":", term),
                     .(spec_id, run = "v2", hr, hr.low, hr.high, p.value)]
v1_int <- v1_coefs[model == "M4" & grepl(":", term),
                     .(spec_id = spec, run = "v1",
                       hr, hr.low = hr_lower, hr.high = hr_upper,
                       p.value)]
forest_dt <- rbindlist(list(v1_int, v2_int, v3_int), fill = TRUE)
SPEC_LBL <- c(spec_A = "A: records",
              spec_B = "B: visits (headline)",
              spec_C = "C: PCA composite",
              spec_D = "D: birding-days")
forest_dt[, spec_lbl := factor(SPEC_LBL[spec_id], levels = rev(SPEC_LBL))]
forest_dt[, run := factor(run, levels = c("v1","v2","v3"))]

theme_pub <- function(s = 9) {
  theme_bw(base_size = s) +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(linewidth = 0.4, colour = "grey20"),
          plot.title = element_text(face = "bold", size = s + 1),
          plot.subtitle = element_text(size = s - 1, colour = "grey30"))
}
COL_RUN <- c(v1 = "#3B4CC0", v2 = "#7F7F7F", v3 = "#B40426")

p_forest <- ggplot(forest_dt,
                    aes(x = hr, y = spec_lbl, colour = run)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  height = 0.18, linewidth = 0.5,
                  position = position_dodge(width = 0.55)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = COL_RUN, name = "Run") +
  scale_x_continuous(trans = "log",
                      breaks = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5)) +
  labs(tag = "a",
        title = "Climate × effort interaction (M4) — 4 effort specs × 3 runs",
        subtitle = "v3 (red): relaxed threshold + event override + 501 modelled species. Bars = 95 % CI.",
        x = "Hazard ratio (log scale)", y = NULL) +
  theme_pub() + theme(legend.position = "top")

# AIC ladder for v3 only
aic_show <- copy(aic_all)
aic_show[, model := factor(model, levels = c("M0","M1","M2","M3","M4"))]
aic_show[, spec_lbl := factor(SPEC_LBL[spec_id], levels = SPEC_LBL)]
p_aic <- ggplot(aic_show, aes(x = dAIC, y = model)) +
  geom_segment(aes(xend = 0, yend = model), linewidth = 0.35,
                 colour = "#B40426") +
  geom_point(size = 2.2, colour = "#B40426") +
  facet_wrap(~ spec_lbl, ncol = 2) +
  labs(tag = "b", title = "v3 AIC ladder by effort spec",
        x = "ΔAIC vs best", y = NULL) +
  theme_pub()

# Sample size table panel
size_dt <- data.table(
  Run = c("v1","v2","v3"),
  `Risk rows` = c(12813, 12813, nrow(risk_full)),
  `Species` = c(333, 333, uniqueN(risk_full$species)),
  `Events` = c(512, 512, sum(risk_full$event)))
size_long <- melt(size_dt, id.vars = "Run", variable.name = "metric",
                    value.name = "n")
p_size <- ggplot(size_long, aes(x = Run, y = n, fill = Run)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = scales::comma(n)), vjust = -0.4,
             size = 2.5, fontface = "bold") +
  scale_y_continuous(trans = "log10", labels = scales::comma_format(),
                      expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_manual(values = COL_RUN, guide = "none") +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  labs(tag = "c", title = "Sample-size scaling across the three runs",
        x = NULL, y = "Count (log)") +
  theme_pub()

fig <- (p_forest / p_size) | p_aic
fig <- fig + plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "Figure 2 v3 — Province headline robustness across 4 effort specs and 3 risk-set runs",
    theme = theme(plot.title = element_text(face = "bold", size = 10)))

ggsave(file.path(V2, "figures", "main",
                  "Figure_2_province_headline_v3_all_specs.pdf"),
       fig, width = 20, height = 14, units = "cm",
       device = grDevices::cairo_pdf)
ggsave(file.path(V2, "figures", "main",
                  "Figure_2_province_headline_v3_all_specs.png"),
       fig, width = 20, height = 14, units = "cm", dpi = 600)
log("wrote Figure_2_province_headline_v3_all_specs.{pdf,png}")

log("")
log("══════════════════════════════════════════════════════════")
log("              v3 × 4 SPEC × M0-M4 COMPLETE")
log("══════════════════════════════════════════════════════════")
print(coefs_all[model == "M4" & grepl(":", term),
                  .(spec_id,
                    HR = round(hr, 3),
                    `95%CI` = sprintf("%.3f, %.3f", hr.low, hr.high),
                    p = signif(p.value, 3))])
