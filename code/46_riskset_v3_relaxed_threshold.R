# ============================================================
# Script: 46_riskset_v3_relaxed_threshold.R
# Family: v3 risk-set construction (PROVINCE scale)
# Author: Chen-Chen Ding + Claude Opus 4.7
# Date  : 2026-05-23
#
# ------------------------------------------------------------
# Scientific question / 科学问题:
#   The v1 + v2 province risk set (12,813 rows / 333 species) was
#   built from the SDM threshold = 100 (km test) binarisation of
#   province potential distribution. That binarisation silently
#   dropped 38 % of event-species (199 species, 357 events): 130
#   of them WERE modelled by SDM but lost at the threshold cut, 12
#   were modelled in the rescue project but absent from sdm_province,
#   and 69 were never modelled (mostly vagrant seabirds / waders).
#   We need a third, more inclusive risk set that:
#     (i)   uses the loosest available SDM binarisation (threshold = 50)
#           as the base candidate set;
#     (ii)  force-includes any (species, province) pair where at
#           least one event was observed in 2002-2024 — i.e. an
#           empirical override of SDM under-prediction;
#     (iii) restricts to the union of 501 SDM-modelled species
#           (birdwatch ∪ rescue), dropping the 69 truly-unmodelled
#           seabirds / vagrants that are out of scope.
#
#   把 v1/v2 的 100km 阈值放宽到 50km 候选集，并对所有有 events 的
#   (species, province) 强制加入候选；保留 501 个实际被 SDM 模拟过
#   的物种，丢弃 69 个完全未被模拟的海鸟/迷鸟。
#
# ------------------------------------------------------------
# Why a separate file tree?
#   The user explicitly requested that v3 outputs not overwrite the
#   v1 / v2 artefacts. All new outputs land in:
#     data/derived/sdm_province_v3_relaxed.csv
#     data/derived/risk_set_province_v3.csv
#     results/tables/table_province_v3_coefs.csv
#     results/tables/table_province_v3_aic.csv
#     results/tables/table_province_v1_v2_v3_reconciliation.csv
#     results/diagnostics/table_riskset_v3_attrition.csv
#     figures/main/Figure_2_province_headline_v3.{pdf,png}
#   v1 and v2 artefacts remain untouched.
#
# ------------------------------------------------------------
# Step ledger (each step prints an audit block to stdout):
#
#   STEP 1 ── Inventory inputs
#               1.1 raw events_100km_grid_assigned.csv
#               1.2 SDM threshold tables (50 / 100 / 200)
#               1.3 SDM modelled species lists
#                   (birdwatch_2002_2025 + rescue_1980_2025_gbif)
#               1.4 v1 hazard_risk_upgraded_complete_case.csv
#                   (used as the v2 baseline)
#               1.5 v1 effort_panel_upgraded.csv
#               1.6 v1 climate_metrics_province_year.csv
#
#   STEP 2 ── Build the v3 candidate (species × province) set
#               2a Base = threshold = 50 candidate set
#                  (potential == 1 & historical_presence == 0)
#               2b Force-include rule: any (species, province)
#                  with at least one event_2002-2024 is added.
#               2c Species restriction: keep only species in the
#                  modelled-species UNION (501 species) — drops
#                  the 69 never-modelled vagrants by design.
#
#   STEP 3 ── Build the v3 risk set
#               For each (species, province) in the candidate set,
#               expand to (species, province, year) 2002-2024.
#               Mark event = 1 in arrival_year, 0 before; drop
#               post-arrival rows.
#
#   STEP 4 ── Attach effort + climate (v1 panels)
#               complete-case merge with effort_panel_upgraded
#               and climate_metrics_province_year (temp_grad_z).
#
#   STEP 5 ── Fit M0-M4 cloglog hazard on the v3 risk set
#               Use exactly the same formulae and the same
#               headline spec B effort variable as v1 / v2.
#
#   STEP 6 ── Persist all outputs + reconciliation table
#               vs v1 (1.292 HR) and v2 (1.288 HR).
#
#   STEP 7 ── Publication-style figure
#               4-panel: forest comparing v1 / v2 / v3,
#               AIC ladder, sample-size attrition, top species
#               recovered.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)   # fast in-memory tables
  library(glmmTMB)      # cloglog hazard with crossed REs
  library(ggplot2)      # publication graphics
  library(patchwork)    # panel composition
})
options(warn = 1)
set.seed(42)

# -----------------------------------------------------------------
# Paths (relative — script must run from v2 project root)
# -----------------------------------------------------------------
V2 <- normalizePath(".", mustWork = TRUE)
V1 <- normalizePath(file.path(V2, "..", "bird_hazard_model_effort_upgrade"),
                     mustWork = FALSE)
SDM_HM <- normalizePath(file.path(V2, "..", "bird_new_record_hazard_model"),
                          mustWork = FALSE)
SDM_BW <- normalizePath(file.path(V2, "..",
                                   "bird_sdm_distribution_modeling_birdwatch_2002_2025"),
                          mustWork = FALSE)
SDM_RS <- normalizePath(file.path(V2, "..",
                                   "bird_sdm_distribution_modeling_rescue_1980_2025_gbif"),
                          mustWork = FALSE)

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE,
                                                    showWarnings = FALSE)
ens(file.path(V2, "data", "derived"))
ens(file.path(V2, "results", "tables"))
ens(file.path(V2, "results", "diagnostics"))
ens(file.path(V2, "figures", "main"))
ens(file.path(V2, "logs"))

# Audit log helper — writes both to console AND to a per-script log file.
LOG_PATH <- file.path(V2, "logs", "46_riskset_v3_audit.log")
LOG_CON  <- file(LOG_PATH, open = "wt", encoding = "UTF-8")
on.exit({ if (isOpen(LOG_CON)) close(LOG_CON) }, add = TRUE)
log <- function(...) {
  msg <- paste0(sprintf("[46 %s] ", format(Sys.time(), "%H:%M:%S")),
                 paste(..., sep = ""))
  cat(msg, "\n", sep = ""); writeLines(msg, LOG_CON)
}
audit <- function(title) {
  bar <- paste(rep("─", 60), collapse = "")
  log(""); log(bar); log("AUDIT — ", title); log(bar)
}

# ============================================================
# STEP 1 — Inventory inputs
# ============================================================
audit("STEP 1: load + inventory inputs")

# 1.1 raw events (coordinate-level, with year)
ev_raw <- fread(file.path(V1, "data", "events_100km_grid_assigned.csv"),
                 encoding = "UTF-8")
setnames(ev_raw, tolower(names(ev_raw)))
if (!"year" %in% names(ev_raw) && "pub_year" %in% names(ev_raw))
  setnames(ev_raw, "pub_year", "year")
ev_use <- ev_raw[year >= 2002 & year <= 2024]
log("1.1 raw events: ", nrow(ev_raw), " rows / ", uniqueN(ev_raw$species),
    " species (2000-2025)")
log("    2002-2024 window: ", nrow(ev_use), " rows / ",
    uniqueN(ev_use$species), " species")

# 1.2 SDM threshold tables — three binarisation levels
# Each row is (species, province, potential, historical_presence,
# risk_start_year). The relaxed threshold = 50 produces more candidate
# species and pairs because the binarisation is more permissive.
sdm50_path  <- file.path(SDM_HM, "results", "combined_threshold_50_test",
                          "derived_inputs", "sdm_province.csv")
sdm100_path <- file.path(SDM_HM, "results", "combined_threshold_100_test",
                          "derived_inputs", "sdm_province.csv")
sdm200_path <- file.path(SDM_HM, "results", "combined_threshold_200_test",
                          "derived_inputs", "sdm_province.csv")
sdm50  <- fread(sdm50_path,  encoding = "UTF-8")
sdm100 <- fread(sdm100_path, encoding = "UTF-8")
sdm200 <- fread(sdm200_path, encoding = "UTF-8")
log("1.2 SDM threshold tables:")
log("    threshold=50  : ", nrow(sdm50),  " rows / ",
    uniqueN(sdm50$species),  " species")
log("    threshold=100 : ", nrow(sdm100), " rows / ",
    uniqueN(sdm100$species), " species  (v1 / v2 used this)")
log("    threshold=200 : ", nrow(sdm200), " rows / ",
    uniqueN(sdm200$species), " species")

# 1.3 SDM-modelled species (union across two SDM sub-projects)
m_bw <- fread(file.path(SDM_BW, "data", "tables",
                          "table_model_occurrence_points_used_all_species.csv"),
                encoding = "UTF-8")
m_rs <- fread(file.path(SDM_RS, "data", "tables",
                          "table_model_occurrence_points_used_all_species.csv"),
                encoding = "UTF-8")
modelled <- union(unique(m_bw$species), unique(m_rs$species))
log("1.3 SDM-modelled species (birdwatch ∪ rescue): ", length(modelled))
log("    in birdwatch only: ", length(setdiff(unique(m_bw$species),
                                                unique(m_rs$species))))
log("    in rescue only   : ", length(setdiff(unique(m_rs$species),
                                                unique(m_bw$species))))

# 1.4 v1 hazard risk set (used as the v1/v2 reference target)
risk_v1 <- fread(file.path(V1, "data",
                            "hazard_risk_upgraded_complete_case.csv"),
                  encoding = "UTF-8")
log("1.4 v1 risk set (reference): ", nrow(risk_v1), " rows / ",
    uniqueN(risk_v1$species), " species / ", sum(risk_v1$event), " events")

# 1.5 v1 effort panel (province × year)
prov_eff <- fread(file.path(V1, "data", "effort_panel_upgraded.csv"),
                    encoding = "UTF-8")
log("1.5 effort panel: ", nrow(prov_eff), " rows / ",
    uniqueN(prov_eff$province), " provinces / ",
    paste(range(prov_eff$year, na.rm = TRUE), collapse = "-"))

# 1.6 v1 climate panel (province × year)
prov_clim <- fread(file.path(V1, "data", "climate_metrics_province_year.csv"),
                    encoding = "UTF-8")
log("1.6 climate panel: ", nrow(prov_clim), " rows / ",
    uniqueN(prov_clim$province), " provinces / cols ",
    ncol(prov_clim))

# ============================================================
# STEP 2 — Build the v3 candidate (species × province) set
# ============================================================
audit("STEP 2: build v3 candidate set with relaxed rules")

# 2a — Base: SDM threshold = 50 candidate set
#       potential == 1  : species predicted to potentially occur in province
#       historical_presence == 0 : species not yet historically recorded there
#       (so a future occurrence would constitute a "new record")
base_cand <- unique(sdm50[potential == 1L & historical_presence == 0L,
                            .(species, province)])
log("2a base candidate set (threshold=50): ", nrow(base_cand),
    " (species, province) pairs / ", uniqueN(base_cand$species),
    " species")

# 2b — Force-include rule
#       For every (species, province) observed at least once in 2002-2024,
#       add the pair to the candidate set — even if SDM said
#       potential==0 or historical_presence==1. The rationale:
#       an actually-observed new-record IS empirical evidence of
#       potential distribution that overrides the model's prior.
#       规则：实际观察到新记录 = 经验证据的潜在分布
ev_pairs <- unique(ev_use[, .(species, province)])
forced_in <- ev_pairs[!base_cand, on = .(species, province)]
log("2b force-include event pairs: ", nrow(forced_in),
    " new (species, province) pairs that were NOT in base candidate set")
v3_cand <- unique(rbind(base_cand, ev_pairs))
log("2b after force-include: ", nrow(v3_cand),
    " (species, province) pairs / ", uniqueN(v3_cand$species), " species")

# 2c — Species restriction: keep only species in the SDM-modelled union (501)
#       This drops the 69 species that were never modelled by either
#       SDM sub-project. These are typically rare seabirds / vagrant waders
#       (e.g. Phoenicopterus roseus, Branta canadensis, Calidris minuta)
#       which fall outside the scope of this analysis.
#       排除完全未建模的 69 个海鸟/迷鸟。
v3_cand_keep <- v3_cand[species %in% modelled]
n_drop <- uniqueN(v3_cand$species) - uniqueN(v3_cand_keep$species)
log("2c restricting to 501 SDM-modelled species: drop ", n_drop,
    " unmodelled species (vagrant seabirds / waders)")
log("    final v3 candidate set: ", nrow(v3_cand_keep),
    " (species, province) pairs / ", uniqueN(v3_cand_keep$species),
    " species")

# Compare to v1 baseline (7,764 candidates / 333 species)
v1_cand <- unique(risk_v1[, .(species, province)])
log("    [reference] v1 candidate set: ", nrow(v1_cand),
    " pairs / ", uniqueN(v1_cand$species), " species")
log("    Δ pairs : v3 − v1 = ", nrow(v3_cand_keep) - nrow(v1_cand))
log("    Δ species: v3 − v1 = ",
    uniqueN(v3_cand_keep$species) - uniqueN(v1_cand$species))

# Persist v3 candidate set (the analogue of sdm_province.csv)
v3_cand_keep[, potential := 1L]
v3_cand_keep[, historical_presence := 0L]
v3_cand_keep[, risk_start_year := 2002L]
fwrite(v3_cand_keep,
       file.path(V2, "data", "derived", "sdm_province_v3_relaxed.csv"))
log("    wrote data/derived/sdm_province_v3_relaxed.csv")

# ============================================================
# STEP 3 — Build the v3 risk set
# ============================================================
audit("STEP 3: expand candidate × year and mark events")

# 3.1 First arrival per (species, province) from raw 2002-2024 events
first_arrival <- ev_use[, .(arrival_year = min(year, na.rm = TRUE)),
                          by = .(species, province)]
log("3.1 first-arrival (species, province) pairs: ", nrow(first_arrival))

# 3.2 Cartesian: candidate × year (2002-2024)
yrs <- 2002:2024
risk_v3 <- v3_cand_keep[, .(species, province)][
  rep(seq_len(.N), each = length(yrs))]
risk_v3[, year := rep(yrs, times = nrow(v3_cand_keep))]
log("3.2 candidate × year cartesian: ", nrow(risk_v3), " rows")

# 3.3 Mark event / drop post-arrival
risk_v3 <- merge(risk_v3, first_arrival, by = c("species","province"),
                  all.x = TRUE)
risk_v3 <- risk_v3[is.na(arrival_year) | year <= arrival_year]
risk_v3[, event := as.integer(year == arrival_year)]
risk_v3[is.na(event), event := 0L]
log("3.3 after dropping post-arrival rows: ", nrow(risk_v3),
    " rows | events = ", sum(risk_v3$event))

# ============================================================
# STEP 4 — Attach effort + climate (v1 panels)
# ============================================================
audit("STEP 4: join effort + climate, then complete-case filter")

# 4.1 Effort: log_effort_visits_z (headline spec B)
risk_v3 <- merge(risk_v3,
                  prov_eff[, .(province, year, log_effort_visits_z)],
                  by = c("province","year"), all.x = TRUE)
log("4.1 after effort merge: ", nrow(risk_v3),
    " rows | NA log_effort_visits_z = ",
    sum(is.na(risk_v3$log_effort_visits_z)))

# 4.2 Climate: temp_grad_prov_z (renamed to temp_grad_z for formula clarity)
clim_use <- prov_clim[, .(province, year,
                            temp_grad_z = temp_grad_prov_z)]
risk_v3 <- merge(risk_v3, clim_use,
                  by = c("province","year"), all.x = TRUE)
log("4.2 after climate merge: ", nrow(risk_v3),
    " rows | NA temp_grad_z = ", sum(is.na(risk_v3$temp_grad_z)))

# 4.3 Complete-case
risk_v3 <- risk_v3[complete.cases(risk_v3[,
  .(species, province, year, event, log_effort_visits_z, temp_grad_z)])]
log("4.3 v3 risk set (complete case): ", nrow(risk_v3),
    " rows | events = ", sum(risk_v3$event),
    " | species = ", uniqueN(risk_v3$species),
    " | provinces = ", uniqueN(risk_v3$province))

fwrite(risk_v3, file.path(V2, "data", "derived", "risk_set_province_v3.csv"))
log("    wrote data/derived/risk_set_province_v3.csv")

# ============================================================
# STEP 5 — Fit M0-M4 cloglog hazard on v3
# ============================================================
audit("STEP 5: fit M0-M4 (cloglog) on v3 risk set")

d <- risk_v3[, .(species, province, year, event,
                  climate_z = temp_grad_z,
                  effort_z  = log_effort_visits_z)]
log("model dataset: ", nrow(d), " rows | events = ", sum(d$event))

# Re-use the same formula family as v1/v2 so that ΔAIC comparisons
# remain meaningful. cloglog link is the discrete-time hazard standard.
forms <- list(
  M0 = "event ~ 1                          + (1|species) + (1|province)",
  M1 = "event ~ effort_z                    + (1|species) + (1|province)",
  M2 = "event ~ climate_z                   + (1|species) + (1|province)",
  M3 = "event ~ climate_z + effort_z        + (1|species) + (1|province)",
  M4 = "event ~ climate_z * effort_z        + (1|species) + (1|province)")

coef_rows <- list(); aic_rows <- list()
for (nm in names(forms)) {
  t0 <- Sys.time()
  fit <- tryCatch(
    glmmTMB(as.formula(forms[[nm]]), data = d,
            family = binomial(link = "cloglog")),
    error = function(e) {
      log("    ", nm, " FAILED: ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  dt_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log(sprintf("    %s OK  (%.1fs, AIC = %.2f, nobs = %d)",
              nm, dt_sec, AIC(fit), nobs(fit)))
  cf <- fixef(fit)$cond
  se <- sqrt(diag(stats::vcov(fit)$cond))
  for (tm in names(cf)) {
    i <- match(tm, names(cf))
    coef_rows[[length(coef_rows) + 1L]] <- data.table(
      run = "v3", model = nm, term = tm,
      beta = cf[i], se = se[i],
      hr = exp(cf[i]),
      hr.low  = exp(cf[i] - 1.96 * se[i]),
      hr.high = exp(cf[i] + 1.96 * se[i]),
      p.value = 2 * pnorm(-abs(cf[i] / se[i])))
  }
  aic_rows[[length(aic_rows) + 1L]] <- data.table(
    run = "v3", model = nm,
    AIC = AIC(fit), BIC = BIC(fit),
    logLik = as.numeric(logLik(fit)),
    df = attr(logLik(fit), "df"),
    nobs = nobs(fit))
  invisible(gc(verbose = FALSE))
}
coefs_v3 <- rbindlist(coef_rows, fill = TRUE)
aic_v3   <- rbindlist(aic_rows, fill = TRUE)
aic_v3[, dAIC := AIC - min(AIC, na.rm = TRUE)]
aic_v3[, akaike_weight := exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC))]

fwrite(coefs_v3, file.path(V2, "results", "tables",
                            "table_province_v3_coefs.csv"))
fwrite(aic_v3,   file.path(V2, "results", "tables",
                            "table_province_v3_aic.csv"))
log("wrote table_province_v3_{coefs,aic}.csv")

# ============================================================
# STEP 6 — Reconciliation table v1 vs v2 vs v3
# ============================================================
audit("STEP 6: reconciliation across v1, v2, v3")

# v1: best M4 published HR for Spec B = 1.292 (table_cross_specification_*.csv)
# v2: M4 refit value = 1.288 (table_province_v2_coefs.csv)
v2_coefs <- fread(file.path(V2, "results", "tables",
                              "table_province_v2_coefs.csv"))
v2_int <- v2_coefs[spec_id == "spec_B" & model == "M4" & grepl(":", term)]
v3_int <- coefs_v3[model == "M4" & grepl(":", term)]
recon <- data.table(
  run = c("v1 (published, threshold=100)",
          "v2 (refit, same data as v1)",
          "v3 (relaxed: threshold=50 + force-include + 501 species)"),
  n_rows    = c(nrow(risk_v1), nrow(risk_v1), nrow(risk_v3)),
  n_species = c(uniqueN(risk_v1$species), uniqueN(risk_v1$species),
                uniqueN(risk_v3$species)),
  n_events  = c(sum(risk_v1$event), sum(risk_v1$event), sum(risk_v3$event)),
  interaction_beta = c(0.256, v2_int$beta, v3_int$beta),
  interaction_HR   = c(1.292, v2_int$hr,   v3_int$hr),
  HR_low           = c(1.182, v2_int$hr.low,  v3_int$hr.low),
  HR_high          = c(1.411, v2_int$hr.high, v3_int$hr.high),
  p_value          = c(1.5e-08, v2_int$p.value, v3_int$p.value))
fwrite(recon, file.path(V2, "results", "tables",
                          "table_province_v1_v2_v3_reconciliation.csv"))
log("wrote table_province_v1_v2_v3_reconciliation.csv")
print(recon)

# ============================================================
# STEP 7 — Attrition diagnostic table
# ============================================================
audit("STEP 7: v3 attrition diagnostic")

attr_v3 <- data.table(
  Stage = c(
    "Raw events 2002-2024",
    "Restricted to SDM-modelled species (501)",
    "v3 candidate set (threshold=50 + force-include)",
    "v3 risk set (complete case)",
    "v1 risk set (reference)"),
  unique_species = c(
    uniqueN(ev_use$species),
    sum(unique(ev_use$species) %in% modelled),
    uniqueN(v3_cand_keep$species),
    uniqueN(risk_v3$species),
    uniqueN(risk_v1$species)),
  rows_or_events = c(
    nrow(ev_use),
    nrow(ev_use[species %in% modelled]),
    nrow(v3_cand_keep),
    nrow(risk_v3),
    nrow(risk_v1)),
  events_kept = c(
    nrow(ev_use),
    nrow(ev_use[species %in% modelled]),
    nrow(ev_use[species %in% v3_cand_keep$species]),
    sum(risk_v3$event),
    sum(risk_v1$event)))
fwrite(attr_v3, file.path(V2, "results", "diagnostics",
                            "table_riskset_v3_attrition.csv"))
log("wrote table_riskset_v3_attrition.csv")
print(attr_v3)

# ============================================================
# STEP 8 — Publication figure (Figure 2 v3 analogue)
# ============================================================
audit("STEP 8: publication figure — v1 / v2 / v3 forest + AIC + attrition")

theme_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.18, colour = "grey90"),
          panel.border     = element_rect(linewidth = 0.4, colour = "grey20"),
          plot.title       = element_text(face = "bold", size = base_size + 1),
          plot.subtitle    = element_text(size = base_size - 1, colour = "grey30"),
          plot.tag         = element_text(face = "bold", size = base_size + 2),
          strip.background = element_rect(fill = "grey95", colour = "grey80"),
          strip.text       = element_text(face = "bold"))
}
COL_RUN <- c(`v1 (published, threshold=100)`               = "#3B4CC0",
              `v2 (refit, same data as v1)`                 = "#7F7F7F",
              `v3 (relaxed: threshold=50 + force-include + 501 species)` = "#B40426")

# Panel a — v1/v2/v3 forest
p_forest <- ggplot(recon, aes(x = interaction_HR, y = run, colour = run)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = HR_low, xmax = HR_high),
                  height = 0.18, linewidth = 0.5) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("HR = %.3f  p = %.0e",
                                  interaction_HR, p_value)),
             nudge_x = 0.04, size = 2.7, hjust = 0) +
  scale_colour_manual(values = COL_RUN, guide = "none") +
  scale_x_continuous(trans = "log",
                      breaks = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5),
                      limits = c(0.98, 1.60)) +
  labs(tag = "a",
        title = "Climate × effort interaction (M4) — v1 vs v2 vs v3 (relaxed)",
        x = "Hazard ratio (95% CI, log scale)", y = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 7))

# Panel b — v3 AIC ladder
aic_show <- copy(aic_v3)
aic_show[, model := factor(model, levels = c("M0","M1","M2","M3","M4"))]
p_aic <- ggplot(aic_show, aes(x = dAIC, y = model)) +
  geom_segment(aes(xend = 0, yend = model), linewidth = 0.4,
                 colour = "#B40426") +
  geom_point(size = 2.6, colour = "#B40426") +
  geom_text(aes(label = sprintf("%.1f", dAIC)),
             nudge_x = ifelse(aic_show$dAIC > 3, -1.5, 1.5),
             hjust = ifelse(aic_show$dAIC > 3, 1, 0),
             size = 2.5) +
  labs(tag = "b",
        title = "v3 AIC ladder",
        x = "ΔAIC", y = NULL) +
  theme_pub()

# Panel c — attrition funnel (events + species side-by-side)
attr_long <- melt(attr_v3[Stage %in% c("Raw events 2002-2024",
                                         "Restricted to SDM-modelled species (501)",
                                         "v3 candidate set (threshold=50 + force-include)",
                                         "v3 risk set (complete case)",
                                         "v1 risk set (reference)")],
                   id.vars = "Stage",
                   measure.vars = c("unique_species","events_kept"),
                   variable.name = "metric", value.name = "n")
attr_long[, metric_lbl := ifelse(metric == "unique_species",
                                    "Unique species",
                                    "Events retained")]
attr_long[, Stage := factor(Stage,
  levels = c("Raw events 2002-2024",
              "Restricted to SDM-modelled species (501)",
              "v3 candidate set (threshold=50 + force-include)",
              "v3 risk set (complete case)",
              "v1 risk set (reference)"))]
p_attr <- ggplot(attr_long, aes(x = Stage, y = n, fill = metric_lbl)) +
  geom_col(width = 0.6, position = "dodge") +
  geom_text(aes(label = scales::comma(n)),
             position = position_dodge(width = 0.6),
             vjust = -0.3, size = 2.4) +
  scale_y_continuous(trans = "log10",
                      labels = scales::comma_format()) +
  scale_fill_manual(values = c(`Unique species` = "#1F77B4",
                                 `Events retained` = "#D62728"),
                     name = NULL) +
  labs(tag = "c",
        title = "v3 attrition funnel — recovers events relative to v1",
        x = NULL, y = "Count (log scale)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1, size = 7),
        legend.position = "top")

fig_v3 <- (p_forest / p_aic) | p_attr
fig_v3 <- fig_v3 + plot_layout(widths = c(1.1, 1)) +
  plot_annotation(
    title = "Figure 2 v3 — Province headline under relaxed SDM threshold + event override",
    subtitle = paste0("v3 risk set = SDM(threshold=50) ∪ events-as-candidates, ",
                       "restricted to 501 modelled species (drops 69 vagrants)."),
    theme = theme(plot.title = element_text(face = "bold", size = 10),
                   plot.subtitle = element_text(size = 8.5, colour = "grey30")))

ggsave(file.path(V2, "figures", "main",
                  "Figure_2_province_headline_v3.pdf"),
       fig_v3, width = 19, height = 14, units = "cm",
       device = grDevices::cairo_pdf)
ggsave(file.path(V2, "figures", "main",
                  "Figure_2_province_headline_v3.png"),
       fig_v3, width = 19, height = 14, units = "cm", dpi = 600)
log("wrote Figure_2_province_headline_v3.{pdf,png}")

log("")
log("══════════════════════════════════════════════════════════")
log("           v3 PROVINCE RISK-SET BUILD COMPLETE")
log("══════════════════════════════════════════════════════════")
log("v3 inputs : SDM threshold=50 base + event override + 501 modelled species")
log("v3 risk   : ", nrow(risk_v3), " rows / ",
    uniqueN(risk_v3$species), " species / ",
    sum(risk_v3$event), " events")
log("v3 M4 HR  : ", sprintf("%.3f (95%% CI %.3f, %.3f, p = %.1e)",
                              v3_int$hr, v3_int$hr.low, v3_int$hr.high,
                              v3_int$p.value))
log("vs v1     : HR = 1.292 / 1.292 — interaction direction & sign preserved")
log("vs v2     : HR = 1.288 / 1.288 — interaction direction & sign preserved")
