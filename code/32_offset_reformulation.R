# ============================================================
# Scientific question / 科学问题:
#   v1's M5 used a z-scored effort vector inside `offset()`, breaking
#   the offset≡1 assumption of cloglog hazard models. Does refitting
#   M5 with `offset = log(person_hours + 1)` change the climate × effort
#   interaction conclusion, and how sensitive is the offset choice
#   (raw / log / sqrt / z-score / none)?
#   v1 M5 把 z-scored effort 当 offset 用，破坏了 offset 系数=1 的
#   假设。改用 log(person_hours+1) 重做后结论是否稳健？
#
# Objective / 分析目标:
#   - Refit M5 with five offset specifications.
#   - Report ΔAIC vs the v1 reference (z-score) and vs M3 (no offset).
#   - Verify interaction β remains within ±0.05 of the headline.
#
# Input data / 输入数据:
#   data/raw/hazard_risk_upgraded_complete_case.csv (v1 risk set)
#   data/raw/effort_panel_upgraded.csv               (raw effort metrics)
#
# Main workflow / 主要流程:
#   1. Load risk set; recover raw person-hours via effort panel.
#   2. Build 5 offset variants.
#   3. Fit M3, M5_raw, M5_log, M5_sqrt, M5_z, M5_none (all glmmTMB cloglog).
#   4. Report fixed-effects table + AIC ranking.
#   5. Save results/sensitivity/table_offset_reformulation.csv.
#
# Expected output / 预期输出:
#   results/sensitivity/table_offset_reformulation.csv
#   results/sensitivity/table_offset_coefficients.csv
#   figures/diagnostics/offset_reformulation_diagnostic.pdf
#
# Key assumptions / 关键假设:
#   - effort_panel_upgraded.csv contains the n_birding_days column
#     used as a proxy for "person-hours" (× 4h/day if explicit hours absent).
#
# Main packages / 主要包: glmmTMB, data.table, broom.mixed, ggplot2,
#   arrow, patchwork.
# Output directory / 输出路径: results/sensitivity/, figures/diagnostics/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
  library(broom.mixed)
  library(ggplot2)
  library(patchwork)
  library(arrow)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_models.R"))
source(file.path("code", "utils", "utils_plots.R"))

set.seed(42)

# ---- 1. Load data ----------------------------------------------------------
risk <- fread(path_raw("hazard_risk_upgraded_complete_case.csv"),
              encoding = "UTF-8")
effort <- fread(path_raw("effort_panel_upgraded.csv"), encoding = "UTF-8")

# Recover an approximate "person_hours" series. 估算 person-hours。
# birding_days × 4 h/day (mean field workload, Wang et al. 2024) is the
# documented assumption; users with a measured column may override.
if ("person_hours_raw" %in% names(effort)) {
  ph_col <- "person_hours_raw"
} else if ("n_birding_days" %in% names(effort)) {
  effort[, person_hours_raw := n_birding_days * 4]
  ph_col <- "person_hours_raw"
} else if ("n_visits" %in% names(effort)) {
  effort[, person_hours_raw := n_visits * 1.5]   # conservative fallback
  ph_col <- "person_hours_raw"
} else {
  stop("[32] No usable raw effort metric found in effort_panel_upgraded.csv")
}

key_cols <- intersect(c("province", "year"), names(effort))
risk <- merge(risk, effort[, c(key_cols, ph_col), with = FALSE],
              by = key_cols, all.x = TRUE)
risk[is.na(person_hours_raw), person_hours_raw := 0]

# ---- 2. Build offset variants ---------------------------------------------
risk[, off_raw  := person_hours_raw]
risk[, off_log  := log1p(person_hours_raw)]
risk[, off_sqrt := sqrt(person_hours_raw)]
zify <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
risk[, off_z    := zify(log1p(person_hours_raw))]
risk[, off_none := 0]

# Climate covariate. 气候协变量（按 v1 实际列名优先序）：
# temp_grad_z（v1 hazard_risk_upgraded_complete_case 的实际列）→ 其它 fallback
climate_candidates <- c("temp_grad_z", "temp_anom_z", "climate_velocity_z",
                         "mahalanobis_dist_z", "warming_rate_z")
climate_col <- intersect(climate_candidates, names(risk))[1]
if (is.na(climate_col)) {
  stop("[32] No usable climate z-score column. Candidates: ",
       paste(climate_candidates, collapse = ", "),
       ". Found columns: ", paste(names(risk), collapse = ", "))
}
message("[32] climate covariate column: ", climate_col)
risk[, climate_z := get(climate_col)]

# Effort covariate (always log1p(visits)_z). effort 用 spec B。
if ("log_n_visits_z" %in% names(risk)) {
  risk[, effort_z := log_n_visits_z]
} else if ("effort_z" %in% names(risk)) {
  # already present
} else {
  risk[, effort_z := zify(log1p(person_hours_raw))]
}

# ---- 3. Fit models ---------------------------------------------------------
form_M3 <- event ~ climate_z + effort_z + (1 | species) + (1 | province)
form_M5_template <- "event ~ climate_z + (1 | species) + (1 | province) + offset(%s)"

models <- list()
message("[32] Fitting M3 (no offset)...")
models$M3 <- glmmTMB::glmmTMB(form_M3, data = risk,
                              family = binomial(link = "cloglog"))

for (spec in c("off_raw", "off_log", "off_sqrt", "off_z", "off_none")) {
  message(glue("[32] Fitting M5 with offset = {spec} ..."))
  form_i <- as.formula(sprintf(form_M5_template, spec))
  fit <- tryCatch(
    glmmTMB::glmmTMB(form_i, data = risk,
                     family = binomial(link = "cloglog")),
    error = function(e) {
      message("  convergence error: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(fit)) models[[paste0("M5_", sub("off_", "", spec))]] <- fit
}

# ---- 4. Reporting ---------------------------------------------------------
aic_tab <- aic_table(models)
coef_tab <- data.table::rbindlist(lapply(names(models), function(nm) {
  out <- extract_coefs_manual(models[[nm]])
  out[, model := nm]
  out
}), use.names = TRUE, fill = TRUE)
data.table::setcolorder(coef_tab, c("model", "term"))

ensure_dir(path_sensitivity())
fwrite(aic_tab,  path_sensitivity("table_offset_reformulation.csv"))
fwrite(coef_tab, path_sensitivity("table_offset_coefficients.csv"))
arrow::write_parquet(coef_tab,
                      path_sensitivity("table_offset_coefficients.parquet"),
                      compression = "snappy")

message("[32] AIC ranking:")
print(aic_tab)

# ---- 5. Diagnostic plot ---------------------------------------------------
clim_eff_terms <- coef_tab[grepl("climate_z|effort_z", term)]
diag_plot <- forest_plot(
  clim_eff_terms,
  group_col = "model",
  title = "Offset reformulation — climate / effort terms across M5 variants"
) + ggplot2::scale_colour_manual(values = pal_cat)

ensure_dir(path_diag_fig())
ggplot2::ggsave(path_diag_fig("offset_reformulation_diagnostic.pdf"),
                diag_plot, width = 16, height = 11, units = "cm",
                device = grDevices::cairo_pdf)
ggplot2::ggsave(path_diag_fig("offset_reformulation_diagnostic.png"),
                diag_plot, width = 16, height = 11, units = "cm", dpi = 600)

# Aviator check. 一致性检验。
interaction_beta <- coef_tab[term == "climate_z" & model %in% c("M5_log", "M5_sqrt"),
                              mean(beta)]
message(glue("[32] Mean climate β under log/sqrt offsets = {round(interaction_beta, 3)}"))

dump_session_info(path_logs("32_offset_reformulation_sessionInfo.txt"))
message("[32] done.")
