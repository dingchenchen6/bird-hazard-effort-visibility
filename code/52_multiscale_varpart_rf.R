# ============================================================
# Scientific question / 科学问题:
#   How much variance in new-bird-record hazard is uniquely
#   explained by climate, effort, and their interaction — and does
#   this decomposition change across spatial scales?
#   在省、市、县、100km 网格四个尺度上，气候、调查努力和交互
#   各自解释多少方差？随机森林变量重要性排序是否跨尺度一致？
#
# Objective / 分析目标:
#   - Variance partitioning via nested glmmTMB model sequence
#     (marginal R² increments) at 4 scales
#   - Random forest (ranger) permutation importance at 4 scales
#   - Cross-scale comparison figures (stacked bar + dot plot)
#
# Input data / 输入数据:
#   - outputs_multiscale/data/derived/risk_set_{province,prefecture,
#     county,grid_100km}.csv  (produced by 51_multiscale_full.R)
#
# Expected output / 预期输出:
#   outputs_multiscale/results/tables/
#     table_multiscale_variance_decomposition.csv
#     table_multiscale_rf_variable_importance.csv
#   outputs_multiscale/figures/main/
#     fig_multiscale_variance_decomposition.{pdf,png}
#     fig_multiscale_rf_variable_importance.{pdf,png}
#   outputs_multiscale/logs/
#     52_multiscale_varpart_rf.log
#
# Key assumptions / 关键假设:
#   - Risk sets from script 51 must exist before running this.
#   - Climate variables are province-inherited at prefecture/county;
#     grid has its own climate_velocity_z.
#   - RF ignores random effects (standard approach); class weights
#     correct for extreme event imbalance.
#
# Main packages / 主要包: data.table, glmmTMB, ranger, ggplot2.
# Output directory / 输出路径: outputs_multiscale/
# ============================================================

# ---- 0. CLI args ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
arg <- list()
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--base-dir")        { arg$base_dir   <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--output-dir") { arg$output_dir <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--skip-county"){ arg$skip_county <- TRUE; i <- i + 1L }
  else { i <- i + 1L }
}
`%||%` <- function(a, b) if (is.null(a)) b else a

BASE_DIR   <- arg$base_dir   %||% "/Users/dingchenchen/Documents/New records/bird-new-distribution-records/tasks/bird_hazard_model_effort_upgrade_v2"
OUTPUT_DIR <- arg$output_dir %||% "outputs_multiscale"
SKIP_COUNTY <- isTRUE(arg$skip_county)

OUT <- file.path(BASE_DIR, OUTPUT_DIR)
ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
for (sd in c("results/tables", "results/diagnostics", "figures/main", "logs"))
  ens(file.path(OUT, sd))

log_file <- file.path(OUT, "logs", "52_multiscale_varpart_rf.log")
log <- function(...) {
  msg <- sprintf("[52 %s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = ""))
  cat(msg); cat(msg, file = log_file, append = TRUE)
}
log("=== 52_multiscale_varpart_rf.R START ===")

# ---- 1. Packages ----------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
  library(ranger)
  library(ggplot2)
})
options(warn = 1)
set.seed(42)

# ---- 2. Helpers -----------------------------------------------------------
zify <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

compute_r2 <- function(fit) {
  tryCatch({
    vc <- VarCorr(fit)
    X <- model.matrix(fit)
    beta <- fixef(fit)$cond
    fix_var <- var(as.numeric(X %*% beta))
    re_var <- sum(sapply(vc$cond, function(x) sum(diag(x))))
    dist_var <- pi^2 / 6   # cloglog distribution-specific variance
    list(marginal = fix_var / (fix_var + re_var + dist_var),
         conditional = (fix_var + re_var) / (fix_var + re_var + dist_var),
         fix_var = fix_var, re_var = re_var)
  }, error = function(e) {
    list(marginal = NA, conditional = NA, fix_var = NA, re_var = NA)
  })
}

# ---- 3. Load risk sets from script 51 ------------------------------------
log("=== Step 3: Loading risk sets ===")

scales <- c("province", "prefecture", "grid_100km")
if (!SKIP_COUNTY) scales <- c(scales, "county")
# 确保 county 在正确位置
scales <- intersect(c("province", "prefecture", "county", "grid_100km"), scales)

risk_sets <- list()
for (sc in scales) {
  f <- file.path(OUT, "data", "derived", paste0("risk_set_", sc, ".csv"))
  if (file.exists(f)) {
    dt <- fread(f, encoding = "UTF-8")
    setnames(dt, tolower(names(dt)))
    log(sc, ": ", nrow(dt), " rows, ", sum(dt$event, na.rm = TRUE), " events")
    risk_sets[[sc]] <- dt
  } else {
    log("WARNING: ", f, " not found, skipping ", sc)
  }
}
if (length(risk_sets) == 0L) stop("No risk sets found. Run 51_multiscale_full.R first.")

# ---- 4. Variance partitioning at each scale -------------------------------
log("=== Step 4: Variance partitioning ===")

pal_contrib <- c(
  "Climate" = "#d94801", "Effort" = "#2171b5",
  "Interaction" = "#6a51a3", "Year" = "#666666",
  "Random effects" = "#aaaaaa")

vd_all <- list()

for (sc in names(risk_sets)) {
  log("  Variance partitioning: ", sc)
  dt <- copy(risk_sets[[sc]])

  # Determine unit column name
  unit_col <- if ("unit_id" %in% names(dt)) "unit_id"
              else if ("grid_id" %in% names(dt)) "grid_id"
              else "province"
  dt[, unit_col_val := get(unit_col)]
  dt[, unit_col_val := factor(unit_col_val)]
  dt[, species := factor(species)]

  # Ensure required columns exist
  if (!"climate_z" %in% names(dt)) dt[, climate_z := zify(climate_velocity_z)]
  if (!"effort_z"  %in% names(dt)) {
    eff_col <- intersect(c("log_n_events_z", "log_effort_visits_z",
                           "effort_pc1_z", "log_n_visits_z"), names(dt))[1]
    if (!is.na(eff_col)) dt[, effort_z := get(eff_col)]
    else { log("  no effort column in ", sc, ", skip"); next }
  }
  if (!"year_c" %in% names(dt)) dt[, year_c := year - 2013]

  re_form <- if (unit_col == "province") "(1|species) + (1|province)"
             else paste0("(1|species) + (1|unit_col_val)")

  nested_forms <- list(
    M0_null    = as.formula(paste("event ~ year_c +", re_form)),
    M1_climate = as.formula(paste("event ~ year_c + climate_z +", re_form)),
    M2_effort  = as.formula(paste("event ~ year_c + effort_z +", re_form)),
    M3_additive = as.formula(paste("event ~ year_c + climate_z + effort_z +", re_form)),
    M4_full    = as.formula(paste("event ~ year_c + climate_z * effort_z +", re_form))
  )

  r2_res <- list()
  for (mname in names(nested_forms)) {
    fit <- tryCatch(
      glmmTMB(nested_forms[[mname]], data = dt,
              family = binomial(link = "cloglog")),
      error = function(e) { log("    ", mname, " failed: ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    r2 <- compute_r2(fit)
    aic_val <- tryCatch(AIC(fit), error = function(e) NA_real_)
    log(sprintf("    %s: mR2=%.4f cR2=%.4f AIC=%.1f",
                mname, r2$marginal, r2$conditional, aic_val))
    r2_res[[mname]] <- data.table(scale = sc, model = mname,
                                    marginal_r2 = r2$marginal,
                                    conditional_r2 = r2$conditional,
                                    aic = aic_val)
  }

  if (length(r2_res) < 5L) {
    log("  incomplete model sequence for ", sc, ", skip decomposition")
    next
  }
  dt_r2 <- rbindlist(r2_res)

  r0 <- dt_r2[model == "M0_null",    marginal_r2]
  r1 <- dt_r2[model == "M1_climate", marginal_r2]
  r2 <- dt_r2[model == "M2_effort",  marginal_r2]
  r3 <- dt_r2[model == "M3_additive",marginal_r2]
  r4 <- dt_r2[model == "M4_full",    marginal_r2]

  vd_all[[sc]] <- data.table(
    scale = sc,
    Year        = r0,
    Climate     = r1 - r0,
    Effort      = r2 - r0,
    Joint       = r3 - r1 - (r2 - r0),
    Interaction = r4 - r3,
    Total_mR2   = r4
  )
  log("  ", sc, " decomposition: Year=", round(r0, 4),
      " Climate=", round(r1 - r0, 4),
      " Effort=", round(r2 - r0, 4),
      " Interaction=", round(r4 - r3, 4))
}

vd_dt <- rbindlist(vd_all, fill = TRUE)
fwrite(vd_dt, file.path(OUT, "results", "tables",
                        "table_multiscale_variance_decomposition.csv"))
log("Variance decomposition table saved: ", nrow(vd_dt), " rows")

# ---- 5. Variance decomposition figure ------------------------------------
if (nrow(vd_dt) > 0L) {
  # 转换为 long format
  vd_long <- melt(vd_dt, id.vars = "scale",
                  measure.vars = c("Year", "Climate", "Effort", "Joint", "Interaction"),
                  variable.name = "component", value.name = "delta_r2")
  vd_long[, pct := 100 * delta_r2 / vd_dt[match(vd_long$scale, scale), Total_mR2]]

  vd_long[, scale := factor(scale,
    levels = c("province", "prefecture", "county", "grid_100km"),
    labels = c("Province", "Prefecture", "County", "100km grid"))]

  p_vd <- ggplot(vd_long, aes(x = scale, y = delta_r2, fill = component)) +
    geom_col(alpha = 0.85, position = "stack") +
    geom_text(data = vd_long[delta_r2 > 0.002],
              aes(label = sprintf("%.1f%%", pct)),
              position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = pal_contrib, name = "Component") +
    labs(x = NULL, y = expression(Delta~marginal~R^2),
         title = "Variance decomposition across spatial scales",
         subtitle = "Climate, effort, joint overlap, and interaction contributions to hazard") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom")

  ggsave(file.path(OUT, "figures", "main", "fig_multiscale_variance_decomposition.pdf"),
         p_vd, width = 16, height = 9, units = "cm", device = grDevices::cairo_pdf)
  ggsave(file.path(OUT, "figures", "main", "fig_multiscale_variance_decomposition.png"),
         p_vd, width = 16, height = 9, units = "cm", dpi = 600)
  log("Variance decomposition figure saved")
}

# ---- 6. Random Forest variable importance at each scale -------------------
log("=== Step 6: Random Forest variable importance ===")

pal_varimp <- c("Climate" = "#d94801", "Effort" = "#2171b5",
                "Year" = "#666666", "Other" = "#aaaaaa")

rf_all <- list()

for (sc in names(risk_sets)) {
  log("  RF importance: ", sc)
  dt <- copy(risk_sets[[sc]])

  # Collect available predictor columns
  climate_cols <- intersect(c("climate_velocity_z", "climate_exposure_z",
                              "warming_rate_z", "mahalanobis_dist_z",
                              "temp_grad_z", "prec_grad_z",
                              "climate_z", "prov_climate_velocity_z",
                              "prov_warming_rate_z", "prov_temp_anom_z",
                              "prov_climate_exposure_z"), names(dt))
  effort_cols  <- intersect(c("effort_z", "log_n_events_z",
                              "log_effort_visits_z", "effort_pc1_z",
                              "log_n_visits_z", "log_effort_days_z",
                              "log_n_observers_z"), names(dt))
  year_col <- intersect(c("year_c", "year"), names(dt))[1]

  keep_cols <- c("event", year_col, climate_cols, effort_cols)
  keep_cols <- keep_cols[!is.na(keep_cols)]
  rf_data <- dt[, ..keep_cols]
  rf_data <- rf_data[complete.cases(rf_data)]

  if (nrow(rf_data) < 500L) {
    log("  too few complete rows (", nrow(rf_data), "), skip"); next
  }

  # year → year_c if not already
  if ("year" %in% names(rf_data) && !"year_c" %in% names(rf_data)) {
    rf_data[, year_c := year - 2013]
    rf_data[, year := NULL]
  }

  rf_data[, event := factor(event, levels = c("0", "1"))]
  n_events <- sum(rf_data$event == "1")
  n_nonev  <- sum(rf_data$event == "0")

  log("  RF data: ", nrow(rf_data), " rows, ", n_events, " events")

  rf_fit <- tryCatch(
    ranger(event ~ ., data = rf_data,
           num.trees = 500L,
           mtry = floor(sqrt(ncol(rf_data) - 1L)),
           min.node.size = 10L,
           class.weights = c(1, n_nonev / pmax(n_events, 1L)),
           importance = "permutation",
           seed = 42L,
           verbose = FALSE),
    error = function(e) { log("  RF failed: ", conditionMessage(e)); NULL })

  if (is.null(rf_fit)) next

  vi <- data.table(
    scale = sc,
    variable = names(rf_fit$variable.importance),
    importance = as.numeric(rf_fit$variable.importance))
  vi <- vi[order(-importance)]
  vi[, category := fcase(
    variable %in% c("climate_z", "climate_velocity_z", "climate_exposure_z",
                    "warming_rate_z", "mahalanobis_dist_z", "temp_grad_z",
                    "prec_grad_z", "prov_climate_velocity_z",
                    "prov_warming_rate_z", "prov_temp_anom_z",
                    "prov_climate_exposure_z"),
    "Climate",
    variable %like% "effort" | variable %like% "log_n_" | variable %like% "log_effort",
    "Effort",
    variable == "year_c", "Year",
    default = "Other")]

  rf_all[[sc]] <- vi
  log("  top 3: ", paste(vi$variable[1:3], collapse = ", "))
}

rf_dt <- rbindlist(rf_all, fill = TRUE)
fwrite(rf_dt, file.path(OUT, "results", "tables",
                        "table_multiscale_rf_variable_importance.csv"))
log("RF importance table saved: ", nrow(rf_dt), " rows")

# ---- 7. RF importance figure ----------------------------------------------
if (nrow(rf_dt) > 0L) {
  rf_dt[, scale := factor(scale,
    levels = c("province", "prefecture", "county", "grid_100km"),
    labels = c("Province", "Prefecture", "County", "100km grid"))]

  # Standardize importance within each scale for comparability
  rf_dt[, imp_std := importance / max(importance, na.rm = TRUE), by = scale]

  p_rf <- ggplot(rf_dt, aes(x = reorder(variable, imp_std), y = imp_std,
                              fill = category)) +
    geom_col(alpha = 0.85, show.legend = TRUE) +
    coord_flip() +
    facet_wrap(~ scale, scales = "free_x", ncol = 2) +
    scale_fill_manual(values = pal_varimp, name = "Category") +
    labs(x = NULL, y = "Standardized permutation importance",
         title = "Random Forest variable importance across spatial scales",
         subtitle = "Permutation importance from ranger (500 trees), standardized per scale") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          strip.text = element_text(face = "bold"))

  ggsave(file.path(OUT, "figures", "main", "fig_multiscale_rf_variable_importance.pdf"),
         p_rf, width = 22, height = 16, units = "cm", device = grDevices::cairo_pdf)
  ggsave(file.path(OUT, "figures", "main", "fig_multiscale_rf_variable_importance.png"),
         p_rf, width = 22, height = 16, units = "cm", dpi = 600)
  log("RF importance figure saved")
}

# ---- 8. Session info ------------------------------------------------------
sink(file.path(OUT, "logs", "52_sessionInfo.txt"))
print(sessionInfo())
sink()

log("=== 52_multiscale_varpart_rf.R DONE ===")
