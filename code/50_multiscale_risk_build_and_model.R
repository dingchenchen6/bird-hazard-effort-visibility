# ============================================================
# 50_multiscale_risk_build_and_model.R (Fixed Version)
# ============================================================
# 修复内容：
# 1. 省级风险集：直接使用事件数据的province列（无需空间连接）
# 2. 市县风险集：暂停（现有shapefile不完整，仅13个市/7个县）
# 3. 100km网格：使用已有的 risk_set_grid_100km_v2.csv
# 4. 修复所有路径和函数调用问题
#
# 科学问题：在省、市、县、100km网格水平构建新纪录风险集并建模
# ============================================================

# ---- Parse command line arguments ----
args <- commandArgs(trailingOnly = TRUE)
arg_list <- list()
i <- 1
while (i <= length(args)) {
  if (args[i] == "--base-dir") {
    arg_list$base_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output-dir") {
    arg_list$output_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--year-min") {
    arg_list$year_min <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--year-max") {
    arg_list$year_max <- as.integer(args[i + 1]); i <- i + 2
  } else { i <- i + 1 }
}

BASE_DIR <- arg_list$base_dir %||%
  "/Users/dingchenchen/Documents/New records/bird-new-distribution-records/tasks/bird_hazard_model_effort_upgrade_v2"
OUTPUT_DIR <- arg_list$output_dir %||% "outputs_multiscale"
YEAR_MIN <- arg_list$year_min %||% 2002L
YEAR_MAX <- arg_list$year_max %||% 2024L

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Setup output directories ----
OUTPUT_PATH <- file.path(BASE_DIR, OUTPUT_DIR)
dir.create(OUTPUT_PATH, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("data/derived", "results/tables", "results/diagnostics",
                  "figures/main", "figures/diagnostics", "logs")) {
  dir.create(file.path(OUTPUT_PATH, subdir), recursive = TRUE, showWarnings = FALSE)
}

# ---- Logging function ----
log_msg <- function(...) {
  msg <- paste0("[50 ", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
  cat(msg)
  log_file <- file.path(OUTPUT_PATH, "logs", "50_multiscale.log")
  cat(msg, file = log_file, append = TRUE)
}

log_msg("=== 50_multiscale_risk_build_and_model.R (Fixed) ===")
log_msg("Base dir: ", BASE_DIR)
log_msg("Output dir: ", OUTPUT_PATH)
log_msg("Year window: ", YEAR_MIN, "-", YEAR_MAX)

# ---- Load packages ----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(sf)
  library(glmmTMB)
  library(ggplot2)
  library(patchwork)
  library(glue)
})

sf::sf_use_s2(FALSE)
set.seed(42)

# ---- Helpers ----
zify <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

path_raw <- function(...) file.path(BASE_DIR, "data", "raw", ...)
path_derived <- function(...) file.path(BASE_DIR, "data", "derived", ...)
path_output <- function(...) file.path(OUTPUT_PATH, ...)

# ============================================================
# STEP 1: Load events and SDM candidates
# ============================================================
log_msg("=== Step 1: Loading events and SDM candidates ===")

events_file <- path_raw("events_100km_grid_assigned.csv")
if (!file.exists(events_file)) {
  stop("[50] Cannot find events file: ", events_file)
}
events <- fread(events_file, encoding = "UTF-8")
setnames(events, tolower(names(events)))

# Handle year column
if (!"year" %in% names(events) && "pub_year" %in% names(events)) {
  setnames(events, "pub_year", "year")
}
events[, year := as.integer(year)]

# Filter to study window
events <- events[year >= YEAR_MIN & year <= YEAR_MAX &
                 !is.na(longitude) & !is.na(latitude)]
log_msg("Events loaded: ", nrow(events), " rows, ", uniqueN(events$species), " species")
log_msg("Provinces in events: ", length(unique(events$province)))

# SDM candidates
sdm_file <- path_raw("grid_100km_risk_set.csv")
if (file.exists(sdm_file)) {
  sdm_province <- unique(fread(sdm_file, select = c("species", "province"),
                              encoding = "UTF-8"))
  log_msg("SDM candidates: ", nrow(sdm_province), " (species x province)")
} else {
  sdm_province <- unique(events[, .(species, province)])
  log_msg("SDM candidates (from events): ", nrow(sdm_province), " pairs")
}

# ============================================================
# STEP 2: Build PROVINCE-level risk set (direct from events)
# ============================================================
log_msg("=== Step 2: Building PROVINCE-level risk set ===")

# Events already have province column - no spatial join needed!
# 事件数据已有province列，无需空间连接！

# First arrival per (species, province)
first_arrival <- events[, .(arrival_year = min(year, na.rm = TRUE)),
                         by = .(species, province)]
log_msg("First arrival pairs: ", nrow(first_arrival))

# Build candidate (species x province) from SDM
candidates <- unique(sdm_province[, .(species, province)])

# Expand to years
year_seq <- YEAR_MIN:YEAR_MAX
# Use CJ (cross join) instead of merge with by=NULL
risk_prov <- CJ(
  species = candidates$species,
  province = candidates$province,
  year = year_seq
)
log_msg("Risk set before filtering: ", nrow(risk_prov), " rows")

# Merge first arrival
risk_prov <- merge(risk_prov, first_arrival,
                    by = c("species", "province"), all.x = TRUE)

# Mark events
risk_prov[, event := as.integer(year == arrival_year)]
risk_prov[is.na(event), event := 0L]

# Censor after arrival
risk_prov <- risk_prov[is.na(arrival_year) | year <= arrival_year]
risk_prov[, arrival_year := NULL]

log_msg("Final PROVINCE risk set: ", nrow(risk_prov), " rows, ",
        sum(risk_prov$event), " events")

# Add province-level climate
clim_file <- path_raw("climate_metrics_province_year.csv")
if (file.exists(clim_file)) {
  clim <- fread(clim_file, encoding = "UTF-8")
  setnames(clim, tolower(names(clim)))
  risk_prov <- merge(risk_prov, clim, by = c("province", "year"), all.x = TRUE)
  log_msg("Climate attached to province risk set")
}

# Add effort (allocate province effort to species proportionally)
effort_file <- path_raw("effort_panel_upgraded.csv")
if (file.exists(effort_file)) {
  effort <- fread(effort_file, encoding = "UTF-8")
  setnames(effort, tolower(names(effort)))
  
  # Province-level effort: sum across species/units
  effort_prov <- effort[, .(effort = sum(effort, na.rm = TRUE)),
                        by = .(province, year)]
  
  # Z-score effort BY PROVINCE (not species - effort_prov has no species column)
  effort_prov[, effort_z := zify(effort), by = province]
  
  risk_prov <- merge(risk_prov, effort_prov[, .(province, year, effort_z)],
                     by = c("province", "year"), all.x = TRUE)
  log_msg("Effort attached to province risk set")
}

# Add unit_id for modeling
risk_prov[, unit_id := paste0("prov_", province)]
risk_prov[, unit_id := factor(unit_id)]

# Save province risk set
prov_out <- path_output("data/derived", "risk_set_province.parquet")
if (requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(risk_prov, prov_out)
  log_msg("Saved: ", basename(prov_out))
} else {
  fwrite(risk_prov, gsub("\\.parquet$", ".csv", prov_out))
  log_msg("Saved as CSV: risk_set_province.csv")
}

# ============================================================
# STEP 3: Try to build PREFECTURE/COUNTY risk sets
# ============================================================
log_msg("=== Step 3: Attempting PREFECTURE/COUNTY risk sets ===")
log_msg("NOTE: Skipped - shapefiles incomplete (13 prefectures, 7 counties only)")
log_msg("Need complete prefecture/county shapefiles for full analysis")

pref_risk <- NULL
county_risk <- NULL

# ============================================================
# STEP 4: Build 100km GRID risk set
# ============================================================
log_msg("=== Step 4: Building 100km GRID risk set ===")

grid_file <- path_derived("risk_set_grid_100km_v2.csv")
if (file.exists(grid_file)) {
  grid_risk <- fread(grid_file, encoding = "UTF-8")
  setnames(grid_risk, tolower(names(grid_risk)))
  log_msg("Grid risk set loaded: ", nrow(grid_risk), " rows")
  
  # Filter to study window
  if ("year" %in% names(grid_risk)) {
    grid_risk <- grid_risk[year >= YEAR_MIN & year <= YEAR_MAX]
    log_msg("After year filter: ", nrow(grid_risk), " rows")
  }
  
  # Save
  grid_out <- path_output("data/derived", "risk_set_grid_100km.parquet")
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(grid_risk, grid_out)
    log_msg("Saved: ", basename(grid_out))
  }
  
} else {
  log_msg("WARNING: grid_100km_v2 not found, skipping grid analysis")
  grid_risk <- NULL
}

# ============================================================
# STEP 5: Fit hazard models (PROVINCE level as example)
# ============================================================
log_msg("=== Step 5: Fitting hazard models (province level) ===")

if (exists("risk_prov") && nrow(risk_prov) > 0) {
  
  # Prepare model data
  mdata <- copy(risk_prov)
  
  # Z-score predictors (only numeric columns that are not already z-scored)
  clim_cols <- names(mdata)[grepl("^(temp|precip|climate|warming|mahalanobis)", names(mdata)) & 
                              !grepl("_z$", names(mdata))]
  
  for (col in clim_cols) {
    if (is.numeric(mdata[[col]])) {
      mdata[, paste0(col, "_z") := zify(get(col)), by = species]
    }
  }
  
  # Effort is already z-scored from earlier step
  if ("effort_z" %in% names(mdata)) {
    log_msg("Effort already z-scored")
  }
  
  # Fit M0: null model (intercept only)
  log_msg("Fitting M0 (null model)...")
  m0 <- try(glmmTMB(event ~ 1 + (1|species) + (1|unit_id),
                  data = mdata, family = "binomial"(link = "cloglog")),
            silent = TRUE)
  
  if (!inherits(m0, "try-error")) {
    log_msg("M0 fitted successfully")
    
    # Extract coefficients
    m0_coef <- as.data.table(summary(m0)$coefficients$cond, keep.rownames = TRUE)
    m0_coef$model <- "M0"
    m0_coef$scale <- "province"
    
    # Save
    coef_file <- path_output("results/tables", "table_multiscale_coefficients.csv")
    fwrite(m0_coef, coef_file)
    log_msg("Coefficients saved: ", basename(coef_file))
    
    # AIC comparison
    aic_df <- data.table(
      model = "M0",
      scale = "province",
      AIC = AIC(m0),
      n_obs = nrow(mdata),
      n_events = sum(mdata$event)
    )
    aic_file <- path_output("results/tables", "table_multiscale_aic.csv")
    fwrite(aic_df, aic_file)
    log_msg("AIC saved: ", basename(aic_file))
    
  } else {
    log_msg("WARNING: M0 fitting failed")
  }
  
} else {
  log_msg("No province risk set available for modeling")
}

# ============================================================
# STEP 6: Session info
# ============================================================
log_msg("=== Step 6: Saving session info ===")

sink(path_output("logs", "50_multiscale_sessionInfo.txt"))
print(sessionInfo())
sink()

log_msg("=== DONE ===")
log_msg("Output directory: ", OUTPUT_PATH)
log_msg("Check logs at: ", file.path(OUTPUT_PATH, "logs", "50_multiscale.log"))
