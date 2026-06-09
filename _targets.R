# ============================================================
# v2 targets pipeline / DAG.
# Run: targets::tar_make()
# Inspect: targets::tar_visnetwork()
# Selective re-run: targets::tar_make(names = c("fig5", "selfcheck"))
# ============================================================

suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
})

# Source all utils so targets sees the function definitions.
tar_option_set(
  packages = c("data.table", "arrow", "fst", "yaml", "glmmTMB", "DHARMa",
               "performance", "MuMIn", "broom.mixed", "blockCV", "xgboost",
               "pROC", "PRROC", "ape", "sf", "terra", "exactextractr",
               "ggplot2", "ggtext", "patchwork", "scales", "viridisLite",
               "future.apply", "glue", "fs", "ggridges"),
  format = "rds",
  storage = "worker",
  retrieval = "worker",
  error = "continue"
)

tar_source(files = c(
  "code/utils/utils_data.R",
  "code/utils/utils_models.R",
  "code/utils/utils_spatial.R",
  "code/utils/utils_plots.R"
))

# --- Helper to convert a script into a tar_target ---------------------------
# Each script returns its primary on-disk artefact path so we can track it.
run_script <- function(script_path) {
  source(script_path, local = TRUE, chdir = FALSE)
  invisible(script_path)
}

list(
  # ---- Phase 0: derived data ------------------------------------------------
  tar_target(grid_native_climate,
             run_script("code/28_grid_native_climate.R"),
             format = "file"),
  tar_target(grid_native_effort,
             run_script("code/28b_grid_native_effort.R"),
             format = "file"),
  tar_target(grid_event_redef,
             {
               deps <- list(grid_native_climate, grid_native_effort)
               run_script("code/33_grid_event_definition_fix.R")
             },
             format = "file"),
  tar_target(prefecture_county_hazard,
             run_script("code/38_prefecture_county_hazard.R"),
             format = "file"),
  tar_target(unified_multi_metric_scale,
             {
               deps <- list(grid_event_redef, prefecture_county_hazard)
               run_script("code/39_unified_multi_metric_multi_scale.R")
             },
             format = "file"),

  # ---- Phase 1: critical fixes ---------------------------------------------
  tar_target(offset_reformulation,
             run_script("code/32_offset_reformulation.R"),
             format = "file"),
  tar_target(forecast_skill,
             run_script("code/30_forecast_skill_decay.R"),
             format = "file"),
  tar_target(spatial_block_cv,
             run_script("code/26_spatial_block_cv.R"),
             format = "file"),
  tar_target(morans_diag,
             run_script("code/27_morans_i_diagnostics.R"),
             format = "file"),

  # ---- Phase 2: ensemble & MAUP --------------------------------------------
  tar_target(cmip6_ensemble,
             run_script("code/29_cmip6_ensemble_prediction.R"),
             format = "file"),
  tar_target(maup_sensitivity,
             run_script("code/31_maup_sensitivity.R"),
             format = "file"),

  # ---- Phase 3: figures ----------------------------------------------------
  tar_target(figures_main,
             {
               deps <- list(offset_reformulation, forecast_skill,
                            spatial_block_cv, morans_diag,
                            cmip6_ensemble, maup_sensitivity,
                            prefecture_county_hazard,
                            unified_multi_metric_scale)
               run_script("code/34_publication_figures_main.R")
             },
             format = "file"),
  tar_target(figures_supp,
             {
               deps <- list(offset_reformulation, forecast_skill,
                            spatial_block_cv, morans_diag,
                            cmip6_ensemble, maup_sensitivity)
               run_script("code/35_publication_figures_supplementary.R")
             },
             format = "file"),

  # ---- Phase 4: data dictionary + self-check -------------------------------
  tar_target(data_dictionary,
             {
               deps <- list(offset_reformulation, spatial_block_cv,
                             morans_diag, forecast_skill,
                             maup_sensitivity, cmip6_ensemble)
               run_script("code/36_data_dictionary_export.R")
             },
             format = "file"),
  tar_target(selfcheck,
             {
               deps <- list(data_dictionary, figures_main, figures_supp)
               run_script("code/37_reproducibility_selfcheck.R")
             },
             format = "file"),

  # ---- Phase 5: manuscript -------------------------------------------------
  tar_render(manuscript_v2,
             "manuscript/manuscript_v2.Rmd",
             output_file = "manuscript_v2.md")
)
