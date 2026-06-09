source("renv/activate.R")
# v2 .Rprofile — activates renv and sets project-wide options.
# Loaded automatically by R when the working directory is this project.

# renv (graceful if lockfile not yet present)
if (file.exists("renv/activate.R")) source("renv/activate.R")

# Defaults
options(
  stringsAsFactors  = FALSE,
  scipen            = 999,
  digits            = 6,
  warn              = 1,
  encoding          = "UTF-8"
)

# Tile / device
if (interactive()) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Source utils only when a non-targets script asks for them. _targets.R
# explicitly sources utils as part of its `tar_source()` call.
.utils_path <- file.path("code", "utils")
if (dir.exists(.utils_path)) {
  attr(.utils_path, "v2_utils_ready") <- TRUE
}
rm(.utils_path)

# Project root convenience
options(v2.project_root = normalizePath(".", mustWork = FALSE))
