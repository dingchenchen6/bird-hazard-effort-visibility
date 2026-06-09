# ============================================================
# Scientific question / 科学问题:
#   Single end-to-end runner. Calls `targets::tar_make()` so that
#   re-running the project respects the DAG (only re-executes nodes
#   whose dependencies changed).
#   v2 入口：调 targets::tar_make()，按 DAG 增量执行。
#
# Usage / 用法:
#   Rscript code/00_run_all.R
#   Rscript code/00_run_all.R --names fig5,selfcheck   # selective rerun
# ============================================================

suppressPackageStartupMessages({
  library(targets)
  library(glue)
})

args <- commandArgs(trailingOnly = TRUE)
names_subset <- NULL
if (length(args) > 0L) {
  if (any(grepl("^--names", args))) {
    arg_idx <- which(grepl("^--names", args))
    if (grepl("=", args[arg_idx])) {
      names_subset <- strsplit(sub("--names=", "", args[arg_idx]), ",")[[1]]
    } else if (length(args) > arg_idx) {
      names_subset <- strsplit(args[arg_idx + 1L], ",")[[1]]
    }
  }
}

t0 <- Sys.time()
message(glue("[00] tar_make start {t0}"))

if (is.null(names_subset)) {
  targets::tar_make()
} else {
  message(glue("[00] selective: {paste(names_subset, collapse=', ')}"))
  targets::tar_make(names = names_subset)
}

t1 <- Sys.time()
message(glue("[00] tar_make end {t1} (elapsed {round(difftime(t1, t0, units='mins'),1)} min)"))

# Persist sessionInfo on every full run. 持久化 sessionInfo。
if (exists("dump_session_info")) dump_session_info()
