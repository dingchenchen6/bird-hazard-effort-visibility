# ============================================================
# Scientific question / 科学问题:
#   Provide a single, validated data-access layer so that all v2
#   scripts share identical assumptions about file paths, encodings,
#   schemas and large-file backends (parquet / fst).
#   提供统一的数据访问层，让所有 v2 脚本共享相同的路径、编码、
#   schema 与大文件后端（parquet / fst）假设。
#
# Objective / 分析目标:
#   - Centralise raw → derived file conversion (CSV → parquet).
#   - Validate columns and types against YAML schemas.
#   - Lazy-load large risk-set tables via arrow.
#
# Input data / 输入数据:
#   data/raw/*.csv (symlinks to v1)
#   data_dictionary/schema_*.yaml
#
# Main workflow / 主要流程:
#   1. Path helpers and project-root anchor.
#   2. Read / write derived parquet (with snappy compression).
#   3. Schema validation (column names, types, nullability).
#   4. Lazy arrow dataset opener for the 5.5 GB grid_50km_risk_set.
#
# Expected output / 预期输出:
#   - functions exported by side-effect (no on-disk artefacts here).
#
# Key assumptions / 关键假设:
#   - Working directory is the v2 project root.
#   - arrow >= 14, fst >= 0.9 installed (see renv.lock).
#
# Main packages / 主要包: arrow, fst, yaml, data.table, glue, fs.
# Output directory / 输出路径: NA (library code).
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(fst)
  library(yaml)
  library(glue)
  library(fs)
})

# ---- 0. Project root anchor -------------------------------------------------
# Use here::here if available, else fall back to getwd(). 项目根目录锚点。
v2_root <- function() {
  if (requireNamespace("here", quietly = TRUE)) {
    return(here::here())
  }
  getwd()
}

v2_path <- function(...) file.path(v2_root(), ...)

# Sub-directory shortcuts. 子目录快捷方式。
path_raw         <- function(...) v2_path("data", "raw", ...)
path_derived     <- function(...) v2_path("data", "derived", ...)
path_spatial     <- function(...) v2_path("data", "spatial", ...)
path_results     <- function(...) v2_path("results", ...)
path_tables      <- function(...) v2_path("results", "tables", ...)
path_diagnostics <- function(...) v2_path("results", "diagnostics", ...)
path_sensitivity <- function(...) v2_path("results", "sensitivity", ...)
path_forecasts   <- function(...) v2_path("results", "forecasts", ...)
path_figures     <- function(...) v2_path("figures", ...)
path_main_fig    <- function(...) v2_path("figures", "main", ...)
path_supp_fig    <- function(...) v2_path("figures", "supplementary", ...)
path_diag_fig    <- function(...) v2_path("figures", "diagnostics", ...)
path_manuscript  <- function(...) v2_path("manuscript", ...)
path_logs        <- function(...) v2_path("logs", ...)

# Ensure directory exists. 确保目录存在。
ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(p)
}

# ---- 1. CSV → parquet migration helper -------------------------------------
# Pulls a v1 CSV (symlinked into data/raw/) and writes it as parquet under
# data/derived/, with optional row-level filtering. 把 v1 CSV 转成 parquet。
csv_to_parquet <- function(csv_name, parquet_name = NULL,
                           where = NULL, compression = "snappy",
                           overwrite = FALSE) {
  src <- path_raw(csv_name)
  if (!file.exists(src)) stop("Missing raw CSV: ", src)
  if (is.null(parquet_name)) {
    parquet_name <- sub("\\.csv$", ".parquet", csv_name)
  }
  dst <- path_derived(parquet_name)
  ensure_dir(dirname(dst))
  if (file.exists(dst) && !overwrite) {
    message("[skip] parquet exists: ", dst)
    return(invisible(dst))
  }

  message("[load] ", src)
  dt <- data.table::fread(src, encoding = "UTF-8", showProgress = FALSE)
  if (!is.null(where)) dt <- dt[eval(where, envir = dt)]
  arrow::write_parquet(dt, dst, compression = compression)
  message("[write] ", dst, "  (", nrow(dt), " rows)")
  invisible(dst)
}

# ---- 2. Lazy load helpers ---------------------------------------------------
# Open a parquet (or fall back to CSV) as a data.table. parquet 优先。
load_derived <- function(name, as_data_table = TRUE) {
  candidates <- c(
    path_derived(name),
    path_derived(paste0(name, ".parquet")),
    path_derived(paste0(name, ".fst")),
    path_derived(paste0(name, ".csv"))
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Cannot find derived file for name: ", name)
  ext <- tolower(tools::file_ext(hit))
  dat <- switch(ext,
    parquet = arrow::read_parquet(hit),
    fst     = fst::read_fst(hit, as.data.table = TRUE),
    csv     = data.table::fread(hit, encoding = "UTF-8"),
    stop("Unsupported extension: ", ext)
  )
  if (as_data_table) data.table::setDT(dat)
  attr(dat, "v2_source") <- hit
  dat
}

# Open a large parquet as an arrow Dataset (no data in RAM). 大文件懒加载。
open_dataset_v2 <- function(name) {
  candidate <- path_derived(paste0(name, ".parquet"))
  if (!file.exists(candidate)) {
    stop("Run csv_to_parquet(\"", name, ".csv\") first; missing ", candidate)
  }
  arrow::open_dataset(candidate)
}

# ---- 3. Schema validation ---------------------------------------------------
# Validate a data.table against a YAML schema. 校验 schema。
validate_schema <- function(dt, schema_path, strict = TRUE) {
  if (!file.exists(schema_path)) stop("Missing schema: ", schema_path)
  schema <- yaml::read_yaml(schema_path)
  vars <- schema$variables
  if (is.null(vars)) stop("Schema YAML lacks $variables: ", schema_path)
  declared <- vapply(vars, function(v) v$name, character(1))
  missing  <- setdiff(declared, names(dt))
  extra    <- setdiff(names(dt), declared)
  ok <- length(missing) == 0L
  if (length(missing)) {
    msg <- glue("Schema mismatch: missing {length(missing)} column(s) — {paste(missing, collapse=', ')}")
    if (strict) stop(msg) else warning(msg)
  }
  if (length(extra)) {
    message("[schema] extra columns: ", paste(extra, collapse = ", "))
  }
  for (v in vars) {
    if (!v$name %in% names(dt)) next
    want <- v$type
    got  <- class(dt[[v$name]])[1]
    if (!is.null(want) && !is.null(got) && !identical(want, got)) {
      msg <- glue("Column {v$name}: expected {want}, got {got}")
      if (isTRUE(v$strict_type %||% FALSE) && strict) stop(msg) else message("[schema] ", msg)
    }
  }
  invisible(ok)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- 4. Convenience: write a CSV + parquet pair under results/tables --------
# Writes a final result table in both formats so collaborators can open
# either; the parquet is canonical. 结果表同时保留 csv + parquet。
write_result_table <- function(dt, name) {
  ensure_dir(path_tables())
  csv <- path_tables(paste0(name, ".csv"))
  parq <- path_tables(paste0(name, ".parquet"))
  data.table::fwrite(dt, csv)
  arrow::write_parquet(dt, parq, compression = "snappy")
  invisible(c(csv = csv, parquet = parq))
}

# ---- 5. sessionInfo dump ----------------------------------------------------
dump_session_info <- function(path = v2_path("sessionInfo.txt")) {
  fp <- file(path, "wt", encoding = "UTF-8")
  on.exit(close(fp))
  writeLines(c(
    paste0("# sessionInfo() captured ", format(Sys.time(), tz = "UTC"), " UTC"),
    capture.output(sessionInfo()),
    "",
    "# Loaded namespaces:",
    capture.output(loadedNamespaces())
  ), fp)
  invisible(path)
}

# ---- 6. Manifest of v1 raw symlinks ----------------------------------------
# Returns a tibble describing what is in data/raw/. 报告 raw/ 内容。
raw_manifest <- function() {
  files <- list.files(path_raw(), full.names = TRUE)
  data.table::data.table(
    file        = basename(files),
    bytes       = file.size(files),
    is_symlink  = file.info(files)$isdir == FALSE & !is.na(Sys.readlink(files)) & nzchar(Sys.readlink(files)),
    real_path   = vapply(files, function(f) {
                    rl <- Sys.readlink(f)
                    if (nzchar(rl)) rl else f
                  }, character(1))
  )
}
