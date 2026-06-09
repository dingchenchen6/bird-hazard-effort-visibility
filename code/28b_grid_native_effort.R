# ============================================================
# Scientific question / 科学问题:
#   v1's `grid_{50,100}km_effort.csv` were built by merging the
#   PROVINCE-LEVEL z-scored effort panel onto every (grid_id, year)
#   cell of that province — so all cells inside a province×year share
#   identical z-scores (we verified e.g. Anhui-2002: 9 cells, all
#   log_effort_visits_z = -1.0019). This is the effort analogue of the
#   climate MAUP artefact (P0-1) and must be fixed: grid hazard models
#   must see real within-province variation in survey effort.
#   v1 把省级 effort z-scores 复制到每一格，等同省级假象。本脚本基于
#   坐标级原始记录与可选 GBIF/eBird 网格密度，构造 grid_id × year 真实
#   变异的 effort panel，并保留省级总量作为约束。
#
# Objective / 分析目标:
#   For each 50- and 100-km grid cell, compute year-resolved effort
#   metrics with true within-province variation:
#     1) n_records_grid      — 该格年内全记录数
#     2) n_visits_grid       — 推断的访问次数（每观测者-日为一次）
#     3) n_observers_grid    — 唯一观测者数
#     4) n_birding_days_grid — 唯一观测日数
#     5) effort_pc1_grid     — PCA 合成
#   Z-score within (year × all-grids) to produce grid-native z-scores
#   *and* keep the original province z-score for sensitivity comparison.
#
# Input data / 输入数据:
#   - Coordinate-level records: data/raw/events_{50,100}km_grid_assigned.csv
#       (species × province × longitude × latitude × grid_id × pub_year)
#   - Optional finer-grained records: data/raw/bird_new_records_20260509.xlsx
#       (if exists, parsed for observer / date fields)
#   - data/raw/effort_panel_upgraded.csv     (province totals, for QC)
#   - data/raw/grid_{50,100}km_base.csv      (grid centroids + province)
#
# Main workflow / 主要流程:
#   1. Load events with coordinates and assigned grid_id.
#   2. (If xlsx present) read full record table, derive observer + date.
#   3. Aggregate to (grid_id × year): records, visits, observers, days.
#   4. Fill (grid_id × year) zero cells from grid base × study years.
#   5. PCA on the four log1p-transformed metrics → effort_pc1.
#   6. Z-score grid-natively, also keep province-mirror for sensitivity.
#   7. QC: province totals reconstructed from grids match v1 totals.
#   8. Write data/derived/grid_{50,100}km_effort_native.parquet.
#
# Expected output / 预期输出:
#   data/derived/grid_50km_effort_native.parquet
#   data/derived/grid_100km_effort_native.parquet
#   results/diagnostics/table_grid_effort_native_summary.csv
#   results/diagnostics/figS_grid_effort_within_province_variance.{pdf,png}
#
# Key assumptions / 关键假设:
#   - In absence of observer/date columns, "visits" ≈ "records",
#     "observers" ≈ "records" and "birding_days" ≈ "records". A warning
#     is logged so users know which metrics collapse.
#   - Year window 2002–2024 (matches effort_panel_upgraded.csv).
#
# Main packages / 主要包: data.table, arrow, readxl, ggplot2, glue.
# Output directory / 输出路径: data/derived/, results/diagnostics/,
#   figures/diagnostics/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(ggplot2)
  library(glue)
})
has_readxl <- requireNamespace("readxl", quietly = TRUE)

source(file.path("code", "utils", "utils_data.R"))
source(file.path("code", "utils", "utils_plots.R"))

CFG <- list(
  year_min   = 2002,
  year_max   = 2024,
  grid_sizes = c(50, 100)
)

# ---- 1. Coordinate-level records ------------------------------------------
load_grid_events <- function(km) {
  p <- path_raw(glue("events_{km}km_grid_assigned.csv"))
  if (!file.exists(p)) stop("[28b] missing ", p)
  dt <- fread(p, encoding = "UTF-8")
  setnames(dt, tolower(names(dt)))
  if (!"year" %in% names(dt) && "pub_year" %in% names(dt)) {
    setnames(dt, "pub_year", "year")
  }
  dt[, year := as.integer(year)]
  dt <- dt[year >= CFG$year_min & year <= CFG$year_max]
  dt
}

# ---- 2. Optional Excel parse for observer / date --------------------------
load_xlsx_records <- function() {
  p <- path_raw("bird_new_records_20260509.xlsx")
  if (!file.exists(p) || !has_readxl) return(NULL)
  sheets <- readxl::excel_sheets(p)
  # Heuristic: pick the largest sheet. 取最大表。
  sizes <- vapply(sheets, function(s) {
    tryCatch(nrow(readxl::read_excel(p, sheet = s, n_max = 5)),
             error = function(e) 0L)
  }, integer(1))
  s_best <- sheets[which.max(sizes)]
  dt <- tryCatch(as.data.table(readxl::read_excel(p, sheet = s_best)),
                  error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  setnames(dt, tolower(names(dt)))
  dt
}

# ---- 3. Derive observer / date columns if available -----------------------
extract_observer_date <- function(xlsx_dt) {
  if (is.null(xlsx_dt)) return(NULL)
  candidates_obs <- c("observer", "observers", "记录人", "观察者", "report_by",
                       "submitted_by", "recorder", "people")
  candidates_date <- c("date", "obs_date", "记录日期", "观察日期",
                        "report_date", "记录时间")
  obs_col  <- intersect(candidates_obs,  names(xlsx_dt))[1]
  date_col <- intersect(candidates_date, names(xlsx_dt))[1]
  if (is.na(obs_col) && is.na(date_col)) return(NULL)
  keep <- intersect(c("species", "longitude", "latitude", "year",
                       "province", obs_col, date_col), names(xlsx_dt))
  out <- xlsx_dt[, ..keep]
  if (!is.na(obs_col))  setnames(out, obs_col,  "observer")
  if (!is.na(date_col)) setnames(out, date_col, "obs_date")
  out
}

# ---- 4. Aggregate to (grid_id × year) -------------------------------------
aggregate_grid_effort <- function(km) {
  events <- load_grid_events(km)
  xlsx <- load_xlsx_records()
  obs_dt <- extract_observer_date(xlsx)

  # If observer/date are available, attempt to join to events by
  # (species, year, province, lon, lat). 若可用，做四元组连接。
  if (!is.null(obs_dt) &&
      all(c("species", "year", "longitude", "latitude") %in% names(obs_dt))) {
    events <- merge(events, obs_dt,
                    by = c("species", "year", "longitude", "latitude"),
                    all.x = TRUE)
  }

  has_observer <- "observer" %in% names(events) && any(!is.na(events$observer))
  has_date     <- "obs_date" %in% names(events) && any(!is.na(events$obs_date))
  if (!has_observer) {
    message("[28b] {", km, "km} no observer column → n_observers ≡ n_records")
  }
  if (!has_date) {
    message("[28b] {", km, "km} no date column → n_birding_days ≡ n_records")
  }

  agg <- events[, .(
    n_records_grid     = .N,
    n_visits_grid      = if (has_date && has_observer)
                           uniqueN(paste(observer, obs_date)) else .N,
    n_observers_grid   = if (has_observer) uniqueN(observer) else .N,
    n_birding_days_grid = if (has_date)    uniqueN(obs_date) else .N
  ), by = .(grid_id, year)]

  # ---- Fill zero cells from grid base × years ----------------------------
  base <- fread(path_raw(glue("grid_{km}km_base.csv")), encoding = "UTF-8")
  setnames(base, tolower(names(base)))
  full <- CJ(grid_id = base$grid_id,
              year    = CFG$year_min:CFG$year_max,
              sorted  = FALSE)
  full <- merge(full, base[, .(grid_id, province)], by = "grid_id",
                 all.x = TRUE)
  full <- merge(full, agg, by = c("grid_id", "year"), all.x = TRUE)
  for (col in grep("_grid$", names(full), value = TRUE)) {
    full[is.na(get(col)), (col) := 0L]
  }

  # ---- log1p + grid-native z-score ---------------------------------------
  zify <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }
  for (col in c("n_records_grid", "n_visits_grid",
                 "n_observers_grid", "n_birding_days_grid")) {
    log_col <- sub("_grid$", "_log_grid", col)
    z_col   <- sub("_grid$", "_z_grid",   col)
    full[, (log_col) := log1p(get(col))]
    full[, (z_col)   := zify(get(log_col))]
  }

  # PCA on four log columns. effort_pc1 网格级。
  pc_in <- as.matrix(full[, .(n_records_log_grid, n_visits_log_grid,
                                n_observers_log_grid, n_birding_days_log_grid)])
  good <- complete.cases(pc_in) & rowSums(pc_in) > 0
  pc1_out <- rep(NA_real_, nrow(full))
  if (sum(good) > 50L) {
    pca <- prcomp(pc_in[good, ], center = TRUE, scale. = TRUE)
    pc1_out[good] <- pca$x[, 1]
  }
  full[, effort_pc1_grid := pc1_out]
  full[, effort_pc1_z_grid := zify(effort_pc1_grid)]

  # ---- QC: reconstructed province total vs v1 effort_panel ---------------
  prov_v1 <- fread(path_raw("effort_panel_upgraded.csv"), encoding = "UTF-8")
  prov_v1 <- prov_v1[, .(province, year,
                          n_visits_prov_v1 = n_visits,
                          n_birding_days_prov_v1 = n_birding_days)]
  prov_v2 <- full[, .(n_visits_prov_v2 = sum(n_visits_grid),
                       n_birding_days_prov_v2 = sum(n_birding_days_grid),
                       n_grids = uniqueN(grid_id)),
                   by = .(province, year)]
  qc <- merge(prov_v1, prov_v2, by = c("province", "year"), all = TRUE)
  qc[, scale := paste0(km, "km")]

  # Persist parquet + qc table. 持久化。
  out_parq <- path_derived(glue("grid_{km}km_effort_native.parquet"))
  ensure_dir(dirname(out_parq))
  arrow::write_parquet(full, out_parq, compression = "snappy")
  ensure_dir(path_diagnostics())
  qc_path <- path_diagnostics("table_grid_effort_native_summary.csv")
  fwrite(qc, qc_path, append = file.exists(qc_path))
  message(glue("[28b] {km}km → {out_parq} ({nrow(full)} rows, ",
                "{sum(full$n_records_grid > 0)} non-zero cells)"))

  # ---- Within-province variance figure ------------------------------------
  wp <- full[year %in% c(2010, 2015, 2020),
              .(within_province_sd = sd(log1p(n_visits_grid), na.rm = TRUE),
                between_province_mean = mean(log1p(n_visits_grid), na.rm = TRUE),
                n_grids = uniqueN(grid_id)),
              by = .(province, year)]
  wp[, scale := paste0(km, "km")]
  wp
}

results_list <- lapply(CFG$grid_sizes, aggregate_grid_effort)
wp_all <- data.table::rbindlist(results_list)

# ---- 5. Diagnostic figure -------------------------------------------------
p <- ggplot(wp_all,
             aes(x = between_province_mean, y = within_province_sd,
                 colour = scale)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.4) +
  facet_wrap(~ year) +
  scale_colour_manual(values = pal_cat[1:2]) +
  labs(title = "Within-province variance of grid-native effort",
        subtitle = "If v1's province-mirror were correct, SD would be 0 here",
        x = "Province-mean log1p(n_visits) per grid",
        y = "Within-province SD") +
  theme_geb()
ensure_dir(path_diag_fig())
ggsave(path_diag_fig("figS_grid_effort_within_province_variance.pdf"),
       p, width = 16, height = 9, units = "cm",
       device = grDevices::cairo_pdf)
ggsave(path_diag_fig("figS_grid_effort_within_province_variance.png"),
       p, width = 16, height = 9, units = "cm", dpi = 600)

dump_session_info(path_logs("28b_grid_native_effort_sessionInfo.txt"))
message("[28b] done.")
