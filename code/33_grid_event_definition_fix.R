# ============================================================
# Scientific question / 科学问题:
#   v1 (code/07_hazard_model_grid.R:177-178) collapsed multi-species
#   first-arrival events within the same (grid, year) into a single row,
#   AND simply exploded every species across every grid cell of every
#   year. The risk set therefore (a) under-counted ties and (b) was
#   not restricted to the SDM-suitable (species, province) candidate
#   set that v1's province-level analysis used. v2 fixes both:
#     (1) event = earliest record per (species, grid).
#     (2) Cartesian is built ONLY for (species, province) pairs that
#         pass the SDM province threshold filter inherited from v1
#         (encoded in `data/raw/grid_*km_risk_set.csv`), then expanded
#         to all grid cells within those provinces. This preserves v1's
#         "省级 SDM 阈值完全风险集" semantics at the grid scale.
#   v1 网格事件被坍缩，且每物种向所有格子展开（忽略 SDM 适宜性）。
#   v2 改为：事件 = (species, grid) 首次到达；笛卡尔积限制在通过
#   SDM 阈值的 (species, province) 候选集合内再展开到对应省份的全部
#   格子，沿用 v1 的"省级 SDM 阈值完全风险集"语义。
#
# Objective / 分析目标:
#   - Re-classify each (species, grid) pair's first detection year.
#   - Construct a species-grid-year risk set:
#       row = (species, grid, year)
#       event = 1 in the arrival year; 0 before; row dropped after.
#   - Restrict to (species, province) candidate pairs from the v1
#     SDM-thresholded risk set.
#   - Attach grid-NATIVE climate (script 28) and grid-NATIVE effort
#     (script 28b) so that within-province variation is preserved.
#
# Input data / 输入数据:
#   data/raw/events_{50,100}km_grid_assigned.csv   (point-level events)
#   data/raw/grid_{50,100}km_risk_set.csv          (v1 SDM-threshold risk set)
#   data/raw/grid_{50,100}km_base.csv              (grid metadata)
#   data/derived/grid_{50,100}km_effort_native.parquet   (script 28b)
#   data/derived/grid_{50,100}km_climate_native.parquet  (script 28)
#
# Expected output / 预期输出:
#   data/derived/events_{50,100}km_grid_assigned_v2.parquet  (first arrival)
#   data/derived/events_{50,100}km_grid_risk_set_v2.parquet  (full risk set)
#   results/diagnostics/table_grid_event_redefinition_summary.csv
#
# Key assumptions / 关键假设:
#   - Study window 2002–2024 (matches effort panel).
#   - When grid-native climate/effort are missing, falls back to v1
#     province-mirror raw CSVs and logs a warning so the analyst can
#     re-run scripts 28 / 28b before producing the final figures.
#
# Main packages / 主要包: data.table, arrow, glue.
# Output directory / 输出路径: data/derived/, results/diagnostics/.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(glue)
})

source(file.path("code", "utils", "utils_data.R"))

CFG <- list(
  grid_sizes_km = c(50, 100),
  year_min      = 2002,
  year_max      = 2024
)

summary_rows <- list()

rebuild_for_grid <- function(km) {
  message(glue("\n[33] ===== Grid {km} km ====="))

  # ---- 1. Load raw events --------------------------------------------------
  evt_path  <- path_raw(glue("events_{km}km_grid_assigned.csv"))
  base_path <- path_raw(glue("grid_{km}km_base.csv"))
  sdm_rs_path <- path_raw(glue("grid_{km}km_risk_set.csv"))
  if (!file.exists(evt_path))  stop("[33] missing ", evt_path)
  if (!file.exists(base_path)) stop("[33] missing ", base_path)
  events <- fread(evt_path, encoding = "UTF-8")
  base   <- fread(base_path, encoding = "UTF-8")
  setnames(events, tolower(names(events)))
  setnames(base,   tolower(names(base)))
  if (!"grid_id" %in% names(events) && "cell_id" %in% names(events))
    setnames(events, "cell_id", "grid_id")
  if (!"species" %in% names(events) && "species_cn" %in% names(events))
    setnames(events, "species_cn", "species")
  if (!"year" %in% names(events) && "pub_year" %in% names(events))
    setnames(events, "pub_year", "year")
  stopifnot(all(c("grid_id", "species", "year") %in% names(events)))
  events[, year := as.integer(year)]
  events <- events[year >= CFG$year_min & year <= CFG$year_max]

  # ---- 2. First arrival per (species, grid) -------------------------------
  first_arrival <- events[, .(arrival_year = min(year, na.rm = TRUE)),
                          by = .(species, grid_id)]
  message(glue("[33] {km}km: {nrow(first_arrival)} (species, grid) first-arrivals"))

  # ---- 3. SDM province threshold candidate set ----------------------------
  # v1 already encoded the SDM threshold filter in grid_*km_risk_set.csv;
  # we extract the candidate (species, province) pairs from it.
  # 从 v1 的 grid 风险集里抽取通过 SDM 阈值的 (species, province) 候选集合。
  if (file.exists(sdm_rs_path)) {
    # The v1 risk set is large; read minimal columns via arrow if csv lacks
    # streaming. Use data.table::fread with select to keep RAM low.
    candidates <- fread(sdm_rs_path, encoding = "UTF-8",
                         select = c("species", "province"))
    candidates <- unique(candidates)
    setnames(candidates, tolower(names(candidates)))
    message(glue("[33] {km}km: SDM threshold filter loaded ",
                  "({nrow(candidates)} species×province candidate pairs)."))
  } else {
    warning(glue("[33] {km}km: SDM province-threshold file missing ({sdm_rs_path}); ",
                  "falling back to unfiltered cartesian."))
    candidates <- NULL
  }

  # ---- 4. Build cartesian within SDM candidates ---------------------------
  # base table provides (grid_id, province) so we can restrict the
  # expansion to grids inside candidate provinces.
  if (!is.null(candidates)) {
    sp_prov <- candidates[, .(species, province)]
    # Each candidate (species, province) -> all grids in that province
    grid_by_prov <- base[, .(grid_id, province)]
    sp_grid <- merge(sp_prov, grid_by_prov, by = "province",
                      allow.cartesian = TRUE)
    risk <- merge(sp_grid[, .(species, grid_id)],
                   data.table::CJ(grid_id = unique(sp_grid$grid_id),
                                   year    = CFG$year_min:CFG$year_max,
                                   sorted  = FALSE),
                   by = "grid_id", allow.cartesian = TRUE)
    risk <- merge(risk, sp_grid, by = c("species", "grid_id"))
  } else {
    risk <- data.table::CJ(species = unique(first_arrival$species),
                            grid_id = base$grid_id,
                            year    = CFG$year_min:CFG$year_max,
                            sorted  = FALSE)
    risk <- merge(risk, base[, .(grid_id, province)], by = "grid_id",
                  all.x = TRUE)
  }
  message(glue("[33] {km}km: SDM-restricted cartesian = {nrow(risk)} rows ",
                "(vs unfiltered {format(uniqueN(risk$species) * nrow(base) * length(CFG$year_min:CFG$year_max), big.mark=',')} rows)"))

  # ---- 5. Mark event / drop post-arrival rows -----------------------------
  risk <- merge(risk, first_arrival, by = c("species", "grid_id"),
                all.x = TRUE)
  risk <- risk[is.na(arrival_year) | year <= arrival_year]
  risk[, event := as.integer(year == arrival_year)]
  risk[is.na(event), event := 0L]

  # ---- 6. Attach grid-native effort (script 28b) --------------------------
  eff_native_path <- path_derived(glue("grid_{km}km_effort_native.parquet"))
  if (file.exists(eff_native_path)) {
    eff_native <- arrow::read_parquet(eff_native_path) |>
      data.table::as.data.table()
    eff_native[, province := NULL]    # avoid duplicate column on merge
    risk <- merge(risk, eff_native, by = c("grid_id", "year"), all.x = TRUE)
    message(glue("[33] {km}km: grid-native effort joined ({ncol(eff_native)} cols)."))
  } else {
    warning(glue("[33] {km}km: grid-native effort missing — falling back ",
                  "to v1 province-mirror in {path_raw(glue('grid_{km}km_effort.csv'))}. ",
                  "Run script 28b to fix."))
    eff_v1 <- fread(path_raw(glue("grid_{km}km_effort.csv")),
                     encoding = "UTF-8")
    setnames(eff_v1, tolower(names(eff_v1)))
    risk <- merge(risk, eff_v1, by = c("grid_id", "year"), all.x = TRUE,
                   suffixes = c("", ".v1eff"))
  }

  # ---- 7. Attach grid-native climate (script 28) --------------------------
  clim_native_path <- path_derived(glue("grid_{km}km_climate_native.parquet"))
  if (file.exists(clim_native_path)) {
    clim_native <- arrow::read_parquet(clim_native_path) |>
      data.table::as.data.table()
    risk <- merge(risk, clim_native, by = "grid_id", all.x = TRUE)
    message(glue("[33] {km}km: grid-native climate joined."))
  } else {
    warning(glue("[33] {km}km: grid-native climate missing — falling back ",
                  "to v1 province-mirror columns in raw risk set. ",
                  "Run script 28 to fix."))
  }

  # ---- 8. Attach grid base metadata --------------------------------------
  risk <- merge(risk, base, by = c("grid_id", "province"),
                all.x = TRUE)

  # ---- 9. Persist --------------------------------------------------------
  out_pairs <- path_derived(glue("events_{km}km_grid_assigned_v2.parquet"))
  out_risk  <- path_derived(glue("events_{km}km_grid_risk_set_v2.parquet"))
  ensure_dir(dirname(out_pairs))
  arrow::write_parquet(first_arrival, out_pairs, compression = "snappy")
  arrow::write_parquet(risk, out_risk, compression = "snappy")
  message(glue("[33] {km}km: wrote {out_risk} ({nrow(risk)} rows, ",
                "{sum(risk$event)} events)"))

  # ---- 10. Diagnostic ---------------------------------------------------
  summary_rows[[as.character(km)]] <<- data.table::data.table(
    grid_km            = km,
    n_candidates_sp_prov = if (!is.null(candidates)) nrow(candidates) else NA_integer_,
    n_species_grid     = nrow(first_arrival),
    n_events           = sum(risk$event),
    n_risk_rows        = nrow(risk),
    n_raw_events       = nrow(events),
    risk_set_uses_sdm_threshold = !is.null(candidates),
    grid_native_climate_used    = file.exists(clim_native_path),
    grid_native_effort_used     = file.exists(eff_native_path)
  )
}

invisible(lapply(CFG$grid_sizes_km, rebuild_for_grid))

# ---- 11. Diagnostic summary --------------------------------------------
diag_dt <- data.table::rbindlist(summary_rows)
ensure_dir(path_diagnostics())
fwrite(diag_dt,
       path_diagnostics("table_grid_event_redefinition_summary.csv"))
message("[33] Summary:")
print(diag_dt)

dump_session_info(path_logs("33_grid_event_definition_fix_sessionInfo.txt"))
message("[33] done.")
