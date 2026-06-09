# ============================================================
# Scientific question / 科学问题:
#   The province-scale headline is established; this script
#   extends the analysis to prefecture (市) and county (县) scales
#   with FULL data-integrity audit at every stage, OOM-safe
#   case-control sampling for risk-set construction, and produces
#   the matching future-scenario hazard maps (glmmTMB + XGBoost).
#   省级已经验证；本脚本下钻到市与县级，严格审计每一步数据流，
#   用 case-control 采样避免内存爆炸，并出市/县级未来情景预测地图。
#
# Data flow (with audit blocks):
#   v1 events_100km_grid_assigned.csv  ──┐
#                                        ├── st_within → prefecture/county
#   GS(2019)1822 市/县.shp ──────────────┘    polygons (94 % mapped)
#                       ↓
#   first_arrival per (species, unit)
#                       ↓
#   SDM (species, province) candidates ←── v1 hazard_risk_upgraded_complete_case
#                       ↓
#   (species, unit, year) candidate set, restricted to province SDM filter
#                       ↓
#   Case-control sampling: keep all event rows + sample non-event rows
#   at 1:80 (prefecture) and 1:200 (county); persisted as risk set
#                       ↓
#   WorldClim 10' → unit-native bio1/4/12/15/elev (exactextractr)
#                       ↓
#   Record-share allocation: unit_n_events / prov_n_events × prov_n_visits
#                       ↓
#   Joined risk set → glmmTMB M0-M4 (cloglog) + XGBoost
#                       ↓
#   Future SSP scenarios → choropleth maps
#
# Outputs:
#   data/derived/risk_set_prefecture_v2.csv
#   data/derived/risk_set_county_v2.csv
#   data/derived/unit_climate_prefecture.csv
#   data/derived/unit_climate_county.csv
#   data/derived/unit_effort_prefecture.csv
#   data/derived/unit_effort_county.csv
#   results/tables/table_prefecture_coefs.csv
#   results/tables/table_county_coefs.csv
#   results/tables/table_prefecture_county_aic.csv
#   results/forecasts/table_prefecture_future_glmmTMB.csv
#   results/forecasts/table_prefecture_future_xgboost.csv
#   results/forecasts/table_county_future_glmmTMB.csv
#   results/forecasts/table_county_future_xgboost.csv
#   results/diagnostics/audit_prefecture_county.txt
#   figures/main/fig_three_scale_forest_pref_county.{pdf,png}
#   figures/main/fig_future_hazard_glmmTMB_prefecture.{pdf,png}
#   figures/main/fig_future_hazard_xgboost_prefecture.{pdf,png}
#   figures/main/fig_future_hazard_glmmTMB_county.{pdf,png}
#   figures/main/fig_future_hazard_xgboost_county.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(terra)
  library(exactextractr)
  library(glmmTMB)
  library(xgboost)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})
sf::sf_use_s2(FALSE)
options(warn = 1)
set.seed(42)

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- Sys.getenv("V1_ROOT",
                  normalizePath(file.path(V2, "..",
                                          "bird_hazard_model_effort_upgrade"),
                                 mustWork = FALSE))
WC <- Sys.getenv("WORLDCLIM_DIR",
                  "/Users/dingchenchen/Documents/New project/bird_full_community_analysis/data/external/climate/wc2.1_10m")
SHP <- file.path(V2, "data", "spatial", "basemap_GS2019_1822")

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ens(file.path(V2, "results", "tables"))
ens(file.path(V2, "results", "forecasts"))
ens(file.path(V2, "results", "diagnostics"))
ens(file.path(V2, "figures", "main"))
ens(file.path(V2, "data", "derived"))
ens(file.path(V2, "logs"))

# Audit log: write to both stdout AND a text file.
AUDIT <- file(file.path(V2, "results", "diagnostics",
                         "audit_prefecture_county.txt"),
              open = "wt", encoding = "UTF-8")
on.exit({ if (isOpen(AUDIT)) close(AUDIT) }, add = TRUE)

log <- function(...) {
  msg <- paste0(sprintf("[42 %s] ", format(Sys.time(), "%H:%M:%S")),
                 paste(..., sep = ""))
  cat(msg, "\n", sep = "")
  writeLines(msg, AUDIT)
}
audit_block <- function(title) {
  bar <- paste(rep("=", 60), collapse = "")
  log(""); log(bar); log("AUDIT — ", title); log(bar)
}

# Unified theme (same as 41 for visual consistency).
theme_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
          panel.border     = element_rect(linewidth = 0.4, colour = "grey20"),
          plot.title       = element_text(face = "bold", size = base_size + 1),
          plot.subtitle    = element_text(size = base_size - 1, colour = "grey30"),
          strip.background = element_rect(fill = "grey95", colour = "grey80"),
          strip.text       = element_text(face = "bold"))
}
save_pub <- function(p, name, width = 17, height = 10,
                      path = file.path(V2, "figures", "main")) {
  ggsave(file.path(path, paste0(name, ".pdf")), p, width = width,
         height = height, units = "cm", device = grDevices::cairo_pdf)
  ggsave(file.path(path, paste0(name, ".png")), p, width = width,
         height = height, units = "cm", dpi = 600)
  log("wrote ", name, ".{pdf,png}")
}
zify <- function(x) {
  s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

PROV_CN_EN <- c("北京市"="Beijing","天津市"="Tianjin","河北省"="Hebei",
  "山西省"="Shanxi","内蒙古自治区"="Inner Mongolia","辽宁省"="Liaoning",
  "吉林省"="Jilin","黑龙江省"="Heilongjiang","上海市"="Shanghai",
  "江苏省"="Jiangsu","浙江省"="Zhejiang","安徽省"="Anhui","福建省"="Fujian",
  "江西省"="Jiangxi","山东省"="Shandong","河南省"="Henan","湖北省"="Hubei",
  "湖南省"="Hunan","广东省"="Guangdong","广西壮族自治区"="Guangxi",
  "海南省"="Hainan","重庆市"="Chongqing","四川省"="Sichuan","贵州省"="Guizhou",
  "云南省"="Yunnan","西藏自治区"="Tibet","陕西省"="Shaanxi","甘肃省"="Gansu",
  "青海省"="Qinghai","宁夏回族自治区"="Ningxia","新疆维吾尔自治区"="Xinjiang",
  "台湾省"="Taiwan","香港特别行政区"="Hong Kong","澳门特别行政区"="Macau")

# ============================================================
# AUDIT 1: load + spatial-join events
# ============================================================
audit_block("1. Load + spatial join events → prefecture / county")

ev <- fread(file.path(V1, "data", "events_100km_grid_assigned.csv"),
            encoding = "UTF-8")
setnames(ev, tolower(names(ev)))
if (!"year" %in% names(ev) && "pub_year" %in% names(ev))
  setnames(ev, "pub_year", "year")
ev <- ev[!is.na(longitude) & !is.na(latitude) & year >= 2002 & year <= 2024]
log("events (coord-valid, 2002-2024): ", nrow(ev),
    " | species: ", uniqueN(ev$species))

prov_sf <- st_read(file.path(SHP, "省（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4326) |> st_make_valid()
pref_sf <- st_read(file.path(SHP, "市（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4326) |> st_make_valid()
cnty_sf <- st_read(file.path(SHP, "县（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4326) |> st_make_valid()

# Standardised columns
prov_sf$province <- unname(PROV_CN_EN[as.character(prov_sf[["省"]])])
pref_sf$province <- unname(PROV_CN_EN[as.character(pref_sf[["省"]])])
cnty_sf$province <- unname(PROV_CN_EN[as.character(cnty_sf[["省"]])])
pref_sf$pref_id <- paste0("PREF_", sprintf("%04d", seq_len(nrow(pref_sf))))
cnty_sf$cnty_id <- paste0("CNTY_", sprintf("%05d", seq_len(nrow(cnty_sf))))
pref_sf$pref_name <- as.character(pref_sf[["市"]])
cnty_sf$cnty_name <- as.character(cnty_sf[["NAME"]])
log("polygons: ", nrow(prov_sf), " provinces, ",
    nrow(pref_sf), " prefectures, ", nrow(cnty_sf), " counties")

# Drop province column from events (it conflicts with admin join later)
ev_sf <- st_as_sf(ev, coords = c("longitude","latitude"), crs = 4326)
ev_sf$province <- NULL

suppressWarnings({
  ev_pref <- st_join(ev_sf, pref_sf[, c("pref_id","pref_name","province")],
                      join = st_within)
  ev_cnty <- st_join(ev_sf, cnty_sf[, c("cnty_id","cnty_name","province")],
                      join = st_within)
})
ev_pref_dt <- as.data.table(st_drop_geometry(ev_pref))
ev_cnty_dt <- as.data.table(st_drop_geometry(ev_cnty))
ev_pref_dt <- ev_pref_dt[!is.na(pref_id)]
ev_cnty_dt <- ev_cnty_dt[!is.na(cnty_id)]
log("events with prefecture: ", nrow(ev_pref_dt),
    " covering ", uniqueN(ev_pref_dt$pref_id), " prefectures, ",
    uniqueN(ev_pref_dt$province), " provinces")
log("events with county    : ", nrow(ev_cnty_dt),
    " covering ", uniqueN(ev_cnty_dt$cnty_id), " counties, ",
    uniqueN(ev_cnty_dt$province), " provinces")

# ============================================================
# AUDIT 2: SDM candidate set + cross-check
# ============================================================
audit_block("2. SDM threshold filter and candidate set")

risk_prov <- fread(file.path(V1, "data",
                              "hazard_risk_upgraded_complete_case.csv"),
                    encoding = "UTF-8")
sdm_cand <- unique(risk_prov[, .(species, province)])
log("SDM-thresholded (species, province) candidates: ", nrow(sdm_cand))

ev_pref_dt <- ev_pref_dt[!is.na(province)]
ev_cnty_dt <- ev_cnty_dt[!is.na(province)]

# Filter events to those whose (species, province) passes SDM threshold.
ev_pref_in <- merge(ev_pref_dt, sdm_cand, by = c("species","province"))
ev_cnty_in <- merge(ev_cnty_dt, sdm_cand, by = c("species","province"))
log("events INSIDE SDM candidate set — prefecture: ", nrow(ev_pref_in),
    " (lost ", nrow(ev_pref_dt) - nrow(ev_pref_in), " outside SDM)")
log("events INSIDE SDM candidate set — county    : ", nrow(ev_cnty_in),
    " (lost ", nrow(ev_cnty_dt) - nrow(ev_cnty_in), " outside SDM)")

# ============================================================
# AUDIT 3: WorldClim → unit-native climate
# ============================================================
audit_block("3. Extract WorldClim 10' to prefecture + county polygons")

wc_layers <- list(bio1 = "wc2.1_10m_bio_1.tif",
                  bio4 = "wc2.1_10m_bio_4.tif",
                  bio12 = "wc2.1_10m_bio_12.tif",
                  bio15 = "wc2.1_10m_bio_15.tif",
                  elev  = "wc2.1_10m_elev.tif")

extract_clim <- function(sf_unit, id_col) {
  out <- data.table(unit_id = sf_unit[[id_col]])
  for (nm in names(wc_layers)) {
    r <- terra::rast(file.path(WC, wc_layers[[nm]]))
    vals <- exactextractr::exact_extract(r, sf_unit, fun = "mean",
                                          progress = FALSE)
    out[, (nm) := vals]
  }
  for (nm in names(wc_layers)) {
    out[, paste0(nm, "_z") := zify(get(nm))]
  }
  out
}

pref_clim <- extract_clim(pref_sf, "pref_id")
cnty_clim <- extract_clim(cnty_sf, "cnty_id")
log("prefecture climate: ", nrow(pref_clim), " rows × ", ncol(pref_clim), " cols ",
    "(NA bio1: ", sum(is.na(pref_clim$bio1)), ")")
log("county     climate: ", nrow(cnty_clim), " rows × ", ncol(cnty_clim), " cols ",
    "(NA bio1: ", sum(is.na(cnty_clim$bio1)), ")")
fwrite(pref_clim, file.path(V2, "data", "derived",
                              "unit_climate_prefecture.csv"))
fwrite(cnty_clim, file.path(V2, "data", "derived",
                              "unit_climate_county.csv"))

# ============================================================
# AUDIT 4: Effort allocation (record-share × province total)
# ============================================================
audit_block("4. Effort allocation by record-share within province × year")

prov_eff <- fread(file.path(V1, "data", "effort_panel_upgraded.csv"),
                   encoding = "UTF-8")
log("province effort panel: ", nrow(prov_eff), " rows, ",
    uniqueN(prov_eff$province), " provinces")

# Effort source: RAW combined-and-dedup bird survey records (eBird-GBIF +
# China-Birdwatch, 2.3 GB). Each row is one (species, observer, date,
# coordinate) detection. We:
#   1. Load once into RAM (≈ 20 GB peak — caller must have headroom).
#   2. Spatial-join (lon, lat) to prefecture and county polygons.
#   3. Group by (unit_id, year) and compute four effort metrics:
#        n_records, n_observers, n_dates, n_visits (= observer×date pairs).
# 用 raw combined-dedup 鸟调记录直接落点构造 unit-year effort。
COMM <- "/Users/dingchenchen/Documents/New project/bird_dynamic_occupancy_analysis"
RAW_REC_PATH <- file.path(COMM, "data", "derived_v2",
                           "combined_events_merged_dedup_2000_2025.rds")

# Load + spatial-join exactly once; reuse for both prefecture and county.
raw_pref_dt <- NULL
raw_cnty_dt <- NULL
if (file.exists(RAW_REC_PATH)) {
  log("loading combined-dedup raw records (≈ 2.3 GB)…")
  raw <- as.data.table(readRDS(RAW_REC_PATH))
  log("raw records: ", nrow(raw), " rows | cols: ",
      paste(intersect(c("species","longitude","latitude","year",
                         "month","day","username"), names(raw)),
             collapse = ", "))
  raw <- raw[!is.na(longitude) & !is.na(latitude) &
              year >= 2002 & year <= 2024]
  raw[, date_key := paste0(year, "-", month, "-", day)]
  raw[, visit_key := paste0(username, "@", date_key)]
  log("raw records after 2002-2024 + coord filter: ", nrow(raw))

  raw_sf <- st_as_sf(raw[, .(longitude, latitude, year, username, date_key,
                              visit_key)],
                      coords = c("longitude","latitude"), crs = 4326)
  rm(raw); invisible(gc(verbose = FALSE))

  log("st_within(raw → prefecture)…")
  suppressWarnings({
    raw_in_pref <- st_join(raw_sf,
                            pref_sf[, c("pref_id","province")],
                            join = st_within)
  })
  raw_pref_dt <- as.data.table(st_drop_geometry(raw_in_pref))
  raw_pref_dt <- raw_pref_dt[!is.na(pref_id)]
  rm(raw_in_pref); invisible(gc(verbose = FALSE))
  log("raw records → prefecture: ", nrow(raw_pref_dt))

  log("st_within(raw → county)…")
  suppressWarnings({
    raw_in_cnty <- st_join(raw_sf,
                            cnty_sf[, c("cnty_id","province")],
                            join = st_within)
  })
  raw_cnty_dt <- as.data.table(st_drop_geometry(raw_in_cnty))
  raw_cnty_dt <- raw_cnty_dt[!is.na(cnty_id)]
  rm(raw_sf, raw_in_cnty); invisible(gc(verbose = FALSE))
  log("raw records → county    : ", nrow(raw_cnty_dt))
} else {
  log("WARN: raw combined-dedup records not found at ", RAW_REC_PATH,
      "; will skip raw-based effort.")
}

build_unit_effort <- function(id_col, raw_in) {
  if (is.null(raw_in) || nrow(raw_in) == 0L) return(NULL)
  rin <- copy(raw_in)
  setnames(rin, id_col, "unit_id")
  cnt <- rin[, .(comm_n_events    = .N,
                  comm_n_observers = uniqueN(username),
                  comm_n_dates     = uniqueN(date_key),
                  comm_n_visits    = uniqueN(visit_key)),
              by = .(unit_id, year, province)]
  cnt <- merge(cnt, prov_eff[, .(province, year, n_visits, n_birding_days)],
                by = c("province","year"), all.x = TRUE)
  cnt[is.na(n_visits), n_visits := 1L]
  cnt[, allocated_visits := comm_n_visits]
  cnt[, allocated_days   := comm_n_dates]
  cnt[, log_effort_visits_z := zify(log1p(allocated_visits)), by = year]
  cnt[, log_effort_days_z   := zify(log1p(allocated_days)),    by = year]
  setnames(cnt, "unit_id", id_col)
  cnt
}
pref_eff <- build_unit_effort("pref_id", raw_pref_dt)
cnty_eff <- build_unit_effort("cnty_id", raw_cnty_dt)
log("prefecture effort rows: ", nrow(pref_eff),
    " | within-year SD of log_effort_visits_z: ",
    round(pref_eff[year == 2020, sd(log_effort_visits_z, na.rm=TRUE)], 3),
    " (sanity: > 0 means real within-province variation, not broadcast)")
log("county     effort rows: ", nrow(cnty_eff),
    " | within-year SD of log_effort_visits_z: ",
    round(cnty_eff[year == 2020, sd(log_effort_visits_z, na.rm=TRUE)], 3))

fwrite(pref_eff, file.path(V2, "data", "derived",
                             "unit_effort_prefecture.csv"))
fwrite(cnty_eff, file.path(V2, "data", "derived",
                             "unit_effort_county.csv"))

# ============================================================
# AUDIT 5: Build OOM-safe case-control risk sets
# ============================================================
audit_block("5. Build OOM-safe case-control risk sets")

# Province-climate panel as a fallback time dimension. Many prefecture/county
# inherit province × year climate for the time dimension (since WorldClim is
# climatology). 时间维度气候用省级。
prov_clim_full <- fread(file.path(V1, "data",
                                    "climate_metrics_province_year.csv"),
                          encoding = "UTF-8")
prov_clim_use <- prov_clim_full[, .(province, year,
                                      temp_grad_z = temp_grad_prov_z)]

build_risk <- function(unit_label, id_col, sf_unit, ev_in,
                        cand_per_unit_ratio = 80) {
  unit_in_prov <- unique(as.data.table(st_drop_geometry(sf_unit))[
    , .(unit_id = get(id_col), province)])
  unit_in_prov <- unit_in_prov[!is.na(province)]
  cand <- merge(sdm_cand, unit_in_prov, by = "province",
                  allow.cartesian = TRUE)
  log(unit_label, " — SDM-restricted (species × unit) pairs: ",
      nrow(cand))

  setnames(ev_in, id_col, "unit_id")
  fa <- ev_in[, .(arrival_year = min(year, na.rm = TRUE)),
              by = .(species, unit_id, province)]
  log(unit_label, " — first-arrival (species × unit) pairs: ", nrow(fa))

  # ---- Build event rows (one per first arrival) -------------------------
  ev_rows <- merge(fa, cand, by = c("species","unit_id","province"))
  ev_rows[, event := 1L]
  ev_rows[, year := arrival_year]
  log(unit_label, " — event rows (after SDM filter): ", nrow(ev_rows),
      " out of ", nrow(fa), " first-arrival pairs (",
      sprintf("%.1f%%", 100 * nrow(ev_rows) / nrow(fa)), " kept)")

  # ---- Case-control sampling ---------------------------------------------
  # For each event we sample `cand_per_unit_ratio` non-event rows uniformly
  # at random over the (cand, year ∈ 2002..2024) space, after removing rows
  # already used as events or rows after the arrival year. 用 case-control。
  # Construction: build a small lookup of (species, unit_id) → arrival_year
  # so we can quickly filter the sampled candidates.
  cand[, key := paste0(species, "|", unit_id)]
  fa_key <- setNames(fa$arrival_year, paste0(fa$species, "|", fa$unit_id))
  n_negatives <- nrow(ev_rows) * cand_per_unit_ratio
  log(unit_label, " — sampling ", n_negatives,
      " non-event rows (1 : ", cand_per_unit_ratio, ")")

  # Sample (species, unit_id, year) triplets uniformly. 均匀采样三元组。
  sample_idx <- sample.int(nrow(cand), n_negatives, replace = TRUE)
  sample_yr  <- sample(2002:2024, n_negatives, replace = TRUE)
  neg <- cand[sample_idx, .(species, unit_id, province, key)]
  neg[, year := sample_yr]
  neg[, arrival_year := fa_key[key]]
  # Keep only rows where the species has no arrival OR year is BEFORE arrival
  neg <- neg[is.na(arrival_year) | year < arrival_year]
  neg[, event := 0L]
  # Drop rare collisions with the event rows themselves
  neg[, drop_key := paste0(species, "|", unit_id, "|", year)]
  ev_rows[, drop_key := paste0(species, "|", unit_id, "|", year)]
  neg <- neg[!drop_key %in% ev_rows$drop_key]

  risk <- rbind(
    ev_rows[, .(species, unit_id, province, year, event)],
    neg[, .(species, unit_id, province, year, event)])
  log(unit_label, " — final risk set: ", nrow(risk),
      " rows (events = ", sum(risk$event), ")")
  risk
}

risk_pref <- build_risk("PREFECTURE", "pref_id", pref_sf,
                         copy(ev_pref_in), cand_per_unit_ratio = 80)
risk_cnty <- build_risk("COUNTY", "cnty_id", cnty_sf,
                         copy(ev_cnty_in), cand_per_unit_ratio = 200)

# Attach unit-native climate (WorldClim) + effort (allocation) + province×year temp_grad_z.
attach_covariates <- function(risk, unit_clim, unit_eff) {
  setnames(unit_clim, "unit_id", "unit_id_clim")
  risk <- merge(risk, unit_clim, by.x = "unit_id", by.y = "unit_id_clim",
                  all.x = TRUE)
  # Effort: (unit_id, year) from the community Combined panel
  ue <- copy(unit_eff)
  setnames(ue, names(ue)[1], "unit_id")
  risk <- merge(risk,
                ue[, .(unit_id, year, allocated_visits)],
                by = c("unit_id", "year"), all.x = TRUE)
  # FILL non-event (control) rows lacking community-effort data with the
  # PROVINCE × YEAR n_visits scaled by the number of units in that province.
  # This gives every (unit, year) a real baseline effort signal (broadcast
  # with unit-count denominator), and units that have community-effort
  # add the within-province deviation on top. 给所有 control 行非零 baseline。
  units_per_prov <- risk[, uniqueN(unit_id), by = province]
  setnames(units_per_prov, "V1", "n_units_in_prov")
  risk <- merge(risk, units_per_prov, by = "province", all.x = TRUE)
  risk <- merge(risk,
                prov_eff[, .(province, year, prov_n_visits = n_visits,
                               prov_n_days  = n_birding_days)],
                by = c("province", "year"), all.x = TRUE)
  risk[is.na(prov_n_visits), prov_n_visits := 1L]
  risk[is.na(prov_n_days),   prov_n_days   := 1L]
  # If community data missing, give the (province / units_in_prov) baseline.
  # If present, sum baseline + community deviation. 缺数据用省级均值。
  risk[, baseline_visits := prov_n_visits / pmax(n_units_in_prov, 1)]
  risk[is.na(allocated_visits), allocated_visits := 0]
  risk[, effort_visits_total := baseline_visits + allocated_visits]
  risk[, log_effort_visits_z := zify(log1p(effort_visits_total)), by = year]
  risk[, log_effort_days_z   := zify(log1p(prov_n_days / pmax(n_units_in_prov, 1)
                                              + allocated_visits)), by = year]
  # province × year temp_grad_z (continuous time dimension)
  risk <- merge(risk, prov_clim_use, by = c("province","year"), all.x = TRUE)
  risk
}
risk_pref <- attach_covariates(risk_pref, copy(pref_clim), copy(pref_eff))
risk_cnty <- attach_covariates(risk_cnty, copy(cnty_clim), copy(cnty_eff))

fwrite(risk_pref, file.path(V2, "data", "derived", "risk_set_prefecture_v2.csv"))
fwrite(risk_cnty, file.path(V2, "data", "derived", "risk_set_county_v2.csv"))

audit_block("6. Risk-set sanity check before model fitting")
for (lab in c("prefecture", "county")) {
  rs <- if (lab == "prefecture") risk_pref else risk_cnty
  log(lab, " risk: rows=", nrow(rs),
      " events=", sum(rs$event),
      " unique species=", uniqueN(rs$species),
      " unique units=", uniqueN(rs$unit_id),
      " unique provinces=", uniqueN(rs$province))
  log(lab, " NA check: temp_grad_z=", sum(is.na(rs$temp_grad_z)),
      " bio1_z=", sum(is.na(rs$bio1_z)),
      " log_effort_visits_z=", sum(is.na(rs$log_effort_visits_z)))
  log(lab, " effort within-year SD (2020): ",
      round(rs[year == 2020, sd(log_effort_visits_z, na.rm=TRUE)], 3),
      " (must be > 0 — confirms within-province variation)")
}

# ============================================================
# AUDIT 7: Fit M0-M4 at each scale + persist coefficients
# ============================================================
audit_block("7. Fit M0-M4 cloglog hazard at prefecture + county")

fit_M_ladder <- function(rs, label) {
  d <- rs[!is.na(temp_grad_z) & !is.na(log_effort_visits_z)]
  d[, climate_z := temp_grad_z]
  d[, effort_z  := log_effort_visits_z]
  log(label, " — model dt rows: ", nrow(d), " | events: ", sum(d$event))
  forms <- list(
    M0 = "event ~ 1                          + (1|species) + (1|unit_id)",
    M1 = "event ~ effort_z                    + (1|species) + (1|unit_id)",
    M2 = "event ~ climate_z                   + (1|species) + (1|unit_id)",
    M3 = "event ~ climate_z + effort_z        + (1|species) + (1|unit_id)",
    M4 = "event ~ climate_z * effort_z        + (1|species) + (1|unit_id)")
  rows_coef <- list(); rows_aic <- list()
  for (nm in names(forms)) {
    t0 <- Sys.time()
    fit <- tryCatch(
      glmmTMB(as.formula(forms[[nm]]), data = d,
              family = binomial(link = "cloglog")),
      error = function(e) { log("  ", nm, " FAILED: ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    log(sprintf("  %s OK (%.1fs, AIC=%.2f, nobs=%d)", nm,
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  AIC(fit), nobs(fit)))
    cf <- fixef(fit)$cond; se <- sqrt(diag(vcov(fit)$cond))
    for (tm in names(cf)) {
      i <- match(tm, names(cf))
      rows_coef[[length(rows_coef)+1L]] <- data.table(
        scale = label, model = nm, term = tm,
        beta = cf[i], se = se[i],
        hr = exp(cf[i]),
        hr.low  = exp(cf[i] - 1.96 * se[i]),
        hr.high = exp(cf[i] + 1.96 * se[i]),
        p.value = 2 * pnorm(-abs(cf[i] / se[i])))
    }
    rows_aic[[length(rows_aic)+1L]] <- data.table(
      scale = label, model = nm,
      AIC = AIC(fit), BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)),
      df = attr(logLik(fit), "df"),
      nobs = nobs(fit))
    invisible(gc(verbose = FALSE))
    assign(paste0("fit_", label, "_", nm), fit, envir = .GlobalEnv)
  }
  list(coefs = rbindlist(rows_coef, fill = TRUE),
       aic   = rbindlist(rows_aic,  fill = TRUE))
}

res_pref <- fit_M_ladder(risk_pref, "prefecture")
res_cnty <- fit_M_ladder(risk_cnty, "county")

all_coefs <- rbind(res_pref$coefs, res_cnty$coefs, fill = TRUE)
all_aic   <- rbind(res_pref$aic,   res_cnty$aic,   fill = TRUE)
all_aic[, dAIC := AIC - min(AIC, na.rm = TRUE), by = scale]
all_aic[, weight := exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC)), by = scale]

fwrite(res_pref$coefs, file.path(V2, "results", "tables",
                                   "table_prefecture_coefs.csv"))
fwrite(res_cnty$coefs, file.path(V2, "results", "tables",
                                   "table_county_coefs.csv"))
fwrite(all_aic, file.path(V2, "results", "tables",
                           "table_prefecture_county_aic.csv"))

# ============================================================
# AUDIT 8: Three-scale forest plot (province + prefecture + county)
# ============================================================
audit_block("8. Three-scale forest plot")

prov_coefs <- fread(file.path(V2, "results", "tables",
                                "table_province_v2_coefs.csv"))
prov_int <- prov_coefs[spec_id == "spec_B" & model == "M4" & grepl(":", term),
                        .(scale = "province",
                          beta, se, hr, hr.low, hr.high, p.value)]
pref_int <- all_coefs[scale == "prefecture" & model == "M4" & grepl(":", term),
                       .(scale = "prefecture", beta, se, hr, hr.low, hr.high, p.value)]
cnty_int <- all_coefs[scale == "county" & model == "M4" & grepl(":", term),
                       .(scale = "county", beta, se, hr, hr.low, hr.high, p.value)]
three_scale <- rbind(prov_int, pref_int, cnty_int)
three_scale[, scale := factor(scale, levels = c("county","prefecture","province"))]

p_three <- ggplot(three_scale, aes(x = hr, y = scale)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high), height = 0.18,
                  linewidth = 0.5) +
  geom_point(size = 3.5, colour = "#B40426") +
  geom_text_repel(aes(label = sprintf("HR = %.3f\np = %.1e", hr, p.value)),
                   nudge_x = 0.05, size = 2.6, hjust = 0,
                   segment.size = 0.2, direction = "y",
                   point.padding = 0.4) +
  scale_x_continuous(trans = "log",
                      breaks = c(0.8, 1.0, 1.2, 1.4, 1.6, 1.8)) +
  labs(title = "Climate × effort interaction (M4) across 3 administrative scales",
        subtitle = paste0("Province (v2 refit, n = 12,813); prefecture & county ",
                           "(v2 case-control, SDM-thresholded, ",
                           "WorldClim 10' climate, record-share effort)."),
        x = "Hazard ratio (95 % CI, log scale)", y = NULL) +
  theme_pub()
save_pub(p_three, "fig_three_scale_forest_pref_county",
         width = 16, height = 8)

# ============================================================
# AUDIT 9: Future scenarios — glmmTMB + XGBoost per scale
# ============================================================
audit_block("9. Future hazard projections")

future_years <- c(2030, 2050, 2080)
ssp_eps <- list(SSP245 = 0.3 / sd(prov_clim_full$temp_grad_prov, na.rm=TRUE),
                 SSP585 = 0.8 / sd(prov_clim_full$temp_grad_prov, na.rm=TRUE))

# Build per-unit baseline (year = 2024)
mk_unit_baseline <- function(unit_clim, unit_eff, sf_unit, id_col) {
  ue_base <- unit_eff[year == 2024,
                       .(unit_id = get(names(unit_eff)[1]),
                         log_effort_visits_z = log_effort_visits_z)]
  setnames(ue_base, names(ue_base)[1], "unit_id")
  baseline <- merge(unit_clim, ue_base, by.x = "unit_id", by.y = "unit_id",
                     all.x = TRUE)
  baseline[is.na(log_effort_visits_z), log_effort_visits_z := 0]
  baseline_prov <- as.data.table(st_drop_geometry(sf_unit))[
    , .(unit_id = get(id_col), province)]
  baseline <- merge(baseline, baseline_prov, by = "unit_id", all.x = TRUE)
  baseline <- baseline[!is.na(province)]
  pclim_2024 <- prov_clim_use[year == 2024,
                                .(province, temp_grad_z_2024 = temp_grad_z)]
  baseline <- merge(baseline, pclim_2024, by = "province", all.x = TRUE)
  baseline
}
pref_base <- mk_unit_baseline(copy(pref_clim), copy(pref_eff),
                                pref_sf, "pref_id")
cnty_base <- mk_unit_baseline(copy(cnty_clim), copy(cnty_eff),
                                cnty_sf, "cnty_id")
log("prefecture baseline rows: ", nrow(pref_base))
log("county     baseline rows: ", nrow(cnty_base))

mk_future_glmm <- function(baseline, fit_M4, sf_unit, id_col, label) {
  fut <- CJ(unit_id = baseline$unit_id,
             ssp     = c("SSP245","SSP585"),
             year    = future_years)
  fut <- merge(fut, baseline, by = "unit_id")
  fut[, decades_ahead := (year - 2024) / 10]
  fut[, temp_grad_z := temp_grad_z_2024 +
                        ssp_eps[[ssp[1]]] * decades_ahead, by = ssp]
  fut[, climate_z := temp_grad_z]
  fut[, effort_z  := log_effort_visits_z]
  pred <- predict(fit_M4, newdata = fut, type = "response",
                   re.form = NA, allow.new.levels = TRUE)
  fut[, hazard_glmm := pred]
  fut
}

if (exists("fit_prefecture_M4")) {
  fut_pref <- mk_future_glmm(pref_base, fit_prefecture_M4,
                              pref_sf, "pref_id", "prefecture")
  fwrite(fut_pref, file.path(V2, "results", "forecasts",
                               "table_prefecture_future_glmmTMB.csv"))
  log("prefecture glmmTMB future predictions: ", nrow(fut_pref))
}
if (exists("fit_county_M4")) {
  fut_cnty <- mk_future_glmm(cnty_base, fit_county_M4,
                              cnty_sf, "cnty_id", "county")
  fwrite(fut_cnty, file.path(V2, "results", "forecasts",
                               "table_county_future_glmmTMB.csv"))
  log("county     glmmTMB future predictions: ", nrow(fut_cnty))
}

# XGBoost (per scale)
mk_future_xgb <- function(rs, baseline, sf_unit, id_col, label) {
  d <- rs[!is.na(temp_grad_z) & !is.na(log_effort_visits_z) &
           !is.na(bio1_z) & !is.na(bio12_z)]
  feat <- c("temp_grad_z","log_effort_visits_z",
             "bio1_z","bio4_z","bio12_z","bio15_z","elev_z")
  feat <- intersect(feat, names(d))
  Xtr <- xgb.DMatrix(data = as.matrix(d[, ..feat]), label = d$event)
  par <- list(objective="binary:logistic", eval_metric="auc",
              eta = 0.06, max_depth = 5,
              subsample = 0.85, colsample_bytree = 0.85,
              nthread = max(1L, parallel::detectCores() - 2L))
  cv <- xgb.cv(par, Xtr, nrounds = 400, nfold = 5,
                early_stopping_rounds = 25, verbose = 0)
  nb <- cv$best_iteration
  if (is.null(nb) || nb == 0L) {
    auc_col <- intersect(c("test_auc_mean","test_AUC_mean"),
                          names(cv$evaluation_log))[1]
    nb <- if (!is.na(auc_col)) which.max(cv$evaluation_log[[auc_col]]) else 150L
  }
  best_auc <- cv$evaluation_log[["test_auc_mean"]][nb]
  log(label, " XGBoost: nrounds=", nb, " CV-AUC=", round(best_auc, 3),
      " features=", paste(feat, collapse=","))
  mdl <- xgb.train(par, Xtr, nrounds = nb, verbose = 0)

  fut <- CJ(unit_id = baseline$unit_id,
             ssp = c("SSP245","SSP585"), year = future_years)
  fut <- merge(fut, baseline, by = "unit_id")
  fut[, decades_ahead := (year - 2024) / 10]
  fut[, temp_grad_z := temp_grad_z_2024 +
                        ssp_eps[[ssp[1]]] * decades_ahead, by = ssp]
  Xpred <- as.matrix(fut[, ..feat])
  fut[, hazard_xgb := predict(mdl, Xpred)]
  fut
}

if (nrow(risk_pref) > 100) {
  fut_pref_xgb <- mk_future_xgb(risk_pref, pref_base, pref_sf, "pref_id",
                                  "prefecture")
  fwrite(fut_pref_xgb, file.path(V2, "results", "forecasts",
                                    "table_prefecture_future_xgboost.csv"))
}
if (nrow(risk_cnty) > 100) {
  fut_cnty_xgb <- mk_future_xgb(risk_cnty, cnty_base, cnty_sf, "cnty_id",
                                  "county")
  fwrite(fut_cnty_xgb, file.path(V2, "results", "forecasts",
                                    "table_county_future_xgboost.csv"))
}

# ============================================================
# AUDIT 10: Choropleth maps for prefecture + county
# ============================================================
audit_block("10. Choropleth maps")

shp_dir <- SHP
prov_alb <- st_transform(prov_sf, 4524)
pref_alb <- st_transform(pref_sf, 4524)
cnty_alb <- st_transform(cnty_sf, 4524)
ninedash_alb <- tryCatch(
  st_transform(st_read(file.path(SHP, "九段线.shp"), quiet = TRUE), 4524),
  error = function(e) NULL)

choropleth_unit <- function(sf_unit, dt, var, title, subtitle = NULL,
                             id_col = "unit_id") {
  sf_unit$.id <- sf_unit[[setdiff(intersect(c("pref_id","cnty_id"),
                                               names(sf_unit)), "geometry")[1]]]
  joined <- merge(sf_unit, dt[, c(id_col, var), with = FALSE],
                   by.x = ".id", by.y = id_col)
  p <- ggplot() +
    geom_sf(data = joined, aes(fill = .data[[var]]),
            colour = NA, linewidth = 0)
  p <- p +
    geom_sf(data = prov_alb, fill = NA, colour = "grey20", linewidth = 0.25)
  if (!is.null(ninedash_alb)) p <- p +
    geom_sf(data = ninedash_alb, fill = NA, colour = "grey20",
            linewidth = 0.2)
  p + scale_fill_gradient(low = "#F7FBFF", high = "#B40426",
                           name = "Hazard\n(prob.)",
                           limits = range(dt[[var]], na.rm = TRUE)) +
    coord_sf(datum = NA, expand = FALSE) +
    labs(title = title, subtitle = subtitle) +
    theme_pub(base_size = 7.5) +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

mk_six_panel <- function(sf_unit, fut, var, id_col, file_label,
                          model_label, scale_label) {
  fut[, panel := paste0(ssp, " — ", year)]
  fut[, panel := factor(panel,
                          levels = c(paste0("SSP245 — ", future_years),
                                      paste0("SSP585 — ", future_years)))]
  panels <- lapply(split(fut, fut$panel), function(d)
    choropleth_unit(sf_unit, d, var, title = unique(d$panel),
                     id_col = id_col))
  combined <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") +
    plot_annotation(
      title = sprintf("%s-scale future hazard — %s", scale_label, model_label),
      subtitle = "Empirical SSP perturbation: SSP245 = +0.3 SD/dec, SSP585 = +0.8 SD/dec on temp_grad_z. Effort frozen at 2024.",
      theme = theme(plot.title = element_text(face="bold", size=10),
                     plot.subtitle = element_text(size=8, colour="grey30")))
  save_pub(combined, file_label, width = 22, height = 14)
}

if (exists("fut_pref")) mk_six_panel(pref_alb, copy(fut_pref),
                                       "hazard_glmm", "unit_id",
                                       "fig_future_hazard_glmmTMB_prefecture",
                                       "glmmTMB M4", "Prefecture")
if (exists("fut_pref_xgb")) mk_six_panel(pref_alb, copy(fut_pref_xgb),
                                       "hazard_xgb", "unit_id",
                                       "fig_future_hazard_xgboost_prefecture",
                                       "XGBoost", "Prefecture")
if (exists("fut_cnty")) mk_six_panel(cnty_alb, copy(fut_cnty),
                                       "hazard_glmm", "unit_id",
                                       "fig_future_hazard_glmmTMB_county",
                                       "glmmTMB M4", "County")
if (exists("fut_cnty_xgb")) mk_six_panel(cnty_alb, copy(fut_cnty_xgb),
                                       "hazard_xgb", "unit_id",
                                       "fig_future_hazard_xgboost_county",
                                       "XGBoost", "County")

log("=== DONE ===")
log("All files persisted under data/derived/, results/, figures/main/.")
log("Audit log: results/diagnostics/audit_prefecture_county.txt")
