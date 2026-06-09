# ============================================================
# Script: 48_v3_prefecture_county_refit.R
# Family: v3 risk set extended to prefecture (市) and county (县)
# Author: Chen-Chen Ding + Claude Opus 4.7
# Date  : 2026-05-23
#
# ------------------------------------------------------------
# Scientific question / 科学问题:
#   The v3 province risk set (script 46) includes 463 species and
#   817 events — substantially larger than v1's 333 / 512. Does
#   the same recovery apply at the finer (prefecture / county)
#   administrative scales? That is, when we propagate the v3
#   (species × province) candidate set to all prefectures /
#   counties within each candidate province, and use raw-records
#   effort + WorldClim unit-native climate, does the
#   climate × effort interaction remain significant and positive?
#   把 v3 候选 (species × province) 扩展到该省内的所有市/县，
#   用 raw-records effort + WorldClim 气候，跑出三尺度 v3 结果。
#
# ------------------------------------------------------------
# Reuses (from earlier scripts):
#   - data/derived/sdm_province_v3_relaxed.csv             (script 46)
#   - data/derived/unit_climate_{prefecture,county}.csv    (script 42)
#   - data/derived/unit_effort_{prefecture,county}.csv     (script 42)
#   - GS(2019)1822 市/县/省 shapefiles
#
# Key implementation note: prefecture / county effort tables were
# built from the raw combined-dedup bird records (7.48 M records)
# in script 42; we reuse them rather than re-aggregating.
#
# OOM-safe case-control sampling preserved:
#   prefecture : 1:80  non-events per event (≈ 75 K rows)
#   county     : 1:200 non-events per event (≈ 175 K rows)
# Both well within local RAM headroom.
#
# Outputs (NEW files, do NOT overwrite v2 prefecture/county):
#   data/derived/risk_set_prefecture_v3.csv
#   data/derived/risk_set_county_v3.csv
#   results/tables/table_v3_prefecture_coefs.csv
#   results/tables/table_v3_county_coefs.csv
#   results/tables/table_v3_prefecture_county_aic.csv
#   results/tables/table_v3_three_scale_summary.csv
#   figures/main/Figure_3_v3_multiscale.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(glmmTMB)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})
sf::sf_use_s2(FALSE)
options(warn = 1); set.seed(42)

V2  <- normalizePath(".", mustWork = TRUE)
V1  <- normalizePath(file.path(V2, "..", "bird_hazard_model_effort_upgrade"),
                     mustWork = FALSE)
SHP <- file.path(V2, "data", "spatial", "basemap_GS2019_1822")

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE,
                                                    showWarnings = FALSE)
ens(file.path(V2, "logs"))
LOG <- file(file.path(V2, "logs", "48_v3_prefecture_county.log"),
            "wt", encoding = "UTF-8")
on.exit({ if (isOpen(LOG)) close(LOG) }, add = TRUE)
log <- function(...) {
  m <- paste0(sprintf("[48 %s] ", format(Sys.time(), "%H:%M:%S")),
              paste(..., sep = ""))
  cat(m, "\n", sep = ""); writeLines(m, LOG)
}
audit <- function(t) {
  bar <- paste(rep("─", 60), collapse = "")
  log(""); log(bar); log("AUDIT — ", t); log(bar)
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
# STEP 1 — Load inputs (v3 candidates + admin shp + climate + effort)
# ============================================================
audit("STEP 1: load inputs")

# v3 candidate (species × province)
v3_cand <- fread(file.path(V2, "data", "derived",
                             "sdm_province_v3_relaxed.csv"),
                   encoding = "UTF-8")
v3_cand <- unique(v3_cand[, .(species, province)])
log("v3 candidate: ", nrow(v3_cand), " pairs / ",
    uniqueN(v3_cand$species), " species")

# Events
ev <- fread(file.path(V1, "data", "events_100km_grid_assigned.csv"),
            encoding = "UTF-8")
setnames(ev, tolower(names(ev)))
if (!"year" %in% names(ev) && "pub_year" %in% names(ev))
  setnames(ev, "pub_year", "year")
ev <- ev[year >= 2002 & year <= 2024 & !is.na(longitude) & !is.na(latitude)]
ev$province <- NULL  # use admin-derived province

# Polygons
prov_sf <- st_read(file.path(SHP, "省（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4326) |> st_make_valid()
prov_sf$province <- unname(PROV_CN_EN[as.character(prov_sf[["省"]])])
pref_sf <- st_read(file.path(SHP, "市（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4326) |> st_make_valid()
pref_sf$pref_id <- paste0("PREF_", sprintf("%04d", seq_len(nrow(pref_sf))))
pref_sf$province <- unname(PROV_CN_EN[as.character(pref_sf[["省"]])])
cnty_sf <- st_read(file.path(SHP, "县（等积投影）.shp"), quiet = TRUE) |>
  st_transform(4326) |> st_make_valid()
cnty_sf$cnty_id <- paste0("CNTY_", sprintf("%05d", seq_len(nrow(cnty_sf))))
cnty_sf$province <- unname(PROV_CN_EN[as.character(cnty_sf[["省"]])])
log("polygons: ", nrow(prov_sf), " prov / ", nrow(pref_sf), " pref / ",
    nrow(cnty_sf), " cnty")

# Reuse v2 unit-climate + unit-effort tables (from script 42)
pref_clim <- fread(file.path(V2, "data", "derived",
                              "unit_climate_prefecture.csv"))
cnty_clim <- fread(file.path(V2, "data", "derived",
                              "unit_climate_county.csv"))
pref_eff  <- fread(file.path(V2, "data", "derived",
                              "unit_effort_prefecture.csv"))
cnty_eff  <- fread(file.path(V2, "data", "derived",
                              "unit_effort_county.csv"))
log("unit climate + effort tables loaded (from script 42)")

# Province climate panel for time-varying temp_grad_z
prov_clim <- fread(file.path(V1, "data", "climate_metrics_province_year.csv"),
                    encoding = "UTF-8")
prov_clim_use <- prov_clim[, .(province, year,
                                 temp_grad_z = temp_grad_prov_z)]

# ============================================================
# STEP 2 — Spatially map events into prefecture / county
# ============================================================
audit("STEP 2: spatial-join events → prefecture / county")

ev_sf <- st_as_sf(ev, coords = c("longitude","latitude"), crs = 4326)
suppressWarnings({
  ev_pref <- st_join(ev_sf, pref_sf[, c("pref_id","province")],
                      join = st_within)
  ev_cnty <- st_join(ev_sf, cnty_sf[, c("cnty_id","province")],
                      join = st_within)
})
ev_pref_dt <- as.data.table(st_drop_geometry(ev_pref))
ev_cnty_dt <- as.data.table(st_drop_geometry(ev_cnty))
ev_pref_dt <- ev_pref_dt[!is.na(pref_id)]
ev_cnty_dt <- ev_cnty_dt[!is.na(cnty_id)]
log("events → prefecture: ", nrow(ev_pref_dt), " (",
    uniqueN(ev_pref_dt$pref_id), " prefectures)")
log("events → county    : ", nrow(ev_cnty_dt), " (",
    uniqueN(ev_cnty_dt$cnty_id), " counties)")

# ============================================================
# STEP 3 — Build v3 unit-level risk sets
# ============================================================
audit("STEP 3: build v3 unit risk sets with case-control sampling")

build_v3_unit_risk <- function(scale_label, id_col, sf_unit,
                                ev_unit, control_ratio) {
  log("=== ", scale_label, " ===")
  unit_meta <- unique(as.data.table(st_drop_geometry(sf_unit))[
    , .(unit_id = get(id_col), province)])
  unit_meta <- unit_meta[!is.na(province)]
  log("  unit polygons with province: ", nrow(unit_meta))

  # v3 (species × province) → expand to (species, unit_id)
  cand_unit <- merge(v3_cand, unit_meta, by = "province",
                      allow.cartesian = TRUE)
  log("  v3 candidate (species, ", scale_label, "): ",
      nrow(cand_unit))

  # First-arrival per (species, unit) from spatially-mapped events
  ev_unit_use <- copy(ev_unit)
  setnames(ev_unit_use, id_col, "unit_id")
  # CRUCIAL: restrict events to those species in v3 candidate
  ev_unit_use <- ev_unit_use[species %in% v3_cand$species]
  fa <- ev_unit_use[, .(arrival_year = min(year, na.rm = TRUE)),
                    by = .(species, unit_id, province)]
  log("  first-arrival (species, ", scale_label, ") within v3 species: ",
      nrow(fa))

  # FORCE-INCLUDE rule (same as v3 at province scale):
  # if an event exists, the (species, unit_id) pair is in the candidate
  # set even if SDM didn't predict it. (At unit scale, SDM resolution is
  # province-only, so all unit-level events in candidate provinces are
  # already included via the (species × province) cartesian. The force-
  # include is implicit at this scale.)

  # Events rows
  ev_rows <- copy(fa)
  ev_rows[, year := arrival_year][, event := 1L]

  # Case-control sampling of non-event rows
  n_neg <- nrow(ev_rows) * control_ratio
  log("  sampling ", n_neg, " non-event rows (1:", control_ratio, ")")
  set.seed(42)
  cand_unit[, key := paste0(species, "|", unit_id)]
  fa_lookup <- setNames(fa$arrival_year,
                          paste0(fa$species, "|", fa$unit_id))
  neg_idx <- sample.int(nrow(cand_unit), n_neg, replace = TRUE)
  neg <- cand_unit[neg_idx, .(species, unit_id, province, key)]
  neg[, year := sample(2002:2024, .N, replace = TRUE)]
  neg[, arrival_year := fa_lookup[key]]
  neg <- neg[is.na(arrival_year) | year < arrival_year]
  neg[, event := 0L]
  neg[, drop_key := paste0(species, "|", unit_id, "|", year)]
  ev_rows[, drop_key := paste0(species, "|", unit_id, "|", year)]
  neg <- neg[!drop_key %in% ev_rows$drop_key]

  risk <- rbind(
    ev_rows[, .(species, unit_id, province, year, event)],
    neg[, .(species, unit_id, province, year, event)])
  log("  v3 ", scale_label, " risk set: ", nrow(risk),
      " rows | events = ", sum(risk$event))
  risk
}

risk_pref <- build_v3_unit_risk("prefecture", "pref_id", pref_sf,
                                  ev_pref_dt, control_ratio = 80)
risk_cnty <- build_v3_unit_risk("county",     "cnty_id", cnty_sf,
                                  ev_cnty_dt, control_ratio = 200)

# ============================================================
# STEP 4 — Attach unit-native climate + raw-records effort + province climate
# ============================================================
audit("STEP 4: attach covariates")

attach_unit_cov <- function(risk, unit_clim, unit_eff) {
  # unit climate (WorldClim 10') is time-invariant
  uc <- copy(unit_clim)
  setnames(uc, names(uc)[1], "unit_id")
  risk <- merge(risk, uc, by = "unit_id", all.x = TRUE)

  # raw-records based effort (from script 42)
  ue <- copy(unit_eff)
  setnames(ue, names(ue)[1], "unit_id")
  risk <- merge(risk,
                ue[, .(unit_id, year,
                       log_effort_visits_z, log_effort_days_z,
                       allocated_visits)],
                by = c("unit_id","year"), all.x = TRUE)

  # baseline: province / unit_count broadcast for missing rows
  prov_eff_v1 <- fread(file.path(V1, "data", "effort_panel_upgraded.csv"),
                       encoding = "UTF-8")
  upp <- risk[, uniqueN(unit_id), by = province]
  setnames(upp, "V1", "n_units_in_prov")
  risk <- merge(risk, upp, by = "province", all.x = TRUE)
  risk <- merge(risk,
                prov_eff_v1[, .(province, year, prov_n_visits = n_visits)],
                by = c("province","year"), all.x = TRUE)
  risk[is.na(prov_n_visits), prov_n_visits := 1L]
  risk[is.na(allocated_visits), allocated_visits := 0]
  risk[, baseline_visits := prov_n_visits / pmax(n_units_in_prov, 1)]
  risk[, effort_visits_total := baseline_visits + allocated_visits]
  zify <- function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }
  risk[, log_effort_visits_z := zify(log1p(effort_visits_total)),
       by = year]

  # province × year temp_grad_z
  risk <- merge(risk, prov_clim_use,
                by = c("province","year"), all.x = TRUE)
  risk
}

risk_pref <- attach_unit_cov(risk_pref, pref_clim, pref_eff)
risk_cnty <- attach_unit_cov(risk_cnty, cnty_clim, cnty_eff)
log("attached covariates OK")

fwrite(risk_pref, file.path(V2, "data", "derived",
                              "risk_set_prefecture_v3.csv"))
fwrite(risk_cnty, file.path(V2, "data", "derived",
                              "risk_set_county_v3.csv"))

# ============================================================
# STEP 5 — Fit M0-M4 at each unit scale
# ============================================================
audit("STEP 5: fit M0-M4 at prefecture + county (v3)")

fit_M_ladder <- function(rs, scale_label) {
  d <- rs[!is.na(temp_grad_z) & !is.na(log_effort_visits_z)]
  d[, climate_z := temp_grad_z]
  d[, effort_z  := log_effort_visits_z]
  log("--- ", scale_label, " | rows: ", nrow(d), " | events ",
      sum(d$event))
  forms <- list(
    M0 = "event ~ 1                          + (1|species) + (1|unit_id)",
    M1 = "event ~ effort_z                    + (1|species) + (1|unit_id)",
    M2 = "event ~ climate_z                   + (1|species) + (1|unit_id)",
    M3 = "event ~ climate_z + effort_z        + (1|species) + (1|unit_id)",
    M4 = "event ~ climate_z * effort_z        + (1|species) + (1|unit_id)")
  cf_rows <- list(); ai_rows <- list()
  for (nm in names(forms)) {
    t0 <- Sys.time()
    fit <- tryCatch(
      glmmTMB(as.formula(forms[[nm]]), data = d,
              family = binomial(link = "cloglog")),
      error = function(e) { log("   ", nm, " FAILED: ",
                                  conditionMessage(e)); NULL })
    if (is.null(fit)) next
    s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    log(sprintf("   %s OK (%.1fs, AIC = %.2f, nobs = %d)",
                nm, s, AIC(fit), nobs(fit)))
    cfx <- fixef(fit)$cond; se <- sqrt(diag(stats::vcov(fit)$cond))
    for (tm in names(cfx)) {
      i <- match(tm, names(cfx))
      cf_rows[[length(cf_rows)+1L]] <- data.table(
        scale = scale_label, model = nm, term = tm,
        beta = cfx[i], se = se[i],
        hr = exp(cfx[i]),
        hr.low = exp(cfx[i] - 1.96 * se[i]),
        hr.high = exp(cfx[i] + 1.96 * se[i]),
        p.value = 2 * pnorm(-abs(cfx[i] / se[i])))
    }
    ai_rows[[length(ai_rows)+1L]] <- data.table(
      scale = scale_label, model = nm,
      AIC = AIC(fit), BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)), nobs = nobs(fit))
    invisible(gc(verbose = FALSE))
  }
  list(coefs = rbindlist(cf_rows, fill = TRUE),
       aic   = rbindlist(ai_rows, fill = TRUE))
}
res_p <- fit_M_ladder(risk_pref, "prefecture")
res_c <- fit_M_ladder(risk_cnty, "county")

fwrite(res_p$coefs, file.path(V2, "results", "tables",
                                "table_v3_prefecture_coefs.csv"))
fwrite(res_c$coefs, file.path(V2, "results", "tables",
                                "table_v3_county_coefs.csv"))
all_aic <- rbind(res_p$aic, res_c$aic)
all_aic[, dAIC := AIC - min(AIC), by = scale]
all_aic[, akaike_weight := exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC)),
        by = scale]
fwrite(all_aic, file.path(V2, "results", "tables",
                            "table_v3_prefecture_county_aic.csv"))
log("wrote v3 prefecture + county tables")

# ============================================================
# STEP 6 — Three-scale (province + prefecture + county) summary
# ============================================================
audit("STEP 6: three-scale summary")

# province v3
v3p_coefs <- fread(file.path(V2, "results", "tables",
                              "table_province_v3_coefs.csv"))
v3p_int <- v3p_coefs[model == "M4" & grepl(":", term),
                       .(scale = "province",
                         hr, hr.low, hr.high, p.value)]
v3pre_int <- res_p$coefs[model == "M4" & grepl(":", term),
                          .(scale = "prefecture",
                            hr, hr.low, hr.high, p.value)]
v3cnt_int <- res_c$coefs[model == "M4" & grepl(":", term),
                          .(scale = "county",
                            hr, hr.low, hr.high, p.value)]
three <- rbind(v3p_int, v3pre_int, v3cnt_int)
three[, n_risk_rows := c(188870, nrow(risk_pref), nrow(risk_cnty))]
three[, n_events := c(817, sum(risk_pref$event), sum(risk_cnty$event))]
fwrite(three, file.path(V2, "results", "tables",
                          "table_v3_three_scale_summary.csv"))
print(three)

# ============================================================
# STEP 7 — Publication figure
# ============================================================
audit("STEP 7: build v3 multi-scale figure")

theme_pub <- function(s = 9) {
  theme_bw(base_size = s) +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(linewidth = 0.4, colour = "grey20"),
          plot.title = element_text(face = "bold", size = s + 1))
}
COL_SC <- c(province = "#B40426", prefecture = "#1F77B4", county = "#7F7F7F")
three[, scale := factor(scale, levels = c("county","prefecture","province"))]

p_forest <- ggplot(three, aes(x = hr, y = scale, colour = scale)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  height = 0.18, linewidth = 0.5) +
  geom_point(size = 3.5) +
  geom_text(aes(label = sprintf("HR = %.3f, p = %.0e", hr, p.value)),
             nudge_x = 0.04, size = 2.5, hjust = 0) +
  scale_colour_manual(values = COL_SC, guide = "none") +
  scale_x_continuous(trans = "log",
                      breaks = c(0.95, 1.0, 1.1, 1.2, 1.3, 1.4),
                      limits = c(0.95, 1.7)) +
  labs(tag = "a",
        title = "v3 three-scale climate × effort (M4)",
        subtitle = paste0("Province + prefecture + county; all use v3 candidate set + raw-records effort + WorldClim climate."),
        x = "Hazard ratio (95 % CI, log scale)", y = NULL) +
  theme_pub()

# Comparison v2 vs v3 at each unit scale
v2_pref <- fread(file.path(V2, "results", "tables",
                             "table_prefecture_coefs.csv"))
v2_cnty <- fread(file.path(V2, "results", "tables",
                             "table_county_coefs.csv"))
v2_pref_int <- v2_pref[model == "M4" & grepl(":", term),
                         .(scale = "prefecture", run = "v2",
                           hr, hr.low, hr.high, p.value)]
v2_cnty_int <- v2_cnty[model == "M4" & grepl(":", term),
                         .(scale = "county", run = "v2",
                           hr, hr.low, hr.high, p.value)]
v3_pref_int_dt <- res_p$coefs[model == "M4" & grepl(":", term),
                                .(scale = "prefecture", run = "v3",
                                  hr, hr.low, hr.high, p.value)]
v3_cnty_int_dt <- res_c$coefs[model == "M4" & grepl(":", term),
                                .(scale = "county", run = "v3",
                                  hr, hr.low, hr.high, p.value)]
cmp <- rbindlist(list(v2_pref_int, v3_pref_int_dt,
                        v2_cnty_int, v3_cnty_int_dt))
p_cmp <- ggplot(cmp, aes(x = hr, y = run, colour = run)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = hr.low, xmax = hr.high),
                  height = 0.15, linewidth = 0.45,
                  position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  facet_wrap(~ scale, scales = "free_x") +
  scale_colour_manual(values = c(v2 = "#7F7F7F", v3 = "#B40426"),
                       name = "Run") +
  scale_x_continuous(trans = "log",
                      breaks = c(0.9, 1.0, 1.1, 1.2, 1.3)) +
  labs(tag = "b",
        title = "v2 vs v3 at prefecture + county",
        x = "Hazard ratio (log)", y = NULL) +
  theme_pub() + theme(legend.position = "top")

fig <- (p_forest / p_cmp)
fig <- fig + plot_annotation(
  title = "Figure 3 v3 — Multi-scale validation with the relaxed risk set",
  subtitle = "v3 propagates the relaxed (species × province) candidate set to all units within each candidate province.",
  theme = theme(plot.title = element_text(face = "bold", size = 10),
                 plot.subtitle = element_text(size = 8.5, colour = "grey30")))

ens(file.path(V2, "figures", "main"))
ggsave(file.path(V2, "figures", "main", "Figure_3_v3_multiscale.pdf"),
       fig, width = 18, height = 14, units = "cm",
       device = grDevices::cairo_pdf)
ggsave(file.path(V2, "figures", "main", "Figure_3_v3_multiscale.png"),
       fig, width = 18, height = 14, units = "cm", dpi = 600)
log("wrote Figure_3_v3_multiscale.{pdf,png}")

log("")
log("══════════════════════════════════════════════════════════")
log("        v3 PREFECTURE + COUNTY REFIT COMPLETE")
log("══════════════════════════════════════════════════════════")
print(three)
