# METHODS — v2

This document mirrors the Methods section of `manuscript/manuscript_v2.Rmd`
and serves as the canonical reference for parameter choices, software
versions and pipeline dependencies. Sections numbered to match the
manuscript.

## 2.1 Data sources / 数据来源

| Source | Version / vintage | Use |
|--------|--------------------|-----|
| Bird new-record dataset | Compiled May 2026 (`鸟类新纪录20260509.xlsx`) | Event records (species × province × year) |
| CHELSA climate | v2.1 (Karger *et al.* 2017, updated 2023) | Monthly tmean / prec rasters, 1981–2010 baseline + 2002–2024 observed |
| WorldClim 2.5 arc-min | v2.1 | Cross-validation against CHELSA at coarser resolution |
| GBIF birding effort proxy | accessed 2025-12 | Backstop for effort panel construction |
| CMIP6 (ESGF) | ACCESS-CM2, EC-Earth3, MPI-ESM1-2-HR, MIROC6, UKESM1-0-LL | SSP245 / SSP585 climate projections (2050, 2080) |
| Administrative boundaries | Mapping Audit GS(2019)1822 | Province / prefecture / county boundaries (`data/spatial/basemap_GS2019_1822/`) |
| Survey effort panel | Custom occupancy-aggregated panel of e-records, 2002–2024 | Person-visits (n_visits), n_observers, n_birding_days, log person-hours |

All rasters are aggregated to grids via area-weighted means; province
metrics are area-weighted means of within-province grid cells.

## 2.2 Spatial units / 空间单元

- **32 mainland provinces** (excludes Hong Kong, Macau, Taiwan).
- **Grids** at 50 km and 100 km nominal resolution, constructed in
  Albers Equal-Area (**EPSG:4524**), then re-projected to WGS84 for
  attribute joins. Grid cells with <50 % land are dropped.
- Raw point coordinates are WGS84 (**EPSG:4326**).
- All area calculations use EPSG:4524.

## 2.3 Effort metric construction / Effort 指标构建

Four effort specifications:

- **Spec A** – `log1p(n_records)_z` (legacy, retained for back-compatibility).
- **Spec B** – `log1p(n_visits)_z` — primary headline metric.
- **Spec C** – `effort_pc1_z` (PCA of 4 raw effort variables).
- **Spec D** – `log1p(n_birding_days)_z`.

For the offset model **M5** the **raw** scale is used:
`offset = log(person_hours + 1)`. v1's use of a z-scored offset is
explicitly retired (see `GAP_RESOLUTION.md` P0-2).

## 2.4 Hazard model family / 危险率模型层次

Discrete-time complementary log-log models fitted with
`glmmTMB::glmmTMB(family = binomial(link="cloglog"))`.

- **M0** Null: `event ~ 1 + (1|species) + (1|province)`
- **M1** Effort additive: `event ~ effort + (1|species) + (1|province)`
- **M2** Climate additive: `event ~ climate + (1|species) + (1|province)`
- **M3** Additive joint: `event ~ effort + climate + (1|species) + (1|province)`
- **M4** Interaction: `event ~ effort * climate + (1|species) + (1|province)`
- **M5** Offset: `event ~ climate + (1|species) + (1|province) + offset(log_person_hours)`

Grid-scale variants replace `(1|province)` with `(1|grid_id)` and use
**grid-native** climate (recomputed in `code/28_grid_native_climate.R`).

## 2.5 Multi-scale and MAUP design

Models are refitted at **province**, **50 km**, **100 km**, and (for
elasticity tests) **200 km** grids. Coefficient elasticity is reported
in `code/31_maup_sensitivity.R` → `results/sensitivity/`. The interaction
β is required to retain sign and significance at all four scales for
the headline conclusion to stand.

## 2.6 Cross-validation and spatial diagnostics

- **Spatial block CV**: 5 folds, 250 km blocks, `blockCV::cv_spatial`
  with stratified random allocation; reported metrics are mean ± SD
  across folds (`code/26_spatial_block_cv.R`).
- **Moran's I**: Pearson residuals of M4 binned at 50 / 100 / 250 / 500 km
  distance classes; complemented by `DHARMa::testSpatialAutocorrelation`
  (`code/27_morans_i_diagnostics.R`).
- **DHARMa global tests**: uniformity, dispersion, outlier, ZI.

## 2.7 Forecasting and uncertainty / 预测与不确定性

- **CMIP6 ensemble**: 5 GCMs × {SSP245, SSP585} × {2050, 2080}; per-cell
  predictions are summarised as median + IQR; **ensemble disagreement**
  is the cell-wise IQR / median.
- **Forecast skill decay**: XGBoost trained on 1980–2010 and validated
  on rolling 2011–2024 windows; skill curve = AUC vs. lead time;
  per-feature Population Stability Index (PSI) flags covariate shift.

## 2.8 Reproducibility / 可重复性

- All package versions pinned in `renv.lock`.
- DAG orchestrated by `_targets.R` (~60 nodes); `tar_visnetwork()` PNG
  saved to `figures/diagnostics/dag.png`.
- `code/37_reproducibility_selfcheck.R` asserts manuscript numbers match
  the regenerated tables to tolerance (records exact; AIC ± 0.5;
  marginal R² ± 0.001; AUC ± 0.005).
- `sessionInfo.txt` is rewritten at the end of `00_run_all.R`.

See `manuscript/supplementary/supp_methods.Rmd` for full equations and
the CMIP6 ESGF DOI manifest.
