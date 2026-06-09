# Effort-moderated visibility of new bird distribution records across China

**A multi-scale discrete-time hazard framework**

This repository contains the analysis code, result tables, figures and manuscript
for a study testing whether **climate exposure and survey effort interact —
rather than add — to drive the hazard of new bird distribution records** across
mainland China (2002–2024), and how that interaction behaves across spatial scales
and under CMIP6 future climate.

> **Headline:** New distribution records are best read as *effort-moderated
> visibility events*. At province scale, climate exposure and survey effort combine
> through a positive interaction (the climate signal becomes recordable only where
> observation effort is sufficient); survey effort is the universal driver of
> detection at every scale.

---

## 1. Scientific question

A new species record is generated only when *ecological exposure* (a species
reaching/persisting in a new place under changing climate) and *observation effort*
(people looking, in the right place and season) coincide. Most biodiversity models
treat effort as an additive nuisance or an offset. We test the alternative: that
effort **moderates** the climatic visibility of redistribution, so climate and
effort enter the record-generating process **multiplicatively**.

Four questions:
1. Do climate exposure and survey effort combine interactively, not additively?
2. Is effort a *moderator* of the climate–hazard slope or a *scaling* factor on
   records? (operationalised as an interaction model vs a raw-effort offset model)
3. Is the answer robust to effort metric, risk-set definition, and spatial scale
   (province → prefecture → county → ecological grid)?
4. What does an effort-moderated reading imply for forecasting and monitoring?

---

## 2. Headline results (all traceable to `results/`)

| Result | Value | Source table |
|---|---|---|
| Province climate × effort interaction (Spec B, observer visits) | **HR = 1.288** (1.179–1.407), p = 2.1×10⁻⁸ | `results/tables/table_province_v2_coefs.csv` |
| Interaction across 4 effort metrics | HR 1.18–1.30, all p ≤ 6×10⁻⁵ | same |
| Interaction (M4) vs offset (M5): M4 preferred | **7 of 8** spec×dataset comparisons | `results/tables/table_m5_offset_summary.csv` |
| Endogeneity check: effort **lagged 1 year** still moderates | **HR = 1.322** (1.191–1.467), p = 1.6×10⁻⁷ | `results/tables/table_effort_lag_refit.csv` |
| Multi-scale (conservative refit): attenuation w/o reversal | province 1.288 → prefecture 1.163 → county 1.114 | `results/tables/table_prefecture_county_aic.csv` |
| Relaxed risk set (×15 rows, +60% events): province preserved | **HR = 1.274**, p = 1.1×10⁻¹⁴ | `results/tables/table_v3_three_scale_summary.csv` |
| Robustness boundary (relaxed, fine scale dilutes) | prefecture 1.043 (p=0.26), county 1.021 (p=0.58) | same |
| Residual spatial autocorrelation (Moran's I) | negligible at 100/250/500 km | `results/diagnostics/table_morans_i_residuals.csv` |
| Spatial-block CV (250 km) vs random | **AUC 0.55** ± 0.05 vs 0.65 ± 0.01 (transfers weakly) | `results/tables/table_spatial_block_cv.csv` |
| Three-run reconciliation (v1/v2/v3) | HR 1.292 / 1.288 / 1.274 | `results/tables/table_province_v1_v2_v3_reconciliation.csv` |
| Grid-native plug-in (100 km): CRU climate × merged effort | downscaled prediction surface | `results/forecasts/table_grid_100km_plugin_cmip6.csv` |
| CMIP6 ensemble warming (4 GCMs) | SSP585 +3.3 °C (2050), +4.8 °C (2080) | `results/diagnostics/table_cmip6_ensemble_delta_summary.csv` |
| Trait stratification: interaction by migratory strategy | **Resident HR = 1.43**, Migrant ≈ 1.25 (all significant) | `results/tables/table_migratory_stratified_interaction.csv` |
| Predictive validation (rich XGBoost) | interpolation 0.73 · temporal forecast 0.63 · spatial 0.56 | `results/tables/table_prediction_accuracy_comparison.csv` |
| **Relaxed-set parallel analysis** (188,870 rows): endogeneity | lagged-effort HR = **1.292**, p = 9.7×10⁻¹¹ | `results/tables/table_v3_endogeneity_lag.csv` |
| Relaxed-set migratory: long-distance migrants strongest | Long-distance **1.33**, Resident 1.28 (all significant) | `results/tables/table_v3_migratory_stratified.csv` |
| Relaxed-set prediction | interpolation 0.76 · spatial 0.62 · temporal 0.59 | `results/tables/table_v3_prediction_accuracy.csv` |

The province-scale interaction is robust across three risk sets, an offset test and
a one-year-lag endogeneity test. It **attenuates with finer spatial grain** and is
not detectable at the ecological grid scale on its own; effort, by contrast,
predicts records strongly at every scale.

---

## 3. Repository structure

```
bird-hazard-effort-visibility/
├── README.md                  # this file
├── METHODS.md                 # canonical methods reference
├── code/                      # all analysis scripts (R + Python)
│   ├── 26–53                  # spatial CV, Moran's I, grid-native climate/effort,
│   │                          #   offset, multi-scale, v3 relaxed, publication figures
│   ├── 54  54b                # spatial-block CV + endogeneity (effort-lag) refit
│   ├── 55 56 57               # stats workbook (xlsx), beeswarm/raincloud, editable PPTX
│   ├── 59 60 63               # corrected grid climate/effort, grid plug-in (100/50 km)
│   ├── 64 66                  # CMIP6 ensemble futures, polished plug-in maps (polygons,
│   │                          #   cool→warm palette, nine-dash-line basemap, effort scenarios)
│   └── utils/                 # shared helpers
├── results/
│   ├── tables/                # model coefficients, AIC, HR/CI/p (all scales, all specs)
│   ├── forecasts/             # province + grid future-hazard tables
│   └── diagnostics/           # Moran's I, spatial CV, CMIP6 deltas, attrition
├── figures/                   # publication PNGs (+ small vector PDFs)
├── manuscript/                # integrated manuscript (md+docx), cover letter, bilingual
│                              #   summary, methodology memo, pipeline audit, stats xlsx, PPTX
├── data_dictionary/           # variable schemas (no raw data)
├── _targets.R, renv.lock      # reproducibility (targets DAG + package lockfile)
└── .gitignore
```

**Data are not included** (licensed occurrence records, large climate rasters,
multi-GB risk sets). See `data_dictionary/` for schemas and §6 for sources.

---

## 4. Methods in brief (see `METHODS.md`)

- **Response & risk set.** Discrete-time survival: rows are species × spatial-unit ×
  year at risk of a *first* record; the event is the first arrival. Headline
  (conservative) province risk set: 12,813 rows / 333 species / 512 events; a relaxed
  event-override set (188,870 / 463 / 817) stress-tests it.
- **Model family.** `glmmTMB` complementary-log-log hazard with crossed
  `(1|species)+(1|unit)` random effects: M0 null, M1 effort, M2 climate, M3 additive,
  **M4 climate × effort interaction**, M5 raw-effort offset. Inferential contrasts:
  M4 vs M3 (interaction vs addition) and M4 vs M5 (moderation vs scaling).
- **Effort.** Four metrics (records, observer visits [headline], PCA composite,
  birding-days), within-year standardised; offset on the raw `log(person-hours+1)`
  scale.
- **Climate.** Temperature anomaly/gradient, climate velocity, Mahalanobis
  displacement (WorldClim 2.1, CHELSA v2.1); grid-native time-varying temperature
  from **CRU TS 4.09** monthly data (2002–2024).
- **Multi-scale.** Province fit; prefecture/county/grid use the same province risk
  set and are **plug-in predictions** with unit-native climate and effort
  (new records are recorded at province resolution, so the model is fitted there).
- **Diagnostics.** Residual Moran's I; **spatial-block cross-validation** (250 km
  blocks, blockCV); one-year **effort-lag** endogeneity refit; random-forest /
  XGBoost importance.
- **Futures.** Survey effort from the merged **eBird/GBIF + China-Birdwatch** panel;
  future climate from a **4-GCM CMIP6 ensemble** (ACCESS-CM2, MPI-ESM1-2-HR, MIROC6,
  UKESM1-0-LL) under SSP245/585; **SSP-differentiated effort-growth scenarios**.
- **Maps** use the official GS(2019)1822 basemap and render the national boundary and
  the South China Sea nine-dash line.

---

## 5. Figures (in `figures/`)

| Figure | Content |
|---|---|
| `Figure_1_concept_and_workflow` | Study domain, sample structure, climate × effort × visibility logic |
| `Figure_2_v4_province_headline_M5_raincloud` | Province headline: 4-spec forest, AIC ladder, M4-vs-M5, raincloud |
| `Figure_3_v4_three_scale_forest_M5` | Province → prefecture → county interaction forest |
| `Figure_4_variable_importance_v3` + `Figure_v3_robustness_panel` | Relaxed-set robustness boundary |
| `Figure_R1_interaction_raincloud`, `Figure_R2_interaction_beeswarm` | Interaction HR across specs × risk sets (incl. lag) |
| `Figure_5_v4_province_future_*` | Province scenario future hazard |
| `Figure_6_prefecture_plugin_unified`, `Figure_7_county_plugin_unified` | Prefecture/county polygon plug-in maps (cool→warm, nine-dash line) |
| `Figure_8_grid100_native_plugin_hazard`, `Figure_8b_grid50_native_plugin_hazard` | 100/50 km grid plug-in: CRU climate × merged effort, CMIP6 + SSP-effort futures |

---

## 6. Reproducibility & data access

- **Software:** R 4.4+ (`glmmTMB`, `xgboost`, `blockCV`, `DHARMa`, `terra`, `sf`,
  `targets`, `arrow`); Python 3.9+ (`pandas`, `matplotlib`, `python-pptx`,
  `openpyxl`) for the stats workbook, beeswarm/raincloud and editable PPTX.
  Package versions pinned in `renv.lock`.
- **Run order (key scripts):** province + multi-scale models → `54`/`54b` (spatial
  CV, endogeneity) → `59`/`60`/`63` (grid climate/effort, grid plug-in) →
  `64` (CMIP6 futures) → `66` (polished maps) → `55`/`56`/`57` (workbook, figures,
  PPTX).
- **Data sources (not redistributed here):** bird new-record dataset; WorldClim 2.1
  & CHELSA v2.1; CRU TS 4.09; CMIP6 (WorldClim downscaled, SSP245/585); merged
  eBird/GBIF + China Bird Watching Record Center survey effort; GS(2019)1822 basemap
  (Ministry of Natural Resources of China). Licensed/sensitive occurrence data are
  shared under source-specific terms; schemas are in `data_dictionary/`.

---

## 7. Manuscript

A submission-ready **Nature Ecology & Evolution** draft is in
`manuscript/manuscript_NEE_final.{md,docx}` (≤150-word abstract, progressive
introduction with explicit hypotheses H1–H4, Results with subheadings, deep
Discussion, Methods-at-end, numbered citations). The fuller working manuscript (target: *Nature Ecology & Evolution* / *Ecology Letters*)
and a bilingual process-and-results synthesis are in `manuscript/`
(`manuscript_v4_integrated_NEE_EcolLett.{md,docx}`,
`中文_过程与结果_系统总结.{md,docx}`), with a cover letter, methodology memo,
peer-review audit trail, a full statistics workbook
(`bird_hazard_all_statistics.xlsx`) and an editable results deck
(`bird_hazard_results_editable.pptx`).

## 8. Citation

> Ding, C.-C. *Survey effort moderates the climatic visibility of new bird
> distribution records across China: a multi-scale hazard framework.* (in prep).

See `CITATION.cff`. Code: MIT (`LICENSE`). Manuscript & figures: CC-BY 4.0 on
acceptance.

## Author

Chen-Chen Ding · Institute of Ecology, Peking University · chenchending1992@gmail.com
