# Stage 6 — Paper Creation Process Record

**Project:** bird_hazard_model_effort_upgrade_v2
**Pipeline:** academic-pipeline (ARS v3.8.2), `verified_only`
**Run date:** 2026-06-06
**Target tier:** *Nature Ecology & Evolution* / *Ecology Letters*

## 1. Pipeline journey

| Stage | Action | Outcome |
|---|---|---|
| Entry detection | Existing draft + persisted results | Entered at Stage 2.5 |
| Figure optimisation | Re-ran `code/53_v4_polished_figures.R` (bypassed empty renv via `--no-init-file`, used system R 4.5.1 library) | 8 v4 figures regenerated 2026-06-06 09:28–09:30 |
| 2.5 INTEGRITY | 5-phase check + 7-mode AI-failure checklist | **PASS**; 2 issues fixed (removed unsourced 80.5% R²; corrected Spec D ΔAIC 24.8→27.7) |
| 3 REVIEW | 5-panel simulated review (EIC+R1+R2+R3+DA) | **Major revision**, ≈73/100; roadmap with P0/P1/P2 |
| 4 REVISE | Applied all revisable-now items (framing + persisted data) | Effect-size translation; climate-metric breadth; SDM-circularity Methods; administrative-scale relabel; forecast reframe; endogeneity promoted; generality subsection |
| 4.5 FINAL INTEGRITY | Re-verified from scratch | **PASS**; no new unsupported number (effect-size 3.4%→4.4% derived from M0 intercept −3.364) |
| 5 FINALIZE | pandoc MD→DOCX | 3 DOCX produced |
| 6 SUMMARY | This record + cover letter | Complete |

## 2. Deliverables (all in `output/doc/`)

**Submission package**
- `manuscript_v4_integrated_NEE_EcolLett.md` / `.docx` — revised top-tier paper
- `cover_letter_v4.md` — cover letter (NEE/Ecol Lett)
- `中文_过程与结果_系统总结.md` / `.docx` — Chinese process+results synthesis

**Pipeline audit trail**
- `pipeline_stage2.5_integrity_report.md` — integrity gate
- `pipeline_stage3_peer_review.md` — 5-panel review + roadmap
- `pipeline_stage4_response_to_reviewers.md` — point-by-point response
- `pipeline_review_package.docx` — bundled audit trail
- `pipeline_stage6_process_record.md` — this file

**Refreshed figures (`figures/main/`, 2026-06-06)**
- `Figure_2_v4_province_headline_M5_raincloud.{pdf,png}`
- `Figure_3_v4_three_scale_forest_M5.{pdf,png}`
- `Figure_5_v4_province_future_glmmTMB.{pdf,png}` / `..._xgboost`
- `Figure_6_v4_prefecture_future_{glmmTMB,xgboost}.{pdf,png}`
- `Figure_7_v4_county_future_{glmmTMB,xgboost}.{pdf,png}`

## 3. Collaboration quality evaluation (6 dimensions, 0–100, honest)

| Dimension | Score | Evidence |
|---|---:|---|
| Scientific clarity of goal | 88 | User specified taxon, results, tier, honesty boundary across two turns |
| Evidence discipline | 95 | `verified_only` enforced; unsourced 80.5% R² removed; no fabricated numbers |
| Critical depth of review | 82 | Endogeneity, SDM circularity, frozen-effort, grid-vs-admin all surfaced (not sycophantic) |
| Revision faithfulness | 84 | 13/13 review points tracked; new-modelling items honestly deferred not faked |
| Reproducibility | 90 | Every claim → named table; figures regenerated from persisted CSVs |
| Output completeness | 86 | Manuscript+cover+Chinese synthesis+figures+DOCX delivered; PDF blocked by missing LaTeX engine |

**Overall: 87/100.** Main deductions: three P0 scientific gaps (endogeneity
identification, ecological-grid grain, spatial-block CV) remain as roadmap, so the
paper is "strong submittable draft pending two analyses," not "ready to submit
as-is."

## 4. AI self-reflection

- **What went right:** strict provenance prevented the most common AI failure
  (hallucinated/over-stated results); the review genuinely stress-tested the work
  and the revision improved honesty rather than inflating claims.
- **What was deferred, not solved:** the endogeneity (effort-lag) refit,
  grid-native 50/100 km refit, and spatial-block CV are the real gates to this
  tier. They need the R pipeline run with the project's renv restored; I declined
  to fabricate their outputs.
- **Highest-value next action:** run the lagged-effort M4 refit and the 250 km
  spatial-block CV (scripts `26`, plus a lag variant), then re-enter the pipeline
  at Stage 4 to convert two P0 ROADMAP items into reported results.

## 4b. Round-2 update — P0 analyses run + editable outputs (same session)

Two of the three P0 gating analyses were executed (R, system library via
`--no-init-file`) and converted from ROADMAP into reported results:

- **Endogeneity (effort lag, `code/54b`):** province-year-level one-year effort
  lag, n = 12,535, 498 events. Interaction PERSISTED and slightly strengthened —
  Spec B HR = 1.322 (1.191–1.467), p = 1.6 × 10⁻⁷; all four specs significant
  (`table_effort_lag_refit.csv`). Now manuscript Table 3; abstract + Discussion
  updated; "reverse causation" downgraded from top threat to addressed.
  (Note: an initial species-level lag was underpowered — 960 rows — and was
  discarded once effort was confirmed to be province-year-level; the corrected
  province-year lag is the reported one.)
- **Spatial-block CV (`code/54`):** blockCV 250 km, 5 folds. Mean AUC 0.55 ± 0.05
  vs random 0.65 ± 0.01 (`table_spatial_block_cv.csv`). Reported candidly:
  random CV overstates skill; framework is for inference, not out-of-region
  forecasting. Now in Diagnostics + abstract + Limitations.

**Final-integrity re-check (4.5, round 2): PASS** — all new numbers trace to the
two new tables; manuscript DOCX regenerated.

**New editable deliverables (per user request):**
- `bird_hazard_all_statistics.xlsx` — 22 sheets, every coefficient/HR/CI/p/AIC/
  weight/diagnostic/forecast parameter + `headline_interaction` summary.
- `bird_hazard_results_editable.pptx` — 20 slides; **9 native editable charts**
  (column/bar/scatter, double-click to edit) + 10 embedded high-res images.
- `figures/main/Figure_R1_interaction_raincloud.{png,pdf,svg}` — raincloud.
- `figures/main/Figure_R2_interaction_beeswarm.{png,pdf,svg}` — beeswarm.

**Remaining single P0:** grid-native 50/100 km refit (ecological grain).

## 5. Pre-submission checklist

- [ ] Run lagged-effort M4 refit (endogeneity) → report or keep as limitation
- [ ] Run spatial-block CV (script 26) → replace random-CV framing
- [ ] Persist grid-native 50/100 km refit (scripts 28/28b/40) → ecological grain
- [ ] `citation-check` all DOIs; add a China citizen-science effort-growth citation
- [ ] Decide venue → adjust abstract structure & length; convert citations to venue style
- [ ] Remove internal "(remove before submission)" notes from manuscript DOCX
- [ ] Mint Zenodo DOI; finalise data/code availability
