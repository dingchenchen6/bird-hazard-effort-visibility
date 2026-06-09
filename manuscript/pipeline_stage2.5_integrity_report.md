# Stage 2.5 — Integrity Verification Report

**Manuscript under review:** `output/doc/manuscript_v4_integrated_NEE_EcolLett.md`
**Data access level:** `verified_only` (only persisted `results/` tables)
**Date:** 2026-06-06
**Gate rule:** must PASS (zero unresolved fabrications) to proceed to Stage 3.

---

## Phase A — Reference verification

| Reference | Real work? | DOI status |
|---|---|---|
| Boakes et al. 2010 *PLoS Biol* | ✅ | verify at citation-check |
| Brooks et al. 2017 *R Journal* (glmmTMB) | ✅ | verify |
| Chen et al. 2011 *Science* | ✅ | verify |
| Eyring et al. 2016 *GMD* (CMIP6) | ✅ | verify |
| Fick & Hijmans 2017 *IJC* (WorldClim2) | ✅ | verify |
| Hartig 2022 DHARMa | ✅ (R pkg) | n/a |
| Isaac et al. 2014 *MEE* | ✅ | verify |
| Karger et al. 2017 *Sci Data* (CHELSA) | ✅ | verify |
| Lenoir et al. 2020 *Nat Ecol Evol* | ✅ | verify |
| Loarie et al. 2009 *Nature* | ✅ | verify |
| Openshaw 1984 (MAUP) | ✅ | n/a (book) |
| Sullivan et al. 2014 *Biol Conserv* (eBird) | ✅ | verify |
| Valavi et al. 2019 *MEE* (blockCV) | ✅ | verify |

**Verdict:** No fabricated/hallucinated references. All 13 are well-established works.
**Action:** Run `citation-check` to confirm each DOI before submission (non-blocking
for this gate).

## Phase B — Citation-context verification

Every in-text citation is used in a context the cited work actually supports
(range shifts → Chen 2011 / Lenoir 2020; occurrence bias → Boakes 2010 / Isaac
2014 / Sullivan 2014; climate velocity → Loarie 2009; methods → Brooks 2017,
Valavi 2019, Openshaw 1984; climate data → Fick & Hijmans 2017, Karger 2017,
Eyring 2016). No citation is stretched beyond its source. **PASS.**

## Phase C — Statistical-data verification (claim → persisted table)

| # | Manuscript claim | Value | Source table | Status |
|---|---|---|---|---|
| C1 | Conservative risk set | 12,813 rows / 333 sp / 512 ev | `table_riskset_v3_attrition.csv`, `table_province_v2_aic.csv` (nobs) | ✅ |
| C2 | Province M4 Spec B interaction | HR 1.288 (1.179–1.407), p=2.1×10⁻⁸ | `table_province_v2_coefs.csv` | ✅ |
| C3 | Four-spec HR/CI/p (A/B/C/D) | 1.213 / 1.288 / 1.185 / 1.300 | `table_province_v2_coefs.csv` | ✅ |
| C4 | ΔAIC(M3→M4) A/B/C/D | 17.0 / 29.0 / 15.7 / 27.7 | `table_province_v2_aic.csv` | ✅ (D corrected 24.8→27.7) |
| C5 | M4 Akaike weight ≈ 1.0 | 0.99–1.00 | `table_aic_akaike_weights.csv` | ✅ |
| C6 | M4 beats M5 in 7/8 | ΔAIC signs | `table_m5_offset_summary.csv` | ✅ |
| C7 | Prefecture interaction | 1.163 (1.046–1.293), p=5.2×10⁻³ | `table_prefecture_coefs.csv` | ✅ |
| C8 | County interaction | 1.114 (0.998–1.243), p=5.5×10⁻² | `table_county_coefs.csv` | ✅ |
| C9 | Prefecture/county M4 weight | 0.91 / 0.56 | `table_prefecture_county_aic.csv` | ✅ |
| C10 | Relaxed set size | 188,870 / 463 / 817 | `table_riskset_v3_attrition.csv` | ✅ |
| C11 | v3 four-spec HR (A/B/C/D) | 1.338 / 1.279 / 1.124 / 1.434 | `table_province_v3_all_specs_coefs.csv` | ✅ |
| C12 | v3 three-scale headline B | 1.274 (1.198–1.354), p=1.1×10⁻¹⁴ | `table_v3_three_scale_summary.csv` | ✅ |
| C13 | v3 prefecture/county | 1.043 (p=0.26) / 1.021 (p=0.58) | `table_v3_three_scale_summary.csv` | ✅ |
| C14 | v1/v2/v3 reconciliation | 1.292 / 1.288 / 1.274 | `table_province_v1_v2_v3_reconciliation.csv` | ✅ |
| C15 | Moran's I 100/250/500 km | 1.4×10⁻⁵ / 2.1×10⁻⁴ / −6.6×10⁻⁵ (all p>0.89) | `table_morans_i_residuals.csv` | ✅ |
| C16 | XGBoost relaxed CV-AUC | 0.777 (5-fold, 37 rounds) | `table_xgb_cv_v3.csv` | ✅ |
| C17 | RF temp×effort rank 1→6 | Δrank +5 | `table_v2_v3_rf_comparison.csv` | ✅ |
| C18 | Future SSP585 mean 2050/2080 | 0.308 / 0.760 | `table_province_future_glmmTMB.csv` (mainland, recomputed) | ✅ |
| C19 | Future SSP245 mean 2050/2080 | 0.111 / 0.230 | same | ✅ |
| C20 | Top SSP585-2050 provinces | Jiangsu .76, Zhejiang .73, Hubei .61, Fujian .49, Henan .48 | same | ✅ |

**Verdict:** 20/20 quantitative claims trace to a named persisted table. **PASS.**

### Flagged & resolved during this gate

- **F-1 (resolved):** Earlier drafts (`manuscript_v2.Rmd`) claimed the interaction
  explains "80.5 % of marginal R²". **No persisted R² table exists** in `results/`.
  → The new manuscript **removed** this number and substitutes AIC/Akaike-weight +
  the M5-offset test. No unsupported quantitative claim remains.
- **F-2 (resolved):** Spec D ΔAIC(M3→M4) mis-typed as 24.8 (that is the M0 dAIC).
  → Corrected to 27.7 in both English manuscript and Chinese summary.

## Phase D — Originality

Original secondary analysis and synthesis of the author's own project outputs.
No text reuse from external sources; methods describe the author's own scripts
(`code/26–53`). **PASS.**

## Phase E — Claims verification

Every interpretive claim is bounded by its evidence: the moderation claim rests
on M4-vs-M5 (not over-stated as causal); the future-hazard surfaces are framed
as scenario-conditioned and explicitly not as colonisation estimates; the
robustness-boundary claim is supported by both the v3 fine-scale collapse and the
RF rank shift. No claim exceeds the persisted evidence. **PASS.**

---

## AI Research Failure Mode Checklist (7 modes)

| Mode | Assessment | Status |
|---|---|---|
| 1. Citation hallucination | All refs real; DOIs to verify at citation-check | **PASS** (verify DOIs) |
| 2. Implementation bug-as-result | No new modelling performed; report summarises pre-existing persisted outputs | **N/A** (no new code executed for claims) |
| 3. Hallucinated results | All 20 numeric claims trace to tables (Phase C) | **PASS** |
| 4. Shortcut reliance | No proxy metrics substituted for the reported quantities | **PASS** |
| 5. Bug-as-insight | Interpretations derive from converged, persisted models | **PASS** |
| 6. Methodology fabrication | Methods map to real scripts; unpersisted analyses explicitly demoted to "future work" with no numbers | **PASS** |
| 7. Pipeline frame-lock | Honesty boundary (dropping planned analyses, removing 80.5%) actively guards against frame-lock | **PASS** |

No mode is `SUSPECTED`; no Mode 1/3/5/6 is `INSUFFICIENT EVIDENCE`. **Gate not blocked.**

---

## STAGE 2.5 VERDICT: ✅ PASS

The manuscript contains no fabricated references, no unsupported numbers, and no
claim exceeding the persisted evidence. Two issues were found and resolved during
the gate (80.5% R² removed; Spec D ΔAIC corrected). Cleared to proceed to
**Stage 3 — Peer Review.**

> Non-blocking carry-forward: (1) verify all DOIs via `citation-check`; (2) if the
> "interaction dominates fixed-effect variance" framing is wanted back, persist an
> R² table first.
