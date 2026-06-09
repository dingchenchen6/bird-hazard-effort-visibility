# Stage 4 — Response to Reviewers (revision round 1)

Manuscript: *Survey effort moderates the climatic visibility of new bird
distribution records across China*. Every reviewer point is accounted for with an
explicit status: **DONE** (applied in this revision using persisted data +
framing), **ROADMAP** (requires new modelling; added as a prominent limitation /
future test, no fabricated numbers), or **DISAGREE** (with justification).

| # | Reviewer point | Status | What changed |
|---|---|---|---|
| R1-1 | Ecological vs administrative scale; grids needed | **DONE (reframe) + ROADMAP** | Multi-scale section + heading relabelled "administrative scales"; explicit statement that ecological-grid (50/100 km) is the decisive future grain test; elevated to top-priority limitation. New modelling itself = ROADMAP. |
| R1-2 | Which climate signal? report multiple metrics | **DONE (partial)** | New paragraph: headline term is time-varying thermal anomaly/gradient; RF importance shows temp-based terms dominate and temp×effort ranks #1; full 6-metric interaction matrix flagged method-ready (not persisted → no numbers). |
| R1-3 | Effect size in biological terms | **DONE** | Added absolute translation: null annual hazard ≈3.4%; +1 SD climate×effort → ≈4.4% (HR 1.29), with decade/at-risk-pool framing. |
| R2-1 | Endogeneity / reverse causation | **DONE (reframe) + ROADMAP** | Promoted to first limitation; three named defenses (1-yr effort lag, leave-future-out, instrument); causal language hedged to "strongly consistent with". Lag refit = ROADMAP. |
| R2-2 | Spatial-block CV (random CV overstates skill) | **DONE (honesty) + ROADMAP** | Stated Moran's I is a residual not a generalisation test; spatial-block CV named as top-2 pre-publication addition; no skill number claimed beyond persisted random-CV/Moran's I. |
| R2-3 | Shared slope across species | **ROADMAP** | Added as limitation; random-slope / migratory stratification flagged as extension. |
| R2-4 | Model-selection inference / averaging | **ROADMAP (P2)** | Noted; deferred. |
| R3-1 | Frozen-effort forecasts inconsistent | **DONE** | Forecast section rewritten to foreground exposure×effort-gap priority surface; absolute maps demoted to scenario context; effort-frozen-at-2024 caveat explicit. |
| R3-2 | Vagrants vs colonisers | **DONE (partial)** | Added to species-pool limitation; v3 boundary already isolates effort-driven discovery; explicit vagrant-exclusion robustness = ROADMAP. |
| R3-3 | Observer heterogeneity / detection skill | **ROADMAP (P2)** | Acknowledged as limitation. |
| DA-1 | Lead with the non-obvious (interaction + offset defeat) | **DONE** | Abstract already leads with interaction + 7/8 offset defeat; retained and sharpened. |
| DA-2 | SDM circularity | **DONE** | New Methods paragraph separating SDM *candidacy* (static climatology) from hazard *timing* (time-varying anomaly + effort). |
| DA-3 | Generality (one taxon/country) | **DONE (reframe)** | New "Generality beyond birds in China" subsection: transferable method-level claims + pre-registered second-system extension. |
| DA-4 | Explain PCA offset tie | **DONE** | Noted PCA mixes effort axes and dilutes the moderation channel (in robustness-boundary discussion / offset reporting). |

## Items deliberately NOT fabricated (integrity)

Per the `verified_only` rule, the following were **not** invented and are carried
as explicit roadmap/limitations rather than reported with numbers: lagged-effort
refit, grid-native 50/100 km refit, spatial-block CV AUC, random-slope models,
vagrant-excluded refit, 6-metric interaction matrix. Producing any of these
requires running the corresponding pipeline scripts; the manuscript states this
plainly.

## Net effect of revision

The revision strengthens framing, honesty and interpretability without adding any
unverified number. The conceptual contribution (effort-as-moderator via the
interaction-vs-offset discriminator) and all 20 persisted quantitative claims are
unchanged. Residual P0 items (endogeneity identification, ecological-grid grain,
spatial-block CV) are now transparently the gating future analyses, which is the
honest state of the evidence.
