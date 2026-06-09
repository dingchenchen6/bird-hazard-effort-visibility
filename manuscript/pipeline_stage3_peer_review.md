# Stage 3 — Simulated Peer Review (Nature Ecology & Evolution / Ecology Letters tier)

**Manuscript:** *Survey effort moderates the climatic visibility of new bird
distribution records across China*
**Panel:** Editor-in-Chief + Reviewer 1 (macroecology/SDM) + Reviewer 2
(statistical ecology) + Reviewer 3 (citizen science / monitoring) + Devil's
Advocate. Scoring 0–100 per dimension; honest, non-sycophantic.

---

## Reviewer 1 — Macroecology & species distributions

**Summary.** A genuinely interesting reframing: effort as a *moderator* of
climatic visibility rather than a nuisance offset, tested with a discrete-time
hazard framework on Chinese bird new-records. The M4-vs-M5 contrast is elegant
and the v2/v3 robustness boundary is a sophisticated touch.

**Major concerns.**
1. **Ecological vs administrative scale.** The multi-scale claim rests entirely on
   administrative units (province/prefecture/county). Administrative grain
   conflates spatial grain with administrative heterogeneity (area, shape, effort
   institutions). The MAUP argument is only convincing with *regular ecological
   grids* (50/100 km), which the manuscript itself flags as not yet persisted.
   For this tier the grid-native refit is not optional — it is the test that
   separates "real grain attenuation" from "administrative artefact." (P0)
2. **Which climate signal?** The headline interaction uses a single climate term
   (`climate_z`, temperature-gradient/anomaly). A reader cannot tell whether the
   moderation is specific to thermal exposure or general. The 6-climate × 4-effort
   matrix is described but not reported. At minimum, report the interaction for
   ≥3 climate metrics (temp anomaly, climate velocity, Mahalanobis displacement)
   in the main text. (P0)
3. **Effect size in biological terms.** HR = 1.29 per SD of the product is hard to
   interpret. Translate into absolute hazard change for representative
   low/median/high effort × exposure cells, and state how many extra records this
   implies per province-decade. (P1)

**Scores:** Originality 84 · Rigor 70 · Evidence 68 · Coherence 82 · Writing 86.

---

## Reviewer 2 — Statistical ecology

**Summary.** Methodologically careful; the offset-vs-interaction discriminator is
the paper's best idea and is correctly specified (offset coefficient fixed at 1).
The integrity of "only persisted numbers" is commendable.

**Major concerns.**
1. **Endogeneity / reverse causation (the central threat).** Effort is not
   exogenous: observers preferentially visit cells with recent notable records,
   and a new record itself attracts effort. The positive climate × effort
   interaction is consistent with moderation *and* with effort responding to
   climate-driven discoveries. The manuscript acknowledges this in one sentence
   but offers no design-based defense. Options: (i) lag effort by ≥1 year relative
   to the event; (ii) a leave-future-out temporal validation; (iii) an
   instrument for effort (e.g., human population/road density). Without at least
   the lag analysis, "moderator" overreaches. (P0)
2. **Spatial-block CV.** Moran's I being ~0 is reassuring but is a residual test,
   not a generalisation test. Random-CV AUC (and the 0.777 figure) overstate
   transferable skill under spatial structure. The 250 km block CV is scripted
   but unrun; for any predictive claim it must be reported. (P0)
3. **Shared slope across 333 species.** A single interaction slope assumes the
   moderation is taxon-invariant. At least test a random slope `(climate ||
   species)` or stratify by migratory strategy / range size, and report whether
   the fixed interaction survives. (P1)
4. **Multiple comparisons / model-selection inference.** Reporting p-values from
   the selected M4 understates selection uncertainty; consider reporting the
   interaction under model averaging across M3/M4. (P2)

**Scores:** Originality 80 · Rigor 64 · Evidence 66 · Coherence 80 · Writing 84.

---

## Reviewer 3 — Citizen science & monitoring design

**Summary.** The monitoring framing (target climate exposure × effort gap, not raw
hazard) is the most policy-useful contribution and is well argued.

**Major concerns.**
1. **Frozen-effort forecasts are internally inconsistent.** The thesis is that
   effort matters, yet the future projections freeze effort at 2024 and vary only
   climate. The eastern concentration of future hazard then largely restates the
   current effort map. Either (a) project effort under explicit growth scenarios,
   or (b) drop absolute future maps and present only the *climate-exposure ×
   effort-gap* prioritisation surface, which is the paper's real recommendation.
   (P0/P1)
2. **Vagrants vs colonisers.** New records mix genuine range expansions with
   vagrant/irruptive occurrences and taxonomic splits. The v3 boundary partly
   captures this, but the main analysis should state how vagrants are treated and
   ideally show the interaction is robust to excluding obvious vagrants. (P1)
3. **Observer heterogeneity.** "Visits" conflate skilled and novice observers;
   detection probability is not modelled. Acknowledge and, if possible, control
   for a proxy of observer skill. (P2)

**Scores:** Originality 82 · Rigor 70 · Evidence 70 · Coherence 84 · Writing 85.

---

## Devil's Advocate

- **"Is this just 'people find things where they look, more so in changing
  climates'?"** State the non-obvious content sharply: the *interaction* (not the
  main effects) and the *offset defeat* are what distinguish moderation from the
  trivial reading. Make sure the abstract leads with that, not with "effort
  matters."
- **Circularity risk.** If the SDM that defines the risk set already used climate,
  then conditioning on the SDM candidate set and then finding a climate
  interaction is partly built-in. Clarify that the SDM defines *candidacy* (where
  a species could occur) while the hazard model tests *timing* (when a record
  appears), and show the interaction is not an artefact of the SDM's own climate
  inputs. (P0)
- **Generality.** One country, one taxon. For this tier, either demonstrate the
  framework on a second system or argue, with mechanism, why birds-in-China is a
  sufficient proof-of-concept and pre-register the generalisation. (P1)
- **Negative space.** The PCA spec tie in the offset test (relaxed set) is the one
  crack; do not bury it — explain why PCA behaves differently (it mixes effort
  axes and dilutes the moderation channel).

**Scores:** Originality 78 · Rigor 62 · Evidence 64 · Coherence 80 · Writing 84.

---

## Editor-in-Chief — Editorial Decision

**Decision: MAJOR REVISION** (borderline reject-and-resubmit at *Nature Ecol
Evol*; clear major revision at *Ecology Letters*).

The core idea — effort as a moderator of climatic visibility, with a clean
interaction-vs-offset discriminator — is novel and well-suited to the tier. The
honesty discipline is a strength. But three issues currently cap the contribution:
**(1) endogeneity of effort is not defended; (2) the multi-scale and predictive
claims rely on analyses that are scripted but unreported (grid-native scale,
spatial-block CV); (3) generality beyond one taxon/country is asserted, not
shown.** None is fatal; all are addressable.

**Aggregate (weighted) score:** Originality 0.20 · Rigor 0.25 · Evidence 0.25 ·
Coherence 0.15 · Writing 0.15 →
**≈ 72/100** (R1 76, R2 71, R3 76, DA 69 → mean ≈ 73).

---

## Revision Roadmap (prioritised)

### P0 — must address (gate to acceptance)
- **P0-A Endogeneity defense.** Add a ≥1-year effort lag refit of M4 and report
  whether the interaction holds; discuss instruments. *(Needs a refit — pipeline
  task; if not run, must be reframed as an explicit, prominent limitation and the
  word "moderator" softened to "is consistent with moderation".)*
- **P0-B Climate-metric breadth.** Report the headline interaction for ≥3 climate
  metrics (use the persisted v2/v3 tables where available; the RF table already
  ranks temp_anom, temp_grad, climate_velocity, Mahalanobis). Add a paragraph and,
  if available, a panel.
- **P0-C SDM-circularity clarification.** Add a Methods paragraph separating SDM
  *candidacy* from hazard *timing*, explicitly stating the hazard model's climate
  term is the time-varying anomaly, not the static SDM climate inputs.
- **P0-D Scale honesty.** Either report grid-native 50/100 km (pipeline task) or
  reframe the multi-scale section as "administrative-scale" and move the ecological
  grid to a clearly-labelled future test — do not let the reader infer ecological
  grain.

### P1 — strongly recommended
- **P1-A Effect-size translation** (absolute hazard at low/median/high cells;
  extra records per province-decade).
- **P1-B Forecast reframing** — lead with the climate-exposure × effort-gap
  prioritisation; demote or scenario-ise the frozen-effort absolute maps.
- **P1-C Vagrant robustness** — show interaction survives excluding obvious
  vagrants (v3 boundary already gestures at this).
- **P1-D Generality framing** — explicit argument + pre-registered generalisation,
  or a second-system illustration.

### P2 — optional polish
- Random-slope / migratory-strategy stratification; model-averaged interaction;
  observer-skill proxy; explain PCA offset-tie.

### Revisable now with persisted data + framing (no new modelling)
P0-B (partial, via RF/spec tables), P0-C, P0-D, P1-B, P1-D, and the abstract
sharpening from the Devil's Advocate. These will be applied in Stage 4.

### Requires new modelling (flag as roadmap / acknowledged limitation)
P0-A (effort lag), P0-D grid-native refit, spatial-block CV, random-slope.
