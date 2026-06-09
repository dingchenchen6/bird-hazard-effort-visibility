# Survey effort moderates the climatic visibility of new species records

*Formatted for* **Nature Ecology & Evolution** *(Article)*

**Chen-Chen Ding**¹

¹ Institute of Ecology, College of Urban and Environmental Sciences, Peking
University, Beijing 100871, China. e-mail: chenchending1992@gmail.com

---

## Abstract

New species records are a primary currency of climate-change biogeography, yet a
record is created only when ecological exposure and human observation coincide.
Whether survey effort merely scales how many records appear, or actively
conditions whether climate-driven change becomes detectable at all, is untested.
Using 12,813 species-by-province-by-year risk-set rows for 333 birds reported as
new distribution records across mainland China (2002–2024), we show with
discrete-time hazard models that climate exposure and effort combine
multiplicatively: their interaction is positive across four effort metrics
(hazard ratio 1.18–1.30) and decisively outperforms both an additive model and a
raw-effort offset, identifying effort as a moderator rather than a scaling term.
The interaction persists when effort is lagged a year, is strongest in resident
species, and survives a fifteen-fold risk-set expansion, but dissolves beyond a
well-modelled species core and transfers poorly across regions. We argue that new
records are effort-moderated visibility events, reframing how range-shift evidence
should be read and prioritised.

---

## Main

The redistribution of life under climate change is now among the most consequential
signals in ecology, and our knowledge of it is assembled overwhelmingly from
occurrence records — the documented appearance of a species where it had not been
seen before¹,². Such records anchor estimates of range-shift velocity, forecasts
of biodiversity turnover, and the spatial targeting of conservation. They are also
deeply ambiguous. A new record is not produced by ecology alone: it requires that
a species become present or detectable, that observers visit the right place in the
right season, and that a reporting system captures the encounter. The ecological
process we wish to study and the observation process that reveals it are entangled
in the very data we use to study it³.

This entanglement is not random noise. Occurrence data are concentrated in
accessible, populous and repeatedly visited landscapes, and the explosive growth of
citizen science can manufacture the appearance of biodiversity change where
ecological change is slight³⁻⁵. The standard remedy is to treat survey effort as an
additive covariate or as an offset that rescales expected counts. Both share a
hidden assumption: that effort changes the *number* of records independently of
ecological exposure. For new distribution records, that assumption is unlikely to
hold. Climatic exposure may raise the probability that a species reaches or
persists in a new place, but the colonisation becomes a *record* only where
observation pressure is sufficient to register it. If so, climate and effort should
act on records multiplicatively, and effort should behave as a *moderator* of the
climate signal rather than as a uniform multiplier.

The distinction is not semantic; it determines what a new-record map means. If
effort is additive, correcting for it recovers the climate signal and predicted
records track climatic exposure. If effort is a moderator, identical exposure
yields records in well-watched landscapes and silence elsewhere, so any map of
"where ranges are shifting" is partly a map of where people look — and any
conservation priority built on it inherits that bias. Despite its importance, the
moderator-versus-scaling question has not been posed as a direct, falsifiable test.

China offers a uniquely demanding arena in which to pose it. The country spans
tropical, temperate, arid, plateau and high-mountain biomes across steep climatic
gradients, and bird recording has grown by more than an order of magnitude since
2002. New bird records here may reflect genuine redistribution, intensifying
observation, or both, and the two are spatially confounded. Resolving them is a
precondition for using new-record data to infer range shifts or to design
monitoring.

We develop a multi-scale, discrete-time hazard framework for new bird records and
use it to test four hypotheses. **H1**: climate exposure and survey effort combine
interactively, not additively. **H2**: effort moderates the climate–hazard
relationship rather than scaling record counts — a contrast we operationalise as an
interaction model versus a raw-effort offset model, and defend against reverse
causation by lagging effort. **H3**: the interaction is general across effort
metrics, life histories, risk-set definitions and spatial scales, while its
strength is trait- and scale-dependent. **H4**: because the process is partly an
observation process, the framework is informative for inference and temporal
projection but bounded for spatial extrapolation, with direct consequences for
monitoring design.

### A discrete-time hazard framework for new records

We cast each potential new record as a survival problem. A row in the risk set is a
species-by-spatial-unit-by-year combination that was at risk of a first record
before the event year; the response is the binary first arrival. Our conservative
(headline) province risk set comprised 12,813 rows for 333 species and 512
first-arrival events across 32 mainland provincial units, 2002–2024, with each
species considered at risk only where a species distribution model judged its
presence plausible. We fitted a nested family of complementary log-log hazard
models with crossed random intercepts for species and spatial unit (Methods): a
null model, effort-only, climate-only, an additive climate-plus-effort model
(M3), a climate-by-effort interaction (M4), and a climate-plus-raw-effort-offset
model (M5). Two contrasts carry the argument — M4 against M3 (interaction versus
addition) and M4 against M5 (moderation versus scaling).

### Climate and effort combine multiplicatively

At province scale the climate-by-effort interaction was positive and significant
under every effort metric we tested (Fig. 1a, Table 1). The headline observer-visit
specification gave a hazard ratio of 1.29 (95% confidence interval 1.18–1.41,
P = 2.1 × 10⁻⁸); record-based, principal-component and birding-day metrics
reproduced the same direction and magnitude (1.18–1.30, all P ≤ 6 × 10⁻⁵). The
interaction model commanded essentially all of the model weight, beating the
additive model by ΔAIC ≈ 29 (Akaike weight ≈ 1). The biological reading is direct:
elevated effort amplified the hazard gained from climatic exposure, and exposure
amplified the hazard gained from effort, so an additive correction understates
record hazard exactly where exposure and observation coincide and overstates it
where only one is present.

### Effort moderates, rather than scales, record hazard

A positive interaction is necessary but not sufficient to establish moderation,
because a multiplicative scaling of effort would also break additivity. We
therefore pitted the interaction model directly against an offset model in which
effort enters on the raw log-rate scale with its coefficient fixed at one — the
formal encoding of "effort scales records." Across eight comparisons (four metrics ×
two risk-set definitions) the interaction model won seven, by AIC margins of several
hundred to over fifteen hundred (Fig. 1b). Effort changes the *slope* of the
climate–hazard relationship, not merely its intercept, which is the statistical
signature of a moderator. New records emerge where exposure and observation jointly
cross a visibility threshold, not where either acts alone.

The most serious alternative is reverse causation: observers may flock to places
with recent notable records, so a concurrent interaction could reflect effort
chasing discoveries rather than enabling them. Because effort is measured at the
province-year level, we replaced each row's current-year effort with the same
province's effort from the previous year and refitted. Effort that predates the
event year cannot respond to it, yet the lagged interaction was, if anything,
stronger (HR = 1.32, 95% CI 1.19–1.47, P = 1.6 × 10⁻⁷; ΔAIC = 26), ordering the two
in the direction moderation predicts and ruling out the simplest reverse-causation
account.

### A general mechanism with trait-dependent strength

Stratifying by migratory strategy relaxed the assumption of a single slope across
species (Fig. 2a). The interaction was positive and significant in every group but
strongest in resident species (HR = 1.43, 95% CI 1.21–1.69, P = 2.9 × 10⁻⁵) and
weaker, though still clear, in migrants (≈ 1.25). The pattern is ecologically
coherent: a resident's range expansion becomes a record only where observation
effort is sufficient, so effort moderation bites hardest for residents, whereas the
movement and stochasticity of migrant and vagrant occurrences dilute, without
erasing, the effort-moderated climate signal. The mechanism is therefore general
across life histories, while its magnitude is a trait.

### A robustness boundary, and the regime beyond it

The conservative risk set, by design, excludes many observed events. To probe its
limits we built a deliberately inclusive, event-override risk set (188,870 rows, 463
species, 817 events; a fifteen-fold expansion recovering 60% of previously excluded
events). The province interaction not only survived but tightened (HR = 1.27,
95% CI 1.20–1.35, P = 1.1 × 10⁻¹⁴), and a full parallel analysis on this set
reproduced every headline result — the lagged-effort defence (HR = 1.29,
P = 9.7 × 10⁻¹¹) and the migratory generality among them. Yet beyond the
well-modelled core the signal changed character: the interaction collapsed at finer
administrative scales (prefecture 1.04, county 1.02; both non-significant), random
forests demoted the temperature-by-effort term from first to sixth, and — reversing
the conservative pattern — the strongest interaction shifted to long-distance
migrants (HR = 1.33), the very taxa whose records track effort-rich landscapes.
This is not a failure of robustness but a *boundary*: a climate-by-effort visibility
mechanism in the modelled core gives way to an effort-density discovery process in
the long tail of borderline and vagrant species.

### Scale dependence and ecological-grain prediction

Propagating the conservative risk-set logic to finer administrative units preserved
the interaction in sign while attenuating it with grain (province 1.29 → prefecture
1.16 → county 1.11; Fig. 2b), the expected behaviour of a signal carried partly by
broad gradients in exposure and survey infrastructure. Residual spatial
autocorrelation was negligible at all tested distances, so the inference is not an
artefact of unmodelled structure. Because new records are documented at province
resolution, we fitted the model there and projected it onto the ecological grain:
applying the fitted relationship to grid-native climate (CRU monthly temperature,
2002–2024) and grid-native effort from the merged eBird/GBIF and China
Bird-Watching panel yields a continuous national risk surface (Fig. 3) in which the
effort-conditioned climate slope turns positive only where observation is
sufficient, rendering the moderation mechanism as geography.

### Predictive scope is itself a result

How far the framework predicts depends on what is being predicted (Fig. 2c). A
feature-rich gradient-boosted model interpolated well within the observed domain
(AUC 0.73) and forecast forward usably (temporal leave-future-out AUC 0.63, the
axis on which future projections rest), but transferred poorly to held-out
biogeographic regions (spatial-block AUC 0.56) and did not improve there with
feature richness. The record-generating process is region-specific: it can be
inferred and projected in time but should not be extrapolated to unobserved regions
on absolute predicted hazard. Carried forward under a four-model CMIP6 ensemble
(median warming +3.3 °C by 2050 and +4.8 °C by 2080 under SSP5-8.5) with
development-linked effort growth, projected record risk rises through the century
and concentrates where exposure and observation jointly intensify, most steeply in
eastern and central China (Fig. 3).

## Discussion

Our central result is that survey effort does not simply add records or scale them
up; it conditions whether climate-driven change is recorded at all. The
climate-by-effort interaction is positive and decisively favoured over the additive
model under every effort metric, survives the offset test that distinguishes
moderation from scaling, and holds when effort is measured before the event. New
distribution records are therefore best understood as effort-moderated visibility
events: ecological information filtered through an observation process that varies
by orders of magnitude across space and time. High exposure without sufficient
effort stays invisible; high effort in exposed places converts redistribution into
evidence. Treating effort as a nuisance to be subtracted, the field's default,
mis-specifies the record-generating process in the direction that biases
range-shift inference toward well-surveyed landscapes.

The interaction-versus-offset comparison is, we suggest, a reusable instrument.
Much of biodiversity informatics corrects for effort with a covariate or offset and
then proceeds as though the climate signal has been cleansed. Our result shows that
this step is a testable hypothesis rather than a settled procedure: where the offset
model loses decisively to the interaction model, effort is doing more than scaling,
and the "corrected" signal remains conditional on observation intensity. The test
is inexpensive and, we argue, should accompany any effort-aware analysis of
occurrence data.

The trait dependence sharpens the mechanism. That effort moderation is strongest in
residents fits a visibility logic: a resident expanding its range generates a record
only where someone is watching, so its detectability is tightly coupled to local
effort. Migrants and vagrants add a movement component that is partly decoupled from
local conditions, weakening the coupling in the conservative core; but when the
species pool is widened to include the borderline and vagrant taxa that dominate the
relaxed set, those same effort-tracked occurrences make the migrant interaction the
strongest of all. The contrast between the two risk sets is thus not a contradiction
but a window onto two coexisting regimes — climate-effort visibility for a modelled
core, effort-density discovery for the long tail — that any synthesis pooling all
species will blur. Reporting the core and the tail separately is the more honest and
more informative practice.

Scale carries its own lesson. The interaction attenuates from province to county
without reversing, the textbook fingerprint of a signal carried by broad gradients,
and our grid-native projection shows that the same moderation, resolved at the
ecological grain, produces a coherent national surface rather than an artefact of
administrative aggregation. The deliberate decision to fit where events are observed
(the province) and to project, rather than re-estimate, at finer grain is what keeps
that surface honest. Equally honest is the predictive boundary: the process
interpolates and forecasts in time but does not transfer across regions, because
what predicts records in one biogeographic setting need not predict them in another.
This is a property of the phenomenon, not a deficiency of the model, and it has a
clear practical corollary.

For monitoring under climate change, raw projected hazard is a poor guide to where
survey investment is most valuable. Our projections concentrate future risk in
eastern and central China — places where the model expects climatic exposure to
become highly *visible* under existing and growing effort, not necessarily where
range change is most novel. Allocating effort there maximises the rate of new
records but adds little to what is already well observed. The decision-relevant
target is the opposite: units where high projected exposure meets a current effort
gap, where modest additional survey converts latent redistribution into evidence and
most reduces uncertainty. We therefore recommend prioritising on the product of
climatic exposure and effort deficit, and reporting forecast uncertainty alongside
point hazard, rather than ranking regions on projected hazard alone.

Several boundaries frame these conclusions. The interaction slope, though shown to
vary by migratory strategy, remains coarsely structured; trait- and
phylogeny-resolved slopes are a natural next step. Effort, despite the lagged-effort
defence, is not fully exogenous, and slow-moving confounders that make a region both
well-watched and climatically dynamic cannot be excluded without an instrument or a
quasi-experiment in effort allocation. New records mix genuine colonisers with
vagrants and taxonomic revisions, a heterogeneity the relaxed analysis exposes but
does not fully model. Future surfaces inherit uncertainty from both climate
projections and assumed effort trajectories, and the ecological-grain map is a
province-fitted projection, not an independently estimated grid model. None of these
overturns the core finding; each marks where the framework can be strengthened.

New species records will continue to accumulate, faster than ever, as observation
expands. Read as effort-moderated visibility events rather than as direct
ecological signals, they become both more interpretable and more useful: a basis for
separating biodiversity redistribution from observation bias, and for designing
monitoring that learns most where it is currently blindest.

## Methods

### Study system and response variable

The study covered 32 mainland provincial units of China (excluding Hong Kong, Macau
and Taiwan), 2002–2024. The response was the first recorded arrival of a species in
a spatial unit and year, encoded as a binary event in a discrete-time risk set; once
a species records its first arrival in a unit, that unit leaves the species' risk set
in later years.

### Conservative and relaxed risk sets

The conservative (headline) risk set retained, for each species, provinces where a
species distribution model indicated plausible presence and no prior record, with
complete climate and effort information: 12,813 rows, 333 species, 512 events. To
stress-test the SDM filter we built a relaxed event-override set that adopted the
loosest binarisation threshold, force-included every species-province pair with an
observed event, and restricted the pool to SDM-modelled species: 188,870 rows, 463
species, 817 events.

### Climate and effort covariates

Climate covariates comprised temperature and precipitation anomalies and gradients,
climate velocity, warming rate and Mahalanobis climate displacement (WorldClim 2.1,
CHELSA v2.1)⁶,⁷. Grid-native time-varying temperature anomaly was computed per cell
from CRU TS 4.09 monthly data (2002–2024)⁸. Survey effort was represented by four
metrics — record counts, observer visits (headline), a principal-component composite
and birding-days — standardised within year; for the offset model, effort entered on
the raw log(person-hours + 1) scale. Grid and fine-scale effort came from the merged
eBird/GBIF + China Bird-Watching "Combined" panel⁵.

### Hazard models

We fitted discrete-time complementary log-log hazard models in glmmTMB⁹ with crossed
random intercepts for species and spatial unit: M0 (null), M1 (effort), M2 (climate),
M3 (climate + effort), M4 (climate × effort) and M5 (climate + raw-effort offset). We
report hazard ratios with 95% confidence intervals, P values, AIC differences and
Akaike weights. Endogeneity was probed by refitting M4 with effort lagged one year at
the province-year level; trait dependence by refitting separately by migratory
strategy. The model was fitted at province scale (where events are recorded) and
applied as a plug-in to finer administrative units and to 50/100 km ecological grids
with unit-native covariates.

### Diagnostics, prediction and projection

Residual spatial autocorrelation was assessed by Moran's I at 100/250/500 km. Predictive
performance was evaluated with a feature-rich gradient-boosted model (XGBoost) under
three cross-validation regimes — random (interpolation), temporal leave-future-out
(forecast) and 250 km spatial-block (extrapolation; blockCV)¹⁰ — with AUC by the
Mann–Whitney statistic. Future climate used a four-GCM CMIP6 ensemble (ACCESS-CM2,
MPI-ESM1-2-HR, MIROC6, UKESM1-0-LL) under SSP2-4.5 and SSP5-8.5⁴; effort futures used
development-linked growth scenarios. Modifiable-areal-unit elasticity, the
six-climate-metric interaction matrix and per-feature covariate-shift indices are
method-ready extensions not reported numerically here. Maps use the official
GS(2019)1822 basemap and render the national boundary and the South China Sea
nine-dash line.

### Reproducibility

Analyses used R (glmmTMB, xgboost, blockCV, terra, sf, targets, arrow) with versions
pinned in a lockfile and a targets pipeline; every reported value is traceable to a
named result table. Code and derived outputs are at
https://github.com/dingchenchen6/bird-hazard-effort-visibility.

## Data availability

Derived covariate tables, model outputs and data dictionaries are in the project
repository; non-sensitive products will be deposited at Zenodo on acceptance, with
licensed occurrence data shared under source-specific terms. Climate products are
public (WorldClim, CHELSA, CRU TS, CMIP6 via WorldClim/ESGF); the GS(2019)1822
basemap is from the Ministry of Natural Resources of China.

## Code availability

Version-controlled R/Python code, the renv lockfile and the targets pipeline are at
the repository above and will be archived at Zenodo on acceptance.

## References

1. Chen, I.-C. et al. *Science* **333**, 1024–1026 (2011).
2. Lenoir, J. et al. *Nat. Ecol. Evol.* **4**, 1044–1059 (2020).
3. Boakes, E. H. et al. *PLoS Biol.* **8**, e1000385 (2010).
4. Isaac, N. J. B. et al. *Methods Ecol. Evol.* **5**, 1052–1060 (2014).
5. Sullivan, B. L. et al. *Biol. Conserv.* **169**, 31–40 (2014).
6. Fick, S. E. & Hijmans, R. J. *Int. J. Climatol.* **37**, 4302–4315 (2017).
7. Karger, D. N. et al. *Sci. Data* **4**, 170122 (2017).
8. Harris, I. et al. *Sci. Data* **7**, 109 (2020).
9. Brooks, M. E. et al. *R Journal* **9**, 378–400 (2017).
10. Valavi, R. et al. *Methods Ecol. Evol.* **10**, 225–232 (2019).
11. Loarie, S. R. et al. *Nature* **462**, 1052–1055 (2009).
12. Eyring, V. et al. *Geosci. Model Dev.* **9**, 1937–1958 (2016).
13. Openshaw, S. *The Modifiable Areal Unit Problem* (Geo Books, 1984).

*(Full reference details with DOIs in `pipeline_finalization_and_citation_audit.md`;
NEE uses superscript numbered citations — renumber on final typesetting. CMIP6,
climate-velocity and MAUP references appear in Methods.)*

## Display items

- **Fig. 1** Province headline: (a) climate × effort interaction across four effort
  metrics; (b) interaction (M4) vs offset (M5) — moderation, not scaling.
  [`Figure_2_v4_province_headline_M5_raincloud`]
- **Fig. 2** (a) interaction by migratory strategy; (b) attenuation across
  administrative scales; (c) predictive AUC across interpolation/forecast/extrapolation.
  [`Figure_9_migratory_and_prediction`, `Figure_3_v4_three_scale_forest_M5`]
- **Fig. 3** Ecological-grain risk surface: grid-native CRU climate × merged effort,
  current and CMIP6 + effort-scenario futures, cool→warm, with nine-dash line.
  [`Figure_8_grid100_native_plugin_hazard`]
- **Table 1** Province interaction estimates (four effort specifications).
- **Extended Data** robustness boundary (relaxed set), v1/v2/v3 reconciliation,
  Moran's I, CMIP6 ensemble deltas, feature importance, 50 km grid, v3 grid surface.

## Author contributions

C.-C.D. conceived the study, designed and implemented the analyses, and wrote the
manuscript.

## Competing interests

The author declares no competing interests.

## Additional information

**Funding** [to be completed]. **AI-usage disclosure** Generative AI assisted with
drafting, code organisation and prose editing under the author's direct supervision;
all analyses, numerical results and interpretations were produced and verified by the
author against persisted result tables, for which the author takes full
responsibility.
