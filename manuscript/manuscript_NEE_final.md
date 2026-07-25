# Survey effort moderates the climatic visibility of new species records

*Formatted for* **Nature Ecology & Evolution** *(Article)*

**Chen-Chen Ding**¹

¹ Institute of Ecology, College of Urban and Environmental Sciences, Peking
University, Beijing 100871, China. e-mail: chenchending1992@gmail.com

---

## Abstract

New species records are a primary currency of climate-change biogeography, yet a
record is created only when ecological exposure and human observation coincide.
Whether survey effort merely scales how many records appear, or actively conditions
whether climate-driven change becomes detectable, is untested. Using complete
discrete-time risk sets for 463 birds newly recorded across mainland China
(2002–2024; 188,870 species-by-province-by-year rows, 817 first records), we show
that climate exposure and effort combine multiplicatively. Their interaction is
positive and invariant across three species-distribution-model thresholds
(hazard ratio 1.274–1.290), and the interaction model outperforms both an additive
model and a raw-effort offset at every threshold — identifying effort as a
moderator rather than a scaling term. Partitioning explained deviance, effort
accounts for ~79 %, the interaction 12–30 % and climate alone under 6 %: climate
enters the record process chiefly *through* effort. The interaction is carried
specifically by thermal anomalies, is present in every migratory strategy and
strongest in long-distance migrants, and survives lagging effort by a year.
Mechanistic and machine-learning projections agree until covariates leave the
training domain and then diverge sharply, bounding credible forecasts near 2050.
New records are best read as effort-moderated visibility events.

---

## Main

The redistribution of life under climate change is among the most consequential
signals in ecology, and our knowledge of it is assembled overwhelmingly from
occurrence records — the documented appearance of a species where it had not been
seen before¹,². Such records anchor estimates of range-shift velocity, forecasts of
biodiversity turnover and the spatial targeting of conservation. They are also
deeply ambiguous. A new record is not produced by ecology alone: it requires that a
species become present or detectable, that observers visit the right place in the
right season, and that a reporting system capture the encounter. The ecological
process we wish to study and the observation process that reveals it are entangled
in the very data we use to study it³.

This entanglement is not random noise. Occurrence data concentrate in accessible,
populous and repeatedly visited landscapes, and the explosive growth of citizen
science can manufacture the appearance of biodiversity change where ecological
change is slight³⁻⁵. The standard remedy is to treat survey effort as an additive
covariate or as an offset that rescales expected counts. Both share a hidden
assumption: that effort changes the *number* of records independently of ecological
exposure. For new distribution records that assumption is unlikely to hold.
Climatic exposure may raise the probability that a species reaches or persists in a
new place, but the colonisation becomes a *record* only where observation pressure
suffices to register it. If so, climate and effort should act on records
multiplicatively, and effort should behave as a *moderator* of the climate signal
rather than as a uniform multiplier.

The distinction determines what a new-record map means. If effort is additive,
correcting for it recovers the climate signal and predicted records track climatic
exposure. If effort is a moderator, identical exposure yields records in
well-watched landscapes and silence elsewhere, so any map of "where ranges are
shifting" is partly a map of where people look — and any conservation priority
built on it inherits that bias. Despite its importance, the moderator-versus-scaling
question has not been posed as a direct, falsifiable test.

China offers a demanding arena in which to pose it. The country spans tropical,
temperate, arid, plateau and high-mountain biomes across steep climatic gradients,
and bird recording has grown by more than an order of magnitude since 2002. New
bird records here may reflect genuine redistribution, intensifying observation, or
both, and the two are spatially confounded. Resolving them is a precondition for
using new-record data to infer range shifts or to design monitoring.

We develop a multi-scale, discrete-time hazard framework and test four hypotheses.
**H1**: climate exposure and survey effort combine interactively, not additively.
**H2**: effort moderates the climate–hazard relationship rather than scaling record
counts — operationalised as an interaction model versus a raw-effort offset model,
and defended against reverse causation by lagging effort. **H3**: the interaction is
general across climate and effort proxies, life histories, risk-set definitions and
spatial scales, while its strength is proxy-, trait- and scale-dependent. **H4**:
because the process is partly an observation process, the framework is informative
for inference and near-term projection but bounded for extrapolation, with direct
consequences for monitoring design.

### Complete risk sets across species-distribution-model thresholds

We cast each potential new record as a survival problem: a row is a
species-by-province-by-year combination at risk of a first record before the event
year, and the response is the binary first arrival. Crucially, we analyse **complete
risk sets** — no complete-case truncation and no case-control subsampling — and we
rebuild them at three thresholds of the species-distribution model that defines
which species are considered at risk where. Every observed event is force-included,
so the threshold changes the candidate (non-event) denominator without altering the
numerator (Table 1). The three sets contain 188,870, 174,859 and 164,050 rows for
463 species and 817 events each; an earlier, complete-case-truncated set of 12,813
rows and 512 events is retained only as a reference.

**Table 1 | Complete risk sets across species-distribution-model thresholds.**

| Threshold | Candidate pairs | Rows | Species | Provinces | Events |
|---|---:|---:|---:|---:|---:|
| 50 (most permissive) | 8,931 | 188,870 | 463 | 33 | 817 |
| 100 | 8,239 | 174,859 | 463 | 33 | 817 |
| 200 (most restrictive) | 7,739 | 164,050 | 463 | 33 | 817 |
| *truncated reference* | *7,764* | *12,813* | *333* | *32* | *512* |

We fitted a nested family of complementary log-log hazard models with crossed random
intercepts for species and province (Methods): a null model (M0), effort-only (M1),
climate-only (M2), additive climate-plus-effort (M3), a climate-by-effort interaction
(M4) and a climate-plus-raw-effort-offset model (M5). Two contrasts carry the
argument — M4 against M3 (interaction versus addition) and M4 against M5 (moderation
versus scaling).

### Climate and effort combine multiplicatively

**M4 is the best-supported model at every threshold** (Fig. 1a). Under the headline
effort proxy it beats the additive model M3 by ΔAIC = 51–56 and the offset model M5
by 79–81, while the null and climate-only models lag by 332–351 and 317–337. The
interaction hazard ratio is essentially unaffected by the threshold: 1.274
(95 % confidence interval 1.20–1.35, *P* = 1.1 × 10⁻¹⁴) at threshold 50, 1.289
(1.21–1.37, *P* = 6.2 × 10⁻¹⁶) at 100 and 1.290 (1.21–1.37, *P* = 4.7 × 10⁻¹⁶) at
200. Across all five effort proxies the interaction spans 1.12–1.43 (all
*P* < 10⁻⁵).

A second feature of the ladder is diagnostic. Effort alone (M1, ΔAIC = 62–68)
explains far more than climate alone (M2, ΔAIC = 317–337): climate contributes
little as a main effect and much as a modifier. This is the paper's central claim in
its most compact form.

### Effort moderates, rather than scales, record hazard

A positive interaction is necessary but not sufficient to establish moderation,
because multiplicative scaling of effort would also break additivity. We therefore
pitted the interaction model against an offset model in which effort enters on the
raw log-rate scale with its coefficient fixed at one — the formal encoding of
"effort scales records". On complete risk sets the interaction model wins at **all
three thresholds** (Fig. 1a), a cleaner result than on the truncated set, where the
offset model occasionally tied. Effort changes the *slope* of the climate–hazard
relationship, not merely its intercept, which is the statistical signature of a
moderator.

The most serious alternative is reverse causation: observers may flock to places
with recent notable records. Because effort is measured at the province-year level,
we replaced each row's current-year effort with the same province's effort from the
previous year and refitted. Effort that predates the event year cannot respond to
it, yet the lagged interaction was, if anything, stronger (hazard ratio 1.29,
*P* = 9.7 × 10⁻¹¹), ordering the two variables in the direction moderation predicts.

### The interaction runs through a thermal channel

To test whether the result is an artefact of one metric choice we fitted the
interaction for every combination of nine climate and five effort proxies
(Fig. 1b). The interaction is **proxy-specific rather than a generic climate
effect**. Temperature anomaly and its gradient interact strongly with every effort
proxy (1.12–1.43, all *P* < 0.001) and precipitation anomaly moderately
(1.05–1.23); Mahalanobis climate displacement is weak (1.01–1.13); and
climate-velocity, climate-exposure and warming-rate proxies show no interaction at
all (0.88–1.06). Moderation therefore operates through year-to-year thermal
exposure, not through composite indices of climate change. We verified and report
that two proxy pairs in the source panel are perfectly collinear after
standardisation (temperature gradient with temperature anomaly; precipitation
gradient with precipitation anomaly, both *r* = 1.00) and that climate exposure and
warming rate are near-redundant (*r* = 0.89), so the nine proxies span seven
independent axes.

### A general mechanism with trait-dependent strength

Stratifying by migratory strategy relaxed the assumption of a single slope across
species (Fig. 2a). The interaction is positive and significant in **every group at
every threshold**, and the ordering is stable: long-distance migrants strongest
(1.325–1.354), then migrants overall (1.291–1.314), residents (1.278–1.310) and
partial migrants (1.220–1.254). Effort moderation is thus a general property of the
record-generating process rather than a feature of one life history. Notably, the
truncated set had implied the opposite ordering, with residents strongest — an
artefact of the species filtering that removed many migrant and vagrant taxa, and a
concrete illustration of why complete risk sets matter.

### Partitioning the roles of effort, climate and their interaction

Decomposing the deviance explained by M4 relative to the null model gives a
remarkably stable partition across thresholds (Fig. 2b): **unique effort ~79 %, the
climate-by-effort interaction 16 %, unique climate ~4 % and shared components ~1 %**.
Across effort proxies the interaction share ranges 12–30 % and unique climate never
exceeds 6 %. Variance components of M4 are marginal *R*² = 0.121–0.124 and
conditional *R*² = 0.596–0.609: fixed effects explain about a eighth of the
variance and the full model, including species and province random effects, about
three fifths.

The quantitative message is that survey effort is the dominant single driver of
whether and when a record appears, climate contributes almost nothing on its own,
and a substantial 12–30 % arises only once climate is allowed to interact with
effort.

Tree-ensemble variable importance tells a superficially contradictory story
(Fig. 2c): random-forest permutation importance ranks climate proxies at 88 % and
effort at −6 %, and gradient boosting at 55 % versus 21 %. The two approaches answer
different questions. The hazard model conditions on species and province random
effects, so its fixed effects describe within-province temporal variation, where
growing effort is highly informative. The tree models have no such structure and
instead use spatially structured climate proxies to encode *where* a row is; effort,
which rises almost everywhere, then carries little discriminative signal. That
effort's contribution emerges only once species and place are conditioned on is
itself consistent with effort acting as a **conditional moderator** rather than a
free-standing predictor.

### Scale dependence and ecological-grain prediction

Propagating the risk-set logic to finer administrative units preserved the
interaction in sign while attenuating it with grain (province → prefecture →
county), the expected behaviour of a signal carried partly by broad gradients in
exposure and survey infrastructure; residual spatial autocorrelation was negligible
at all tested distances. Because new records are documented at province resolution,
we fitted the model there and projected it onto the ecological grain: applying the
fitted relationship to grid-native climate (CRU monthly temperature, 2002–2024) and
grid-native effort from the merged eBird/GBIF and China Bird-Watching panel yields a
continuous national risk surface (Fig. 3) in which the effort-conditioned climate
slope turns positive only where observation is sufficient, rendering the moderation
mechanism as geography.

### Prediction accuracy depends on the task, not the algorithm

We validated both model classes on the complete risk set under three regimes
(Fig. 4c). For within-domain interpolation the machine-learning model is clearly
better (area under the curve 0.774 versus 0.694; precision-recall 0.014 versus
0.008). For temporal forecasting and for spatial extrapolation to held-out province
blocks, however, the two classes converge (0.615–0.630 and 0.613–0.616): additional
flexibility buys nothing once the prediction target leaves the observed domain.
Calibration slopes tell the same story, falling from 0.91–0.97 under random folds to
0.45–0.65 under temporal and spatial validation, so absolute probabilities — not only
rankings — should be treated cautiously outside interpolation. The difficulty of
forecasting new records is therefore a property of the data-generating process rather
than of any algorithm.

### Two model classes bound the credible forecast horizon

We projected forward with both a mechanistic hazard model and a gradient-boosted
machine-learning model, trained on the same complete risk sets, under a four-model
CMIP6 ensemble (median warming +3.3 °C by 2050 and +4.8 °C by 2080 under SSP5-8.5)
with development-linked effort growth (Fig. 4a). The two classes agree early and
diverge sharply later: relative to 2024 the mechanistic model rises 2.9-fold by
2050 and 13-fold by 2080 under SSP2-4.5 (10.4- and 131-fold under SSP5-8.5), whereas
the machine-learning model saturates near 2.4-fold in every scenario.

The divergence is explained precisely (Fig. 4b). The share of provinces whose
projected climate leaves the training range grows from 0 % in 2030 to 18 % by 2050
and 55 % by 2080 under SSP2-4.5, and reaches 48 % and 91 % under SSP5-8.5, with
effort leaving its training range in 100 % of provinces by 2080. The mechanistic
model extrapolates linearly on the complementary log-log scale and keeps rising — a
131-fold increase is almost certainly an over-extrapolation — while the tree
ensemble cannot extrapolate past its training splits and therefore flattens.
Neither behaviour is correct in isolation; together they **bound the credible
forecast horizon**, which on this evidence reaches roughly 2050 under SSP2-4.5,
where fewer than a fifth of units have left the observed domain, and degrades
sharply thereafter. Rank agreement between the classes is moderate
(Spearman ρ ≈ 0.3–0.6).

### New records are displaced polewards, but the direction is taxon-structured

If effort-moderated visibility tracks a real redistribution, new records should not
fall randomly around historical ranges. Referencing 851 records (564 species) to
each species' BirdLife historical range, they are strongly directionally biased
(Fig. 5). Relative to range centroids the mean bearing is 55.7° (north-east;
concentration *R* = 0.270, Rayleigh *P* < 0.001) with median displacement 1,273 km;
relative to the nearest range edge the mean bearing is 21.1° (nearly due north;
*R* = 0.355, *P* < 0.001) with median displacement 525 km. The edge-referenced
distribution is both tighter and more sharply northward — the cleaner signature of
range extension — with more than a fifth of all records falling in the single north
sector. Ninety-five per cent of the displacement signal is therefore consistent with
poleward redistribution, and 81 records (9.5 %) fall inside the historical range and
represent within-range detection rather than extension.

The aggregate pattern conceals systematic taxonomic structure. Passerines
(*n* = 465) shift east-north-east (66.0°, *R* = 0.380, *P* < 0.001) and raptors
north-east (36.5°, *P* = 0.021), whereas shorebirds trend north-west (323.5°,
*P* = 0.040) and waterfowl are displaced *south-west* (196.6°, *P* = 0.022) and by
the largest distances of any order (median 2,169 km from the centroid, versus
1,190 km for passerines). A single "poleward" summary would therefore misdescribe
almost half the orders: land birds expand along breeding-range margins while
waterbird records accumulate along flyways and toward wintering grounds.

## Discussion

Survey effort does not simply add records or scale them up; it conditions whether
climate-driven change is recorded at all. On complete risk sets the
climate-by-effort interaction is positive, threshold-invariant and decisively
favoured over both the additive and the offset formulation, and it survives
measuring effort a year before the event. New distribution records are therefore
best understood as effort-moderated visibility events: ecological information
filtered through an observation process that varies by orders of magnitude across
space and time. High exposure without sufficient effort stays invisible; high effort
in exposed places converts redistribution into evidence. Treating effort as a
nuisance to be subtracted — the field's default — mis-specifies the
record-generating process in the direction that biases range-shift inference toward
well-surveyed landscapes.

The deviance partition sharpens this into a quantitative claim: effort alone
accounts for about four fifths of what the model explains, climate alone for under a
twentieth, and the interaction for one sixth to nearly one third depending on the
effort proxy. Read carelessly, the dominance of effort might suggest that new
records carry little ecological signal. Read correctly, it says that the ecological
signal is *conditional*: climate becomes visible in the data only in combination
with observation, which is precisely why the interaction — not the climate main
effect — is where the biology appears. The contrast with tree-ensemble importance
reinforces the point, because effort's contribution materialises only once species
identity and place are conditioned on.

That the interaction runs specifically through thermal anomalies, and not through
climate velocity, exposure or warming rate, is informative in two ways. Empirically
it identifies the channel: interannual thermal departures, not long-run
climate-change indices, are what interact with observation to generate records.
Methodologically it warns that "climate" is not interchangeable across proxies in
detection-aware models — a study choosing velocity or exposure as its single
climate axis would have concluded, wrongly, that no interaction exists. The
collinearity we document within the proxy set makes the same point from the other
side: apparent replication across metric names can be a single axis in disguise.

The interaction-versus-offset comparison is a reusable instrument. Much of
biodiversity informatics corrects for effort with a covariate or offset and then
proceeds as though the climate signal has been cleansed. Our result shows this step
is a testable hypothesis rather than settled procedure: where the offset model loses
decisively, effort is doing more than scaling, and the "corrected" signal remains
conditional on observation intensity. The test is inexpensive and should accompany
any effort-aware analysis of occurrence data.

Life-history structure adds nuance without undermining generality. Effort moderation
appears in residents, partial migrants and long-distance migrants alike, and is
strongest in long-distance migrants, whose occurrences accumulate where observation
is dense. That the truncated risk set reversed this ordering is a caution about
species filtering: conclusions about *which* taxa show a detection-mediated climate
signal are far more sensitive to risk-set construction than the existence of the
signal itself.

The two model classes together define where forecasting stops being defensible.
Mechanistic extrapolation and machine-learning saturation are not competing
estimates of the same quantity; they are the two failure modes that bracket the
truth once covariates leave the observed domain. Reporting both, alongside the
fraction of units outside the training range, converts an unfalsifiable
century-scale projection into a bounded near-term one. On this evidence, projections
to mid-century under moderate emissions are supportable; projections to 2080,
especially under high emissions and rapid monitoring growth, are not.

The directional analysis supplies independent, effort-free corroboration. Nothing
in the bearing of a record relative to its species' historical range depends on the
hazard model, yet records fall north-east of range centroids and almost due north of
range edges, with highly significant angular concentration. That the effort-moderated
hazard framework and the geometry of the records point the same way strengthens the
claim that new records carry a genuine redistribution signal rather than only an
observation signal. The taxonomic structure adds mechanism: passerines and raptors
behave as expected of breeding-range expansion, whereas waterfowl move south-west and
farthest, consistent with flyway geometry and wintering-ground redistribution rather
than poleward breeding expansion. Directional analyses that pool orders will average
these opposing signals toward zero.

For monitoring under climate change, raw projected hazard is a poor guide to where
survey investment is most valuable. Projections concentrate future risk where
climatic exposure will become highly *visible* under existing and growing effort,
not necessarily where range change is most novel. Allocating effort there maximises
the rate of new records but adds little to what is already well observed. The
decision-relevant target is the opposite: units where high projected exposure meets
a current effort gap, where modest additional survey converts latent redistribution
into evidence and most reduces uncertainty. We therefore recommend prioritising on
the product of climatic exposure and effort deficit, and reporting forecast
uncertainty alongside point hazard.

Boundaries remain. The interaction slope, though shown to vary by migratory
strategy, is coarsely structured; trait- and phylogeny-resolved slopes are a natural
next step. Effort, despite the lagged-effort defence, is not fully exogenous, and
slow-moving confounders that make a region both well-watched and climatically
dynamic cannot be excluded without an instrument or a quasi-experiment in effort
allocation. New records mix genuine colonisers with vagrants and taxonomic
revisions, a heterogeneity the threshold and migratory analyses expose but do not
fully model. Future surfaces inherit uncertainty from both climate projections and
assumed effort trajectories, and the ecological-grain map is a province-fitted
projection, not an independently estimated grid model.

New species records will continue to accumulate, faster than ever, as observation
expands. Read as effort-moderated visibility events rather than as direct ecological
signals, they become both more interpretable and more useful: a basis for separating
biodiversity redistribution from observation bias, and for designing monitoring that
learns most where it is currently blindest.

## Methods

### Study system and response variable

The study covered 33 provincial units of mainland China, 2002–2024. The response was
the first recorded arrival of a species in a province and year, encoded as a binary
event in a discrete-time risk set; once a species records its first arrival in a
unit, that unit leaves the species' risk set in later years.

### Complete risk sets and threshold sensitivity

Risk sets were rebuilt from source at three species-distribution-model binarisation
thresholds (50, 100, 200). At each threshold the candidate set comprised
species-province pairs with modelled potential presence and no prior record, plus
forced inclusion of every pair with an observed 2002–2024 event (an observed record
being empirical evidence that overrides the model prior), restricted to the union of
species modelled by the two source distribution-model projects. Candidates were
expanded over years and censored after first arrival. No complete-case truncation or
case-control subsampling was applied at any stage; the only exclusions were rows
lacking the headline climate or effort covariate.

### Climate and effort proxies

Nine climate proxies were evaluated: temperature and precipitation anomalies and
gradients, climate velocity, precipitation velocity, climate exposure, warming rate
and Mahalanobis climate displacement (derived from WorldClim 2.1 and CHELSA
v2.1)⁶,⁷. Pairwise collinearity was audited and is reported. Grid-native
time-varying temperature anomaly was computed per cell from CRU TS 4.09 monthly data
(2002–2024)⁸. Five effort proxies were evaluated: record counts, observer visits
(headline), observers, birding-days and a principal-component composite,
standardised within year; for the offset model effort entered on the raw log-rate
scale. Grid and fine-scale effort came from the merged eBird/GBIF and China
Bird-Watching panel⁵.

### Hazard models and decomposition

Discrete-time complementary log-log hazard models were fitted in glmmTMB⁹ with
crossed random intercepts for species and province: M0 (null), M1 (effort),
M2 (climate), M3 (climate + effort), M4 (climate × effort) and M5 (climate +
raw-effort offset). We report hazard ratios with 95 % confidence intervals,
*P* values and AIC differences. Endogeneity was probed by refitting M4 with effort
lagged one year at the province-year level; trait dependence by refitting separately
by migratory strategy. Explained deviance was partitioned from the ladder
log-likelihoods into unique effort (M2 → M3), unique climate (M1 → M3), shared and
interaction increment (M3 → M4) components, with McFadden *R*² relative to M0;
marginal and conditional *R*² followed a Nakagawa-style decomposition of fixed,
random and distribution-specific variance. Relative importance used random-forest
permutation importance (probability forest with class weights) and gradient-boosting
gain on the same complete risk sets.

### Prediction and projection

Predictive performance was evaluated under three cross-validation regimes — random
(interpolation), temporal leave-future-out (forecast) and 250 km spatial-block
(extrapolation)¹⁰ — with area under the curve computed from the Mann–Whitney
statistic. Future climate used a four-GCM CMIP6 ensemble (ACCESS-CM2,
MPI-ESM1-2-HR, MIROC6, UKESM1-0-LL) under SSP2-4.5 and SSP5-8.5 for 2030, 2050 and
2080¹²; effort futures used development-linked growth scenarios (+0.3 and +0.6
standard deviations per decade). Projections were made with both the mechanistic
hazard model and a gradient-boosted model, and we report, for every scenario, the
share of units whose covariates leave the training range. Maps use the official
GS(2019)1822 basemap and render the national boundary and the South China Sea
nine-dash line.

### Directional and displacement analysis

New-record coordinates were referenced to each species' BirdLife (Birds of the World
2024) historical range polygon. For every record we computed the bearing and
great-circle distance to the range centroid and to the nearest point on the range
edge, retaining 851 records (564 species) with usable polygons. Bearings were binned
into 16 compass sectors for display and summarised with circular statistics: the mean
resultant direction, the concentration *R* and a Rayleigh test of angular uniformity.
Records falling inside the historical range have an edge distance of zero and were
retained and reported separately rather than discarded. Displacement distributions
were compared across taxonomic orders with at least 25 records.

### Reproducibility

Analyses used R (glmmTMB, xgboost, ranger, blockCV, terra, sf, arrow) with versions
pinned in a lockfile; every reported value is traceable to a named result table.
Code and derived outputs are at
<https://github.com/dingchenchen6/bird-hazard-complete-riskset> (complete-risk-set
analysis suite) and <https://github.com/dingchenchen6/bird-hazard-effort-visibility>
(manuscript and main analysis).

## Data availability

Derived covariate tables, model outputs and data dictionaries are in the project
repositories; non-sensitive products will be deposited at Zenodo on acceptance, with
licensed occurrence data shared under source-specific terms. Climate products are
public (WorldClim, CHELSA, CRU TS, CMIP6); the GS(2019)1822 basemap is from the
Ministry of Natural Resources of China.

## Code availability

Version-controlled R code, the renv lockfile and the analysis pipeline are at the
repositories above and will be archived at Zenodo on acceptance.

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

## Figures

**Fig. 1 | The interaction is threshold-invariant and runs through a thermal
channel.**
**a**, Model ladder M0–M5 on the three complete risk sets under the headline effort
proxy; bars show ΔAIC relative to the best model. M4 (climate × effort) is best at
every threshold, beating the additive model M3 by 51–56 and the effort-as-offset
model M5 by 79–81; effort-only (M1, 62–68) far outperforms climate-only (M2,
317–337). **b**, Interaction hazard ratio for all nine climate × five effort proxy
combinations at threshold 50, with threshold robustness for the three leading
climate proxies. Thermal anomaly/gradient (1.12–1.43) and precipitation anomaly
(1.05–1.23) interact with effort; climate-velocity, exposure and warming-rate
proxies do not (0.88–1.06). Collinear proxy pairs are annotated.
*Source: `Fig_A_model_ladder`, `Fig_B_proxy_matrix`.*

**Fig. 2 | Generality across life histories, and the roles of effort, climate and
their interaction.**
**a**, Interaction hazard ratio by migratory strategy at all three thresholds; every
group is positive and significant, with long-distance migrants strongest
(1.325–1.354). **b**, Partition of the deviance explained by M4: unique effort
~79 %, interaction 16 %, unique climate ~4 %, shared ~1 %, stable across thresholds.
**c**, Random-forest permutation and gradient-boosting importance on the same data,
which rank climate above effort — the opposite of **b** — because tree ensembles
lack the species and province random effects that condition the hazard model.
*Source: `Fig_C_migratory`, `Fig_D_variance_decomposition`,
`Fig_E_relative_importance`.*

**Fig. 3 | Ecological-grain risk surface.** Province-fitted hazard relationship
applied to grid-native CRU temperature anomaly and merged survey effort at 100 km
(and 50 km), for the present and under CMIP6 warming with SSP-differentiated effort
growth; cool-to-warm palette, official basemap with the nine-dash line. Projected
risk rises where exposure and observation intensify jointly.
*Source: `Figure_8_grid100_native_plugin_hazard`,
`Figure_8b_grid50_native_plugin_hazard`.*

**Fig. 4 | Two model classes bound the credible forecast horizon.**
**a**, Mean province hazard projected by the mechanistic hazard model and by
gradient boosting, for both emissions pathways and all three thresholds (log scale).
The mechanistic model rises to 13-fold (SSP2-4.5) and 131-fold (SSP5-8.5) by 2080;
the machine-learning model saturates near 2.4-fold. **b**, Percentage of provinces
whose projected climate leaves the training range — 18 % by 2050 and 55 % by 2080
under SSP2-4.5, 48 % and 91 % under SSP5-8.5 — explaining where and why the two
classes diverge. **c**, Model accuracy on the complete risk set by validation regime:
machine learning wins at interpolation (0.774 vs 0.694) but the classes converge for
temporal forecasting (0.615–0.630) and spatial extrapolation (0.613–0.616).
*Source: `Fig_F_future_mech_vs_ml`, `Fig_H_model_accuracy`; province projection maps
in `Fig_G_future_maps_mech` and `Fig_G_future_maps_ml`.*

**Fig. 5 | New records are displaced polewards, with taxon-specific directions.**
**a**, Sixteen-sector wind roses of record bearing relative to the historical range
centroid (mean 55.7°, *R* = 0.270) and the nearest range edge (mean 21.1°,
*R* = 0.355); both Rayleigh *P* < 0.001, *n* = 851 records. **b**, Composite roses by
taxonomic order: passerines east-north-east, raptors north-east, shorebirds
north-west and waterfowl south-west. **c**, Polar bearing × distance view with mean
vectors. **d**, Displacement distributions by order and cumulative distributions;
median 1,273 km from the centroid and 525 km beyond the nearest edge, with 81 of 851
records (9.5 %) inside the historical range.
*Source: `Fig_I_windrose_overall`, `Fig_I_windrose_by_order`,
`Fig_J_polar_bearing_distance`, `Fig_K_distance_distributions`.*

### Extended Data

| Item | Content |
|---|---|
| ED Table 1 | Risk-set construction summary across thresholds (`table_riskset_construction_summary.csv`) |
| ED Table 2 | Full model ladder, all thresholds × effort proxies (`table_A_model_ladder.csv`) |
| ED Table 3 | Climate × effort proxy matrix (`table_B_proxy_matrix.csv`) |
| ED Table 4 | Migratory stratification (`table_C_migratory.csv`) |
| ED Table 5 | Deviance decomposition and marginal/conditional *R*² (`table_D_*.csv`) |
| ED Table 6 | Relative importance, random forest and gradient boosting (`table_E_relative_importance.csv`) |
| ED Table 7 | Future projections and mechanistic–ML concordance (`table_F_*.csv`) |
| ED Table 8 | Climate-proxy collinearity audit (`table_G_climate_proxy_collinearity.csv`) |
| ED Fig. 1 | Endogeneity: lagged-effort interaction |
| ED Fig. 2 | Predictive validation across interpolation, temporal and spatial regimes |
| ED Fig. 3 | Multi-scale attenuation (province → prefecture → county) |
| ED Table 9 | Model accuracy by validation regime, including folds (`table_H_model_accuracy*.csv`) |
| ED Table 10 | Directional summary: mean bearings, concentration and Rayleigh tests (`table_I_directional_summary.csv`) |
| ED Table 11 | Displacement distances overall and by order (`table_J_distance_summary.csv`) |

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
