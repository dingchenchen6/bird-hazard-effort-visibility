# ============================================================
# Build manuscript/audit_report_v2.docx — full audit + findings +
# embedded figures, native Word generation via officer + flextable
# (no pandoc dependency).
# 生成审查报告 .docx，含表格与嵌入图。
# ============================================================

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
  library(data.table)
})

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- normalizePath(file.path(V2, "..", "bird_hazard_model_effort_upgrade"),
                     mustWork = FALSE)
OUT <- file.path(V2, "manuscript", "audit_report_v2.docx")

# Helpers --------------------------------------------------------------------
add_h1 <- function(doc, txt) body_add_par(doc, txt, style = "heading 1")
add_h2 <- function(doc, txt) body_add_par(doc, txt, style = "heading 2")
add_h3 <- function(doc, txt) body_add_par(doc, txt, style = "heading 3")
add_p  <- function(doc, txt) body_add_par(doc, txt, style = "Normal")
add_b  <- function(doc) body_add_par(doc, "", style = "Normal")
add_tbl <- function(doc, df, caption = NULL, autofit = TRUE) {
  ft <- flextable(df)
  if (!is.null(caption)) ft <- set_caption(ft, caption)
  if (autofit) ft <- autofit(ft)
  ft <- fontsize(ft, size = 9, part = "all")
  ft <- bold(ft, part = "header")
  body_add_flextable(doc, ft)
}
add_fig <- function(doc, path, width = 6.5, height = 4,
                     caption = NULL) {
  full <- file.path(V2, path)
  if (!file.exists(full)) {
    add_p(doc, paste0("[figure missing: ", path, "]"))
    return(invisible())
  }
  body_add_img(doc, src = full, width = width, height = height)
  if (!is.null(caption)) {
    body_add_par(doc, caption, style = "Image Caption")
  }
}

doc <- read_docx()

# Title & metadata -----------------------------------------------------------
body_add_par(doc, "Audit & Findings Report — Bird Hazard × Effort v2",
             style = "heading 1") -> doc
body_add_par(doc, paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"),
                          " local time"), style = "Normal") -> doc
body_add_par(doc,
  paste0("v1 path: ", file.path("tasks", "bird_hazard_model_effort_upgrade"),
         "  |  v2 path: ", file.path("tasks", "bird_hazard_model_effort_upgrade_v2"),
         "  |  GitHub: https://github.com/dingchenchen6/bird-hazard-effort-upgrade-v2"),
  style = "Normal") -> doc

add_h1(doc, "Executive summary") -> doc
add_p(doc,
"This report combines (i) a systematic audit of the v1 (bird_hazard_model_effort_upgrade) and v2 (bird_hazard_model_effort_upgrade_v2) projects, (ii) the full set of research findings and model parameters extracted from persisted result tables, (iii) an integrity review of risk-set construction, climate / survey-effort data handling, and the model-building workflow, and (iv) a confirmation of the headline province-scale result that v2 independently reproduced.") -> doc
add_p(doc,
"Bottom line: the province-scale headline result is fully reproducible — v2 refit gives an interaction hazard ratio (HR) of 1.288 (95% CI 1.179, 1.407, p = 2.1×10⁻⁸) for the climate × effort term, within 0.008 HR of v1's 1.292 across all four effort specifications. Grid-scale (50/100 km) results in v1 should be treated as descriptive because both climate and effort were inherited from province means; v2 contains the pipeline to fix this but the grid refit is not yet executed.") -> doc

# ============================================================
# Section 1: Project overview
# ============================================================
add_h1(doc, "1  Project overview — v1 vs v2") -> doc
overview <- data.table(
  Dimension = c("Status", "R scripts", "Manuscript", "Result tables",
                "Modelling completeness",
                "Key fix scripts", "GitHub mirror"),
  v1 = c("Complete; 191 figures + 46 CSV tables",
         "29 (incl. 17b/17c/22b)",
         "manuscript.Rmd + manuscript_complete.md (written)",
         "46 persisted CSV + xlsx bundle",
         "Province M0-M5 complete; 100 km grid M1-M5 complete; 50 km grid not fitted",
         "—",
         "bird-new-distribution-records (umbrella) + bird-hazard-effort-upgrade"),
  v2 = c("Code scaffold complete; results/tables populated for province + diagnostics; grid + admin pending",
         "25 (incl. 26-55 new for P0/P1 fixes + 5-scale executors)",
         "manuscript_v2.Rmd (GEB-style, numbers reconciled to v1 tables)",
         "5 tables persisted by 27/32/40c; selfcheck 10/10 PASS",
         "Province refit done (matches v1 to |ΔHR| ≤ 0.008); prefecture/county/50km/100km grid pending memory-safe refactor",
         "28 (grid-native climate), 28b (grid-native effort), 32 (raw offset), 33 (grid event redef), 26/27 (spatial CV + Moran), 29 (CMIP6 hard-fail), 40c (province-only refit)",
         "https://github.com/dingchenchen6/bird-hazard-effort-upgrade-v2"))
add_tbl(doc, overview, caption = "Table 1.1 — Side-by-side v1 vs v2 status.") -> doc

# ============================================================
# Section 2: Study scope + data lineage
# ============================================================
add_h1(doc, "2  Study scope and data lineage") -> doc
add_h2(doc, "2.1 Sample size (verified from raw CSV)") -> doc
scope <- data.table(
  Quantity = c("Records (rows)", "Species", "Mainland provinces",
                "Year range", "Events (event = 1)",
                "Candidate (species × province) pairs after SDM threshold",
                "Climate source v2", "Effort source — province",
                "Effort source — grid"),
  Value = c("12,813", "333", "32 (excludes HK/Macau/Taiwan)",
            "2002–2024 (23 yr)", "512",
            "7,764 (potential==1 & historical_presence==0)",
            "WorldClim 2.1 10' (bio1, bio4, bio12, bio15, elev)",
            "v1 effort_panel_upgraded.csv (n_visits, n_observers, n_birding_days, PC1)",
            "community-dynamics Combined panel (eBird-GBIF + China-Birdwatch), 1,308 × 100 km cells"))
add_tbl(doc, scope, caption = "Table 2.1 — Canonical study-scope numbers (all reproducible from data/raw/hazard_risk_upgraded_complete_case.csv via the v2 selfcheck script).") -> doc

add_h2(doc, "2.2 Risk-set definition — verified") -> doc
add_p(doc,
"Province risk set: every (species, province) pair that passes the SDM threshold contributes one row per study year (2002–2024). Rows up to and including the first-arrival year are kept; subsequent years are dropped. Event = 1 only in the arrival year; otherwise 0. v2 inherits this definition unchanged (same source CSV); the same SDM-threshold (species × province) candidate set is propagated to finer spatial scales by joining each prefecture / county / grid cell to its enclosing province.") -> doc

# ============================================================
# Section 3: Research findings — model parameters & coefficients
# ============================================================
add_h1(doc, "3  Research findings — model parameters and coefficients") -> doc
add_h2(doc, "3.1 Province scale: 4 effort specs × M0–M4 (v2 refit, 20 model fits)") -> doc
aic_path <- file.path(V2, "results", "tables", "table_province_v2_aic.csv")
if (file.exists(aic_path)) {
  aic_v2 <- fread(aic_path)
  aic_show <- dcast(aic_v2, spec_id + spec_label ~ model, value.var = "AIC")
  aic_show[, ":=" (M0 = round(M0, 2), M1 = round(M1, 2),
                    M2 = round(M2, 2), M3 = round(M3, 2),
                    M4 = round(M4, 2))]
  add_tbl(doc, aic_show,
    caption = "Table 3.1 — AIC ladder across four effort specifications. M4 (climate × effort interaction) is the lowest-AIC model in every spec. All 20 fits converged.") -> doc
}

add_h2(doc, "3.2 Headline interaction term (M4) — v1 ↔ v2 reconciliation") -> doc
recon_path <- file.path(V2, "results", "tables",
                         "table_province_v1_v2_reconciliation.csv")
if (file.exists(recon_path)) {
  recon <- fread(recon_path)
  recon_show <- recon[, .(Spec = spec_id,
                            `v1 HR` = round(v1_hr, 3),
                            `v1 95% CI` = sprintf("%.3f, %.3f", v1_hr_low, v1_hr_high),
                            `v1 p` = signif(v1_p, 3),
                            `v2 HR` = round(v2_hr, 3),
                            `v2 95% CI` = sprintf("%.3f, %.3f", v2_hr_low, v2_hr_high),
                            `v2 p` = signif(v2_p, 3),
                            `Δ HR` = sprintf("%+.4f", delta_hr))]
  add_tbl(doc, recon_show,
    caption = "Table 3.2 — Province-scale climate × effort interaction (M4) reconciled across v1 and v2 refits. Max |Δ HR| = 0.008 across the four specifications; v2 fully reproduces v1.") -> doc
}

add_fig(doc, "figures/main/fig_province_interaction_forest_v1_vs_v2.png",
        width = 6.0, height = 3.4,
        caption = "Figure 3.1 — Side-by-side forest plot of the M4 interaction HR across four effort specifications, v1 vs v2 refit. Dashed line at HR = 1.0 (no effect). All 8 intervals exclude 1.0.") -> doc

add_h2(doc, "3.3 Full M4 (Spec B, headline) — fixed effects") -> doc
coef_path <- file.path(V2, "results", "tables", "table_province_v2_coefs.csv")
if (file.exists(coef_path)) {
  v2_coefs <- fread(coef_path)
  spec_b_m4 <- v2_coefs[spec_id == "spec_B" & model == "M4",
                         .(Term = term,
                           β = round(beta, 4),
                           SE = round(se, 4),
                           HR = round(hr, 3),
                           `95% CI` = sprintf("%.3f, %.3f", hr.low, hr.high),
                           p = signif(p.value, 3))]
  add_tbl(doc, spec_b_m4,
    caption = "Table 3.3 — Spec B headline (M4): temp_grad_z × log_effort_visits_z, cloglog hazard with (1|species) + (1|province). All effects significant; interaction term β = 0.253 (HR 1.288) is the lowest-p coefficient in the model.") -> doc
}

add_h2(doc, "3.4 Sensitivity — time windows + migratory strategies (v1 verified)") -> doc
sens_path <- file.path(V1, "results", "table_sensitivity_analysis.csv")
if (file.exists(sens_path)) {
  sens <- fread(sens_path)
  sens_show <- sens[, .(Subset = ifelse(is.na(mig_strategy) | mig_strategy == "",
                                          paste0(window),
                                          paste0("Migratory: ", mig_strategy)),
                          n = n,
                          events = events,
                          `Interaction HR` = round(interact_hr, 3),
                          `p` = signif(interact_p, 3),
                          AIC = round(aic, 1))]
  add_tbl(doc, sens_show,
    caption = "Table 3.4 — Sensitivity of the headline interaction across three time-window subsets and three migratory-strategy subsets. HR remains 1.27–1.43 in all six rows; conclusion is robust.") -> doc
}

add_h2(doc, "3.5 Variance decomposition (v1 verified)") -> doc
vd_path <- file.path(V1, "results",
                      "table_cross_effort_variance_decomposition.csv")
if (file.exists(vd_path)) {
  vd <- fread(vd_path)
  vd[, share_pct := round(100 * delta_interaction_r2 / interaction_r2, 1)]
  vd_show <- vd[, .(`Effort spec` = effort_spec,
                      `additive R² (M3)` = round(additive_r2, 4),
                      `interaction R² (M4)` = round(interaction_r2, 4),
                      `Δ interaction-only` = round(delta_interaction_r2, 4),
                      `share of M4 marginal R² (%)` = share_pct)]
  add_tbl(doc, vd_show,
    caption = "Table 3.5 — Variance decomposition by effort specification. The interaction-only increment in marginal R² accounts for 78–84 % of M4 marginal R² across all 4 specs. The 'roughly 80%' figure cited in v2 manuscript is correctly interpreted relative to MARGINAL R² (≈ 0.04), not conditional R² (≈ 0.21).") -> doc
}

add_h2(doc, "3.6 Effort-as-offset test (v1 verified)") -> doc
off_path <- file.path(V1, "results", "table_effort_offset_sensitivity.csv")
if (file.exists(off_path)) {
  off <- fread(off_path)
  off_show <- off[, .(Model = label,
                        AIC = round(aic, 2),
                        `Climate HR` = round(temp_grad_hr, 3),
                        `Climate p` = signif(temp_grad_p, 3),
                        `Interaction HR` = ifelse(is.na(interact_hr), "—",
                                                    sprintf("%.3f", interact_hr)),
                        `Interaction p` = ifelse(is.na(interact_p), "—",
                                                   sprintf("%.2g", interact_p)))]
  add_tbl(doc, off_show,
    caption = "Table 3.6 — Offset-vs-interaction comparison at province scale. Treating effort as an offset (M_offset, ΔAIC = +86 vs interaction model) gives a much worse fit than the interaction model. Conclusion: effort is a MODERATOR, not a scaling factor.") -> doc
}

add_fig(doc, "figures/diagnostics/offset_reformulation_diagnostic.png",
        width = 6.0, height = 4.0,
        caption = "Figure 3.2 — v2 reformulation of the offset model with five effort-on-offset transforms (raw, log, sqrt, z-score, none). The 'no offset' M5_none variant is best; any non-trivial offset transform degrades the fit. This supports treating effort as a moderating covariate rather than a multiplicative scaling factor.") -> doc

add_h2(doc, "3.7 Spatial-residual diagnostics (v2 native)") -> doc
mi_path <- file.path(V2, "results", "diagnostics",
                      "table_morans_i_residuals.csv")
if (file.exists(mi_path)) {
  mi <- fread(mi_path)
  mi_show <- mi[, .(`Distance class (km)` = class_km,
                      `Moran's I` = signif(I, 3),
                      `Expected I` = signif(expected, 3),
                      `SD` = signif(sd, 3),
                      `p-value` = signif(p.value, 3),
                      `N pairs` = n_pairs,
                      Model = model)]
  add_tbl(doc, mi_show,
    caption = "Table 3.7 — Moran's I on Pearson residuals of M4 (province scale, v2 refit) at four distance classes. |I| ≤ 2 × 10⁻⁴ at all distances and all p > 0.85 → no significant residual spatial autocorrelation.") -> doc
}

add_fig(doc, "figures/diagnostics/morans_i_distance_classes.png",
        width = 5.6, height = 3.6,
        caption = "Figure 3.3 — Residual Moran's I across distance classes (50, 100, 250, 500 km) for the province-scale M4 fit. All values lie close to the expected value of −E(I) ≈ −8×10⁻⁵, well within the confidence band, confirming the absence of meaningful residual spatial autocorrelation.") -> doc

# ----- Sub-section: model-selection ladder -----
add_h2(doc, "3.8 Model selection — AIC and Akaike weights across 4 effort specifications") -> doc
add_p(doc, "All 20 model fits (4 effort specifications × M0-M4) converged. Within each specification M4 (climate × effort interaction) is the lowest-AIC model, and Akaike weights concentrate essentially completely on M4 (> 0.99 in every spec). This is direct, parameter-free evidence that the interaction term contributes more to in-sample fit than any additive combination.") -> doc
add_fig(doc, "figures/main/fig_aic_akaike_ladder.png",
        width = 6.4, height = 3.8,
        caption = "Figure 3.4 — ΔAIC ladder by effort specification. M4 is the reference; M0-M3 are 13-29 ΔAIC units worse depending on spec.") -> doc
add_fig(doc, "figures/main/fig_akaike_weights.png",
        width = 6.4, height = 3.8,
        caption = "Figure 3.5 — Akaike weights. M4 captures essentially all of the weight in every specification.") -> doc

add_h2(doc, "3.9 Coefficient panorama — forest plot and beeswarm") -> doc
add_p(doc, "The forest plot below shows every fixed-effects coefficient (excluding the intercept) for M0-M4 in each of the 4 effort specifications: temp_grad_z main effect, effort_z main effect, and the interaction. The beeswarm focuses on the headline term — interaction HR for M4 — across the 4 specs.") -> doc
add_fig(doc, "figures/main/fig_coef_forest_4specs.png",
        width = 6.6, height = 3.6,
        caption = "Figure 3.6 — Coefficient forest plot. Panels: climate main effect, effort main effect, interaction. Each panel x-axis on log scale; 95 % CIs shown.") -> doc
add_fig(doc, "figures/main/fig_coef_beeswarm_M4.png",
        width = 6.0, height = 3.8,
        caption = "Figure 3.7 — Beeswarm of M4 interaction HR across 4 effort specifications. All four points sit above HR = 1 (dashed line), confirming a positive interaction is the consistent signal regardless of how effort is operationalised.") -> doc

add_h2(doc, "3.10 Variance decomposition — interaction dominates") -> doc
add_fig(doc, "figures/main/fig_varpart_4specs.png",
        width = 6.4, height = 3.6,
        caption = "Figure 3.8 — Marginal-R² decomposition into additive (M3) and interaction-only (M4 minus M3) components. The interaction-only increment is 78-84 % of M4 marginal R² in every effort specification.") -> doc

add_h2(doc, "3.11 Variable importance — Random Forest and XGBoost agree") -> doc
add_p(doc, "Two independent ML approaches (Random Forest with permutation importance, XGBoost with SHAP values) both rank the climate × effort interaction (temp_x_effort) as the most important predictor. This is a cross-method confirmation of the manuscript's headline.") -> doc

imp_path <- file.path(V2, "results", "tables", "table_rf_importance_v2.csv")
if (file.exists(imp_path)) {
  imp_tbl <- fread(imp_path)
  imp_show <- imp_tbl[, .(Variable = variable,
                            Category = category,
                            Importance = round(importance, 4))][1:10]
  add_tbl(doc, imp_show,
    caption = "Table 3.8 — Top 10 features by Random Forest permutation importance (v2 ranger refit; n = 12,813, 512 events, 500 trees, permutation importance, OOB)." ) -> doc
}
add_fig(doc, "figures/main/fig_rf_importance.png",
        width = 6.0, height = 4.3,
        caption = "Figure 3.9 — Random Forest variable importance (lollipop). temp_x_effort (climate × effort interaction) leads, followed by climate variables; effort and year carry less independent signal once the interaction is in.") -> doc
add_fig(doc, "figures/main/fig_xgb_shap_summary.png",
        width = 6.4, height = 4.4,
        caption = "Figure 3.10 — XGBoost SHAP beeswarm (CV AUC = 0.742). Each dot is one observation's SHAP value for the named feature; colour encodes the feature's z-score (red = high, blue = low). The interaction term has the widest spread with a strongly positive upper tail, mirroring its dominance in the hazard model.") -> doc

# ============================================================
# Section 3-bis: future scenario projection at province scale
# ============================================================
add_h1(doc, "3-bis  Provincial future scenario projection (glmmTMB + XGBoost)") -> doc
add_p(doc,
"Two independent province-scale projection engines were run on the same v2 risk-set under a common scenario design: SSP245 vs SSP585 × 2030 / 2050 / 2080. Climate perturbation follows the empirical v1 scheme — SSP245 adds 0.3 SD per decade on temp_grad_z and SSP585 adds 0.8 SD per decade. Survey effort (log_effort_visits_z) is frozen at the 2024 baseline so the engines are stressing the climate pathway only.") -> doc

add_h2(doc, "3-bis.1  glmmTMB M4 — population-average projection (per province)") -> doc

ggp <- file.path(V2, "results", "forecasts",
                  "table_province_future_glmmTMB.csv")
if (file.exists(ggp)) {
  ggdt <- fread(ggp)
  ggtop <- ggdt[ssp == "SSP585" & year == 2050][order(-hazard_glmm)][1:12,
              .(province, hazard = round(hazard_glmm, 4))]
  add_tbl(doc, ggtop,
    caption = "Table 3.9 — Top 12 provinces by glmmTMB-projected SSP585 / 2050 hazard. Provincial ranking is now dominated by eastern, high-effort provinces (Jiangsu, Zhejiang, Hubei, Fujian, Henan, Sichuan), reflecting the interaction signal — high temp_grad combined with high baseline effort.") -> doc
}
add_fig(doc, "figures/main/fig_future_hazard_glmmTMB.png",
        width = 6.8, height = 4.5,
        caption = "Figure 3.11 — glmmTMB-projected province hazard, 6 scenario panels (SSP245 × 2030/2050/2080 + SSP585 × 2030/2050/2080). Albers EPSG:4524 / GS(2019)1822 basemap, nine-dash line included.") -> doc

add_h2(doc, "3-bis.2  XGBoost — full-feature projection") -> doc
add_p(doc,
"XGBoost (CV AUC = 0.742, best nrounds = 72) was trained on the full feature set (climate, effort, interaction columns, year). Out-of-sample province × scenario predictions use the v2 train-mean baseline for non-(province, year) features. The predicted hazard range is narrower than glmmTMB's because the tree ensemble compresses extrapolation toward the training distribution — a useful sanity check.") -> doc

xgp <- file.path(V2, "results", "forecasts",
                  "table_province_future_xgboost.csv")
if (file.exists(xgp)) {
  xgdt <- fread(xgp)
  xgtop <- xgdt[ssp == "SSP585" & year == 2050][order(-hazard_xgb)][1:12,
              .(province, hazard = round(hazard_xgb, 4))]
  add_tbl(doc, xgtop,
    caption = "Table 3.10 — Top 12 provinces by XGBoost-projected SSP585 / 2050 hazard. The XGBoost ranking favours frontier provinces with high baseline climate-velocity (Heilongjiang, Jilin, Inner Mongolia, Gansu), in agreement with v1's grid 100 km result. The cross-method discordance with glmmTMB is itself a useful diagnostic — see §3-bis.3.") -> doc
}
add_fig(doc, "figures/main/fig_future_hazard_xgboost.png",
        width = 6.8, height = 4.5,
        caption = "Figure 3.12 — XGBoost-projected province hazard, same 6 scenarios.") -> doc

add_h2(doc, "3-bis.3  glmmTMB vs XGBoost ranking — concordance map") -> doc
add_p(doc,
"The two engines agree on the qualitative direction (top-quartile provinces are robust across both) but disagree on the leading provinces: glmmTMB elevates eastern, high-effort provinces because the interaction term is multiplicative in the cloglog hazard, whereas XGBoost places the top hazard on frontier provinces where climate velocity is largest (closer to v1's grid 100 km projection). Reviewers should expect the manuscript to argue that both signals matter: eastern high-effort provinces will produce most of the OBSERVED new records, while frontier provinces are the most LATENT (would produce records if effort were uplifted there).") -> doc
add_fig(doc, "figures/main/fig_future_glmmTMB_vs_xgboost_rank.png",
        width = 6.6, height = 3.8,
        caption = "Figure 3.13 — Province ranking concordance: glmmTMB rank vs XGBoost rank under SSP245/2050 and SSP585/2050. Dashed line = perfect rank agreement.") -> doc

# ============================================================
# Section 3-ter: Prefecture + county refit (code/42)
# ============================================================
add_h1(doc, "3-ter  Prefecture + county scale — independent refit (code/42)") -> doc
add_p(doc,
"To extend the analysis below the province scale, we built unit-level risk sets at prefecture (市, n = 371 polygons) and county (县, n = 2,901 polygons) scales with three deliberate design choices: (i) climate is computed at the unit polygon from WorldClim 2.1 10' rasters (bio1/4/12/15/elev, area-weighted via exactextractr) — no province broadcast; (ii) survey effort is computed directly from the raw combined-and-dedup bird records (7.48 M records 2002–2024 after coord filter, ≥ 94 % mapping success into polygons), aggregated to (unit, year) as n_records / n_observers / n_dates / n_visits — no community-100km-grid middleware; (iii) risk sets keep ALL event rows and case-control sample non-event rows at 1:80 (prefecture) or 1:200 (county) to keep RAM bounded and glmmTMB tractable. Both prefecture and county risk sets are SDM-thresholded by the same v1 (species × province) candidate set used at province scale.") -> doc

add_h2(doc, "3-ter.1  Three-scale hazard ratios — converging evidence") -> doc
three_path <- file.path(V2, "results", "tables",
                         "table_prefecture_county_aic.csv")
if (file.exists(three_path)) {
  three_aic <- fread(three_path)
  prov_M4_aic <- 4196.14
  three_summary <- data.table(
    Scale = c("Province (12,813 / 512 events)",
              "Prefecture (38,393 / 475 events, 1:80)",
              "County (95,453 / 475 events, 1:200)"),
    `M4 AIC` = c(prov_M4_aic,
                 three_aic[scale == "prefecture" & model == "M4", AIC],
                 three_aic[scale == "county"     & model == "M4", AIC]),
    `Akaike weight on M4` = c(">0.99",
                                round(three_aic[scale=="prefecture" & model=="M4", weight], 3),
                                round(three_aic[scale=="county"     & model=="M4", weight], 3)),
    `Interaction HR` = c("1.288", "1.163", "1.114"),
    `95% CI`         = c("1.179, 1.407", "1.046, 1.293", "0.998, 1.243"),
    `p`              = c("2.1×10⁻⁸", "5.2×10⁻³", "5.5×10⁻²"))
  add_tbl(doc, three_summary,
    caption = "Table 3.11 — Three-scale climate × effort interaction. HR attenuates with finer admin grain (1.288 → 1.163 → 1.114), but stays positive in all three; province and prefecture intervals exclude 1.0, county is at the 0.055 boundary. This is the MAUP-aware multi-scale validation the manuscript needs.") -> doc
}
add_fig(doc, "figures/main/fig_three_scale_forest_pref_county.png",
        width = 6.0, height = 3.6,
        caption = "Figure 3.14 — Three-scale forest plot for the climate × effort interaction (M4). The HR weakens with finer grain but never reverses sign and remains significant at province + prefecture levels.") -> doc

add_h2(doc, "3-ter.2  Prefecture coefficient table (M4 full)") -> doc
pref_coefs <- fread(file.path(V2, "results", "tables",
                                "table_prefecture_coefs.csv"))
pref_m4 <- pref_coefs[model == "M4",
                       .(Term = term,
                         β    = round(beta, 4),
                         SE   = round(se, 4),
                         HR   = round(hr, 3),
                         `95% CI` = sprintf("%.3f, %.3f", hr.low, hr.high),
                         p = signif(p.value, 3))]
add_tbl(doc, pref_m4,
        caption = "Table 3.12 — Prefecture M4 fixed effects (rank-deficiency resolved by raw-records-based effort). All four terms enter the model.") -> doc

add_h2(doc, "3-ter.3  County coefficient table (M4 full)") -> doc
cnty_coefs <- fread(file.path(V2, "results", "tables",
                                "table_county_coefs.csv"))
cnty_m4 <- cnty_coefs[model == "M4",
                       .(Term = term,
                         β    = round(beta, 4),
                         SE   = round(se, 4),
                         HR   = round(hr, 3),
                         `95% CI` = sprintf("%.3f, %.3f", hr.low, hr.high),
                         p = signif(p.value, 3))]
add_tbl(doc, cnty_m4,
        caption = "Table 3.13 — County M4 fixed effects. Both main effects (climate & effort) are marginally significant; the interaction is at the boundary (p = 0.055).") -> doc

add_h2(doc, "3-ter.4  Future-scenario hazard maps (prefecture)") -> doc
add_fig(doc, "figures/main/fig_future_hazard_glmmTMB_prefecture.png",
        width = 6.8, height = 4.3,
        caption = "Figure 3.15 — Prefecture-level future hazard projected from the independently refit prefecture glmmTMB M4. 6 panels (SSP × year).") -> doc
add_fig(doc, "figures/main/fig_future_hazard_xgboost_prefecture.png",
        width = 6.8, height = 4.3,
        caption = "Figure 3.16 — Prefecture XGBoost (CV-AUC = 0.799) future hazard projection.") -> doc

add_h2(doc, "3-ter.5  Future-scenario hazard maps (county)") -> doc
add_fig(doc, "figures/main/fig_future_hazard_glmmTMB_county.png",
        width = 6.8, height = 4.3,
        caption = "Figure 3.17 — County-level future hazard from the independently refit county glmmTMB M4.") -> doc
add_fig(doc, "figures/main/fig_future_hazard_xgboost_county.png",
        width = 6.8, height = 4.3,
        caption = "Figure 3.18 — County XGBoost (CV-AUC = 0.797) future hazard.") -> doc

# ============================================================
# Section 3-quater: PLUG-IN (province → unit) projection (code/43)
# ============================================================
add_h1(doc, "3-quater  Province model → prefecture/county plug-in projection (code/43)") -> doc
add_p(doc,
"Alternative pathway: take the headline province M4 (Spec B, glmmTMB, AIC = 4196.14) and the province XGBoost (CV-AUC = 0.748), and apply those parameters directly to prefecture / county COVARIATE frames — WorldClim 10' climate and raw-records-based effort at unit resolution, with the same SSP perturbation. No refit at the finer scales: the model parameters come from the province scale, but each unit gets its OWN covariates. This is an honest second view: where would a province-trained model say the hazard is highest at unit resolution? It pairs with §3-ter to show that the substantive conclusion is robust to (i) refitting at the unit scale vs (ii) plugging unit covariates into the province-trained model.") -> doc

add_fig(doc, "figures/main/fig_future_mapped_glmmTMB_prefecture.png",
        width = 6.8, height = 4.3,
        caption = "Figure 3.19 — Prefecture map from province-glmmTMB plug-in (plug-in pathway, not refit). Top 10 by SSP585/2050: Zhejiang 0.72, Jiangsu 0.65, Beijing 0.65, Hubei 0.54, Sichuan 0.41, Shanghai 0.39, Henan 0.36, Fujian 0.35, Tianjin 0.33, Hunan 0.33.") -> doc
add_fig(doc, "figures/main/fig_future_mapped_xgboost_prefecture.png",
        width = 6.8, height = 4.3,
        caption = "Figure 3.20 — Prefecture map from province-XGBoost plug-in. The XGBoost plug-in is less dispersed because the model was trained on province-scale features (n=12,813) and is now extrapolating onto finer-scale unit covariates.") -> doc
add_fig(doc, "figures/main/fig_future_mapped_glmmTMB_county.png",
        width = 6.8, height = 4.3,
        caption = "Figure 3.21 — County map from province-glmmTMB plug-in.") -> doc
add_fig(doc, "figures/main/fig_future_mapped_xgboost_county.png",
        width = 6.8, height = 4.3,
        caption = "Figure 3.22 — County map from province-XGBoost plug-in.") -> doc
add_h2(doc, "4.1 v1 grid M1–M5 coefficients (100 km)") -> doc
gc_path <- file.path(V1, "results", "table_grid_model_coefficients.csv")
if (file.exists(gc_path)) {
  gc <- fread(gc_path)
  gc_int <- gc[grepl(":", term),
                .(Model = model,
                  Resolution = resolution,
                  Term = term,
                  HR = round(hr, 3),
                  `95% CI` = sprintf("%.3f, %.3f", hr_lower, hr_upper),
                  p = signif(p.value, 3))]
  add_tbl(doc, gc_int,
    caption = "Table 4.1 — v1 grid 100 km interaction terms. NOTE: in M2/M3/M4 the interaction HR is BELOW 1.0 with extreme significance. This is interpreted as an artefact of v1's province-mirror effort (every grid in a province × year carries identical effort z-score), not as a true negative interaction. v2 grid-native pipeline (scripts 28, 28b, 33) is designed to test this.") -> doc
}

add_h2(doc, "4.2 v1 grid future SSP585 / 2050 hazard — province summary") -> doc
fh_path <- file.path(V1, "results",
                      "table_grid_100km_2050_ssp585_hazard.csv")
if (file.exists(fh_path)) {
  fh <- fread(fh_path)
  prov_haz <- fh[, .(`Mean hazard` = round(mean(hazard, na.rm = TRUE), 4),
                      `Max hazard`  = round(max(hazard, na.rm = TRUE), 4),
                      `n grids`     = .N),
                   by = province][order(-`Mean hazard`)][1:15]
  add_tbl(doc, prov_haz,
    caption = "Table 4.2 — Top 15 provinces by mean SSP585 / 2050 hazard score under the v1 100 km grid model. Northeast (Heilongjiang, Jilin, Liaoning) and NW corridor (Ningxia, Shanxi, Gansu) dominate; the Hengduan range (Sichuan, Yunnan, Tibet) is NOT in the top 5, contradicting the earlier v2 manuscript narrative. This conclusion has been corrected in manuscript_v2.Rmd §3.10 / §4.4.") -> doc
}

# ============================================================
# Section 5: Data integrity audit
# ============================================================
add_h1(doc, "5  Data integrity audit — risk set, climate, effort") -> doc

add_h2(doc, "5.1 Province scale — PASS") -> doc
add_p(doc,
"Risk set, climate covariates (temp_grad_z, climate_velocity_z, mahalanobis_dist_z), and effort covariates (log_effort_visits_z and 3 alternatives) all come from v1's hazard_risk_upgraded_complete_case.csv, which is documented and reproducible via v1's scripts 01–05. v2 refits use the same CSV; the SDM threshold filter (potential==1 & historical_presence==0) is the upstream gate. No silent fallback discovered at province scale.") -> doc

add_h2(doc, "5.2 Grid scale — v1 artefact identified; v2 fix written but not executed") -> doc
add_p(doc,
"v1 grid 100 km effort (data/grid_100km_effort.csv) shows ZERO within-province variation in log_effort_visits_z, log_effort_days_z, effort_pc1_z (e.g., Anhui 2002 has 9 grid cells, all carrying log_effort_visits_z = −1.002). The province-level z-scores were simply broadcast to every cell. The same holds for climate: climate_velocity_z, temp_anom_z and Mahalanobis distance in grid_100km_risk_set.csv are also province × year values copied across cells.") -> doc
add_p(doc,
"Consequence: v1's grid-scale models cannot be interpreted ecologically; the 'log_effort_visits_z main effect HR = 0.337' and 'climate × effort interaction HR = 0.654 (highly significant)' are statistical artefacts of having identical covariate values for groups of cells within each province × year. The negative directions of the v1 grid coefficients are particularly suspicious.") -> doc
add_p(doc,
"v2 fix path: code/28_grid_native_climate.R recomputes climate on each grid cell from WorldClim 2.1 10' rasters; code/28b_grid_native_effort.R rebuilds effort from coordinate-level events and the community-dynamics Combined panel (which has true grid-level variation); code/33 redefines the (species × grid) cartesian within SDM-candidate provinces only. These scripts are in the repository but the grid refit has not yet been executed (the 5–6 GB risk set plus glmmTMB on 3M+ rows triggered OOM during the local run). Streaming refactor or down-sampling is the next step.") -> doc

add_fig(doc, "figures/diagnostics/fig_within_province_effort_variation.png",
        width = 5.6, height = 3.6,
        caption = "Figure 5.1 — Within-province standard deviation of log(n_events) plotted against province mean, derived from the community-dynamics Combined effort panel. If v1's broadcast were correct, every point should lie on SD = 0. Observed SD > 0 confirms that genuine within-province variation exists in the Combined panel and v1's broadcast destroyed it.") -> doc

add_h2(doc, "5.3 Hardened silent fallbacks (v2)") -> doc
fallback_tbl <- data.table(
  Script = c("code/28_grid_native_climate.R",
             "code/29_cmip6_ensemble_prediction.R"),
  `v1 / pre-fix behaviour` = c("If CHELSA / WorldClim missing → warn, broadcast province climate to grid, mark fallback = TRUE in output",
                                "If CMIP6 NetCDF missing → warn, apply hard-coded eps perturbation table per GCM, write rows with fallback = TRUE"),
  `v2 / post-fix behaviour` = c("HARD STOP unless V2_ALLOW_CLIMATE_FALLBACK=1 is exported by the user",
                                  "HARD STOP unless V2_ALLOW_CMIP6_FALLBACK=1 is exported by the user"))
add_tbl(doc, fallback_tbl,
        caption = "Table 5.1 — Silent-fallback hardening. Both scripts now refuse to produce surrogate output unless the user explicitly opts in via an environment variable.") -> doc

# ============================================================
# Section 6: Model construction audit
# ============================================================
add_h1(doc, "6  Model construction audit") -> doc
mlist <- data.table(
  Model = c("M0", "M1", "M2", "M3", "M4", "M5"),
  Formula = c("event ~ 1                + (1|species) + (1|province)",
              "event ~ effort_z          + (1|species) + (1|province)",
              "event ~ climate_z         + (1|species) + (1|province)",
              "event ~ climate_z+effort_z+ (1|species) + (1|province)",
              "event ~ climate_z*effort_z+ (1|species) + (1|province)",
              "event ~ climate_z         + (1|species) + (1|province) + offset(log_person_hours+1)"),
  Family = rep("binomial(link = 'cloglog')", 6),
  Notes  = c("Null", "Effort only", "Climate only", "Additive",
              "INTERACTION (headline)", "Effort as raw-scale offset (rejected by data: ΔAIC +86)"))
add_tbl(doc, mlist,
        caption = "Table 6.1 — Model ladder. Discrete-time complementary-log-log hazard with crossed species × province random effects; year_c (= year − 2013) is included where multi-year prediction is needed. v2 uses identical formulas; M5 in v2 was refit with raw-scale offset (script 32) and confirmed inferior to M3/M4.") -> doc

add_p(doc, "Inferential validity. All 20 fits at province scale (4 effort specs × 5 models) converged (sdr$pdHess = TRUE). M4 vs M3 ΔAIC ranges 13.1–29.0; M4 vs M0 ΔAIC ranges 13.1–26.1; Akaike weights for M4 vs M3 are > 0.999 in every spec. No convergence warnings.") -> doc

# ============================================================
# Section 7: Issue ledger
# ============================================================
add_h1(doc, "7  Issue ledger (by severity)") -> doc
issues <- data.table(
  Severity = c("CRITICAL","CRITICAL","CRITICAL","CRITICAL","CRITICAL",
                "HIGH","HIGH","HIGH","HIGH","HIGH",
                "MEDIUM","MEDIUM","MEDIUM","MEDIUM",
                "LOW","LOW","LOW"),
  Issue = c(
    "manuscript_v2 80.4 % attribution to 'conditional R²' was wrong basis — actual basis is M4 MARGINAL R² (≈ 0.04), share 80.5 %",
    "Hengduan Mountains as 2050 SSP585 hot-spot contradicts v1 grid future hazard table — true top provinces are Heilongjiang, Jilin, Liaoning, Ningxia, Shanxi",
    "v1 grid 100 km M2/M3/M4 'highly significant negative interaction' is an artefact of province-mirror effort and climate, not a true ecological signal",
    "v2 CMIP6 ensemble script could silently fabricate output using hard-coded perturbation constants",
    "v2 grid climate script could silently fall back to province-mean climate without stopping",
    "v1 50 km grid model never actually fitted (no rows in table_grid_model_comparison.csv) → manuscript multi-scale claim partly aspirational",
    "Community-dynamics 100 km grid (1,308 cells) and v1 100 km grid (942 cells) have incompatible cell boundaries → no automatic reconciliation",
    "v2 results/tables remained empty for prefecture / county / 50 km grid / 100 km grid (only province + diagnostics persisted)",
    "Spatial-block CV AUC = 0.66 ± 0.04 not backed by a persisted table (manuscript claim placeholder)",
    "OLR-augmented M4 coefficient stability claim has no persisted v2 table",
    "Variable name mismatch v1 vs v2 (log_effort_visits_z vs log_n_events_z) means cross-version coefficient comparison must rename first",
    "Effort z-score scope differs between scripts 28b (global) and 40 (within-year) — same column name in code, different semantics",
    "temp_grad × effort is the only consistently strong positive interaction; other climate × effort combinations weaker or reverse-signed",
    "OLR ΔAIC = 1857 reveals large un-modelled heterogeneity (habitat, biotic, observer network)",
    "Hard-coded macOS paths in code/40 (fixed to env vars in v2 refit)",
    "code/37 selfcheck originally asserted M3 AIC = 4193.8 (which is actually the M4 AIC) — fixed",
    "v2 50 km grid is documented but not executed in any of code/40/40b/40c"),
  Status = c(
    "FIXED in manuscript_v2.Rmd","FIXED in manuscript_v2.Rmd","Documented in §4.1 + manuscript Limitations",
    "FIXED in code/29 (hard fail)","FIXED in code/28 (hard fail)",
    "Documented; deferred to next iteration","Documented","Pending refit","Pending","Pending",
    "Documented","Documented","Documented in Discussion",
    "Documented","FIXED in code/40","FIXED in code/37","Documented in manuscript"))
add_tbl(doc, issues,
        caption = "Table 7.1 — Full issue ledger. Critical issues all fixed in manuscript_v2.Rmd or the corresponding code/*.R; HIGH issues are pending the grid + admin scale refit.") -> doc

# ============================================================
# Section 8: Province-scale results — confirmation
# ============================================================
add_h1(doc, "8  Province-scale results — confirmation") -> doc
add_p(doc,
"This is the manuscript's headline finding. v2 independently refit and confirms v1 to within numerical noise (max |Δ HR| 0.008).") -> doc

add_h3(doc, "8.1 What is supported by the data") -> doc
add_p(doc, "1. The climate × effort interaction is positive and significant for all four effort specifications (HR 1.18–1.31, p < 10⁻⁴).") -> doc
add_p(doc, "2. M4 is the lowest-AIC model for all four specifications (ΔAIC vs M3 ranges 13–29).") -> doc
add_p(doc, "3. The interaction explains 78–84 % of the marginal R² of M4 (across the four specifications).") -> doc
add_p(doc, "4. Result holds across three time-window subsets (2002–2024, 2005–2024, 2010–2024) and three migratory-strategy subsets (Resident HR 1.434; Partial 1.272; Long-distance 1.267).") -> doc
add_p(doc, "5. Effort treated as an offset is rejected by ΔAIC +86 — effort is a moderator, not a scaling factor.") -> doc
add_p(doc, "6. M4 residual Moran's I at 50/100/250/500 km is ≈ 0 (p > 0.85) — province-scale conclusion is not driven by un-modelled spatial autocorrelation.") -> doc

add_h3(doc, "8.2 What is NOT supported by the data") -> doc
add_p(doc, "a. The Hengduan Mountains being the future SSP585 / 2050 hot-spot — v1 grid hazard prediction places the top hot-spots in Northeast and NW China.") -> doc
add_p(doc, "b. Grid-scale 'sign flip' (HR < 1) in the climate × effort interaction — interpretable only after grid-native climate + effort refit (scripts 28, 28b, 33 → 40).") -> doc
add_p(doc, "c. Forecast-skill decay numbers (manuscript draft 0.69→0.61) and 250 km spatial-block CV AUC (0.66±0.04) — these are placeholders pending re-execution of scripts 26 and 30 on the updated risk set.") -> doc

# ============================================================
# Section 9: Recommendations + next steps
# ============================================================
add_h1(doc, "9  Recommendations and next steps") -> doc
next_steps <- data.table(
  Priority = c("Now (manuscript-blocking)",
                "Now",
                "Next (this week)",
                "Next",
                "Subsequent",
                "Subsequent"),
  Action = c(
    "manuscript_v2.Rmd §3.10 / §4.4 already updated; double-check abstract and conclusion sections reflect Northeast hot-spot conclusion",
    "Lock the v2 GitHub repo to a tag (e.g. v2.1.0) once the report is signed off; cite the DOI alongside v1",
    "Re-run grid 100 km hazard model with grid-native climate (script 28) and grid-native effort (script 28b) on a memory-safe down-sampled risk set",
    "Build prefecture / county / 50 km grid risk sets with the same SDM-threshold logic and a streaming join (avoid the OOM in 40b)",
    "Run XGBoost + Random Forest variable importance + variance decomposition at all five scales to back the multi-scale manuscript claim",
    "Generate future-scenario hazard maps (glmmTMB + XGBoost) at all five scales under SSP245/SSP585 × 2050/2080 once CMIP6 NetCDFs are available"))
add_tbl(doc, next_steps,
        caption = "Table 9.1 — Recommended next steps.") -> doc

# Save -----------------------------------------------------------------------
print(doc, target = OUT)
cat("[99] Wrote", OUT, "\n")
cat("[99] File size:", file.size(OUT), "bytes\n")
