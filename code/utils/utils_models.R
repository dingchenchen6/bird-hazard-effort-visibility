# ============================================================
# Scientific question / 科学问题:
#   Consolidate the duplicated modelling helpers that were scattered
#   across v1 scripts (extract_coefs_manual, aic_table, dharma_full,
#   variance_decomp) into a single audited module.
#   把 v1 中重复定义的建模辅助函数收编到统一模块，避免不一致。
#
# Objective / 分析目标:
#   - Provide a single `extract_coefs_manual()` honouring conditional
#     vs zero-inflated components.
#   - Standardise AIC / BIC tables and Δ ranking.
#   - Wrap DHARMa global tests + a marginal/conditional R² helper.
#   - Implement the variance-decomposition formula used in Methods 2.3.
#
# Input data / 输入数据:
#   Fitted glmmTMB / lme4 / glm model objects.
#
# Main workflow / 主要流程:
#   1. Coefficient extraction with hazard-ratio + 95 % CI on cloglog.
#   2. AIC table with ΔAIC, Akaike weights, evidence ratios.
#   3. DHARMa simulation + global tests.
#   4. Variance decomposition (additive vs interaction, joint overlap).
#
# Expected output / 预期输出: NA — library functions.
# Key assumptions / 关键假设:
#   - glmmTMB 1.1.9, lme4 1.1, DHARMa 0.4.6, MuMIn 1.47.
# Main packages / 主要包: glmmTMB, lme4, DHARMa, MuMIn, broom.mixed,
#   performance, data.table, glue.
# Output directory / 输出路径: NA.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(glmmTMB)
  library(DHARMa)
  library(performance)
  library(MuMIn)
  library(broom.mixed)
  library(glue)
})

# ---- 1. Coefficient extraction ---------------------------------------------
# Returns a data.table with HR=exp(β) and 95 % Wald CI for cloglog links.
# 对 cloglog 链接函数，HR = exp(β) 自然对应离散时间风险率比。
extract_coefs_manual <- function(fit,
                                 component = c("cond", "zi"),
                                 fdr = TRUE) {
  component <- match.arg(component)
  is_glmmTMB <- inherits(fit, "glmmTMB")
  if (is_glmmTMB) {
    tab <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE,
                             component = component)
  } else {
    tab <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE)
  }
  data.table::setDT(tab)
  cols <- intersect(c("term", "estimate", "std.error", "statistic",
                      "p.value", "conf.low", "conf.high"), names(tab))
  tab <- tab[, ..cols]
  data.table::setnames(tab, "estimate", "beta")
  tab[, hr        := exp(beta)]
  tab[, hr.low    := exp(conf.low)]
  tab[, hr.high   := exp(conf.high)]
  if (fdr && "p.value" %in% names(tab)) {
    tab[, p.fdr := stats::p.adjust(p.value, method = "BH")]
  }
  tab[]
}

# ---- 2. AIC / BIC ranking table --------------------------------------------
aic_table <- function(model_list, names_ = NULL) {
  if (is.null(names_)) names_ <- names(model_list)
  if (is.null(names_)) names_ <- paste0("M", seq_along(model_list))
  res <- data.table::data.table(
    model = names_,
    df    = vapply(model_list, function(m) attr(stats::logLik(m), "df"),
                   numeric(1)),
    logLik = vapply(model_list, function(m) as.numeric(stats::logLik(m)),
                   numeric(1)),
    AIC = vapply(model_list, stats::AIC, numeric(1)),
    BIC = vapply(model_list, stats::BIC, numeric(1))
  )
  res[, dAIC := AIC - min(AIC)]
  res[, aic_w := exp(-0.5 * dAIC) / sum(exp(-0.5 * dAIC))]
  res[, evidence_ratio := max(aic_w) / aic_w]
  data.table::setorder(res, AIC)
  res[]
}

# ---- 3. DHARMa wrapper -----------------------------------------------------
dharma_full <- function(fit, n = 1000, seed = 42, plot = FALSE) {
  set.seed(seed)
  sim <- DHARMa::simulateResiduals(fit, n = n, plot = plot)
  list(
    sim       = sim,
    uniform   = DHARMa::testUniformity(sim, plot = FALSE),
    dispersion = DHARMa::testDispersion(sim, plot = FALSE),
    outlier   = DHARMa::testOutliers(sim, plot = FALSE),
    zi        = tryCatch(DHARMa::testZeroInflation(sim, plot = FALSE),
                         error = function(e) NULL)
  )
}

# ---- 4. VIF for glmmTMB ----------------------------------------------------
vif_glmmTMB <- function(fit) {
  X <- stats::model.matrix(fit)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(X) < 2) return(data.table::data.table(term = colnames(X), vif = NA_real_))
  R <- stats::cor(X)
  vifs <- diag(solve(R))
  data.table::data.table(term = names(vifs), vif = unname(vifs))
}

# ---- 5. Marginal / conditional R² ------------------------------------------
marginal_R2 <- function(fit) {
  out <- tryCatch(performance::r2(fit), error = function(e) NULL)
  if (is.null(out)) {
    return(list(R2_marginal = NA_real_, R2_conditional = NA_real_))
  }
  list(
    R2_marginal    = as.numeric(out$R2_marginal),
    R2_conditional = as.numeric(out$R2_conditional)
  )
}

# ---- 6. Variance decomposition (additive vs interaction) -------------------
# Splits R² into:
#   R²(climate)  — added when climate is dropped
#   R²(effort)   — added when effort is dropped
#   R²(interaction) — added when the interaction is dropped
#   joint overlap — R²(full) - sum of marginal components
# 这是 Methods 2.3 中报告 80.4 % 的口径。
variance_decomp <- function(fit_full, fit_no_climate,
                            fit_no_effort, fit_no_interaction) {
  r2 <- function(m) marginal_R2(m)$R2_marginal
  R_full   <- r2(fit_full)
  R_noC    <- r2(fit_no_climate)
  R_noE    <- r2(fit_no_effort)
  R_noI    <- r2(fit_no_interaction)
  contrib_climate     <- R_full - R_noC
  contrib_effort      <- R_full - R_noE
  contrib_interaction <- R_full - R_noI
  joint_overlap <- R_full - (contrib_climate + contrib_effort + contrib_interaction)
  data.table::data.table(
    component = c("climate", "effort", "interaction", "joint_overlap", "total_marginal_R2"),
    R2        = c(contrib_climate, contrib_effort, contrib_interaction,
                  joint_overlap, R_full),
    share     = c(contrib_climate, contrib_effort, contrib_interaction,
                  joint_overlap, R_full) / R_full
  )
}

# ---- 7. Pretty-print summary -----------------------------------------------
print_fit_summary <- function(fit, label = "") {
  cat(glue::glue("\n── {label} {Sys.time()} ──────────────────────────────\n\n"))
  print(summary(fit))
  cat("\nFixed effects (HR scale):\n")
  print(extract_coefs_manual(fit))
  invisible(NULL)
}
