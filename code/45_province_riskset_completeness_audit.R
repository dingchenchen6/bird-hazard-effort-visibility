# ============================================================
# Scientific question / 科学问题:
#   The v1 province risk set has 12,813 rows / 333 species / 512
#   events, but raw events 2002-2024 contain 519 species / 930
#   records. This script audits the attrition chain in full,
#   identifies the precise mechanism of loss (SDM modelling vs SDM
#   threshold vs complete-case merge), and produces a reconciliation
#   table for transparency.
#   省级风险集完整性审计：从原始 1026 records / 565 species 一路追到
#   12,813 rows / 333 species / 512 events，定位每一步损失原因。
#
# Outputs:
#   results/diagnostics/table_province_riskset_completeness.csv
#   results/diagnostics/table_missing_event_species_detail.csv
#   results/diagnostics/figure_riskset_attrition_funnel.{pdf,png}
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
options(warn = 1)

V2 <- normalizePath(".", mustWork = TRUE)
V1 <- normalizePath(file.path(V2, "..", "bird_hazard_model_effort_upgrade"),
                     mustWork = FALSE)
SDM_PROV <- normalizePath(file.path(V2, "..", "bird_new_record_hazard_model",
                                     "results", "combined_threshold_100_test",
                                     "derived_inputs", "sdm_province.csv"),
                           mustWork = FALSE)
SDM_BW <- normalizePath(file.path(V2, "..",
                                   "bird_sdm_distribution_modeling_birdwatch_2002_2025"),
                         mustWork = FALSE)
SDM_RS <- normalizePath(file.path(V2, "..",
                                   "bird_sdm_distribution_modeling_rescue_1980_2025_gbif"),
                         mustWork = FALSE)

ens <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE,
                                                    showWarnings = FALSE)
ens(file.path(V2, "results", "diagnostics"))

log <- function(...) cat(sprintf("[45 %s] ", format(Sys.time(), "%H:%M:%S")),
                          ..., "\n", sep = "")

# ---- 1. Stage counts ----------------------------------------------------
ev <- fread(file.path(V1, "data", "events_100km_grid_assigned.csv"),
            encoding = "UTF-8")
setnames(ev, tolower(names(ev)))
if (!"year" %in% names(ev) && "pub_year" %in% names(ev))
  setnames(ev, "pub_year", "year")

ev_raw    <- copy(ev)
ev_window <- ev[year >= 2002 & year <= 2024]

m_bw <- fread(file.path(SDM_BW, "data", "tables",
                          "table_model_occurrence_points_used_all_species.csv"),
                encoding = "UTF-8")
m_rs <- fread(file.path(SDM_RS, "data", "tables",
                          "table_model_occurrence_points_used_all_species.csv"),
                encoding = "UTF-8")
sdm_modelled <- union(unique(m_bw$species), unique(m_rs$species))

sdm  <- fread(SDM_PROV, encoding = "UTF-8")
sdm_threshold_pass <- unique(sdm[potential == 1L & historical_presence == 0L,
                                  species])

risk <- fread(file.path(V1, "data",
                          "hazard_risk_upgraded_complete_case.csv"),
                encoding = "UTF-8")
risk_sp <- unique(risk$species)

# ---- 2. Build attrition table -------------------------------------------
stages <- data.table(
  Stage = c(
    "1. Raw events (xlsx → events_100km_grid_assigned.csv, 2000-2025)",
    "2. Filter to study window 2002-2024",
    "3. SDM modelling ran (birdwatch + rescue projects, union)",
    "4. Pass SDM threshold (potential=1 & historical_presence=0)",
    "5. Complete-case merge with effort + climate (final risk set)"),
  Records = c(
    nrow(ev_raw),
    nrow(ev_window),
    NA,  # SDM modelling is per-species, not per-event
    NA,
    nrow(risk)),
  `Unique species` = c(
    uniqueN(ev_raw$species),
    uniqueN(ev_window$species),
    length(sdm_modelled),
    length(sdm_threshold_pass),
    uniqueN(risk$species)),
  `Events kept` = c(
    nrow(ev_raw),
    nrow(ev_window),
    sum(ev_window$species %in% sdm_modelled),
    sum(ev_window$species %in% sdm_threshold_pass),
    sum(risk$event)))
stages[, `Event loss vs window` := nrow(ev_window) -
        as.numeric(ifelse(is.na(`Events kept`), nrow(ev_window),
                            `Events kept`))]
stages[Stage == "1. Raw events (xlsx → events_100km_grid_assigned.csv, 2000-2025)",
       `Event loss vs window` := NA]
fwrite(stages, file.path(V2, "results", "diagnostics",
                          "table_province_riskset_completeness.csv"))
log("wrote table_province_riskset_completeness.csv")
print(stages)

# ---- 3. Missing event-species detail ------------------------------------
not_in_sdm <- setdiff(unique(ev_window$species), unique(sdm$species))
missing_detail <- ev_window[species %in% not_in_sdm,
                              .(events_lost = .N,
                                provinces = paste(unique(province), collapse=";"),
                                first_year = min(year, na.rm=TRUE),
                                last_year  = max(year, na.rm=TRUE)),
                              by = species]
missing_detail[, in_birdwatch_modelled := species %in% unique(m_bw$species)]
missing_detail[, in_rescue_modelled    := species %in% unique(m_rs$species)]
missing_detail[, modelling_status := fcase(
  in_birdwatch_modelled, "modelled_in_birdwatch_but_dropped_by_threshold",
  in_rescue_modelled,    "modelled_in_rescue_but_not_in_sdm_province",
  default = "not_modelled_at_all")]
setorder(missing_detail, -events_lost)
fwrite(missing_detail,
       file.path(V2, "results", "diagnostics",
                  "table_missing_event_species_detail.csv"))
log("wrote table_missing_event_species_detail.csv (",
    nrow(missing_detail), " species)")

cat("\n=== Missing event-species by modelling status ===\n")
print(missing_detail[, .(n_species = .N,
                          total_events_lost = sum(events_lost)),
                       by = modelling_status])

# ---- 4. Attrition funnel figure -----------------------------------------
fig_dt <- data.table(
  Stage = factor(c("1. Raw\nevents\n2000-25",
                    "2. 2002-2024\nwindow",
                    "3. SDM\nmodelling\nran",
                    "4. SDM\nthreshold\npass",
                    "5. Final\nrisk set\n(complete case)"),
                   levels = c("1. Raw\nevents\n2000-25",
                               "2. 2002-2024\nwindow",
                               "3. SDM\nmodelling\nran",
                               "4. SDM\nthreshold\npass",
                               "5. Final\nrisk set\n(complete case)")),
  Species = stages$`Unique species`,
  Events  = c(stages$Records[1:2], stages$`Events kept`[3:5]))

p_sp <- ggplot(fig_dt, aes(x = Stage, y = Species, fill = Stage)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::comma(Species)),
             vjust = -0.4, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("#1F77B4","#3B6FB4","#2CA02C","#D17F0E","#B40426"),
                     guide = "none") +
  labs(title = "(a) Unique species at each attrition stage",
        y = "Unique species", x = NULL) +
  expand_limits(y = max(fig_dt$Species) * 1.18) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", size=10))

p_ev <- ggplot(fig_dt, aes(x = Stage, y = Events, fill = Stage)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::comma(Events)),
             vjust = -0.4, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("#1F77B4","#3B6FB4","#2CA02C","#D17F0E","#B40426"),
                     guide = "none") +
  scale_y_continuous(trans = "log10",
                      labels = scales::comma_format()) +
  labs(title = "(b) Event-records retained (log scale)",
        y = "Events / records", x = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", size=10))

if (!requireNamespace("patchwork", quietly = TRUE))
  install.packages("patchwork", repos="https://cloud.r-project.org")
suppressPackageStartupMessages(library(patchwork))
fig <- p_sp / p_ev + plot_annotation(
  title = "Province risk-set attrition funnel",
  subtitle = "Raw 1,026 records / 565 species → final 12,813 rows / 333 species / 512 events",
  theme = theme(plot.title=element_text(face="bold", size=10),
                 plot.subtitle=element_text(size=8.5, colour="grey30")))

ens(file.path(V2, "figures", "diagnostics"))
ggsave(file.path(V2, "figures", "diagnostics",
                  "figure_riskset_attrition_funnel.pdf"),
       fig, width = 17, height = 14, units = "cm",
       device = grDevices::cairo_pdf)
ggsave(file.path(V2, "figures", "diagnostics",
                  "figure_riskset_attrition_funnel.png"),
       fig, width = 17, height = 14, units = "cm", dpi = 600)
log("wrote figure_riskset_attrition_funnel.{pdf,png}")

log("=== DONE ===")
