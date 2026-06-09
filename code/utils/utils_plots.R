# ============================================================
# Scientific question / 科学问题:
#   Enforce a single, GEB-compliant aesthetic across all v2 figures
#   (theme, colour palettes, basemap, panel labels) so reviewers see a
#   coherent visual story rather than 25 mismatched scripts.
#   保证 v2 全套图件视觉风格统一、符合 GEB 制图规范。
#
# Objective / 分析目标:
#   - theme_geb() — base ggplot2 theme, 8 pt font, Helvetica fallback.
#   - China basemap via GS(2019)1822 shapefile family (Albers).
#   - Colour-blind palettes for diverging risk maps.
#   - Reusable helpers: forest_plot(), pdp_panel(), ridge_panel(),
#     bivariate_legend(), panel_letters().
#
# Input data / 输入数据:
#   data.frames / data.tables from results/ + sf basemap.
#
# Main workflow / 主要流程: NA — library functions.
# Expected output / 预期输出: NA.
# Key assumptions / 关键假设: ggplot2 >= 3.4, ggtext >= 0.1.2,
#   patchwork >= 1.2, ggridges >= 0.5.6, ggspatial >= 1.1, scales.
# Main packages / 主要包: ggplot2, ggtext, patchwork, ggridges,
#   ggspatial, scales, viridisLite, RColorBrewer.
# Output directory / 输出路径: NA.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtext)
  library(patchwork)
  library(scales)
  library(viridisLite)
})

# ---- 1. Master theme -------------------------------------------------------
# GEB column widths: 8.5 cm single, 17.5 cm double. 8 pt body font.
# 期刊单栏 8.5 / 双栏 17.5 cm；正文字体 8 pt。
theme_geb <- function(base_size = 8, base_family = "") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major   = ggplot2::element_line(linewidth = 0.2,
                                                 colour = "grey90"),
      panel.border       = ggplot2::element_rect(linewidth = 0.4,
                                                 colour = "grey20"),
      axis.text          = ggplot2::element_text(colour = "grey20"),
      axis.title         = ggplot2::element_text(colour = "grey10"),
      strip.background   = ggplot2::element_rect(fill = "grey95",
                                                 colour = "grey80"),
      strip.text         = ggplot2::element_text(face = "bold"),
      legend.position    = "right",
      legend.key.size    = ggplot2::unit(0.3, "cm"),
      legend.title       = ggplot2::element_text(face = "bold"),
      plot.title         = ggplot2::element_text(face = "bold",
                                                 size = base_size + 1),
      plot.tag           = ggplot2::element_text(face = "bold",
                                                 size = base_size + 1)
    )
}

# ---- 2. Colour-blind palettes ---------------------------------------------
# Diverging palette anchored at zero — Okabe-Ito + viridis hybrid.
pal_div <- function(n = 11) {
  scales::col_numeric(
    palette = c("#3B4CC0", "#7AA0E0", "#B8D8F2",
                "#F5F5F5",
                "#F2C8B8", "#E08C7A", "#B40426"),
    domain  = c(-1, 1)
  )(seq(-1, 1, length.out = n))
}

# Sequential viridis-derived. 顺序色谱。
pal_seq  <- viridisLite::viridis
pal_seq_inferno <- viridisLite::inferno

# Categorical (≤ 8 levels). 分类色板。
pal_cat <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
             "#9467BD", "#8C564B", "#E377C2", "#7F7F7F")

scale_fill_div_zero  <- function(...) ggplot2::scale_fill_gradient2(low = "#3B4CC0",
                                                                    mid = "#F5F5F5",
                                                                    high = "#B40426",
                                                                    midpoint = 0, ...)
scale_colour_div_zero <- function(...) ggplot2::scale_colour_gradient2(low = "#3B4CC0",
                                                                       mid = "#F5F5F5",
                                                                       high = "#B40426",
                                                                       midpoint = 0, ...)

# ---- 3. Basemap helper ----------------------------------------------------
china_basemap_gs2019_1822 <- function(crs = "EPSG:4524") {
  # Lazy-load province + ninedash + national. 仅在调用时读取，避免污染。
  src <- function(layer) {
    if (exists("read_gs2019_basemap", mode = "function")) {
      return(read_gs2019_basemap(layer))
    }
    stop("utils_spatial.R must be sourced before china_basemap_gs2019_1822().")
  }
  prov     <- sf::st_transform(src("province"),   crs)
  national <- tryCatch(sf::st_transform(src("national"), crs), error = function(e) NULL)
  ninedash <- tryCatch(sf::st_transform(src("ninedash"), crs), error = function(e) NULL)
  list(province = prov, national = national, ninedash = ninedash, crs = crs)
}

add_basemap_layers <- function(p, basemap,
                                fill_var = NULL, fill_data = NULL,
                                show_ninedash = TRUE) {
  if (!is.null(fill_data) && !is.null(fill_var)) {
    fill_sf <- merge(basemap$province, fill_data, by.x = "name", by.y = "province",
                     all.x = TRUE)
    p <- p + ggplot2::geom_sf(data = fill_sf,
                              ggplot2::aes(fill = .data[[fill_var]]),
                              colour = "grey50", linewidth = 0.1)
  } else {
    p <- p + ggplot2::geom_sf(data = basemap$province,
                              fill = "grey98", colour = "grey60",
                              linewidth = 0.1)
  }
  if (!is.null(basemap$national)) {
    p <- p + ggplot2::geom_sf(data = basemap$national, fill = NA,
                              colour = "grey20", linewidth = 0.4)
  }
  if (show_ninedash && !is.null(basemap$ninedash)) {
    p <- p + ggplot2::geom_sf(data = basemap$ninedash, fill = NA,
                              colour = "grey20", linewidth = 0.3)
  }
  p + ggplot2::coord_sf(crs = basemap$crs, datum = NA, expand = FALSE)
}

# ---- 4. Reusable plot helpers ---------------------------------------------
# Forest plot with HR + 95 % CI. 森林图（HR + 95% CI）。
forest_plot <- function(dt, x_col = "hr", lo_col = "hr.low",
                        hi_col = "hr.high", term_col = "term",
                        group_col = NULL, title = NULL, vline = 1) {
  dt <- data.table::as.data.table(dt)
  dt[[term_col]] <- factor(dt[[term_col]], levels = unique(dt[[term_col]]))
  aes_call <- if (is.null(group_col)) {
    ggplot2::aes(x = .data[[x_col]], y = .data[[term_col]])
  } else {
    ggplot2::aes(x = .data[[x_col]], y = .data[[term_col]],
                 colour = .data[[group_col]])
  }
  ggplot2::ggplot(dt, aes_call) +
    ggplot2::geom_vline(xintercept = vline, linetype = "dashed",
                         colour = "grey50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[lo_col]],
                                          xmax = .data[[hi_col]]),
                            height = 0.15, linewidth = 0.4) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_x_continuous(trans = "log",
                                 breaks = c(0.5, 1, 2, 4)) +
    ggplot2::labs(title = title, x = "Hazard ratio (log scale)",
                  y = NULL) +
    theme_geb()
}

# Partial-dependence panel from a feature x prediction grid. PDP 子图。
pdp_panel <- function(pdp_dt, x_col, y_col = "yhat", group_col = NULL,
                       facet_col = NULL, title = NULL) {
  aes_call <- if (is.null(group_col)) {
    ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])
  } else {
    ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]],
                  colour = .data[[group_col]])
  }
  p <- ggplot2::ggplot(pdp_dt, aes_call) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::labs(title = title) +
    theme_geb()
  if (!is.null(facet_col)) p <- p + ggplot2::facet_wrap(facet_col, scales = "free")
  p
}

# Panel letter helper (a, b, c, …). 子图字母。
panel_letters <- function(plots, prefix = "(", suffix = ")") {
  tags <- paste0(prefix, letters[seq_along(plots)], suffix)
  Map(function(p, t) p + ggplot2::labs(tag = t), plots, tags)
}

# Save publication-grade PDF + PNG pair at 600 dpi.
save_pub <- function(plot, name, width = 17.5, height = 12,
                     units = "cm", path = NULL) {
  if (is.null(path)) path <- v2_path("figures", "main")
  ensure_dir(path)
  pdf_path <- file.path(path, paste0(name, ".pdf"))
  png_path <- file.path(path, paste0(name, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height,
                   units = units, device = grDevices::cairo_pdf)
  ggplot2::ggsave(png_path, plot, width = width, height = height,
                   units = units, dpi = 600)
  invisible(c(pdf = pdf_path, png = png_path))
}
