#' Publication-Ready Visualization Theme and Palettes
#'
#' A clean, publication-ready style for scientific figures.
#' Designed for Cell, Nature, and similar high-impact journals.
#'
#' @name publication_theme
NULL

#' Publication Color Palettes
#'
#' Color palettes following scientific publication standards.
#'
#' @param n Number of colors needed
#' @param palette Palette type: "categorical", "sequential", "diverging", "utility"
#' @param reverse Reverse the palette order
#'
#' @return Vector of hex color codes
#' @export
#'
#' @examples
#' pub_colors(4)
#' pub_colors(6, "categorical")
#' pub_colors(9, "sequential")
pub_colors <- function(n = 4, palette = c("categorical", "sequential", "diverging", "utility"),
                       reverse = FALSE) {
  palette <- match.arg(palette)

  colors <- switch(palette,
    categorical = c(
      "#E41A1C",  # red

"#377EB8",  # blue
      "#4DAF4A",  # green
      "#984EA3",  # purple
      "#FF7F00",  # orange
      "#FFFF33",  # yellow
      "#A65628",  # brown
      "#F781BF"   # pink
    ),
    sequential = c(
      "#DEEBF7",  # lightest
      "#C6DBEF",
      "#9ECAE1",
      "#6BAED6",
      "#4292C6",
      "#2171B5",
      "#08519C",
      "#08306B"   # darkest
    ),
    diverging = c(
      "#D73027",  # negative (red)
      "#F46D43",
      "#FDAE61",
      "#FEE090",
      "#FFFFBF",  # neutral
      "#E0F3F8",
      "#ABD9E9",
      "#74ADD1",
      "#4575B4",  # positive (blue)
      "#1A9850"   # alternative positive (green)
    ),
    utility = c(
      "#4393C3",  # significant (steel blue)
      "#BDBDBD",  # non-significant (grey)
      "#757575",  # reference line (dark grey)
      "#FF7F00",  # highlight (orange)
      "#2166AC",  # cold
      "#B2182B"   # hot
    )
  )

  if (n > length(colors)) {
    colors <- colorRampPalette(colors)(n)
  } else {
    colors <- colors[1:n]
  }

  if (reverse) colors <- rev(colors)
  colors
}

#' Publication Theme for ggplot2
#'
#' A clean, minimal theme suitable for high-impact scientific publications.
#'
#' @param base_size Base font size (default 11pt for standalone figures)
#' @param base_family Font family
#' @param grid Show major grid lines
#'
#' @return A ggplot2 theme object
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   theme_publication()
#' }
theme_publication <- function(base_size = 11, base_family = "", grid = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  grid_color <- if (grid) "grey92" else "transparent"

  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Panel
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = ggplot2::element_line(color = grid_color, linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),

      # Axes
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.4),
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1),
      axis.title = ggplot2::element_text(color = "black", size = base_size),
      axis.line = ggplot2::element_blank(),

      # Text
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "grey30", hjust = 0),
      plot.caption = ggplot2::element_text(size = base_size - 2, color = "grey50", hjust = 1),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 2),

      # Legend
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      legend.title = ggplot2::element_text(face = "bold", size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 1),

      # Margins
      plot.margin = ggplot2::margin(10, 12, 10, 10),

      # Strip (for facets)
      strip.background = ggplot2::element_rect(fill = "grey95", color = "black", linewidth = 0.4),
      strip.text = ggplot2::element_text(face = "bold", size = base_size)
    )
}

#' Publication Scale for Continuous Fill
#'
#' Continuous color scale for heatmaps and density plots.
#'
#' @param palette "sequential", "diverging", or "viridis"
#' @param midpoint For diverging scales, the midpoint value
#' @param limits Scale limits
#' @param ... Additional arguments to scale function
#'
#' @return A ggplot2 scale object
#' @export
scale_fill_publication <- function(palette = c("sequential", "diverging", "viridis"),
                                   midpoint = 0, limits = NULL, ...) {
  palette <- match.arg(palette)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  if (palette == "sequential") {
    ggplot2::scale_fill_gradientn(
      colors = c("#DEEBF7", "#9ECAE1", "#3182BD", "#08306B"),
      limits = limits,
      ...
    )
  } else if (palette == "diverging") {
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = midpoint,
      limits = limits,
      ...
    )
  } else {
    ggplot2::scale_fill_viridis_c(option = "viridis", limits = limits, ...)
  }
}

#' Publication Scale for Discrete Colors
#'
#' Discrete color scale for categorical data.
#'
#' @param n Number of categories (auto-detected if NULL)
#' @param ... Additional arguments to scale function
#'
#' @return A ggplot2 scale object
#' @export
scale_color_publication <- function(n = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  if (is.null(n)) n <- 8
  ggplot2::scale_color_manual(values = pub_colors(n, "categorical"), ...)
}

#' @rdname scale_color_publication
#' @export
scale_colour_publication <- scale_color_publication

#' @rdname scale_color_publication
#' @export
scale_fill_publication_d <- function(n = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  if (is.null(n)) n <- 8
  ggplot2::scale_fill_manual(values = pub_colors(n, "categorical"), ...)
}

#' Add Panel Label to Plot
#'
#' Add a panel label (A, B, C, etc.) to a ggplot.
#'
#' @param label The label text (e.g., "A", "B")
#' @param x X position (0-1, default top-left)
#' @param y Y position (0-1, default top-left)
#' @param size Font size
#'
#' @return A ggplot2 annotation layer
#' @export
add_panel_label <- function(label, x = 0.02, y = 0.98, size = 14) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  ggplot2::annotation_custom(
    grob = grid::textGrob(
      label,
      x = grid::unit(x, "npc"),
      y = grid::unit(y, "npc"),
      hjust = 0,
      vjust = 1,
      gp = grid::gpar(fontface = "bold", fontsize = size)
    )
  )
}

#' Significance Formatting
#'
#' Convert p-values to publication-standard significance notation.
#'
#' @param p P-value or vector of p-values
#' @param style "stars" for asterisks, "text" for ns/*, "exact" for p = 0.023
#'
#' @return Character string or vector of significance labels
#' @export
format_pvalue <- function(p, style = c("stars", "text", "exact")) {
  style <- match.arg(style)

  format_single <- function(pval) {
    if (is.na(pval)) return("NA")

    if (style == "exact") {
      if (pval < 0.001) return("p < 0.001")
      return(sprintf("p = %.3f", pval))
    }

    if (pval >= 0.05) {
      if (style == "stars") return("") else return("ns")
    } else if (pval < 0.001) {
      return("***")
    } else if (pval < 0.01) {
      return("**")
    } else {
      return("*")
    }
  }

  sapply(p, format_single)
}

#' Save Publication Figure
#'
#' Save a figure with publication-ready settings.
#'
#' @param plot A ggplot object
#' @param filename Output filename
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution (default 300)
#' @param format "png", "pdf", or "tiff"
#'
#' @export
save_pub_figure <- function(plot, filename, width = 6, height = 5,
                            dpi = 300, format = c("png", "pdf", "tiff")) {
  format <- match.arg(format)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  ext <- paste0(".", format)
  if (!grepl(paste0(ext, "$"), filename)) {
    filename <- paste0(filename, ext)
  }

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )

  message(sprintf("Saved: %s (%dx%d inches, %d DPI)", filename, width, height, dpi))
}

#' Create Multi-Panel Figure
#'
#' Combine multiple plots into a publication-ready multi-panel figure.
#'
#' @param ... ggplot objects to combine
#' @param ncol Number of columns
#' @param nrow Number of rows
#' @param labels Panel labels ("AUTO" for A, B, C..., or vector)
#' @param label_size Size of panel labels
#' @param guides How to handle legends: "collect", "keep", "auto"
#'
#' @return A patchwork object
#' @export
combine_panels <- function(..., ncol = NULL, nrow = NULL, labels = "AUTO",
                           label_size = 14, guides = "collect") {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package required. Install with: install.packages('patchwork')")
  }

  plots <- list(...)

  # Auto-generate labels
  if (identical(labels, "AUTO")) {
    labels <- LETTERS[seq_along(plots)]
  }

  # Add labels to plots
  if (!is.null(labels) && length(labels) == length(plots)) {
    plots <- lapply(seq_along(plots), function(i) {
      plots[[i]] + ggplot2::labs(tag = labels[i])
    })
  }

  # Combine with patchwork
  combined <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow)

  if (guides == "collect") {
    combined <- combined + patchwork::plot_layout(guides = "collect")
  }

  combined & ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = label_size))
}
