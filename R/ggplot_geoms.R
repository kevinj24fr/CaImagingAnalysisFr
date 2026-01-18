#' Custom ggplot2 Geoms for Calcium Imaging
#'
#' Specialized geoms and scales for visualizing calcium imaging data.
#' These provide intuitive plotting interfaces for common visualization tasks.
#'
#' @name ggplot_geoms
#' @keywords internal
NULL

#' Geom for calcium trace heatmaps
#'
#' Creates a raster/heatmap visualization of calcium traces optimized
#' for fluorescence data.
#'
#' @param mapping Aesthetic mapping
#' @param data Data to display
#' @param stat Statistical transformation
#' @param position Position adjustment
#' @param ... Additional arguments
#' @param interpolate Interpolate raster
#' @param na.rm Remove NAs
#'
#' @return ggplot2 layer
#' @export
#'
#' @examples
#' \dontrun
#' library(ggplot2)
#' ggplot(trace_data, aes(x = time, y = cell_id, fill = dff)) +
#'   geom_calcium_heatmap() +
#'   scale_fill_calcium()
#' }
geom_calcium_heatmap <- function(mapping = NULL, data = NULL,
                                  stat = "identity", position = "identity",
                                  ..., interpolate = TRUE, na.rm = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = ggplot2::GeomRaster,
    position = position,
    params = list(
      interpolate = interpolate,
      na.rm = na.rm,
      ...
    )
  )
}

#' Geom for spike train raster plots
#'
#' Creates a raster plot of spike trains.
#'
#' @param mapping Aesthetic mapping (requires x for time, y for cell)
#' @param data Data to display
#' @param stat Statistical transformation
#' @param position Position adjustment
#' @param spike_height Height of spike marks (0-1)
#' @param ... Additional arguments
#'
#' @return ggplot2 layer
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(spike_data, aes(x = time, y = cell_id)) +
#'   geom_spike_train()
#' }
geom_spike_train <- function(mapping = NULL, data = NULL,
                              stat = "identity", position = "identity",
                              spike_height = 0.8, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Use geom_segment for spike representation
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSpikeTrain,
    position = position,
    params = list(
      spike_height = spike_height,
      na.rm = FALSE,
      ...
    )
  )
}

#' GeomSpikeTrain ggproto object
#'
#' @rdname geom_spike_train
#' @format NULL
#' @usage NULL
#' @keywords internal
#' @export
GeomSpikeTrain <- NULL

.init_geom_spike_train <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(NULL)
  }

  GeomSpikeTrain <<- ggplot2::ggproto(
    "GeomSpikeTrain",
    ggplot2::Geom,
    required_aes = c("x", "y"),
    default_aes = ggplot2::aes(
      colour = "black",
      linewidth = 0.5,
      alpha = NA
    ),

    draw_panel = function(data, panel_params, coord, spike_height = 0.8) {
      coords <- coord$transform(data, panel_params)

      # Create segment data for spikes
      half_height <- spike_height / 2

      grid::segmentsGrob(
        x0 = coords$x,
        y0 = coords$y - half_height * 0.1,
        x1 = coords$x,
        y1 = coords$y + half_height * 0.1,
        gp = grid::gpar(
          col = coords$colour,
          lwd = coords$linewidth * ggplot2::.pt,
          alpha = coords$alpha
        )
      )
    }
  )
}

#' Geom for event-aligned traces
#'
#' Plot traces aligned to events with confidence bands.
#'
#' @param mapping Aesthetic mapping
#' @param data Data to display
#' @param stat Statistical transformation
#' @param position Position adjustment
#' @param show_ci Show confidence interval
#' @param ci_alpha Alpha for CI band
#' @param ... Additional arguments
#'
#' @return ggplot2 layer
#' @export
geom_aligned_trace <- function(mapping = NULL, data = NULL,
                                stat = "summary", position = "identity",
                                show_ci = TRUE, ci_alpha = 0.3, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  layers <- list(
    # Mean line
    ggplot2::stat_summary(
      mapping = mapping,
      data = data,
      geom = "line",
      fun = mean,
      position = position,
      ...
    )
  )

  if (show_ci) {
    # Add confidence ribbon
    layers <- c(
      list(
        ggplot2::stat_summary(
          mapping = mapping,
          data = data,
          geom = "ribbon",
          fun.data = ggplot2::mean_se,
          alpha = ci_alpha,
          position = position
        )
      ),
      layers
    )
  }

  layers
}

#' Color scale for calcium imaging
#'
#' Optimized color scale for fluorescence data (dF/F).
#'
#' @param ... Additional arguments passed to scale_fill_gradientn
#' @param palette Color palette: "viridis", "magma", "inferno", "plasma", "calcium"
#' @param midpoint Midpoint for diverging scale (default: 0)
#' @param limits Scale limits
#'
#' @return ggplot2 scale
#' @export
#'
#' @examples
#' \dontrun{
#' ggplot(data, aes(x, y, fill = dff)) +
#'   geom_raster() +
#'   scale_fill_calcium()
#' }
scale_fill_calcium <- function(..., palette = "viridis", midpoint = 0, limits = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  if (palette == "calcium") {
    # Custom calcium imaging palette (blue-black-yellow-red)
    colors <- c("#0000AA", "#000044", "#000000", "#444400", "#AAAA00", "#FF4400", "#FF0000")
    return(ggplot2::scale_fill_gradientn(colors = colors, limits = limits, ...))
  }

  # Use viridis palettes
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    colors <- switch(palette,
      viridis = viridisLite::viridis(256),
      magma = viridisLite::magma(256),
      inferno = viridisLite::inferno(256),
      plasma = viridisLite::plasma(256),
      viridisLite::viridis(256)
    )
    ggplot2::scale_fill_gradientn(colors = colors, limits = limits, ...)
  } else {
    ggplot2::scale_fill_viridis_c(option = substr(palette, 1, 1), limits = limits, ...)
  }
}

#' Color scale for correlation matrices
#'
#' Diverging scale optimized for correlation values.
#'
#' @param ... Additional arguments
#' @param limits Scale limits (default: -1 to 1)
#' @param midpoint Center of scale
#'
#' @return ggplot2 scale
#' @export
scale_fill_correlation <- function(..., limits = c(-1, 1), midpoint = 0) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = midpoint,
    limits = limits,
    ...
  )
}

#' Facet by cell
#'
#' Convenience function to facet plots by cell ID.
#'
#' @param ncol Number of columns
#' @param nrow Number of rows
#' @param scales Scale parameter for facet_wrap
#' @param ... Additional arguments to facet_wrap
#'
#' @return ggplot2 facet specification
#' @export
#'
#' @examples
#' \dontrun{
#' ggplot(trace_data, aes(x = time, y = dff)) +
#'   geom_line() +
#'   facet_cell(ncol = 5)
#' }
facet_cell <- function(ncol = NULL, nrow = NULL, scales = "free_y", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  ggplot2::facet_wrap(
    ~ cell_id,
    ncol = ncol,
    nrow = nrow,
    scales = scales,
    ...
  )
}

#' Plot calcium traces
#'
#' Convenience function for plotting multiple calcium traces.
#'
#' @param data data.frame with columns: time, value/dff, cell_id
#' @param y_var Variable to plot on y-axis
#' @param color_var Variable for coloring
#' @param facet Facet by cell_id
#' @param title Plot title
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#' \dontrun{
#' dt <- traces_to_dt(traces)
#' plot_traces(dt)
#' }
plot_traces <- function(data, y_var = "value", color_var = NULL,
                        facet = TRUE, title = "Calcium Traces") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$time, y = .data[[y_var]]))

  if (!is.null(color_var)) {
    p <- p + ggplot2::aes(color = .data[[color_var]])
  }

  p <- p +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::labs(
      title = title,
      x = "Time",
      y = y_var
    ) +
    ggplot2::theme_minimal()

  if (facet && "cell_id" %in% names(data)) {
    p <- p + facet_cell(scales = "free_y")
  }

  p
}

#' Plot correlation matrix
#'
#' @param cor_matrix Correlation matrix
#' @param title Plot title
#' @param cluster Hierarchically cluster cells
#'
#' @return ggplot2 object
#' @export
plot_correlation_matrix <- function(cor_matrix, title = "Pairwise Correlations",
                                    cluster = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  if (cluster) {
    # Hierarchical clustering for ordering
    hc <- hclust(as.dist(1 - cor_matrix))
    ord <- hc$order
    cor_matrix <- cor_matrix[ord, ord]
  }

  # Convert to long format
  n <- nrow(cor_matrix)
  cell_ids <- rownames(cor_matrix) %||% paste0("cell_", seq_len(n))

  df <- expand.grid(
    cell_1 = factor(cell_ids, levels = cell_ids),
    cell_2 = factor(cell_ids, levels = cell_ids)
  )
  df$correlation <- as.vector(cor_matrix)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$cell_1, y = .data$cell_2, fill = .data$correlation)) +
    ggplot2::geom_tile() +
    scale_fill_correlation() +
    ggplot2::labs(title = title, x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed()
}

#' Plot spike raster
#'
#' @param spikes Matrix of spikes (cells x time) or data.frame
#' @param time_vector Time vector
#' @param title Plot title
#'
#' @return ggplot2 object
#' @export
plot_spike_raster <- function(spikes, time_vector = NULL, title = "Spike Raster") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  if (is.matrix(spikes)) {
    n_cells <- nrow(spikes)
    n_frames <- ncol(spikes)

    if (is.null(time_vector)) {
      time_vector <- seq_len(n_frames)
    }

    # Convert to long format with only spike times
    spike_data <- data.frame(
      cell_id = integer(),
      time = numeric()
    )

    for (i in seq_len(n_cells)) {
      spike_times <- which(spikes[i, ] > 0)
      if (length(spike_times) > 0) {
        spike_data <- rbind(spike_data, data.frame(
          cell_id = i,
          time = time_vector[spike_times]
        ))
      }
    }
  } else {
    spike_data <- spikes
  }

  ggplot2::ggplot(spike_data, ggplot2::aes(x = .data$time, y = .data$cell_id)) +
    ggplot2::geom_point(shape = "|", size = 2) +
    ggplot2::labs(
      title = title,
      x = "Time",
      y = "Cell"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Plot aligned responses (PETH-style)
#'
#' @param aligned_data List from align_to_events or data.frame
#' @param show_trials Show individual trial traces
#' @param show_mean Show mean trace
#' @param show_sem Show standard error
#' @param title Plot title
#'
#' @return ggplot2 object
#' @export
plot_aligned_responses <- function(aligned_data, show_trials = FALSE,
                                   show_mean = TRUE, show_sem = TRUE,
                                   title = "Event-Aligned Response") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Convert aligned array to long format
  if (is.list(aligned_data) && "aligned" %in% names(aligned_data)) {
    arr <- aligned_data$aligned
    time_vec <- aligned_data$time_vector

    n_cells <- dim(arr)[1]
    n_time <- dim(arr)[2]
    n_trials <- dim(arr)[3]

    df <- data.frame(
      cell_id = rep(rep(seq_len(n_cells), each = n_time), times = n_trials),
      time = rep(time_vec, times = n_cells * n_trials),
      trial = rep(seq_len(n_trials), each = n_cells * n_time),
      value = as.vector(aperm(arr, c(2, 1, 3)))
    )
  } else {
    df <- aligned_data
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$value))

  if (show_trials) {
    p <- p + ggplot2::geom_line(
      ggplot2::aes(group = interaction(.data$cell_id, .data$trial)),
      alpha = 0.2,
      linewidth = 0.3
    )
  }

  if (show_mean || show_sem) {
    if (show_sem) {
      p <- p + ggplot2::stat_summary(
        fun.data = ggplot2::mean_se,
        geom = "ribbon",
        alpha = 0.3
      )
    }
    if (show_mean) {
      p <- p + ggplot2::stat_summary(
        fun = mean,
        geom = "line",
        linewidth = 1
      )
    }
  }

  p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    ggplot2::labs(
      title = title,
      x = "Time from event",
      y = "Response"
    ) +
    ggplot2::theme_minimal()

  if ("cell_id" %in% names(df) && length(unique(df$cell_id)) > 1) {
    p <- p + ggplot2::facet_wrap(~ cell_id, scales = "free_y")
  }

  p
}

# Note: %||% operator defined in aaa_utils.R
