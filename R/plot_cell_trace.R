#' Plot Corrected Trace with Detected Spikes
#' 
#' Creates a comprehensive visualization of calcium traces with spike detection
#' results, including the original signal, deconvolved activity, and detected spikes.
#' 
#' @param corrected_df Output of [calcium_correction()]
#' @param cell Name of a cell column
#' @param method Spike inference method (default: "oasis")
#' @param show_spikes Whether to show detected spikes (default: TRUE)
#' @param show_deconvolved Whether to show deconvolved signal (default: TRUE)
#' @param colors Color scheme for the plot
#' @param verbose Whether to show progress messages (default: FALSE)
#' 
#' @return ggplot object
#' 
#' @examples
#' # Basic usage
#' raw <- generate_synthetic_data(3, 500)
#' corrected <- calcium_correction(raw)
#' p <- plot_cell_trace(corrected, "Cell_1")
#' print(p)
#' 
#' # Custom options
#' p <- plot_cell_trace(corrected, "Cell_1", method = "caiman", show_spikes = FALSE)
#' 
#' @export
plot_cell_trace <- function(corrected_df, 
                           cell,
                           method = "oasis",
                           show_spikes = TRUE,
                           show_deconvolved = TRUE,
                           colors = NULL,
                           verbose = FALSE) {
  
  # Validate inputs
  if (!is.data.frame(corrected_df)) {
    stop("corrected_df must be a data frame")
  }
  
  if (!cell %in% names(corrected_df)) {
    stop("Cell '", cell, "' not found in corrected_df")
  }
  
  if (!"Time" %in% names(corrected_df)) {
    stop("corrected_df must contain a 'Time' column")
  }
  
  # Set default colors if not provided
  if (is.null(colors)) {
    colors <- list(
      signal = "steelblue",
      deconvolved = "orange",
      spikes = "red",
      background = "lightgray"
    )
  }
  
  if (verbose) {
    message("Creating plot for cell: ", cell)
  }
  
  # Prepare data
  dat <- data.frame(
    Time = corrected_df$Time,
    Signal = corrected_df[[cell]]
  )
  
  # Perform spike inference if requested
  if (show_spikes || show_deconvolved) {
    if (verbose) message("  Performing spike inference...")
    
    spikes_result <- tryCatch({
      infer_spikes(dat$Signal, method = method, verbose = verbose)
    }, error = function(e) {
      warning("Spike inference failed: ", e$message)
      data.frame(
        fit = rep(0, length(dat$Signal)),
        spike = rep(0, length(dat$Signal))
      )
    })
    
    dat$Deconvolved <- spikes_result$fit
    dat$Spike <- spikes_result$spike > 0
  }
  
  # Create base plot
  p <- ggplot2::ggplot(data = dat, ggplot2::aes(x = Time)) +
    ggplot2::geom_line(ggplot2::aes(y = Signal), color = colors$signal, linewidth = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Calcium Trace:", cell),
      x = "Time",
      y = "Fluorescence"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )
  
  # Add deconvolved signal
  if (show_deconvolved && "Deconvolved" %in% names(dat)) {
    p <- p + ggplot2::geom_line(
      ggplot2::aes(y = Deconvolved), 
      color = colors$deconvolved, 
      alpha = 0.7,
      linewidth = 0.6
    )
  }
  
  # Add detected spikes
  if (show_spikes && "Spike" %in% names(dat)) {
    spike_data <- dat[dat$Spike, ]
    if (nrow(spike_data) > 0) {
      p <- p + ggplot2::geom_point(
        data = spike_data,
        ggplot2::aes(y = Signal),
        color = colors$spikes,
        size = 1.5,
        alpha = 0.8
      )
    }
  }
  
  # Add legend
  if (show_deconvolved || show_spikes) {
    p <- p + ggplot2::theme(legend.position = "bottom")
  }
  
  return(p)
}

#' Plot Multiple Cell Traces
#' 
#' Creates a faceted plot showing multiple cell traces simultaneously.
#' 
#' @param corrected_df Output of [calcium_correction()]
#' @param cells Vector of cell names to plot (default: all cells)
#' @param ncol Number of columns in the faceted plot (default: 3)
#' @param method Spike inference method
#' @param show_spikes Whether to show detected spikes
#' @param show_deconvolved Whether to show deconvolved signal
#' @param verbose Whether to show progress messages
#' 
#' @return ggplot object
#' 
#' @examples
#' raw <- generate_synthetic_data(6, 500)
#' corrected <- calcium_correction(raw)
#' p <- plot_multiple_cells(corrected, ncol = 2)
#' 
#' @export
plot_multiple_cells <- function(corrected_df,
                               cells = NULL,
                               ncol = 3,
                               method = "oasis",
                               show_spikes = TRUE,
                               show_deconvolved = TRUE,
                               verbose = FALSE) {
  
  # Validate ncol parameter
  if (!is.numeric(ncol) || length(ncol) != 1 || ncol <= 0) {
    stop("ncol must be a positive numeric value")
  }
  
  # Get cell columns if not specified
  if (is.null(cells)) {
    config <- get_config()
    cells <- names(corrected_df)[grepl(config$cell_pattern, names(corrected_df))]
  }
  
  # Validate cells
  missing_cells <- setdiff(cells, names(corrected_df))
  if (length(missing_cells) > 0) {
    stop("Cells not found in data: ", paste(missing_cells, collapse = ", "))
  }
  
  if (verbose) {
    message("Creating multi-cell plot for ", length(cells), " cells")
  }
  
  # Create plots for each cell
  plot_list <- lapply(cells, function(cell) {
    plot_cell_trace(
      corrected_df, 
      cell, 
      method = method,
      show_spikes = show_spikes,
      show_deconvolved = show_deconvolved,
      verbose = verbose
    ) + ggplot2::labs(title = cell)
  })
  
  # Combine plots using patchwork if available
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- Reduce(`+`, plot_list) + 
      patchwork::plot_layout(ncol = ncol)
    return(combined_plot)
  } else {
    # Fallback to grid.arrange
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      return(gridExtra::grid.arrange(grobs = plot_list, ncol = ncol))
    } else {
      warning("Neither patchwork nor gridExtra available. Returning list of plots.")
      return(plot_list)
    }
  }
}
