#' Calcium Imaging Analysis Pipeline using targets
#' 
#' Defines a reproducible workflow for calcium imaging analysis using
#' the `targets` package. This pipeline includes data generation,
#' correction, spike inference, and visualization.
#' 
#' @param n_cells Number of cells for synthetic data (default: 5)
#' @param n_time Number of time points (default: 1000)
#' @param spike_prob Spike probability (default: 0.02)
#' @param correction_method Correction method (default: "modern")
#' @param spike_method Spike inference method (default: "oasis")
#' @param span Loess span parameter (default: 0.45)
#' @param normalize Whether to normalize traces (default: TRUE)
#' 
#' @return List of targets for the pipeline
#' 
#' @examples
#' # Basic pipeline
#' pipeline <- calcium_pipeline()
#' 
#' # Custom parameters
#' pipeline <- calcium_pipeline(
#'   n_cells = 10,
#'   n_time = 2000,
#'   correction_method = "legacy",
#'   spike_method = "caiman"
#' )
#' 
#' @export
calcium_pipeline <- function(n_cells = NULL,
                            n_time = NULL,
                            spike_prob = NULL,
                            correction_method = "modern",
                            spike_method = "oasis",
                            span = 0.45,
                            normalize = TRUE) {
  
  # Check if targets is available
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required for the pipeline. Install with: install.packages('targets')")
  }
  
  # Get defaults from config
  config <- get_config()
  if (is.null(n_cells)) n_cells <- config$default_n_cells
  if (is.null(n_time)) n_time <- config$default_n_time
  if (is.null(spike_prob)) spike_prob <- config$default_spike_prob
  
  # Validate parameters
  validate_numeric_param(n_cells, 1, 1000, "n_cells")
  validate_numeric_param(n_time, 10, 100000, "n_time")
  validate_numeric_param(spike_prob, config$min_spike_prob, config$max_spike_prob, "spike_prob")
  validate_numeric_param(span, config$min_span, config$max_span, "span")
  
  if (!correction_method %in% c("modern", "legacy")) {
    stop("correction_method must be 'modern' or 'legacy'")
  }
  
  if (!spike_method %in% config$supported_methods) {
    stop("spike_method must be one of: ", paste(config$supported_methods, collapse = ", "))
  }
  
  # Define pipeline targets
  list(
    # Data generation
    targets::tar_target(
      raw_data,
      generate_synthetic_data(
        n_cells = n_cells,
        n_time = n_time,
        spike_prob = spike_prob,
        verbose = FALSE
      ),
      description = "Generate synthetic calcium imaging data"
    ),
    
    # Data correction
    targets::tar_target(
      corrected_data,
      calcium_correction(
        raw_data,
        method = correction_method,
        span = span,
        normalize = normalize,
        verbose = FALSE
      ),
      description = "Apply calcium correction to raw data"
    ),
    
    # Spike inference for all cells
    targets::tar_target(
      spike_results,
      {
        config <- get_config()
        cell_cols <- names(corrected_data)[grepl(config$cell_pattern, names(corrected_data))]
        
        results <- lapply(cell_cols, function(cell) {
          tryCatch({
            trace <- corrected_data[[cell]]
            spikes <- infer_spikes(trace, method = spike_method, verbose = FALSE)
            list(
              cell = cell,
              spikes = spikes,
              n_spikes = sum(spikes$spike > 0, na.rm = TRUE),
              spike_rate = sum(spikes$spike > 0, na.rm = TRUE) / length(trace) * 100
            )
          }, error = function(e) {
            list(
              cell = cell,
              spikes = data.frame(fit = rep(0, length(corrected_data[[cell]])), 
                                spike = rep(0, length(corrected_data[[cell]]))),
              n_spikes = 0,
              spike_rate = 0,
              error = e$message
            )
          })
        })
        
        names(results) <- cell_cols
        results
      },
      description = "Perform spike inference on all cells"
    ),
    
    # Summary statistics
    targets::tar_target(
      summary_stats,
      {
        # Extract statistics from spike results
        stats <- data.frame(
          Cell = names(spike_results),
          N_Spikes = sapply(spike_results, function(x) x$n_spikes),
          Spike_Rate = sapply(spike_results, function(x) x$spike_rate),
          stringsAsFactors = FALSE
        )
        
        # Add trace statistics
        config <- get_config()
        cell_cols <- names(corrected_data)[grepl(config$cell_pattern, names(corrected_data))]
        
        trace_stats <- sapply(cell_cols, function(cell) {
          trace <- corrected_data[[cell]]
          c(
            Mean = mean(trace, na.rm = TRUE),
            SD = sd(trace, na.rm = TRUE),
            Min = min(trace, na.rm = TRUE),
            Max = max(trace, na.rm = TRUE)
          )
        })
        
        stats$Mean_Amplitude <- trace_stats["Mean", ]
        stats$SD_Amplitude <- trace_stats["SD", ]
        stats$Min_Amplitude <- trace_stats["Min", ]
        stats$Max_Amplitude <- trace_stats["Max", ]
        
        stats
      },
      description = "Calculate summary statistics"
    ),
    
    # Visualization plots
    targets::tar_target(
      trace_plots,
      {
        config <- get_config()
        cell_cols <- names(corrected_data)[grepl(config$cell_pattern, names(corrected_data))]
        
        plots <- lapply(cell_cols, function(cell) {
          tryCatch({
            plot_cell_trace(
              corrected_data,
              cell,
              method = spike_method,
              show_spikes = TRUE,
              show_deconvolved = TRUE,
              verbose = FALSE
            )
          }, error = function(e) {
            warning("Failed to create plot for ", cell, ": ", e$message)
            NULL
          })
        })
        
        names(plots) <- cell_cols
        plots
      },
      description = "Generate trace plots for all cells"
    ),
    
    # Multi-cell plot
    targets::tar_target(
      multi_cell_plot,
      {
        tryCatch({
          plot_multiple_cells(
            corrected_data,
            ncol = min(3, length(names(corrected_data)[grepl(get_config()$cell_pattern, names(corrected_data))])),
            method = spike_method,
            show_spikes = TRUE,
            show_deconvolved = TRUE,
            verbose = FALSE
          )
        }, error = function(e) {
          warning("Failed to create multi-cell plot: ", e$message)
          NULL
        })
      },
      description = "Generate multi-cell overview plot"
    ),
    
    # Export results
    targets::tar_target(
      export_results,
      {
        # Create output directory
        output_dir <- "pipeline_results"
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        
        # Export corrected data
        utils::write.csv(corrected_data, 
                        file.path(output_dir, "corrected_data.csv"), 
                        row.names = FALSE)
        
        # Export summary statistics
        utils::write.csv(summary_stats, 
                        file.path(output_dir, "summary_statistics.csv"), 
                        row.names = FALSE)
        
        # Export spike results
        spike_data <- do.call(rbind, lapply(names(spike_results), function(cell) {
          spikes <- spike_results[[cell]]$spikes
          spikes$cell <- cell
          spikes$time <- seq_len(nrow(spikes))
          spikes
        }))
        utils::write.csv(spike_data, 
                        file.path(output_dir, "spike_results.csv"), 
                        row.names = FALSE)
        
        # Save plots
        for (cell in names(trace_plots)) {
          if (!is.null(trace_plots[[cell]])) {
            ggplot2::ggsave(
              file.path(output_dir, paste0("trace_", cell, ".png")),
              trace_plots[[cell]],
              width = 10, height = 6, dpi = 300
            )
          }
        }
        
        if (!is.null(multi_cell_plot)) {
          ggplot2::ggsave(
            file.path(output_dir, "multi_cell_overview.png"),
            multi_cell_plot,
            width = 15, height = 10, dpi = 300
          )
        }
        
        list(
          output_directory = output_dir,
          files_created = list.files(output_dir, full.names = TRUE)
        )
      },
      description = "Export all results to files"
    )
  )
}
