#' Calcium Imaging Analysis Pipeline using targets
#' 
#' Defines a reproducible workflow for calcium imaging analysis using
#' the `targets` package. This pipeline includes data generation,
#' correction, spike inference, and visualization.
#' 
#' @param raw_data_path Path to raw data file
#' @param output_dir Output directory for results
#' @param config Configuration parameters
#' @param ... Additional arguments
#' @return targets pipeline object
#' 
#' @examples
#' # Basic pipeline
#' pipeline <- calcium_pipeline()
#' 
#' # Custom parameters
#' pipeline <- calcium_pipeline(
#'   raw_data_path = "path_to_custom_data.csv",
#'   config = list(
#'     correction_method = "legacy",
#'     spike_method = "caiman"
#'   )
#' )
#' 
#' @export
calcium_pipeline <- function(raw_data_path = NULL, 
                           output_dir = "results", 
                           config = NULL, 
                           ...) {
  
  # Check if targets package is available
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required for the pipeline. Install with: install.packages('targets')")
  }
  
  # Set default configuration
  if (is.null(config)) {
    config <- get_config()
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define pipeline targets
  pipeline <- list(
    # Load and validate data
    raw_data = targets::tar_target(
      name = raw_data,
      command = {
        if (is.null(raw_data_path)) {
          # Generate synthetic data if no path provided
          generate_synthetic_data(
            n_cells = 10,
            n_time = 500,
            spike_prob = 0.1
          )
        } else {
          # Load from file
          validate_calcium_data(raw_data_path)
        }
      }
    ),
    
    # Calcium correction
    corrected_data = targets::tar_target(
      name = corrected_data,
      command = {
        calcium_correction(
          raw_data$traces, 
          method = config$correction_method,
          span = config$default_span
        )
      }
    ),
    
    # Spike inference
    spike_results = targets::tar_target(
      name = spike_results,
      command = {
        infer_spikes(
          corrected_data, 
          method = config$spike_method
        )
      }
    ),
    
    # Quality control
    qc_report = targets::tar_target(
      name = qc_report,
      command = {
        generate_qc_report(
          data = raw_data$traces
        )
      }
    ),
    
    # Summary statistics
    summary_stats = targets::tar_target(
      name = summary_stats,
      command = {
        list(
          n_cells = nrow(corrected_data),
          n_timepoints = ncol(corrected_data),
          mean_spike_rate = mean(spike_results$spike_rates, na.rm = TRUE),
          correction_quality = calculate_correction_quality(raw_data$traces, corrected_data)
        )
      }
    ),
    
    # Generate plots
    trace_plots = targets::tar_target(
      name = trace_plots,
      command = {
        lapply(1:min(5, nrow(corrected_data)), function(i) {
          plot_cell_trace(
            corrected_df = data.frame(
              Time = 1:ncol(corrected_data),
              Cell = corrected_data[i, ]
            ),
            cell = "Cell"
          )
        })
      }
    ),
    
    multi_cell_plot = targets::tar_target(
      name = multi_cell_plot,
      command = {
        multi_panel_calcium_plot(
          traces = corrected_data[1:min(9, nrow(corrected_data)), ],
          spike_results = spike_results$spikes[1:min(9, nrow(corrected_data)), ]
        )
      }
    ),
    
    # Export results
    export_results = targets::tar_target(
      name = export_results,
      command = {
        export_processed_data(
          corrected_data = corrected_data,
          spike_results = spike_results,
          output_dir = output_dir
        )
      }
    )
  )
  
  return(pipeline)
}
