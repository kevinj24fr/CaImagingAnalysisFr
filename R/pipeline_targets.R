#' Calcium Imaging Analysis Pipeline using targets
#'
#' Defines a reproducible workflow for calcium imaging analysis using
#' the `targets` package. This pipeline includes data generation,
#' correction, spike inference, network analysis, and visualization.
#'
#' @param raw_data_path Path to raw data file (CSV or RDS). If NULL, synthetic data is generated.
#' @param output_dir Output directory for results (default: "results")
#' @param n_cells Number of cells for synthetic data (default: 10)
#' @param n_time Number of time points for synthetic data (default: 500)
#' @param correction_method Correction method: "modern" or "legacy" (default: "modern")
#' @param spike_method Spike inference method (default: "oasis")
#' @param include_network Whether to include network analysis (default: TRUE)
#' @param ... Additional arguments passed to pipeline steps
#' @return List of targets for use with targets::tar_make()
#'
#' @details
#' This function generates a list of targets that can be used with the targets

#' package for reproducible analysis. To run the pipeline:
#'
#' \enumerate{
#'   \item Create a _targets.R file in your project root
#'   \item Call \code{calcium_pipeline()} to get the target list
#'   \item Run \code{targets::tar_make()} to execute the pipeline
#' }
#'
#' @examples
#' \dontrun{
#' # In your _targets.R file:
#' library(targets)
#' library(CaImagingAnalysisFr)
#'
#' # Basic pipeline with synthetic data
#' list(
#'   calcium_pipeline()
#' )
#'
#' # Pipeline with custom data
#' list(
#'   calcium_pipeline(
#'     raw_data_path = "data/my_calcium_data.csv",
#'     correction_method = "modern",
#'     spike_method = "oasis"
#'   )
#' )
#' }
#'
#' @export
calcium_pipeline <- function(raw_data_path = NULL,
                             output_dir = "results",
                             n_cells = 10,
                             n_time = 500,
                             correction_method = "modern",
                             spike_method = "oasis",
                             include_network = TRUE,
                             ...) {

  # Check if targets package is available
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required for the pipeline. Install with: install.packages('targets')")
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Store pipeline parameters for use in targets
  # These are captured as values, not references
  .raw_data_path <- raw_data_path
  .output_dir <- output_dir
  .n_cells <- n_cells
  .n_time <- n_time
  .correction_method <- correction_method
  .spike_method <- spike_method

  # Build the pipeline target list
  pipeline_targets <- list(

    # Target 1: Load or generate raw data
    targets::tar_target(
      name = raw_data,
      command = {
        if (is.null(.raw_data_path)) {
          # Generate synthetic data
          CaImagingAnalysisFr::generate_synthetic_data(
            n_cells = .n_cells,
            n_time = .n_time,
            spike_prob = 0.05
          )
        } else if (grepl("\\.csv$", .raw_data_path, ignore.case = TRUE)) {
          utils::read.csv(.raw_data_path, stringsAsFactors = FALSE)
        } else if (grepl("\\.rds$", .raw_data_path, ignore.case = TRUE)) {
          readRDS(.raw_data_path)
        } else {
          stop("Unsupported file format. Use CSV or RDS.")
        }
      }
    ),

    # Target 2: Validate raw data
    targets::tar_target(
      name = validation_results,
      command = {
        CaImagingAnalysisFr::validate_calcium_data(
          raw_data,
          data_type = "calcium_traces",
          output_format = "summary"
        )
      }
    ),

    # Target 3: Calcium correction
    targets::tar_target(
      name = corrected_data,
      command = {
        CaImagingAnalysisFr::calcium_correction(
          raw_data,
          method = .correction_method,
          verbose = TRUE
        )
      }
    ),

    # Target 4: Spike inference for each cell
    targets::tar_target(
      name = spike_results,
      command = {
        config <- CaImagingAnalysisFr::get_config()
        cell_cols <- names(corrected_data)[grepl(config$cell_pattern, names(corrected_data))]

        # Infer spikes for each cell
        results <- lapply(cell_cols, function(cell) {
          tryCatch({
            CaImagingAnalysisFr::infer_spikes(
              corrected_data[[cell]],
              method = .spike_method
            )
          }, error = function(e) {
            list(
              spikes = rep(0, nrow(corrected_data)),
              fit = corrected_data[[cell]],
              error = e$message
            )
          })
        })
        names(results) <- cell_cols

        # Combine into summary
        list(
          by_cell = results,
          total_spikes = sum(sapply(results, function(x) sum(x$spikes > 0, na.rm = TRUE))),
          spike_rates = sapply(results, function(x) mean(x$spikes > 0, na.rm = TRUE))
        )
      }
    ),

    # Target 5: Quality metrics
    targets::tar_target(
      name = quality_metrics,
      command = {
        CaImagingAnalysisFr::calculate_correction_quality(raw_data, corrected_data)
      }
    ),

    # Target 6: Summary statistics
    targets::tar_target(
      name = summary_stats,
      command = {
        config <- CaImagingAnalysisFr::get_config()
        cell_cols <- names(corrected_data)[grepl(config$cell_pattern, names(corrected_data))]

        list(
          n_cells = length(cell_cols),
          n_timepoints = nrow(corrected_data),
          total_spikes = spike_results$total_spikes,
          mean_spike_rate = mean(spike_results$spike_rates, na.rm = TRUE),
          quality_score = if (!is.null(quality_metrics$overall_quality))
            quality_metrics$overall_quality else NA,
          validation_passed = validation_results$validation_passed
        )
      }
    ),

    # Target 7: Generate trace plots
    targets::tar_target(
      name = trace_plots,
      command = {
        config <- CaImagingAnalysisFr::get_config()
        cell_cols <- names(corrected_data)[grepl(config$cell_pattern, names(corrected_data))]
        cells_to_plot <- head(cell_cols, 5)

        lapply(cells_to_plot, function(cell) {
          CaImagingAnalysisFr::plot_cell_trace(corrected_data, cell)
        })
      }
    ),

    # Target 8: Export results
    targets::tar_target(
      name = exported_files,
      command = {
        # Create output directory if needed
        if (!dir.exists(.output_dir)) {
          dir.create(.output_dir, recursive = TRUE)
        }

        # Export corrected data
        corrected_file <- file.path(.output_dir, "corrected_data.csv")
        utils::write.csv(corrected_data, corrected_file, row.names = FALSE)

        # Export spike results
        spike_file <- file.path(.output_dir, "spike_results.rds")
        saveRDS(spike_results, spike_file)

        # Export summary
        summary_file <- file.path(.output_dir, "analysis_summary.rds")
        saveRDS(summary_stats, summary_file)

        list(
          corrected_data = corrected_file,
          spike_results = spike_file,
          summary = summary_file
        )
      }
    )
  )

  # Add network analysis targets if requested
  if (include_network) {
    network_targets <- list(

      # Target N1: Functional connectivity
      targets::tar_target(
        name = network_connectivity,
        command = {
          config <- CaImagingAnalysisFr::get_config()
          cell_cols <- names(corrected_data)[grepl(config$cell_pattern, names(corrected_data))]

          # Extract cell traces as matrix
          trace_matrix <- as.matrix(corrected_data[, cell_cols])

          CaImagingAnalysisFr::functional_connectivity(
            trace_matrix,
            method = "correlation",
            threshold = 0.3
          )
        }
      ),

      # Target N2: Network metrics
      targets::tar_target(
        name = network_metrics,
        command = {
          if (!is.null(network_connectivity$adjacency_matrix)) {
            CaImagingAnalysisFr::graph_metrics(
              network_connectivity$adjacency_matrix
            )
          } else {
            list(error = "No adjacency matrix available")
          }
        }
      ),

      # Target N3: Export network results
      targets::tar_target(
        name = network_export,
        command = {
          network_file <- file.path(.output_dir, "network_results.rds")
          saveRDS(list(
            connectivity = network_connectivity,
            metrics = network_metrics
          ), network_file)
          network_file
        }
      )
    )

    pipeline_targets <- c(pipeline_targets, network_targets)
  }

  return(pipeline_targets)
}


#' Create a Standalone _targets.R File
#'
#' Generates a _targets.R file for running the calcium imaging pipeline.
#'
#' @param output_file Path to write the _targets.R file (default: "_targets.R")
#' @param raw_data_path Path to raw data file (optional)
#' @param ... Additional parameters passed to calcium_pipeline
#' @return Path to the created file (invisibly)
#'
#' @examples
#' \dontrun{
#' # Create a basic _targets.R file
#' create_targets_file()
#'
#' # Create with custom data path
#' create_targets_file(raw_data_path = "data/experiment1.csv")
#' }
#'
#' @export
create_targets_file <- function(output_file = "_targets.R",
                                raw_data_path = NULL,
                                ...) {

  # Capture additional arguments
  extra_args <- list(...)

  # Build the parameter string
  param_lines <- c()
  if (!is.null(raw_data_path)) {
    param_lines <- c(param_lines, sprintf('  raw_data_path = "%s"', raw_data_path))
  }
  for (name in names(extra_args)) {
    val <- extra_args[[name]]
    if (is.character(val)) {
      param_lines <- c(param_lines, sprintf('  %s = "%s"', name, val))
    } else if (is.logical(val)) {
      param_lines <- c(param_lines, sprintf('  %s = %s', name, toupper(as.character(val))))
    } else {
      param_lines <- c(param_lines, sprintf('  %s = %s', name, val))
    }
  }

  params_str <- if (length(param_lines) > 0) {
    paste0("(\n", paste(param_lines, collapse = ",\n"), "\n)")
  } else {
    "()"
  }

  content <- c(
    "# _targets.R file for Calcium Imaging Analysis Pipeline",
    "# Generated by CaImagingAnalysisFr::create_targets_file()",
    sprintf("# Created: %s", Sys.time()),
    "",
    "# Load required packages",
    "library(targets)",
    "library(CaImagingAnalysisFr)",
    "",
    "# Set target options",
    "tar_option_set(",
    '  packages = c("CaImagingAnalysisFr", "ggplot2"),',
    "  format = \"rds\"",
    ")",
    "",
    "# Define the pipeline",
    paste0("calcium_pipeline", params_str)
  )

  writeLines(content, output_file)
  message("Created targets file: ", output_file)
  message("Run the pipeline with: targets::tar_make()")
  invisible(output_file)
}


#' Run Calcium Imaging Pipeline (Standalone Mode)
#'
#' Run a complete calcium imaging analysis pipeline without using targets.
#' This is useful for quick analyses or when reproducibility tracking isn't needed.
#'
#' @param raw_data Data frame with calcium traces, or path to data file
#' @param correction_method Correction method: "modern" or "legacy"
#' @param spike_method Spike inference method
#' @param include_network Whether to include network analysis
#' @param include_plots Whether to generate plots
#' @param verbose Whether to show progress messages
#' @param ... Additional arguments
#' @return List containing all analysis results
#'
#' @examples
#' # Quick analysis with synthetic data
#' results <- run_calcium_analysis()
#'
#' # Analysis with custom data
#' data <- generate_synthetic_data(n_cells = 20, n_time = 1000)
#' results <- run_calcium_analysis(data, spike_method = "oasis")
#'
#' @export
run_calcium_analysis <- function(raw_data = NULL,
                                  correction_method = "modern",
                                  spike_method = "oasis",
                                  include_network = TRUE,
                                  include_plots = TRUE,
                                  verbose = TRUE,
                                  ...) {

  # Step 1: Load or generate data
  if (verbose) message("Step 1/6: Loading data...")
  if (is.null(raw_data)) {
    raw_data <- generate_synthetic_data(n_cells = 10, n_time = 500)
  } else if (is.character(raw_data)) {
    if (grepl("\\.csv$", raw_data, ignore.case = TRUE)) {
      raw_data <- utils::read.csv(raw_data, stringsAsFactors = FALSE)
    } else if (grepl("\\.rds$", raw_data, ignore.case = TRUE)) {
      raw_data <- readRDS(raw_data)
    }
  }

  # Step 2: Validate data
  if (verbose) message("Step 2/6: Validating data...")
  validation <- validate_calcium_data(raw_data, output_format = "summary")

  # Step 3: Calcium correction
  if (verbose) message("Step 3/6: Correcting calcium traces...")
  corrected <- calcium_correction(raw_data, method = correction_method, verbose = verbose)

  # Step 4: Spike inference
  if (verbose) message("Step 4/6: Inferring spikes...")
  config <- get_config()
  cell_cols <- names(corrected)[grepl(config$cell_pattern, names(corrected))]

  spike_results <- lapply(cell_cols, function(cell) {
    tryCatch({
      infer_spikes(corrected[[cell]], method = spike_method)
    }, error = function(e) {
      list(spikes = rep(0, nrow(corrected)), fit = corrected[[cell]], error = e$message)
    })
  })
  names(spike_results) <- cell_cols

  # Step 5: Network analysis (optional)
  network_results <- NULL
  if (include_network && length(cell_cols) > 2) {
    if (verbose) message("Step 5/6: Computing network connectivity...")
    trace_matrix <- as.matrix(corrected[, cell_cols])
    network_results <- tryCatch({
      conn <- functional_connectivity(trace_matrix, method = "correlation", threshold = 0.3)
      metrics <- if (!is.null(conn$adjacency_matrix)) {
        graph_metrics(conn$adjacency_matrix)
      } else {
        NULL
      }
      list(connectivity = conn, metrics = metrics)
    }, error = function(e) {
      list(error = e$message)
    })
  } else {
    if (verbose) message("Step 5/6: Skipping network analysis...")
  }

  # Step 6: Generate plots (optional)
  plots <- NULL
  if (include_plots) {
    if (verbose) message("Step 6/6: Generating plots...")
    plots <- tryCatch({
      cells_to_plot <- head(cell_cols, 3)
      lapply(cells_to_plot, function(cell) {
        plot_cell_trace(corrected, cell)
      })
    }, error = function(e) {
      list(error = e$message)
    })
  } else {
    if (verbose) message("Step 6/6: Skipping plot generation...")
  }

  # Compile summary
  summary_stats <- list(
    n_cells = length(cell_cols),
    n_timepoints = nrow(corrected),
    total_spikes = sum(sapply(spike_results, function(x) sum(x$spikes > 0, na.rm = TRUE))),
    mean_spike_rate = mean(sapply(spike_results, function(x) mean(x$spikes > 0, na.rm = TRUE))),
    validation_passed = validation$validation_passed
  )

  if (verbose) message("Analysis complete!")

  list(
    raw_data = raw_data,
    validation = validation,
    corrected_data = corrected,
    spike_results = spike_results,
    network_results = network_results,
    plots = plots,
    summary = summary_stats,
    parameters = list(
      correction_method = correction_method,
      spike_method = spike_method,
      include_network = include_network
    )
  )
}
