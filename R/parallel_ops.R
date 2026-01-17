#' Parallel Processing Operations
#'
#' Parallelization wrappers using future/furrr for calcium imaging analysis.
#' These functions enable efficient multi-core processing of large datasets.
#'
#' @name parallel_ops
#' @keywords internal
NULL

#' Configure parallel processing backend
#'
#' Sets up the parallel processing strategy using the future package.
#'
#' @param strategy Parallel strategy: "sequential", "multisession", "multicore", "cluster"
#' @param workers Number of workers (default: availableCores() - 1)
#' @param memory_limit Per-worker memory limit in MB (optional)
#'
#' @return Invisible previous plan
#' @export
#'
#' @examples
#' \dontrun{
#' setup_parallel("multisession", workers = 4)
#' # ... do parallel work ...
#' setup_parallel("sequential")  # Reset to sequential
#' }
setup_parallel <- function(strategy = "multisession", workers = NULL,
                           memory_limit = NULL) {
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("future package required. Install with: install.packages('future')")
  }

  if (is.null(workers)) {
    workers <- max(1, future::availableCores() - 1)
  }

  old_plan <- future::plan()

  if (strategy == "sequential") {
    future::plan(future::sequential)
  } else if (strategy == "multisession") {
    future::plan(future::multisession, workers = workers)
  } else if (strategy == "multicore") {
    if (.Platform$OS.type == "windows") {
      warning("multicore not available on Windows, using multisession")
      future::plan(future::multisession, workers = workers)
    } else {
      future::plan(future::multicore, workers = workers)
    }
  } else if (strategy == "cluster") {
    future::plan(future::cluster, workers = workers)
  } else {
    stop("Unknown strategy: ", strategy)
  }

  message(sprintf("Parallel backend: %s with %d workers", strategy, workers))
  invisible(old_plan)
}

#' Check current parallel configuration
#'
#' @return List with current parallel settings
#' @export
parallel_info <- function() {
  if (!requireNamespace("future", quietly = TRUE)) {
    return(list(available = FALSE))
  }

  plan <- future::plan()
  workers <- future::nbrOfWorkers()

  list(
    available = TRUE,
    plan = class(plan)[1],
    workers = workers,
    cores = future::availableCores()
  )
}

#' Apply function to traces in parallel
#'
#' Process each cell's trace in parallel using furrr.
#'
#' @param traces Matrix of traces (cells x time)
#' @param fun Function to apply to each row
#' @param ... Additional arguments to fun
#' @param .progress Show progress bar
#'
#' @return List or matrix of results
#' @export
#'
#' @examples
#' \dontrun{
#' traces <- matrix(rnorm(100 * 1000), nrow = 100)
#' setup_parallel("multisession", workers = 4)
#' results <- parallel_apply_traces(traces, function(x) smooth(x, 5))
#' }
parallel_apply_traces <- function(traces, fun, ..., .progress = TRUE) {
  if (!requireNamespace("furrr", quietly = TRUE)) {
    warning("furrr package not available, using sequential processing")
    return(lapply(seq_len(nrow(traces)), function(i) fun(traces[i, ], ...)))
  }

  # Convert to list of rows
  trace_list <- lapply(seq_len(nrow(traces)), function(i) traces[i, ])

  # Apply in parallel
  furrr::future_map(
    trace_list,
    fun,
    ...,
    .progress = .progress,
    .options = furrr::furrr_options(seed = TRUE)
  )
}

#' Parallel motion correction
#'
#' Correct motion in image stacks using parallel processing.
#'
#' @param movie 3D array (height x width x frames)
#' @param template Reference frame or "mean"
#' @param max_shift Maximum shift in pixels
#' @param batch_size Number of frames per batch
#' @param .progress Show progress bar
#'
#' @return List with corrected movie and shifts
#' @export
parallel_motion_correct <- function(movie, template = "mean", max_shift = 20,
                                    batch_size = 100, .progress = TRUE) {
  if (!requireNamespace("furrr", quietly = TRUE)) {
    warning("furrr package not available, using sequential processing")
    return(motion_correct(movie, template = template, max_shift = max_shift))
  }

  dims <- dim(movie)
  n_frames <- dims[3]

  # Create template
  if (identical(template, "mean")) {
    template <- apply(movie[, , 1:min(100, n_frames)], c(1, 2), mean)
  }

  # Split into batches
  batch_indices <- split(seq_len(n_frames), ceiling(seq_len(n_frames) / batch_size))

  # Process batches in parallel
  batch_results <- furrr::future_map(
    batch_indices,
    function(idx) {
      batch <- movie[, , idx, drop = FALSE]
      shifts <- matrix(0, nrow = length(idx), ncol = 2)

      corrected <- array(0, dim = c(dims[1], dims[2], length(idx)))

      for (i in seq_along(idx)) {
        frame <- batch[, , i]
        shift <- phase_correlation(template, frame, max_shift = max_shift)
        corrected[, , i] <- shift_image(frame, shift)
        shifts[i, ] <- shift
      }

      list(corrected = corrected, shifts = shifts)
    },
    .progress = .progress,
    .options = furrr::furrr_options(seed = TRUE)
  )

  # Combine results
  corrected <- array(0, dim = dims)
  all_shifts <- matrix(0, nrow = n_frames, ncol = 2)

  for (i in seq_along(batch_indices)) {
    idx <- batch_indices[[i]]
    corrected[, , idx] <- batch_results[[i]]$corrected
    all_shifts[idx, ] <- batch_results[[i]]$shifts
  }

  list(
    corrected = corrected,
    shifts = all_shifts,
    template = template
  )
}

#' Parallel spike detection
#'
#' Detect spikes in traces using parallel processing.
#'
#' @param traces Matrix of traces (cells x time)
#' @param method Detection method
#' @param ... Additional arguments to detection function
#' @param .progress Show progress bar
#'
#' @return Matrix of spike predictions
#' @export
parallel_spike_detection <- function(traces, method = "threshold", ...,
                                     .progress = TRUE) {
  if (!requireNamespace("furrr", quietly = TRUE)) {
    warning("furrr package not available, using sequential processing")
    return(threshold_spike_detection(traces, ...)$spike_predictions)
  }

  # Convert to list of rows
  trace_list <- lapply(seq_len(nrow(traces)), function(i) traces[i, ])

  # Process in parallel
  results <- furrr::future_map(
    trace_list,
    function(trace) {
      if (method == "threshold") {
        result <- threshold_spike_detection(matrix(trace, nrow = 1), ...)
      } else if (method == "adaptive") {
        result <- adaptive_threshold_detection(matrix(trace, nrow = 1), ...)
      } else {
        result <- statistical_spike_inference(trace, method = method, ...)
      }

      if (is.matrix(result$spike_predictions)) {
        result$spike_predictions[1, ]
      } else {
        result$spikes
      }
    },
    .progress = .progress,
    .options = furrr::furrr_options(seed = TRUE)
  )

  # Combine into matrix
  do.call(rbind, results)
}

#' Parallel trace extraction
#'
#' Extract traces from ROIs in parallel.
#'
#' @param movie 3D array (height x width x frames)
#' @param rois List of ROI masks
#' @param method Extraction method: "mean", "weighted", "median"
#' @param .progress Show progress bar
#'
#' @return Matrix of traces (cells x time)
#' @export
parallel_extract_traces <- function(movie, rois, method = "mean",
                                    .progress = TRUE) {
  if (!requireNamespace("furrr", quietly = TRUE)) {
    warning("furrr package not available, using sequential processing")
    return(extract_traces(movie, rois, method = method))
  }

  n_frames <- dim(movie)[3]

  # Process ROIs in parallel
  traces <- furrr::future_map(
    rois,
    function(roi) {
      if (is.list(roi) && !is.null(roi$mask)) {
        mask <- roi$mask
      } else if (is.matrix(roi)) {
        mask <- roi
      } else {
        stop("Invalid ROI format")
      }

      mask_idx <- which(mask > 0)

      trace <- vapply(seq_len(n_frames), function(t) {
        vals <- movie[, , t][mask_idx]
        if (method == "mean") {
          mean(vals, na.rm = TRUE)
        } else if (method == "median") {
          median(vals, na.rm = TRUE)
        } else if (method == "weighted") {
          weights <- mask[mask_idx]
          sum(vals * weights) / sum(weights)
        }
      }, numeric(1))

      trace
    },
    .progress = .progress,
    .options = furrr::furrr_options(seed = TRUE)
  )

  do.call(rbind, traces)
}

#' Parallel cross-correlation
#'
#' Compute pairwise cross-correlations in parallel.
#'
#' @param traces Matrix of traces (cells x time)
#' @param max_lag Maximum lag for cross-correlation
#' @param .progress Show progress bar
#'
#' @return Array of cross-correlations (cells x cells x lags)
#' @export
parallel_cross_correlation <- function(traces, max_lag = 10, .progress = TRUE) {
  if (!requireNamespace("furrr", quietly = TRUE)) {
    warning("furrr package not available, using sequential processing")
    # Fall back to sequential
    n_cells <- nrow(traces)
    lags <- seq(-max_lag, max_lag)
    result <- array(0, dim = c(n_cells, n_cells, length(lags)))

    for (i in seq_len(n_cells)) {
      for (j in seq_len(n_cells)) {
        cc <- ccf(traces[i, ], traces[j, ], lag.max = max_lag, plot = FALSE)
        result[i, j, ] <- cc$acf
      }
    }
    return(result)
  }

  n_cells <- nrow(traces)
  lags <- seq(-max_lag, max_lag)

  # Create pairs
  pairs <- expand.grid(i = seq_len(n_cells), j = seq_len(n_cells))
  pairs <- pairs[pairs$i <= pairs$j, ]  # Only upper triangle

  # Compute in parallel
  cc_results <- furrr::future_map2(
    pairs$i, pairs$j,
    function(i, j) {
      cc <- ccf(traces[i, ], traces[j, ], lag.max = max_lag, plot = FALSE)
      list(i = i, j = j, cc = as.numeric(cc$acf))
    },
    .progress = .progress,
    .options = furrr::furrr_options(seed = TRUE)
  )

  # Combine into array
  result <- array(0, dim = c(n_cells, n_cells, length(lags)))

  for (res in cc_results) {
    result[res$i, res$j, ] <- res$cc
    result[res$j, res$i, ] <- rev(res$cc)  # Symmetric
  }

  result
}

#' Parallel batch processing
#'
#' Process multiple files or sessions in parallel.
#'
#' @param file_paths Character vector of file paths
#' @param process_fun Function to apply to each file
#' @param ... Additional arguments to process_fun
#' @param .progress Show progress bar
#'
#' @return List of results
#' @export
parallel_batch_process <- function(file_paths, process_fun, ...,
                                   .progress = TRUE) {
  if (!requireNamespace("furrr", quietly = TRUE)) {
    warning("furrr package not available, using sequential processing")
    return(lapply(file_paths, process_fun, ...))
  }

  furrr::future_map(
    file_paths,
    process_fun,
    ...,
    .progress = .progress,
    .options = furrr::furrr_options(seed = TRUE)
  )
}

#' Shut down parallel workers
#'
#' Returns to sequential processing and frees resources.
#'
#' @export
shutdown_parallel <- function() {
  if (requireNamespace("future", quietly = TRUE)) {
    future::plan(future::sequential)
    message("Parallel workers shut down, returned to sequential processing")
  }
  invisible(NULL)
}
