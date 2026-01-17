#' Neuropil Correction Functions
#'
#' Functions for correcting neuropil contamination in calcium imaging data.
#'
#' @name neuropil_correction
NULL

#' Neuropil Correction
#'
#' Correct for neuropil (background) contamination in calcium traces.
#' This is essential for accurate signal extraction in densely labeled tissue.
#'
#' @param cell_traces Matrix or data frame of cell traces (time x cells)
#' @param neuropil_traces Matrix or data frame of neuropil traces (time x cells)
#' @param method Correction method ("coefficient", "regression", "robust", "adaptive")
#' @param coefficient Fixed neuropil coefficient (for "coefficient" method)
#' @param min_coefficient Minimum allowed coefficient (for "regression"/"adaptive")
#' @param max_coefficient Maximum allowed coefficient (for "regression"/"adaptive")
#' @param baseline_percentile Percentile for baseline estimation
#' @param verbose Show progress messages
#' @return List with corrected traces and coefficients
#'
#' @details
#' Neuropil contamination occurs because each ROI captures signal not just from
#' the target cell but also from surrounding neuropil (dendrites, axons, etc.).
#'
#' Methods:
#' \describe{
#'   \item{coefficient}{Simple subtraction: F_corrected = F_cell - coef * F_neuropil}
#'   \item{regression}{Estimate coefficient by regression on baseline periods}
#'   \item{robust}{Robust regression to handle outliers}
#'   \item{adaptive}{Time-varying coefficient estimation}
#' }
#'
#' Typical coefficient values range from 0.5-0.9, with 0.7 being common.
#'
#' @examples
#' # Simple correction with fixed coefficient
#' result <- neuropil_correct(cell_traces, neuropil_traces, coefficient = 0.7)
#'
#' # Regression-based coefficient estimation
#' result <- neuropil_correct(cell_traces, neuropil_traces, method = "regression")
#'
#' # Access results
#' corrected <- result$corrected
#' coefficients <- result$coefficients
#'
#' @export
neuropil_correct <- function(cell_traces,
                             neuropil_traces,
                             method = c("coefficient", "regression", "robust", "adaptive"),
                             coefficient = 0.7,
                             min_coefficient = 0.3,
                             max_coefficient = 1.0,
                             baseline_percentile = 8,
                             verbose = TRUE) {

  method <- match.arg(method)

  # Convert to matrices if needed
  if (is.data.frame(cell_traces)) {
    cell_mat <- as.matrix(cell_traces[, sapply(cell_traces, is.numeric)])
  } else {
    cell_mat <- as.matrix(cell_traces)
  }

  if (is.data.frame(neuropil_traces)) {
    np_mat <- as.matrix(neuropil_traces[, sapply(neuropil_traces, is.numeric)])
  } else {
    np_mat <- as.matrix(neuropil_traces)
  }

  # Validate dimensions
  if (nrow(cell_mat) != nrow(np_mat)) {
    stop("Cell and neuropil traces must have the same number of time points")
  }

  if (ncol(cell_mat) != ncol(np_mat)) {
    stop("Cell and neuropil traces must have the same number of cells")
  }

  n_time <- nrow(cell_mat)
  n_cells <- ncol(cell_mat)

  if (verbose) {
    message("Neuropil correction: ", n_cells, " cells, ", n_time, " time points")
    message("Method: ", method)
  }

  # Apply correction based on method
  result <- switch(method,
    "coefficient" = neuropil_correct_coefficient(cell_mat, np_mat, coefficient, verbose),
    "regression" = neuropil_correct_regression(cell_mat, np_mat, min_coefficient,
                                               max_coefficient, baseline_percentile, verbose),
    "robust" = neuropil_correct_robust(cell_mat, np_mat, min_coefficient,
                                       max_coefficient, baseline_percentile, verbose),
    "adaptive" = neuropil_correct_adaptive(cell_mat, np_mat, min_coefficient,
                                           max_coefficient, verbose)
  )

  # Convert back to data frame if input was data frame
  if (is.data.frame(cell_traces)) {
    result$corrected <- as.data.frame(result$corrected)
    colnames(result$corrected) <- colnames(cell_mat)
  }

  result$method <- method
  result$n_cells <- n_cells
  result$n_time <- n_time

  return(result)
}

#' Fixed Coefficient Neuropil Correction
#' @keywords internal
neuropil_correct_coefficient <- function(cell_mat, np_mat, coefficient, verbose) {
  if (verbose) message("  Using fixed coefficient: ", coefficient)

  corrected <- cell_mat - coefficient * np_mat

  list(
    corrected = corrected,
    coefficients = rep(coefficient, ncol(cell_mat))
  )
}

#' Regression-based Neuropil Correction
#' @keywords internal
neuropil_correct_regression <- function(cell_mat, np_mat, min_coef, max_coef,
                                        baseline_percentile, verbose) {
  n_cells <- ncol(cell_mat)
  coefficients <- numeric(n_cells)
  corrected <- matrix(0, nrow(cell_mat), ncol(cell_mat))

  for (i in seq_len(n_cells)) {
    cell <- cell_mat[, i]
    np <- np_mat[, i]

    # Identify baseline periods (lowest percentile of cell activity)
    threshold <- quantile(cell, baseline_percentile / 100, na.rm = TRUE)
    baseline_idx <- which(cell <= threshold)

    if (length(baseline_idx) < 10) {
      # Not enough baseline points, use all data
      baseline_idx <- seq_along(cell)
    }

    # Regress cell on neuropil during baseline
    fit <- lm(cell[baseline_idx] ~ np[baseline_idx])
    coef <- coef(fit)[2]

    # Constrain coefficient
    coef <- max(min_coef, min(max_coef, coef))
    coefficients[i] <- coef

    # Apply correction
    corrected[, i] <- cell - coef * np

    if (verbose && i %% 50 == 0) message("  Processed ", i, "/", n_cells, " cells")
  }

  if (verbose) {
    message("  Coefficient range: ", round(min(coefficients), 3), " - ",
            round(max(coefficients), 3))
    message("  Mean coefficient: ", round(mean(coefficients), 3))
  }

  list(
    corrected = corrected,
    coefficients = coefficients
  )
}

#' Robust Neuropil Correction
#' @keywords internal
neuropil_correct_robust <- function(cell_mat, np_mat, min_coef, max_coef,
                                    baseline_percentile, verbose) {
  n_cells <- ncol(cell_mat)
  coefficients <- numeric(n_cells)
  corrected <- matrix(0, nrow(cell_mat), ncol(cell_mat))

  # Check for MASS package
  use_robust <- requireNamespace("MASS", quietly = TRUE)
  if (!use_robust) {
    warning("MASS package not available, using standard regression")
  }

  for (i in seq_len(n_cells)) {
    cell <- cell_mat[, i]
    np <- np_mat[, i]

    # Identify baseline periods
    threshold <- quantile(cell, baseline_percentile / 100, na.rm = TRUE)
    baseline_idx <- which(cell <= threshold)

    if (length(baseline_idx) < 10) {
      baseline_idx <- seq_along(cell)
    }

    # Robust regression
    if (use_robust) {
      fit <- MASS::rlm(cell[baseline_idx] ~ np[baseline_idx])
      coef <- coef(fit)[2]
    } else {
      fit <- lm(cell[baseline_idx] ~ np[baseline_idx])
      coef <- coef(fit)[2]
    }

    # Constrain coefficient
    coef <- max(min_coef, min(max_coef, coef))
    coefficients[i] <- coef

    # Apply correction
    corrected[, i] <- cell - coef * np

    if (verbose && i %% 50 == 0) message("  Processed ", i, "/", n_cells, " cells")
  }

  if (verbose) {
    message("  Coefficient range: ", round(min(coefficients), 3), " - ",
            round(max(coefficients), 3))
  }

  list(
    corrected = corrected,
    coefficients = coefficients
  )
}

#' Adaptive Neuropil Correction
#' @keywords internal
neuropil_correct_adaptive <- function(cell_mat, np_mat, min_coef, max_coef, verbose) {
  n_time <- nrow(cell_mat)
  n_cells <- ncol(cell_mat)

  # Window size for adaptive estimation (e.g., 30 seconds at 30 Hz = 900 frames)
  window_size <- min(900, floor(n_time / 4))
  step_size <- floor(window_size / 2)

  coefficients_matrix <- matrix(0, n_time, n_cells)
  corrected <- matrix(0, n_time, n_cells)

  if (verbose) message("  Adaptive correction with window size ", window_size)

  for (i in seq_len(n_cells)) {
    cell <- cell_mat[, i]
    np <- np_mat[, i]

    # Estimate coefficient in sliding windows
    window_coefs <- numeric(0)
    window_centers <- numeric(0)

    for (start in seq(1, n_time - window_size + 1, by = step_size)) {
      end <- start + window_size - 1
      idx <- start:end

      # Use lowest 20% of cell activity in window
      threshold <- quantile(cell[idx], 0.2, na.rm = TRUE)
      baseline_idx <- idx[cell[idx] <= threshold]

      if (length(baseline_idx) >= 10) {
        fit <- lm(cell[baseline_idx] ~ np[baseline_idx])
        coef <- coef(fit)[2]
        coef <- max(min_coef, min(max_coef, coef))
      } else {
        coef <- 0.7  # Default
      }

      window_coefs <- c(window_coefs, coef)
      window_centers <- c(window_centers, (start + end) / 2)
    }

    # Interpolate coefficients to all time points
    if (length(window_coefs) > 1) {
      time_varying_coef <- approx(window_centers, window_coefs, xout = 1:n_time,
                                  rule = 2)$y
    } else {
      time_varying_coef <- rep(window_coefs[1], n_time)
    }

    coefficients_matrix[, i] <- time_varying_coef
    corrected[, i] <- cell - time_varying_coef * np

    if (verbose && i %% 50 == 0) message("  Processed ", i, "/", n_cells, " cells")
  }

  # Mean coefficient per cell
  mean_coefficients <- colMeans(coefficients_matrix)

  if (verbose) {
    message("  Mean coefficient range: ", round(min(mean_coefficients), 3), " - ",
            round(max(mean_coefficients), 3))
  }

  list(
    corrected = corrected,
    coefficients = mean_coefficients,
    coefficients_time_varying = coefficients_matrix
  )
}

#' Extract Neuropil Traces from Movie
#'
#' Extract neuropil (surround) traces for each ROI from a movie.
#'
#' @param movie 3D array [height x width x frames]
#' @param rois List of ROI masks or labeled matrix
#' @param inner_radius Radius to exclude around cell (pixels)
#' @param outer_radius Outer radius of neuropil annulus (pixels)
#' @param verbose Show progress
#' @return Matrix of neuropil traces (time x cells)
#'
#' @export
extract_neuropil_traces <- function(movie,
                                     rois,
                                     inner_radius = 3,
                                     outer_radius = 15,
                                     verbose = TRUE) {
  n_frames <- dim(movie)[3]
  height <- dim(movie)[1]
  width <- dim(movie)[2]

  # Convert labeled matrix to list if needed
  if (is.matrix(rois) && !is.logical(rois)) {
    labels <- unique(as.vector(rois))
    labels <- labels[labels > 0]
    roi_list <- lapply(labels, function(l) rois == l)
    names(roi_list) <- paste0("NP_", seq_along(roi_list))
    rois <- roi_list
  }

  n_rois <- length(rois)
  if (verbose) message("Extracting neuropil traces for ", n_rois, " ROIs")

  # Initialize output
  np_traces <- matrix(0, nrow = n_frames, ncol = n_rois)

  for (i in seq_along(rois)) {
    cell_mask <- rois[[i]]

    # Create neuropil annulus
    dilated_inner <- dilate_mask(cell_mask, inner_radius)
    dilated_outer <- dilate_mask(cell_mask, outer_radius)
    neuropil_mask <- dilated_outer & !dilated_inner

    # Exclude other ROIs from neuropil mask
    for (j in seq_along(rois)) {
      if (j != i) {
        other_roi <- dilate_mask(rois[[j]], inner_radius)
        neuropil_mask <- neuropil_mask & !other_roi
      }
    }

    # Extract neuropil trace
    if (sum(neuropil_mask) > 0) {
      for (t in seq_len(n_frames)) {
        np_traces[t, i] <- mean(movie[, , t][neuropil_mask], na.rm = TRUE)
      }
    } else {
      # Fallback: use broader annulus
      neuropil_mask <- dilate_outer & !cell_mask
      for (t in seq_len(n_frames)) {
        np_traces[t, i] <- mean(movie[, , t][neuropil_mask], na.rm = TRUE)
      }
    }

    if (verbose && i %% 20 == 0) message("  Processed ", i, "/", n_rois, " ROIs")
  }

  colnames(np_traces) <- if (!is.null(names(rois))) {
    paste0("NP_", gsub("Cell_", "", names(rois)))
  } else {
    paste0("NP_", seq_len(n_rois))
  }

  return(np_traces)
}

#' Estimate Optimal Neuropil Coefficient
#'
#' Estimate the optimal neuropil contamination coefficient using multiple methods.
#'
#' @param cell_trace Single cell trace
#' @param neuropil_trace Corresponding neuropil trace
#' @param methods Methods to use for estimation
#' @return List with coefficient estimates from each method
#'
#' @export
estimate_neuropil_coefficient <- function(cell_trace,
                                          neuropil_trace,
                                          methods = c("regression", "percentile", "correlation")) {
  results <- list()

  # Method 1: Regression on baseline
  if ("regression" %in% methods) {
    threshold <- quantile(cell_trace, 0.08, na.rm = TRUE)
    baseline_idx <- which(cell_trace <= threshold)

    if (length(baseline_idx) >= 10) {
      fit <- lm(cell_trace[baseline_idx] ~ neuropil_trace[baseline_idx])
      results$regression <- max(0, min(1, coef(fit)[2]))
    } else {
      results$regression <- NA
    }
  }

  # Method 2: Match percentiles
  if ("percentile" %in% methods) {
    # Find coefficient that makes corrected trace non-negative at low percentiles
    test_coefs <- seq(0.3, 1.0, by = 0.05)
    best_coef <- 0.7
    best_score <- Inf

    for (coef in test_coefs) {
      corrected <- cell_trace - coef * neuropil_trace
      pct_5 <- quantile(corrected, 0.05, na.rm = TRUE)

      # Score: how far is 5th percentile from zero
      score <- abs(pct_5)
      if (pct_5 >= 0 && score < best_score) {
        best_score <- score
        best_coef <- coef
      }
    }
    results$percentile <- best_coef
  }

  # Method 3: Minimize correlation in corrected trace
  if ("correlation" %in% methods) {
    # Good correction should remove correlation with neuropil
    test_coefs <- seq(0.3, 1.0, by = 0.05)
    best_coef <- 0.7
    best_cor <- 1

    for (coef in test_coefs) {
      corrected <- cell_trace - coef * neuropil_trace
      cor_val <- abs(cor(corrected, neuropil_trace, use = "complete.obs"))

      if (cor_val < best_cor) {
        best_cor <- cor_val
        best_coef <- coef
      }
    }
    results$correlation <- best_coef
  }

  # Summary
  estimates <- unlist(results[sapply(results, function(x) !is.na(x) && is.numeric(x))])
  results$mean <- mean(estimates, na.rm = TRUE)
  results$median <- median(estimates, na.rm = TRUE)

  return(results)
}

#' Plot Neuropil Correction Diagnostics
#'
#' Visualize the effect of neuropil correction.
#'
#' @param cell_trace Original cell trace
#' @param neuropil_trace Neuropil trace
#' @param corrected_trace Corrected cell trace
#' @param coefficient Neuropil coefficient used
#' @param time_range Optional time range to plot
#' @return ggplot object
#'
#' @export
plot_neuropil_correction <- function(cell_trace,
                                     neuropil_trace,
                                     corrected_trace,
                                     coefficient,
                                     time_range = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }

  n_time <- length(cell_trace)
  time <- 1:n_time

  if (!is.null(time_range)) {
    idx <- time_range[1]:min(time_range[2], n_time)
  } else {
    idx <- 1:min(2000, n_time)  # Show first 2000 points by default
  }

  df <- data.frame(
    time = rep(time[idx], 3),
    value = c(cell_trace[idx], neuropil_trace[idx], corrected_trace[idx]),
    trace = factor(rep(c("Cell (raw)", "Neuropil", "Cell (corrected)"), each = length(idx)),
                   levels = c("Cell (raw)", "Neuropil", "Cell (corrected)"))
  )

  ggplot2::ggplot(df, ggplot2::aes(x = time, y = value, color = trace)) +
    ggplot2::geom_line(alpha = 0.8, linewidth = 0.5) +
    ggplot2::facet_wrap(~trace, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = paste("Neuropil Correction (coefficient =", round(coefficient, 3), ")"),
      x = "Time",
      y = "Fluorescence",
      color = NULL
    ) +
    ggplot2::scale_color_manual(values = c("Cell (raw)" = "blue",
                                           "Neuropil" = "red",
                                           "Cell (corrected)" = "darkgreen")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}
