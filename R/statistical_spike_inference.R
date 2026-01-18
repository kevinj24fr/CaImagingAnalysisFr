#' Statistical Spike Inference
#'
#' Statistical methods for spike inference from calcium traces.
#' These functions use classical signal processing approaches (convolution,
#' rolling statistics, correlation) rather than deep learning.
#'
#' @name statistical_spike_inference
NULL

#' Statistical Spike Inference Interface
#'
#' Unified interface for different statistical spike inference methods.
#' These methods use signal processing techniques rather than neural networks.
#'
#' @param trace Calcium trace (numeric vector)
#' @param method Method to use: "rolling" (rolling statistics), "convolution"
#'   (1D convolution filters), or "correlation" (correlation-weighted)
#' @param ... Additional arguments passed to the specific method
#' @return List containing spike times and inferred calcium
#' @keywords internal
#' @examples
#' \dontrun{
#' trace <- rnorm(1000) + cumsum(rnorm(1000, sd = 0.1))
#' result <- statistical_spike_inference(trace, method = "rolling")
#' plot(trace, type = "l")
#' points(result$spike_times, trace[result$spike_times], col = "red", pch = 19)
#' }
statistical_spike_inference <- function(trace, method = c("rolling", "convolution", "correlation"), ...) {
  method <- match.arg(method)
  message("Running statistical spike inference (", method, " method)")

  if (method == "rolling") {
    return(rolling_stats_inference(trace, ...))
  } else if (method == "convolution") {
    return(convolutional_inference(trace, ...))
  } else if (method == "correlation") {
    return(correlation_weighted_inference(trace, ...))
  } else {
    stop("Unknown method: ", method)
  }
}

#' Rolling Statistics Spike Inference
#'
#' Infer spikes using rolling window statistics. Detects events that
#' exceed a locally-adaptive threshold based on rolling mean and SD.
#'
#' @param trace Calcium trace (numeric vector)
#' @param window_size Window size for rolling statistics (default: 20)
#' @param threshold_sd Number of SDs above rolling mean (default: 2)
#' @param ... Additional arguments
#' @return List containing spike times, spike vector, and filtered trace
#' @keywords internal
#' @examples
#' \dontrun{
#' trace <- rnorm(500)
#' result <- rolling_stats_inference(trace, window_size = 30)
#' }
rolling_stats_inference <- function(trace, window_size = 20, threshold_sd = 2, ...) {
  message("Running rolling statistics spike inference")

  n <- length(trace)

  # Normalize trace
  trace_norm <- (trace - mean(trace, na.rm = TRUE)) / sd(trace, na.rm = TRUE)

  # Calculate rolling statistics
  ws <- min(window_size, n %/% 4)
  if (ws < 3) ws <- 3

  rolling_mean <- stats::filter(trace_norm, rep(1/ws, ws), sides = 2)
  rolling_mean[is.na(rolling_mean)] <- mean(trace_norm, na.rm = TRUE)

  residuals <- trace_norm - rolling_mean
  rolling_var <- stats::filter(residuals^2, rep(1/ws, ws), sides = 2)
  rolling_var[is.na(rolling_var)] <- var(trace_norm, na.rm = TRUE)
  rolling_sd <- sqrt(pmax(rolling_var, 1e-10))

  # Filtered trace (residuals scaled by local variability)
  filtered_trace <- residuals / rolling_sd

  # Detect spikes
  threshold <- threshold_sd
  spike_candidates <- which(filtered_trace > threshold)

  # Apply simple refractory period (minimum 2 frames between spikes)
  spike_indices <- numeric(0)
  if (length(spike_candidates) > 0) {
    spike_indices <- spike_candidates[1]
    for (i in 2:length(spike_candidates)) {
      if (spike_candidates[i] - spike_indices[length(spike_indices)] >= 2) {
        spike_indices <- c(spike_indices, spike_candidates[i])
      }
    }
  }

  # Create spike vector
  spikes <- rep(0, n)
  spikes[spike_indices] <- 1

  # Estimate calcium using simple exponential model
  calcium_est <- estimate_calcium(trace, spikes)

  return(list(
    spikes = spikes,
    spike_times = spike_indices,
    calcium_est = calcium_est,
    filtered_trace = as.numeric(filtered_trace),
    method = "rolling_stats",
    parameters = list(window_size = ws, threshold_sd = threshold_sd)
  ))
}

#' Convolutional Spike Inference
#'
#' Infer spikes using 1D convolution filters. Applies edge detection
#' and smoothing kernels to detect transient events.
#'
#' @param trace Calcium trace (numeric vector)
#' @param kernel_size Kernel size for convolution filters (default: 5)
#' @param threshold_quantile Quantile threshold for spike detection (default: 0.9)
#' @param ... Additional arguments
#' @return List containing spike times, spike vector, and feature maps
#' @keywords internal
#' @examples
#' \dontrun{
#' trace <- rnorm(500)
#' result <- convolutional_inference(trace, kernel_size = 7)
#' }
convolutional_inference <- function(trace, kernel_size = 5, threshold_quantile = 0.9, ...) {
  message("Running convolutional spike inference")

  n <- length(trace)

  # Normalize trace
  trace_norm <- (trace - mean(trace, na.rm = TRUE)) / sd(trace, na.rm = TRUE)

  # Define convolution kernels
  # Edge detection kernel (first derivative approximation)
  half_k <- kernel_size %/% 2
  edge_kernel <- c(rep(-1, half_k), 0, rep(1, half_k))
  edge_kernel <- edge_kernel / sum(abs(edge_kernel))

  # Smoothing kernel (moving average)
  smooth_kernel <- rep(1/kernel_size, kernel_size)

  # Peak detection kernel (second derivative approximation)
  peak_kernel <- c(rep(-1, half_k), kernel_size, rep(-1, half_k))
  peak_kernel <- peak_kernel / sum(abs(peak_kernel))

  # Apply convolutions
  edge_features <- convolve_1d(trace_norm, edge_kernel)
  smooth_features <- convolve_1d(trace_norm, smooth_kernel)
  peak_features <- convolve_1d(trace_norm, peak_kernel)

  # Combine features (weighted sum)
  combined_features <- 0.4 * pmax(edge_features, 0) +
                       0.3 * smooth_features +
                       0.3 * pmax(peak_features, 0)

  # Apply ReLU-like activation
  activated_features <- pmax(combined_features, 0)

  # Detect spikes using threshold
  threshold <- quantile(activated_features, threshold_quantile, na.rm = TRUE)
  spike_candidates <- which(activated_features > threshold)

  # Apply refractory period
  spike_indices <- apply_refractory(spike_candidates, min_interval = 3)

  # Create spike vector
  spikes <- rep(0, n)
  spikes[spike_indices] <- 1

  # Estimate calcium
  calcium_est <- estimate_calcium(trace, spikes)

  return(list(
    spikes = spikes,
    spike_times = spike_indices,
    calcium_est = calcium_est,
    features = list(
      edge = edge_features,
      smooth = smooth_features,
      peak = peak_features,
      combined = combined_features,
      activated = activated_features
    ),
    method = "convolution",
    parameters = list(kernel_size = kernel_size, threshold_quantile = threshold_quantile)
  ))
}

#' Correlation-Weighted Spike Inference
#'
#' Infer spikes using correlation-based weighting. Weights each time point
#' by its similarity to neighboring points (similar to attention mechanisms).
#'
#' @param trace Calcium trace (numeric vector)
#' @param window_size Window size for computing local correlations (default: 30)
#' @param threshold_quantile Quantile threshold for spike detection (default: 0.9)
#' @param ... Additional arguments
#' @return List containing spike times, spike vector, and correlation weights
#' @keywords internal
#' @examples
#' \dontrun{
#' trace <- rnorm(500)
#' result <- correlation_weighted_inference(trace, window_size = 50)
#' }
correlation_weighted_inference <- function(trace, window_size = 30, threshold_quantile = 0.9, ...) {
  message("Running correlation-weighted spike inference")

  n <- length(trace)

  # Normalize trace
  trace_norm <- (trace - mean(trace, na.rm = TRUE)) / sd(trace, na.rm = TRUE)

  # Compute local correlation weights
  weights <- compute_local_weights(trace_norm, window_size)

  # Weight the trace by local similarity
  weighted_trace <- trace_norm * weights

  # Also compute multi-scale features
  local_features <- compute_local_features(trace_norm, window_size %/% 2)
  global_features <- compute_global_features(trace_norm)
  temporal_features <- compute_temporal_features(trace_norm, window_size)

  # Combine features
  combined_features <- 0.4 * weighted_trace +
                       0.2 * local_features +
                       0.2 * global_features +
                       0.2 * temporal_features

  # Detect spikes
  feature_scores <- combined_features * weights
  threshold <- quantile(feature_scores, threshold_quantile, na.rm = TRUE)
  spike_candidates <- which(feature_scores > threshold)

  # Apply refractory period
  spike_indices <- apply_refractory(spike_candidates, min_interval = 3)

  # Create spike vector
  spikes <- rep(0, n)
  spikes[spike_indices] <- 1

  # Estimate calcium
  calcium_est <- estimate_calcium(trace, spikes)

  return(list(
    spikes = spikes,
    spike_times = spike_indices,
    calcium_est = calcium_est,
    correlation_weights = weights,
    weighted_trace = weighted_trace,
    combined_features = combined_features,
    method = "correlation_weighted",
    parameters = list(window_size = window_size, threshold_quantile = threshold_quantile)
  ))
}

# Helper functions

#' @keywords internal
convolve_1d <- function(signal, kernel) {
  # 1D convolution with edge handling
  n <- length(signal)
  k <- length(kernel)
  result <- numeric(n)
  half_k <- k %/% 2

  for (i in 1:n) {
    start_idx <- max(1, i - half_k)
    end_idx <- min(n, i + half_k)
    local_signal <- signal[start_idx:end_idx]
    local_kernel <- kernel[1:length(local_signal)]
    if (length(local_kernel) > length(local_signal)) {
      local_kernel <- local_kernel[1:length(local_signal)]
    }
    result[i] <- sum(local_signal * local_kernel, na.rm = TRUE)
  }

  return(result)
}

#' @keywords internal
apply_refractory <- function(spike_indices, min_interval = 3) {
  if (length(spike_indices) <= 1) {
    return(spike_indices)
  }

  filtered <- spike_indices[1]
  last_spike <- spike_indices[1]

  for (i in 2:length(spike_indices)) {
    if (spike_indices[i] - last_spike >= min_interval) {
      filtered <- c(filtered, spike_indices[i])
      last_spike <- spike_indices[i]
    }
  }

  return(filtered)
}

#' @keywords internal
compute_local_weights <- function(trace, window_size) {
  # Compute weights based on local similarity (exponential decay)
  n <- length(trace)
  weights <- numeric(n)

  for (i in 1:n) {
    start_idx <- max(1, i - window_size)
    end_idx <- min(n, i + window_size)

    # Distance-based weights
    distances <- abs(seq(start_idx, end_idx) - i)
    dist_weights <- exp(-distances / (window_size / 2))

    # Similarity-based weights
    local_trace <- trace[start_idx:end_idx]
    similarities <- 1 / (1 + abs(local_trace - trace[i]))

    # Combined weight
    combined <- dist_weights * similarities
    weights[i] <- mean(combined, na.rm = TRUE)
  }

  # Normalize
  weights <- weights / max(weights, na.rm = TRUE)

  return(weights)
}

#' @keywords internal
compute_local_features <- function(trace, window_size) {
  n <- length(trace)
  result <- numeric(n)

  for (i in 1:n) {
    start_idx <- max(1, i - window_size)
    end_idx <- min(n, i + window_size)
    local_trace <- trace[start_idx:end_idx]

    # Weight by distance
    distances <- abs(seq(start_idx, end_idx) - i)
    weights <- exp(-distances / window_size)
    weights <- weights / sum(weights)

    result[i] <- sum(local_trace * weights, na.rm = TRUE)
  }

  return(result)
}

#' @keywords internal
compute_global_features <- function(trace) {
  # Weight by deviation from global statistics
  global_mean <- mean(trace, na.rm = TRUE)
  global_sd <- sd(trace, na.rm = TRUE)

  if (global_sd == 0) global_sd <- 1

  # Z-score weighted features
  z_scores <- (trace - global_mean) / global_sd
  weights <- pmax(z_scores, 0)  # Only positive deviations
  weights <- weights / max(weights + 1e-10, na.rm = TRUE)

  return(weights * trace)
}

#' @keywords internal
compute_temporal_features <- function(trace, window_size) {
  # Emphasize recent history (causal weighting)
  n <- length(trace)
  result <- numeric(n)

  for (i in 1:n) {
    start_idx <- max(1, i - window_size)
    recent_trace <- trace[start_idx:i]

    # Recency weighting
    recency_weights <- seq_along(recent_trace) / length(recent_trace)
    recency_weights <- recency_weights / sum(recency_weights)

    result[i] <- sum(recent_trace * recency_weights, na.rm = TRUE)
  }

  return(result)
}

#' @keywords internal
estimate_calcium <- function(trace, spikes) {
  # Simple exponential calcium model
  n <- length(trace)
  calcium_est <- numeric(n)

  tau_decay <- 10  # Decay time constant in frames

  calcium_est[1] <- 0

  for (i in 2:n) {
    # Decay from previous
    calcium_est[i] <- calcium_est[i-1] * exp(-1/tau_decay)

    # Add spike contribution
    if (spikes[i-1] == 1) {
      calcium_est[i] <- calcium_est[i] + 1
    }
  }

  # Scale to match trace range
  if (max(calcium_est) > 0) {
    calcium_est <- calcium_est / max(calcium_est) * (max(trace) - min(trace)) + min(trace)
  }

  return(calcium_est)
}

# Deprecated aliases for backwards compatibility

#' @rdname statistical_spike_inference
#' @keywords internal
deep_spike_inference <- function(trace, model_type = "rolling", ...) {
  .Deprecated("statistical_spike_inference",
              msg = "deep_spike_inference() is deprecated. Use statistical_spike_inference() instead - this uses statistical methods, not deep learning.")

  # Map old model_type names to new method names
  method <- switch(model_type,
    "lstm" = "rolling",
    "cnn" = "convolution",
    "transformer" = "correlation",
    "rolling"
  )

  statistical_spike_inference(trace, method = method, ...)
}
