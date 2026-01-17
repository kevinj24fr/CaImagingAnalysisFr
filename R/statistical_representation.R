#' Statistical Representation Learning
#'
#' Statistical methods for learning representations and clustering calcium traces.
#' These functions use classical statistical approaches (PCA, k-means) rather than
#' deep learning, providing fast and interpretable results.
#'
#' @name statistical_representation
NULL

#' PCA-Based Representation Learning
#'
#' Learn low-dimensional representations of calcium traces using PCA.
#' This provides an interpretable embedding based on principal components.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param embedding_dim Dimension of learned embeddings (default: 64)
#' @param scale Whether to scale features (default: TRUE)
#' @param center Whether to center features (default: TRUE)
#' @param ... Additional arguments
#' @return List containing embeddings and PCA results
#' @export
#' @examples
#' \dontrun{
#' traces <- matrix(rnorm(100 * 500), nrow = 100)
#' result <- pca_representation(traces, embedding_dim = 32)
#' print(result$explained_variance[1:5])
#' }
pca_representation <- function(traces, embedding_dim = 64, scale = TRUE, center = TRUE, ...) {
  message("Learning PCA-based representation")

  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # PCA for dimensionality reduction

  pca_result <- prcomp(t(traces), scale. = scale, center = center)

  # Take the first embedding_dim components
  n_components <- min(embedding_dim, ncol(pca_result$x))
  embeddings <- pca_result$x[, 1:n_components, drop = FALSE]

  # Calculate explained variance
  explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cumulative_variance <- cumsum(explained_variance)

  return(list(
    embeddings = embeddings,
    pca_result = pca_result,
    explained_variance = explained_variance,
    cumulative_variance = cumulative_variance,
    n_components = n_components,
    method = "pca"
  ))
}

#' Threshold-Based Spike Detection
#'
#' Detect spikes using statistical thresholding (mean + k*SD).
#' Simple and fast method for spike inference.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param threshold_sd Number of standard deviations above mean for threshold (default: 2)
#' @param min_duration Minimum spike duration in frames (default: 1)
#' @param ... Additional arguments
#' @return List containing spike predictions and thresholds used
#' @export
#' @examples
#' \dontrun{
#' traces <- matrix(rnorm(10 * 1000), nrow = 10)
#' result <- threshold_spike_detection(traces, threshold_sd = 2.5)
#' }
threshold_spike_detection <- function(traces, threshold_sd = 2, min_duration = 1, ...) {
  message("Running threshold-based spike detection")

  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  spike_predictions <- matrix(0, n_cells, n_time)
  thresholds_used <- numeric(n_cells)

  for (i in 1:n_cells) {
    trace <- traces[i, ]
    threshold <- mean(trace, na.rm = TRUE) + threshold_sd * sd(trace, na.rm = TRUE)
    thresholds_used[i] <- threshold
    spike_predictions[i, ] <- as.numeric(trace > threshold)
  }

  return(list(
    spike_predictions = spike_predictions,
    thresholds = thresholds_used,
    threshold_sd = threshold_sd,
    method = "threshold"
  ))
}

#' Wavelet Denoising
#'
#' Denoise calcium traces using wavelet thresholding or moving average smoothing.
#' Falls back to simple smoothing if waveslim package is not available.
#'
#' @param traces Matrix of noisy calcium traces (cells x time)
#' @param noise_level Estimated noise level for wavelet threshold (default: 0.1)
#' @param method Denoising method: "wavelet" or "smooth" (default: "wavelet")
#' @param smooth_window Window size for moving average if using "smooth" (default: 5)
#' @param ... Additional arguments
#' @return List containing denoised traces and metrics
#' @export
#' @examples
#' \dontrun{
#' traces <- matrix(rnorm(10 * 1000), nrow = 10)
#' result <- wavelet_denoise(traces, method = "smooth")
#' }
wavelet_denoise <- function(traces, noise_level = 0.1, method = "wavelet", smooth_window = 5, ...) {
  message("Running ", method, " denoising")

  denoised_traces <- traces

  if (method == "wavelet" && requireNamespace("waveslim", quietly = TRUE)) {
    # Wavelet denoising
    for (i in 1:nrow(traces)) {
      trace <- traces[i, ]
      if (length(trace) >= 8) {  # Need minimum length for wavelets
        tryCatch({
          wd <- waveslim::dwt(trace, wf = "haar", n.levels = min(3, floor(log2(length(trace))) - 1))
          threshold <- noise_level * sqrt(2 * log(length(trace)))
          # Soft thresholding
          for (level in paste0("d", 1:min(3, length(wd) - 1))) {
            if (!is.null(wd[[level]])) {
              wd[[level]][abs(wd[[level]]) < threshold] <- 0
            }
          }
          denoised_traces[i, ] <- waveslim::idwt(wd)
        }, error = function(e) {
          # Fall back to smoothing on error
          denoised_traces[i, ] <- stats::filter(trace, rep(1/smooth_window, smooth_window), sides = 2)
          denoised_traces[i, is.na(denoised_traces[i, ])] <- trace[is.na(denoised_traces[i, ])]
        })
      }
    }
  } else {
    # Simple moving average smoothing
    for (i in 1:nrow(traces)) {
      smoothed <- stats::filter(traces[i, ], rep(1/smooth_window, smooth_window), sides = 2)
      denoised_traces[i, ] <- ifelse(is.na(smoothed), traces[i, ], smoothed)
    }
  }

  # Calculate metrics
  snr_before <- mean(abs(traces), na.rm = TRUE) / sd(traces, na.rm = TRUE)
  snr_after <- mean(abs(denoised_traces), na.rm = TRUE) / sd(denoised_traces, na.rm = TRUE)

  return(list(
    denoised_traces = denoised_traces,
    original_traces = traces,
    metrics = list(
      snr_before = snr_before,
      snr_after = snr_after,
      snr_improvement = snr_after - snr_before
    ),
    method = method
  ))
}

#' Adaptive Threshold Spike Detection
#'
#' Detect spikes using locally adaptive thresholding with rolling statistics.
#' More robust than fixed threshold for traces with varying baseline.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param window_size Window size for rolling statistics (default: 50)
#' @param threshold_sd Number of local SDs above local mean (default: 2)
#' @param ... Additional arguments
#' @return List containing spike predictions
#' @export
#' @examples
#' \dontrun{
#' traces <- matrix(rnorm(10 * 1000), nrow = 10)
#' result <- adaptive_threshold_detection(traces, window_size = 100)
#' }
adaptive_threshold_detection <- function(traces, window_size = 50, threshold_sd = 2, ...) {
  message("Running adaptive threshold spike detection")

  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  spike_predictions <- matrix(0, n_cells, n_time)

  for (i in 1:n_cells) {
    trace <- traces[i, ]

    # Adjust window size if needed
    ws <- min(window_size, n_time %/% 4)
    if (ws < 5) ws <- 5

    # Calculate rolling mean and SD
    rolling_mean <- stats::filter(trace, rep(1/ws, ws), sides = 2)
    rolling_mean[is.na(rolling_mean)] <- mean(trace, na.rm = TRUE)

    residuals <- trace - rolling_mean
    rolling_var <- stats::filter(residuals^2, rep(1/ws, ws), sides = 2)
    rolling_var[is.na(rolling_var)] <- var(trace, na.rm = TRUE)
    rolling_sd <- sqrt(pmax(rolling_var, 0))

    # Adaptive threshold
    threshold <- rolling_mean + threshold_sd * rolling_sd
    spike_predictions[i, ] <- as.numeric(trace > threshold)
  }

  return(list(
    spike_predictions = spike_predictions,
    window_size = window_size,
    threshold_sd = threshold_sd,
    method = "adaptive_threshold"
  ))
}

#' PCA + K-means Clustering
#'
#' Cluster calcium traces using PCA for dimensionality reduction followed by k-means.
#' Provides interpretable clustering based on temporal activity patterns.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param n_clusters Number of clusters (default: 5)
#' @param n_components Number of PCA components to use (default: 32)
#' @param nstart Number of random starts for k-means (default: 10)
#' @param ... Additional arguments
#' @return List containing cluster assignments and embeddings
#' @export
#' @examples
#' \dontrun{
#' traces <- matrix(rnorm(100 * 500), nrow = 100)
#' result <- pca_kmeans_clustering(traces, n_clusters = 5)
#' table(result$cluster_assignments)
#' }
pca_kmeans_clustering <- function(traces, n_clusters = 5, n_components = 32, nstart = 10, ...) {
  message("Performing PCA + k-means clustering")

  n_cells <- nrow(traces)

  # PCA for dimensionality reduction
  pca_result <- prcomp(t(traces), scale. = TRUE, center = TRUE)
  n_comp <- min(n_components, ncol(pca_result$x))
  embeddings <- pca_result$x[, 1:n_comp, drop = FALSE]

  # Adjust n_clusters if needed
  if (n_clusters >= n_cells) {
    n_clusters <- max(2, n_cells - 1)
    message("Adjusted n_clusters to ", n_clusters)
  }

  # K-means clustering
  kmeans_result <- kmeans(embeddings, centers = n_clusters, nstart = nstart)
  cluster_assignments <- kmeans_result$cluster

  # Calculate silhouette-like score (simplified)
  within_ss <- kmeans_result$tot.withinss
  between_ss <- kmeans_result$betweenss
  cluster_quality <- between_ss / (within_ss + between_ss)

  return(list(
    embeddings = embeddings,
    cluster_assignments = cluster_assignments,
    n_clusters = n_clusters,
    cluster_centers = kmeans_result$centers,
    cluster_sizes = kmeans_result$size,
    cluster_quality = cluster_quality,
    kmeans_result = kmeans_result,
    pca_result = pca_result,
    method = "pca_kmeans"
  ))
}

#' Temporal Correlation Heatmap
#'
#' Visualize temporal correlation structure of calcium traces.
#' Useful for understanding temporal dependencies in the data.
#'
#' @param traces Matrix of calcium traces (cells x time) or single trace
#' @param cell_indices Indices of cells to visualize (default: 1:5)
#' @param max_lag Maximum lag for correlation (default: 100)
#' @param ... Additional arguments
#' @return ggplot object with correlation heatmaps or list of correlation matrices
#' @export
#' @examples
#' \dontrun{
#' traces <- matrix(rnorm(10 * 500), nrow = 10)
#' result <- temporal_correlation_heatmap(traces, cell_indices = 1:3)
#' }
temporal_correlation_heatmap <- function(traces, cell_indices = 1:5, max_lag = 100, ...) {
  message("Computing temporal correlations")

  if (is.null(dim(traces))) {
    traces <- matrix(traces, nrow = 1)
  }

  n_cells <- min(length(cell_indices), nrow(traces))
  cell_indices <- cell_indices[1:n_cells]
  n_time <- ncol(traces)
  max_lag <- min(max_lag, n_time - 1)

  # Compute autocorrelation for each cell
  correlations <- list()
  for (idx in seq_along(cell_indices)) {
    i <- cell_indices[idx]
    trace <- traces[i, ]
    acf_result <- acf(trace, lag.max = max_lag, plot = FALSE)
    correlations[[idx]] <- list(
      cell = i,
      lags = acf_result$lag,
      acf = acf_result$acf
    )
  }

  # Try to create ggplot visualization if available
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_data <- do.call(rbind, lapply(seq_along(correlations), function(idx) {
      data.frame(
        cell = correlations[[idx]]$cell,
        lag = as.numeric(correlations[[idx]]$lags),
        correlation = as.numeric(correlations[[idx]]$acf)
      )
    }))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = lag, y = correlation, color = factor(cell))) +
      ggplot2::geom_line() +
      ggplot2::labs(title = "Temporal Autocorrelation",
                    x = "Lag (frames)", y = "Autocorrelation",
                    color = "Cell") +
      ggplot2::theme_minimal() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

    return(p)
  }

  return(correlations)
}

#' Benchmark Spike Detection Methods
#'
#' Compare different spike detection methods on the same data.
#' If ground truth is provided, computes accuracy metrics.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param ground_truth Optional ground truth spike matrix (same dimensions as traces)
#' @param methods Methods to compare (default: c("threshold", "adaptive", "mad"))
#' @param ... Additional arguments passed to detection methods
#' @return Data frame with performance comparison
#' @export
#' @examples
#' \dontrun{
#' traces <- matrix(rnorm(10 * 1000), nrow = 10)
#' comparison <- benchmark_spike_methods(traces)
#' print(comparison)
#' }
benchmark_spike_methods <- function(traces, ground_truth = NULL, methods = c("threshold", "adaptive", "mad"), ...) {
  message("Benchmarking spike detection methods")

  results <- list()

  for (method in methods) {
    start_time <- Sys.time()

    if (method == "threshold") {
      pred <- threshold_spike_detection(traces, ...)$spike_predictions
    } else if (method == "adaptive") {
      pred <- adaptive_threshold_detection(traces, ...)$spike_predictions
    } else if (method == "mad") {
      # MAD-based detection
      pred <- matrix(0, nrow(traces), ncol(traces))
      for (i in 1:nrow(traces)) {
        trace <- traces[i, ]
        med <- median(trace, na.rm = TRUE)
        mad_val <- median(abs(trace - med), na.rm = TRUE) * 1.4826
        threshold <- med + 3 * mad_val
        pred[i, ] <- as.numeric(trace > threshold)
      }
    } else {
      next
    }

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    # Calculate metrics
    n_spikes <- sum(pred)
    spike_rate <- n_spikes / (nrow(traces) * ncol(traces))

    result_row <- data.frame(
      method = method,
      n_spikes_detected = n_spikes,
      spike_rate = spike_rate,
      runtime_sec = elapsed,
      stringsAsFactors = FALSE
    )

    # If ground truth provided, calculate accuracy metrics
    if (!is.null(ground_truth) && all(dim(ground_truth) == dim(pred))) {
      tp <- sum(pred == 1 & ground_truth == 1)
      fp <- sum(pred == 1 & ground_truth == 0)
      fn <- sum(pred == 0 & ground_truth == 1)
      tn <- sum(pred == 0 & ground_truth == 0)

      precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
      recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
      f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else 0
      accuracy <- (tp + tn) / length(pred)

      result_row$accuracy <- accuracy
      result_row$precision <- precision
      result_row$recall <- recall
      result_row$f1 <- f1
    }

    results[[method]] <- result_row
  }

  return(do.call(rbind, results))
}

# Deprecated function aliases for backwards compatibility
# These will warn users about the name change

#' @rdname pca_representation
#' @export
self_supervised_learning <- function(...) {
  .Deprecated("pca_representation", msg = "self_supervised_learning() is deprecated. Use pca_representation() instead - this function uses PCA, not deep learning.")
  pca_representation(...)
}

#' @rdname threshold_spike_detection
#' @export
transformer_spike_inference <- function(...) {
  .Deprecated("threshold_spike_detection", msg = "transformer_spike_inference() is deprecated. Use threshold_spike_detection() instead - this function uses thresholding, not transformers.")
  threshold_spike_detection(...)
}

#' @rdname wavelet_denoise
#' @export
transformer_denoising <- function(...) {
  .Deprecated("wavelet_denoise", msg = "transformer_denoising() is deprecated. Use wavelet_denoise() instead - this function uses wavelets/smoothing, not transformers.")
  wavelet_denoise(...)
}

#' @rdname adaptive_threshold_detection
#' @export
transfer_learning_spike_inference <- function(...) {
  .Deprecated("adaptive_threshold_detection", msg = "transfer_learning_spike_inference() is deprecated. Use adaptive_threshold_detection() instead - this function uses adaptive thresholding, not transfer learning.")
  adaptive_threshold_detection(...)
}

#' @rdname pca_kmeans_clustering
#' @export
contrastive_clustering <- function(...) {
  .Deprecated("pca_kmeans_clustering", msg = "contrastive_clustering() is deprecated. Use pca_kmeans_clustering() instead - this function uses PCA + k-means, not contrastive learning.")
  pca_kmeans_clustering(...)
}

#' @rdname temporal_correlation_heatmap
#' @export
visualize_attention <- function(...) {
  .Deprecated("temporal_correlation_heatmap", msg = "visualize_attention() is deprecated. Use temporal_correlation_heatmap() instead - this function visualizes correlations, not attention weights.")
  temporal_correlation_heatmap(...)
}

#' @rdname benchmark_spike_methods
#' @export
compare_deep_models <- function(...) {
  .Deprecated("benchmark_spike_methods", msg = "compare_deep_models() is deprecated. Use benchmark_spike_methods() instead - this function compares statistical methods, not deep learning models.")
  benchmark_spike_methods(...)
}
