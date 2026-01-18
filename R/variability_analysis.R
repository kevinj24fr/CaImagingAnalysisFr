#' Neural Variability Analysis
#'
#' Analyze trial-to-trial variability and noise correlations.
#' Quantify reliability and precision of neural responses.
#'
#' @name variability_analysis
#' @keywords internal
NULL

#' Compute Fano factor
#'
#' Ratio of variance to mean (measure of over/under-dispersion).
#'
#' @param traces Matrix of traces (cells x time) or counts
#' @param window_size Window for computing local Fano factor
#' @param step Step size for sliding window
#'
#' @return List with Fano factors per cell
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' ff <- compute_fano_factor(traces)
#' }
compute_fano_factor <- function(traces, window_size = NULL, step = NULL) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  if (is.null(window_size)) {
    # Global Fano factor
    means <- rowMeans(traces)
    vars <- apply(traces, 1, var)
    fano <- vars / (means + 1e-10)

    return(list(
      fano_factor = fano,
      mean = means,
      variance = vars,
      type = "global"
    ))
  }

  # Sliding window Fano factor
  if (is.null(step)) step <- window_size %/% 2

  n_windows <- floor((n_time - window_size) / step) + 1
  fano_over_time <- matrix(NA, n_cells, n_windows)
  time_points <- numeric(n_windows)

  for (w in seq_len(n_windows)) {
    start <- (w - 1) * step + 1
    end <- start + window_size - 1
    window_data <- traces[, start:end, drop = FALSE]

    means <- rowMeans(window_data)
    vars <- apply(window_data, 1, var)
    fano_over_time[, w] <- vars / (means + 1e-10)
    time_points[w] <- (start + end) / 2
  }

  list(
    fano_factor = rowMeans(fano_over_time, na.rm = TRUE),
    fano_over_time = fano_over_time,
    time_points = time_points,
    type = "sliding"
  )
}

#' Compute noise correlations
#'
#' Trial-to-trial correlation in residual activity.
#'
#' @param traces_list List of trace matrices (one per trial)
#' @param subtract_mean Subtract trial-averaged response
#'
#' @return Noise correlation matrix
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # List of traces from repeated trials
#' noise_corr <- compute_noise_correlations(trial_traces)
#' }
compute_noise_correlations <- function(traces_list, subtract_mean = TRUE) {
  n_trials <- length(traces_list)
  n_cells <- nrow(traces_list[[1]])
  n_time <- ncol(traces_list[[1]])

  # Stack trials
  all_traces <- array(0, dim = c(n_cells, n_time, n_trials))
  for (t in seq_len(n_trials)) {
    all_traces[, , t] <- traces_list[[t]]
  }

  if (subtract_mean) {
    # Compute trial-averaged response
    mean_response <- apply(all_traces, c(1, 2), mean)

    # Compute residuals
    residuals <- array(0, dim = c(n_cells, n_time, n_trials))
    for (t in seq_len(n_trials)) {
      residuals[, , t] <- all_traces[, , t] - mean_response
    }
  } else {
    residuals <- all_traces
  }

  # Flatten residuals (cells x (time * trials))
  flat_residuals <- matrix(0, n_cells, n_time * n_trials)
  for (t in seq_len(n_trials)) {
    cols <- ((t - 1) * n_time + 1):(t * n_time)
    flat_residuals[, cols] <- residuals[, , t]
  }

  # Compute correlation matrix
  noise_corr <- cor(t(flat_residuals))

  list(
    noise_correlations = noise_corr,
    mean_noise_correlation = mean(noise_corr[upper.tri(noise_corr)]),
    n_trials = n_trials
  )
}

#' Compute signal correlations
#'
#' Correlation in trial-averaged responses.
#'
#' @param traces_list List of trace matrices (one per trial)
#'
#' @return Signal correlation matrix
#' @keywords internal
compute_signal_correlations <- function(traces_list) {
  n_trials <- length(traces_list)
  n_cells <- nrow(traces_list[[1]])
  n_time <- ncol(traces_list[[1]])

  # Compute trial-averaged response
  mean_response <- Reduce(`+`, traces_list) / n_trials

  # Signal correlation
  signal_corr <- cor(t(mean_response))

  list(
    signal_correlations = signal_corr,
    mean_signal_correlation = mean(signal_corr[upper.tri(signal_corr)])
  )
}

#' Compute coefficient of variation
#'
#' CV = SD / mean for each cell.
#'
#' @param traces Matrix of traces (cells x time)
#' @param by_trial If list of trials, compute CV across trials
#'
#' @return CV values per cell
#' @keywords internal
compute_cv <- function(traces, by_trial = FALSE) {
  if (is.list(traces) && by_trial) {
    # CV across trials at each time point
    n_trials <- length(traces)
    n_cells <- nrow(traces[[1]])
    n_time <- ncol(traces[[1]])

    # Stack
    stacked <- array(0, dim = c(n_cells, n_time, n_trials))
    for (t in seq_len(n_trials)) {
      stacked[, , t] <- traces[[t]]
    }

    # CV at each time point
    cv_matrix <- apply(stacked, c(1, 2), function(x) sd(x) / (mean(x) + 1e-10))

    return(list(
      cv_matrix = cv_matrix,
      mean_cv = rowMeans(cv_matrix),
      type = "across_trials"
    ))
  }

  # CV across time
  means <- rowMeans(traces)
  sds <- apply(traces, 1, sd)
  cv <- sds / (means + 1e-10)

  list(
    cv = cv,
    mean = means,
    sd = sds,
    type = "across_time"
  )
}

#' Compute response reliability
#'
#' Correlation between repeated presentations of same stimulus.
#'
#' @param traces_list List of trace matrices for repeated trials
#' @param method Reliability method: "correlation", "split_half"
#'
#' @return Reliability measure per cell
#' @keywords internal
compute_reliability <- function(traces_list, method = "correlation") {
  n_trials <- length(traces_list)
  n_cells <- nrow(traces_list[[1]])

  if (n_trials < 2) {
    warning("Need at least 2 trials for reliability")
    return(rep(NA, n_cells))
  }

  if (method == "correlation") {
    # Mean pairwise correlation across trials
    reliability <- sapply(seq_len(n_cells), function(cell) {
      correlations <- c()
      for (i in 1:(n_trials - 1)) {
        for (j in (i + 1):n_trials) {
          correlations <- c(correlations,
                            cor(traces_list[[i]][cell, ],
                                traces_list[[j]][cell, ]))
        }
      }
      mean(correlations, na.rm = TRUE)
    })

  } else if (method == "split_half") {
    # Split-half reliability with Spearman-Brown correction
    half1 <- seq(1, n_trials, by = 2)
    half2 <- seq(2, n_trials, by = 2)

    mean1 <- Reduce(`+`, traces_list[half1]) / length(half1)
    mean2 <- Reduce(`+`, traces_list[half2]) / length(half2)

    r <- sapply(seq_len(n_cells), function(cell) {
      cor(mean1[cell, ], mean2[cell, ])
    })

    # Spearman-Brown correction
    reliability <- 2 * r / (1 + r)
  }

  list(
    reliability = reliability,
    method = method,
    n_trials = n_trials
  )
}

#' Compute response precision
#'
#' Timing precision of responses.
#'
#' @param traces Matrix of traces (cells x time)
#' @param events Event times or detected transients
#' @param frame_rate Frame rate
#'
#' @return Precision measures
#' @keywords internal
compute_response_precision <- function(traces, events = NULL, frame_rate = 10) {
  n_cells <- nrow(traces)

  if (is.null(events)) {
    # Detect events
    events <- detect_transients(traces, frame_rate = frame_rate)
  }

  # Precision as inverse of jitter in peak times
  precision <- sapply(seq_len(n_cells), function(cell) {
    cell_events <- events$events[events$events$cell_id == cell, ]
    if (nrow(cell_events) < 2) return(NA)

    # Jitter in inter-event intervals
    intervals <- diff(cell_events$peak_time)
    if (length(intervals) < 2) return(NA)

    cv_iei <- sd(intervals) / mean(intervals)
    1 / (cv_iei + 0.1)  # Precision = 1/CV
  })

  list(
    precision = precision,
    mean_precision = mean(precision, na.rm = TRUE)
  )
}

#' Compute population variability
#'
#' Shared vs private variance decomposition.
#'
#' @param traces Matrix of traces (cells x time)
#' @param n_components Number of components for factor analysis
#'
#' @return List with shared and private variance
#' @keywords internal
decompose_variability <- function(traces, n_components = 3) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Factor analysis
  fa_result <- tryCatch({
    factanal(t(traces), factors = n_components, scores = "regression")
  }, error = function(e) {
    # Fall back to PCA-based decomposition
    NULL
  })

  if (!is.null(fa_result)) {
    # Shared variance from factor loadings
    loadings <- fa_result$loadings[, 1:n_components]
    shared_var <- rowSums(loadings^2)

    # Private variance (uniqueness)
    private_var <- fa_result$uniquenesses

  } else {
    # PCA-based decomposition
    pca <- prcomp(t(traces), center = TRUE, scale. = TRUE)

    # Shared variance from top PCs
    var_explained <- pca$sdev^2 / sum(pca$sdev^2)
    n_shared <- min(n_components, sum(var_explained > 0.05))

    if (n_shared > 0) {
      shared_var <- rowSums(pca$rotation[, 1:n_shared, drop = FALSE]^2) *
        sum(var_explained[1:n_shared])
    } else {
      shared_var <- rep(0, n_cells)
    }

    private_var <- 1 - shared_var
  }

  total_var <- apply(traces, 1, var)

  list(
    shared_variance = shared_var * total_var,
    private_variance = private_var * total_var,
    shared_fraction = shared_var,
    private_fraction = private_var,
    total_variance = total_var
  )
}

#' Compare variability between conditions
#'
#' @param traces1 Traces for condition 1
#' @param traces2 Traces for condition 2
#' @param metric Variability metric: "fano", "cv", "variance"
#'
#' @return Statistical comparison results
#' @keywords internal
compare_variability <- function(traces1, traces2, metric = "fano") {
  n_cells <- nrow(traces1)

  # Compute metric for each condition
  if (metric == "fano") {
    var1 <- compute_fano_factor(traces1)$fano_factor
    var2 <- compute_fano_factor(traces2)$fano_factor
  } else if (metric == "cv") {
    var1 <- compute_cv(traces1)$cv
    var2 <- compute_cv(traces2)$cv
  } else {
    var1 <- apply(traces1, 1, var)
    var2 <- apply(traces2, 1, var)
  }

  # Paired test (same cells)
  test_result <- wilcox.test(var1, var2, paired = TRUE)

  # Effect size
  diff <- var2 - var1
  effect_size <- mean(diff) / sd(diff)

  list(
    metric = metric,
    condition1 = var1,
    condition2 = var2,
    difference = diff,
    mean_difference = mean(diff),
    p_value = test_result$p.value,
    effect_size = effect_size,
    test = "wilcox_paired"
  )
}

#' Plot variability metrics
#'
#' @param variability Variability analysis result
#' @param type Plot type: "histogram", "scatter", "comparison"
#'
#' @return ggplot object
#' @keywords internal
plot_variability <- function(variability, type = "histogram") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  if (type == "histogram") {
    if ("fano_factor" %in% names(variability)) {
      values <- variability$fano_factor
      xlabel <- "Fano Factor"
    } else if ("cv" %in% names(variability)) {
      values <- variability$cv
      xlabel <- "Coefficient of Variation"
    } else if ("reliability" %in% names(variability)) {
      values <- variability$reliability
      xlabel <- "Reliability"
    } else {
      stop("Unknown variability type")
    }

    df <- data.frame(value = values)

    ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white") +
      ggplot2::geom_vline(xintercept = mean(values, na.rm = TRUE),
                          linetype = "dashed", color = "red") +
      ggplot2::labs(x = xlabel, y = "Count",
                    title = sprintf("Distribution of %s", xlabel)) +
      ggplot2::theme_minimal()

  } else if (type == "comparison" && "condition1" %in% names(variability)) {
    df <- data.frame(
      condition1 = variability$condition1,
      condition2 = variability$condition2
    )

    ggplot2::ggplot(df, ggplot2::aes(x = condition1, y = condition2)) +
      ggplot2::geom_point(alpha = 0.5) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      ggplot2::labs(
        x = "Condition 1",
        y = "Condition 2",
        title = sprintf("%s Comparison (p = %.3g)",
                        variability$metric, variability$p_value)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::coord_fixed()
  }
}
