#' Changepoint Detection for Neural Time Series
#'
#' Detect abrupt changes in neural activity patterns,
#' useful for identifying drug onset, behavioral transitions,
#' and network state changes.
#'
#' @name changepoint
NULL

# ============================================================================
# Changepoint Detection Methods
# ============================================================================

#' Detect Changepoints in Neural Activity
#'
#' Find abrupt changes in mean, variance, or distribution of neural signals.
#'
#' @param x Numeric vector or matrix (cells x time for population analysis)
#' @param method Detection method ("pelt", "binseg", "segneigh", "cusum")
#' @param type Type of change ("mean", "var", "meanvar")
#' @param penalty Penalty type ("BIC", "AIC", "MBIC", "Manual")
#' @param pen_value Manual penalty value (if penalty = "Manual")
#' @param min_segment Minimum segment length
#' @param max_changepoints Maximum number of changepoints (for binseg)
#'
#' @return List with:
#'   - changepoints: Detected changepoint locations
#'   - n_changepoints: Number of changepoints
#'   - segments: Segment statistics
#'   - cost: Total cost
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Detect when drug takes effect
#' cp <- detect_changepoints(mean_activity, method = "pelt", type = "mean")
#' abline(v = cp$changepoints, col = "red")
#' }
detect_changepoints <- function(x,
                                 method = c("pelt", "binseg", "cusum", "segneigh"),
                                 type = c("mean", "var", "meanvar"),
                                 penalty = c("BIC", "AIC", "MBIC", "Manual"),
                                 pen_value = NULL,
                                 min_segment = 2,
                                 max_changepoints = NULL) {

  method <- match.arg(method)
  type <- match.arg(type)
  penalty <- match.arg(penalty)

  # Handle matrix input (use first PC or mean)
  if (is.matrix(x)) {
    if (nrow(x) > 1) {
      pca <- prcomp(t(x), center = TRUE, scale. = FALSE)
      x <- pca$x[, 1]
    } else {
      x <- as.numeric(x)
    }
  }

  x <- as.numeric(x)
  n <- length(x)

  # Compute penalty
  if (is.null(pen_value)) {
    pen_value <- switch(penalty,
      "BIC" = log(n),
      "AIC" = 2,
      "MBIC" = 3 * log(n),
      0
    )
  }

  # Call appropriate method
  result <- switch(method,
    "pelt" = pelt_changepoint(x, type, pen_value, min_segment),
    "binseg" = binseg_changepoint(x, type, pen_value, max_changepoints),
    "cusum" = cusum_changepoint(x, type, pen_value),
    "segneigh" = segneigh_changepoint(x, type, pen_value, min_segment, max_changepoints)
  )

  # Compute segment statistics
  cps <- c(0, result$changepoints, n)
  segments <- lapply(seq_len(length(cps) - 1), function(i) {
    seg_data <- x[(cps[i] + 1):cps[i + 1]]
    list(
      start = cps[i] + 1,
      end = cps[i + 1],
      mean = mean(seg_data),
      var = var(seg_data),
      n = length(seg_data)
    )
  })

  structure(
    list(
      changepoints = result$changepoints,
      n_changepoints = length(result$changepoints),
      segments = segments,
      cost = result$cost,
      method = method,
      type = type,
      n = n
    ),
    class = "changepoint_result"
  )
}

#' PELT Algorithm (Pruned Exact Linear Time)
#' @keywords internal
pelt_changepoint <- function(x, type, pen, min_seg) {
  n <- length(x)

  # Cost function
  cost_fn <- get_cost_function(type)

  # Dynamic programming with pruning
  F <- rep(Inf, n + 1)
  F[1] <- -pen  # No cost at start
  cp <- rep(0L, n + 1)
  R <- list(0L)  # Candidate changepoints

  for (t in seq_len(n)) {
    candidates <- R[[t]]
    costs <- sapply(candidates, function(s) {
      if (t - s < min_seg) return(Inf)
      F[s + 1] + cost_fn(x[(s + 1):t]) + pen
    })

    best_idx <- which.min(costs)
    F[t + 1] <- costs[best_idx]
    cp[t + 1] <- candidates[best_idx]

    # Pruning step
    R[[t + 1]] <- c(candidates[costs <= F[t + 1]], t)
  }

  # Backtrack to find changepoints
  changepoints <- integer(0)
  t <- n
  while (cp[t + 1] > 0) {
    changepoints <- c(cp[t + 1], changepoints)
    t <- cp[t + 1]
  }

  list(changepoints = changepoints, cost = F[n + 1])
}

#' Binary Segmentation Algorithm
#' @keywords internal
binseg_changepoint <- function(x, type, pen, max_cp) {
  n <- length(x)
  if (is.null(max_cp)) max_cp <- n / 2

  cost_fn <- get_cost_function(type)

  # Initial cost
  total_cost <- cost_fn(x)
  changepoints <- integer(0)

  for (k in seq_len(max_cp)) {
    # Find best split in each segment
    cps <- c(0, sort(changepoints), n)
    best_gain <- 0
    best_cp <- 0

    for (i in seq_len(length(cps) - 1)) {
      segment <- x[(cps[i] + 1):cps[i + 1]]
      seg_len <- length(segment)

      if (seg_len < 4) next

      for (t in 2:(seg_len - 1)) {
        left <- segment[1:t]
        right <- segment[(t + 1):seg_len]

        current_cost <- cost_fn(segment)
        split_cost <- cost_fn(left) + cost_fn(right)
        gain <- current_cost - split_cost - pen

        if (gain > best_gain) {
          best_gain <- gain
          best_cp <- cps[i] + t
        }
      }
    }

    if (best_gain <= 0) break

    changepoints <- c(changepoints, best_cp)
    total_cost <- total_cost - best_gain
  }

  list(changepoints = sort(changepoints), cost = total_cost)
}

#' CUSUM Changepoint Detection
#' @keywords internal
cusum_changepoint <- function(x, type, pen) {
  n <- length(x)

  # Cumulative sum statistics
  if (type == "mean") {
    x_centered <- x - mean(x)
    S <- cumsum(x_centered)
    S_max <- max(abs(S))
    cp_idx <- which.max(abs(S))
  } else if (type == "var") {
    x_sq <- (x - mean(x))^2
    S <- cumsum(x_sq - mean(x_sq))
    S_max <- max(abs(S))
    cp_idx <- which.max(abs(S))
  } else {
    # Combined
    x_centered <- x - mean(x)
    S <- cumsum(x_centered)
    S_max <- max(abs(S))
    cp_idx <- which.max(abs(S))
  }

  # Statistical test
  threshold <- sqrt(2 * log(n)) + pen / sqrt(n)

  if (S_max / sqrt(n * var(x)) > threshold) {
    changepoints <- cp_idx
  } else {
    changepoints <- integer(0)
  }

  list(changepoints = changepoints, cost = 0)
}

#' Segment Neighborhood Method
#' @keywords internal
segneigh_changepoint <- function(x, type, pen, min_seg, max_cp) {
  n <- length(x)
  if (is.null(max_cp)) max_cp <- 10

  cost_fn <- get_cost_function(type)

  # Precompute costs for all segments
  cost_matrix <- matrix(Inf, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      if (j - i + 1 >= min_seg) {
        cost_matrix[i, j] <- cost_fn(x[i:j])
      }
    }
  }

  # DP for k changepoints
  best_cp <- NULL
  best_cost <- Inf

  for (k in 0:max_cp) {
    result <- dp_k_changepoints(cost_matrix, n, k, pen)
    if (result$cost < best_cost) {
      best_cost <- result$cost
      best_cp <- result$changepoints
    }
  }

  list(changepoints = best_cp, cost = best_cost)
}

#' Dynamic Programming for k Changepoints
#' @keywords internal
dp_k_changepoints <- function(cost_matrix, n, k, pen) {
  if (k == 0) {
    return(list(changepoints = integer(0), cost = cost_matrix[1, n]))
  }

  # F[t, m] = min cost for first t points with m changepoints
  F <- matrix(Inf, n + 1, k + 2)
  F[1, 1] <- 0
  cp <- array(0L, dim = c(n + 1, k + 2))

  for (t in 1:n) {
    for (m in 1:(k + 1)) {
      for (s in max(1, t - n + 1):(t - 1)) {
        if (m > 1 && F[s, m] < Inf) {
          cost <- F[s, m] + cost_matrix[s + 1, t] + pen
          if (cost < F[t + 1, m + 1]) {
            F[t + 1, m + 1] <- cost
            cp[t + 1, m + 1] <- s
          }
        }
      }
      # No previous changepoint
      if (m == 1) {
        F[t + 1, 2] <- min(F[t + 1, 2], cost_matrix[1, t] + pen)
      }
    }
  }

  # Find best number of segments
  final_cost <- F[n + 1, k + 2]

  # Backtrack
  changepoints <- integer(0)
  m <- k + 2
  t <- n + 1
  while (cp[t, m] > 0) {
    changepoints <- c(cp[t, m], changepoints)
    t <- cp[t, m]
    m <- m - 1
  }

  list(changepoints = changepoints, cost = final_cost)
}

#' Get Cost Function for Changepoint Type
#' @keywords internal
get_cost_function <- function(type) {
  if (type == "mean") {
    function(x) {
      n <- length(x)
      if (n < 2) return(0)
      n * var(x) * (n - 1) / n
    }
  } else if (type == "var") {
    function(x) {
      n <- length(x)
      if (n < 2) return(0)
      s2 <- var(x)
      if (s2 <= 0) return(0)
      n * (log(s2) + 1)
    }
  } else {  # meanvar
    function(x) {
      n <- length(x)
      if (n < 2) return(0)
      s2 <- var(x)
      if (s2 <= 0) return(0)
      n * log(s2)
    }
  }
}

#' @export
print.changepoint_result <- function(x, ...) {
  cat("Changepoint Detection Result\n")
  cat("============================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Type: %s\n", x$type))
  cat(sprintf("Observations: %d\n", x$n))
  cat(sprintf("Changepoints: %d\n", x$n_changepoints))
  if (x$n_changepoints > 0) {
    cat(sprintf("Locations: %s\n", paste(x$changepoints, collapse = ", ")))
  }
  invisible(x)
}

#' Plot Changepoint Results
#'
#' @param x Changepoint result
#' @param data Original data
#' @param ... Additional arguments
#'
#' @export
plot_changepoints <- function(x, data, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  if (is.matrix(data)) {
    if (nrow(data) > 1) {
      pca <- prcomp(t(data), center = TRUE, scale. = FALSE)
      data <- pca$x[, 1]
    } else {
      data <- as.numeric(data)
    }
  }

  df <- data.frame(time = seq_along(data), value = data)

  # Add segment means
  df$segment_mean <- NA
  for (seg in x$segments) {
    df$segment_mean[seg$start:seg$end] <- seg$mean
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = segment_mean), color = "red", linewidth = 1)

  if (x$n_changepoints > 0) {
    p <- p + ggplot2::geom_vline(xintercept = x$changepoints,
                                  linetype = "dashed", color = "blue")
  }

  p + ggplot2::labs(title = "Changepoint Detection",
                    x = "Time", y = "Value") +
    ggplot2::theme_minimal()
}

# ============================================================================
# Online Changepoint Detection
# ============================================================================

#' Online Changepoint Detection (BOCPD)
#'
#' Bayesian Online Changepoint Detection for streaming data.
#'
#' @param x Observations (processed sequentially)
#' @param hazard Hazard function value (constant hazard)
#' @param prior_mean Prior mean for Gaussian
#' @param prior_var Prior variance
#'
#' @return List with run length posteriors and detected changepoints
#' @export
bocpd <- function(x, hazard = 1/100, prior_mean = NULL, prior_var = NULL) {
  n <- length(x)

  if (is.null(prior_mean)) prior_mean <- mean(x)
  if (is.null(prior_var)) prior_var <- var(x)

  # Initialize
  R <- matrix(0, n + 1, n + 1)  # Run length distribution
  R[1, 1] <- 1  # Start with run length 0

  # Sufficient statistics for Gaussian
  means <- rep(prior_mean, n + 1)
  vars <- rep(prior_var, n + 1)
  ns <- rep(0, n + 1)

  detected_cps <- integer(0)
  max_run <- rep(0, n)

  for (t in seq_len(n)) {
    # Predictive probability for each run length
    pred_probs <- numeric(t)
    for (r in 0:(t-1)) {
      if (R[r + 1, t] > 1e-10) {
        # Gaussian predictive
        pred_probs[r + 1] <- dnorm(x[t], means[r + 1], sqrt(vars[r + 1]))
      }
    }

    # Growth probabilities
    R[2:(t+1), t + 1] <- R[1:t, t] * pred_probs * (1 - hazard)

    # Changepoint probability
    R[1, t + 1] <- sum(R[1:t, t] * pred_probs * hazard)

    # Normalize
    R[, t + 1] <- R[, t + 1] / sum(R[, t + 1])

    # Update sufficient statistics
    for (r in 0:t) {
      if (r == 0) {
        means[1] <- prior_mean
        vars[1] <- prior_var
        ns[1] <- 0
      } else if (R[r + 1, t + 1] > 1e-10) {
        # Update posterior
        old_n <- ns[r]
        ns[r + 1] <- old_n + 1
        means[r + 1] <- (old_n * means[r] + x[t]) / ns[r + 1]
        # Variance update (simplified)
        vars[r + 1] <- prior_var / (1 + ns[r + 1])
      }
    }

    # Most likely run length
    max_run[t] <- which.max(R[, t + 1]) - 1

    # Detect changepoint if run length drops
    if (t > 1 && max_run[t] < max_run[t-1] && max_run[t] < 5) {
      detected_cps <- c(detected_cps, t)
    }
  }

  structure(
    list(
      run_length_dist = R,
      max_run_length = max_run,
      changepoints = detected_cps,
      n = n
    ),
    class = "bocpd_result"
  )
}

#' Multiple Changepoint Detection Across Cells
#'
#' Detect changepoints in population activity.
#'
#' @param traces Matrix (cells x time)
#' @param method Detection method
#' @param aggregate How to aggregate across cells ("mean", "pc1", "vote")
#'
#' @return Changepoint result
#' @export
detect_population_changepoints <- function(traces,
                                            method = "pelt",
                                            aggregate = c("mean", "pc1", "vote")) {
  aggregate <- match.arg(aggregate)

  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  if (aggregate == "mean") {
    summary <- colMeans(traces)
    detect_changepoints(summary, method = method)

  } else if (aggregate == "pc1") {
    pca <- prcomp(t(traces), center = TRUE, scale. = FALSE)
    summary <- pca$x[, 1]
    detect_changepoints(summary, method = method)

  } else {
    # Vote: detect changepoints per cell, then cluster
    all_cps <- list()
    for (i in seq_len(n_cells)) {
      cp <- detect_changepoints(traces[i, ], method = method)
      all_cps[[i]] <- cp$changepoints
    }

    # Cluster nearby changepoints
    all_cps_vec <- unlist(all_cps)
    if (length(all_cps_vec) == 0) {
      return(structure(list(changepoints = integer(0), n_changepoints = 0,
                            method = "vote", n = n_time),
                       class = "changepoint_result"))
    }

    # Density-based clustering of changepoints
    counts <- table(cut(all_cps_vec, breaks = seq(1, n_time, by = 10)))
    threshold <- n_cells * 0.3  # 30% of cells agree

    consensus_cps <- as.numeric(names(which(counts >= threshold)))

    # Refine to exact locations
    refined_cps <- sapply(consensus_cps, function(cp) {
      nearby <- all_cps_vec[abs(all_cps_vec - cp) < 10]
      if (length(nearby) > 0) round(median(nearby)) else cp
    })

    structure(
      list(
        changepoints = unique(refined_cps),
        n_changepoints = length(unique(refined_cps)),
        all_changepoints = all_cps,
        method = "vote",
        n = n_time
      ),
      class = "changepoint_result"
    )
  }
}
