#' Sequence Detection
#'
#' Methods for detecting temporal sequences of neural activity.
#' Identifies ordered patterns of cell activation across time.
#'
#' @name sequence_detection
#' @keywords internal
NULL

#' Detect activity sequences
#'
#' Find recurring temporal sequences of cell activations.
#'
#' @param spikes Binary spike matrix (cells x time)
#' @param method Detection method: "template", "bayesian", "correlation"
#' @param min_cells Minimum cells in a sequence
#' @param max_duration Maximum sequence duration in frames
#' @param n_shuffles Number of shuffles for significance testing
#'
#' @return List with detected sequences
#' @export
#'
#' @examples
#' \dontrun{
#' sequences <- detect_sequences(spikes, method = "template")
#' plot_sequences(sequences, spikes)
#' }
detect_sequences <- function(spikes, method = "template", min_cells = 5,
                              max_duration = 50, n_shuffles = 100) {
  n_cells <- nrow(spikes)
  n_time <- ncol(spikes)

  result <- switch(method,
    template = .detect_sequences_template(spikes, min_cells, max_duration),
    bayesian = .detect_sequences_bayesian(spikes, min_cells, max_duration),
    correlation = .detect_sequences_correlation(spikes, min_cells, max_duration),
    stop("Unknown method: ", method)
  )

  # Test significance with shuffles
  if (n_shuffles > 0 && length(result$sequences) > 0) {
    null_scores <- replicate(n_shuffles, {
      shuffled <- .shuffle_spikes_circular(spikes)
      null_result <- switch(method,
        template = .detect_sequences_template(shuffled, min_cells, max_duration),
        bayesian = .detect_sequences_bayesian(shuffled, min_cells, max_duration),
        correlation = .detect_sequences_correlation(shuffled, min_cells, max_duration)
      )
      if (length(null_result$sequences) > 0) {
        max(sapply(null_result$sequences, `[[`, "score"))
      } else {
        0
      }
    })

    # Mark significant sequences
    threshold <- quantile(null_scores, 0.95)
    for (i in seq_along(result$sequences)) {
      result$sequences[[i]]$significant <- result$sequences[[i]]$score > threshold
      result$sequences[[i]]$p_value <- mean(null_scores >= result$sequences[[i]]$score)
    }
  }

  structure(
    c(result, list(
      method = method,
      min_cells = min_cells,
      max_duration = max_duration,
      n_cells = n_cells
    )),
    class = "neural_sequences"
  )
}

.detect_sequences_template <- function(spikes, min_cells, max_duration) {
  n_cells <- nrow(spikes)
  n_time <- ncol(spikes)

  # Find candidate sequence starts (high activity periods)
  activity <- colSums(spikes)
  threshold <- mean(activity) + 2 * sd(activity)
  candidates <- which(activity > threshold)

  sequences <- list()
  used_times <- logical(n_time)

  for (start in candidates) {
    if (used_times[start]) next

    end <- min(start + max_duration - 1, n_time)
    window <- spikes[, start:end, drop = FALSE]

    # Find cells with spikes in window
    active_cells <- which(rowSums(window) > 0)
    if (length(active_cells) < min_cells) next

    # Get first spike time for each active cell
    first_spikes <- sapply(active_cells, function(c) {
      which(window[c, ] > 0)[1]
    })

    # Check if there's a clear temporal order
    order_idx <- order(first_spikes)
    ordered_cells <- active_cells[order_idx]
    ordered_times <- first_spikes[order_idx]

    # Score: correlation between rank and spike time
    score <- cor(seq_along(ordered_times), ordered_times)

    if (!is.na(score) && score > 0.5) {
      sequences[[length(sequences) + 1]] <- list(
        cells = ordered_cells,
        times = ordered_times + start - 1,
        relative_times = ordered_times,
        start_frame = start,
        duration = max(ordered_times),
        score = score
      )

      # Mark times as used
      used_times[start:(start + max(ordered_times))] <- TRUE
    }
  }

  list(sequences = sequences)
}

.detect_sequences_bayesian <- function(spikes, min_cells, max_duration) {
  # Simplified Bayesian sequence detection
  n_cells <- nrow(spikes)
  n_time <- ncol(spikes)

  # Use rank-order correlation across sliding windows
  window_size <- max_duration
  step <- window_size %/% 2

  sequences <- list()

  for (start in seq(1, n_time - window_size, by = step)) {
    window <- spikes[, start:(start + window_size - 1), drop = FALSE]

    # Get activation order
    active_cells <- which(rowSums(window) > 0)
    if (length(active_cells) < min_cells) next

    first_spikes <- sapply(active_cells, function(c) {
      spk <- which(window[c, ] > 0)
      if (length(spk) > 0) spk[1] else NA
    })

    valid <- !is.na(first_spikes)
    if (sum(valid) < min_cells) next

    active_cells <- active_cells[valid]
    first_spikes <- first_spikes[valid]

    # Compute sequence score (how ordered is the activation?)
    order_idx <- order(first_spikes)
    ordered_times <- first_spikes[order_idx]

    # Score based on temporal regularity
    if (length(ordered_times) > 2) {
      intervals <- diff(ordered_times)
      cv <- sd(intervals) / mean(intervals)
      score <- 1 / (1 + cv)  # Higher for regular intervals
    } else {
      score <- 0.5
    }

    if (score > 0.3) {
      sequences[[length(sequences) + 1]] <- list(
        cells = active_cells[order_idx],
        times = first_spikes[order_idx] + start - 1,
        relative_times = first_spikes[order_idx],
        start_frame = start,
        duration = max(first_spikes) - min(first_spikes),
        score = score
      )
    }
  }

  list(sequences = sequences)
}

.detect_sequences_correlation <- function(spikes, min_cells, max_duration) {
  n_cells <- nrow(spikes)
  n_time <- ncol(spikes)

  # Compute time-lagged cross-correlations
  max_lag <- max_duration %/% 2
  lag_corr <- array(0, dim = c(n_cells, n_cells, 2 * max_lag + 1))

  for (i in seq_len(n_cells)) {
    for (j in seq_len(n_cells)) {
      if (i == j) next
      cc <- ccf(spikes[i, ], spikes[j, ], lag.max = max_lag, plot = FALSE)
      lag_corr[i, j, ] <- cc$acf
    }
  }

  # Find cell pairs with strong lagged correlations
  sequences <- list()

  # Build sequences from strongly connected pairs
  best_lags <- apply(lag_corr, c(1, 2), which.max) - max_lag - 1
  best_corrs <- apply(lag_corr, c(1, 2), max)

  # Threshold correlations
  significant <- best_corrs > quantile(best_corrs[best_corrs > 0], 0.9)

  # Build chains
  for (seed_cell in seq_len(n_cells)) {
    chain <- seed_cell
    current <- seed_cell
    total_lag <- 0

    while (TRUE) {
      # Find best follower
      followers <- which(significant[current, ] & best_lags[current, ] > 0)
      followers <- setdiff(followers, chain)

      if (length(followers) == 0) break

      best_follower <- followers[which.max(best_corrs[current, followers])]
      chain <- c(chain, best_follower)
      total_lag <- total_lag + best_lags[current, best_follower]
      current <- best_follower

      if (length(chain) > 20) break  # Prevent infinite loops
    }

    if (length(chain) >= min_cells) {
      sequences[[length(sequences) + 1]] <- list(
        cells = chain,
        relative_times = cumsum(c(0, best_lags[cbind(chain[-length(chain)], chain[-1])])),
        duration = total_lag,
        score = mean(best_corrs[cbind(chain[-length(chain)], chain[-1])])
      )
    }
  }

  # Remove duplicate/overlapping sequences
  sequences <- .remove_duplicate_sequences(sequences)

  list(sequences = sequences)
}

.shuffle_spikes_circular <- function(spikes) {
  # Circular shuffle each cell's spike train independently
  n_cells <- nrow(spikes)
  n_time <- ncol(spikes)

  shuffled <- spikes
  for (i in seq_len(n_cells)) {
    shift <- sample(n_time, 1)
    shuffled[i, ] <- c(spikes[i, (shift + 1):n_time], spikes[i, 1:shift])
  }

  shuffled
}

.remove_duplicate_sequences <- function(sequences) {
  if (length(sequences) < 2) return(sequences)

  # Sort by score
  scores <- sapply(sequences, `[[`, "score")
  order_idx <- order(scores, decreasing = TRUE)
  sequences <- sequences[order_idx]

  # Remove sequences with >50% overlap with higher-scoring sequences
  keep <- rep(TRUE, length(sequences))

  for (i in 2:length(sequences)) {
    cells_i <- sequences[[i]]$cells

    for (j in 1:(i-1)) {
      if (!keep[j]) next

      cells_j <- sequences[[j]]$cells
      overlap <- length(intersect(cells_i, cells_j)) / length(cells_i)

      if (overlap > 0.5) {
        keep[i] <- FALSE
        break
      }
    }
  }

  sequences[keep]
}

#' Detect replay events
#'
#' Find compressed replay of activity sequences.
#'
#' @param spikes Binary spike matrix (cells x time)
#' @param template Sequence template (ordered cell indices or sequence object)
#' @param compression_range Range of compression factors to test
#' @param threshold Correlation threshold for replay detection
#'
#' @return List of replay events
#' @export
detect_replay <- function(spikes, template, compression_range = c(1, 20),
                           threshold = 0.5) {
  n_cells <- nrow(spikes)
  n_time <- ncol(spikes)

  # Extract template order
  if (is.list(template)) {
    template_cells <- template$cells
    template_times <- template$relative_times
  } else {
    template_cells <- template
    template_times <- seq_along(template)
  }

  template_order <- rank(template_times)
  n_template <- length(template_cells)

  replays <- list()

  # Scan through time with various window sizes
  for (compression in seq(compression_range[1], compression_range[2])) {
    window_size <- ceiling(max(template_times) / compression)
    if (window_size < 3) next

    for (start in seq_len(n_time - window_size)) {
      window <- spikes[template_cells, start:(start + window_size - 1), drop = FALSE]

      # Get first spike times in window
      first_spikes <- sapply(seq_len(n_template), function(i) {
        spk <- which(window[i, ] > 0)
        if (length(spk) > 0) spk[1] else NA
      })

      # Skip if not enough cells active
      if (sum(!is.na(first_spikes)) < n_template * 0.5) next

      # Compute rank correlation with template
      valid <- !is.na(first_spikes)
      observed_order <- rank(first_spikes[valid])

      corr <- cor(template_order[valid], observed_order, method = "spearman")

      if (!is.na(corr) && abs(corr) > threshold) {
        replays[[length(replays) + 1]] <- list(
          start_frame = start,
          end_frame = start + window_size - 1,
          compression = compression,
          correlation = corr,
          direction = ifelse(corr > 0, "forward", "reverse"),
          active_cells = template_cells[valid],
          spike_times = first_spikes[valid]
        )
      }
    }
  }

  # Remove overlapping replays (keep highest correlation)
  replays <- .remove_overlapping_replays(replays)

  structure(
    list(
      replays = replays,
      template = list(cells = template_cells, times = template_times),
      n_replays = length(replays)
    ),
    class = "replay_events"
  )
}

.remove_overlapping_replays <- function(replays) {
  if (length(replays) < 2) return(replays)

  # Sort by correlation strength
  corrs <- sapply(replays, function(r) abs(r$correlation))
  order_idx <- order(corrs, decreasing = TRUE)
  replays <- replays[order_idx]

  keep <- rep(TRUE, length(replays))

  for (i in 2:length(replays)) {
    for (j in 1:(i-1)) {
      if (!keep[j]) next

      # Check time overlap
      overlap <- min(replays[[i]]$end_frame, replays[[j]]$end_frame) -
        max(replays[[i]]$start_frame, replays[[j]]$start_frame)

      if (overlap > 0) {
        keep[i] <- FALSE
        break
      }
    }
  }

  replays[keep]
}

#' Compute sequence reactivation strength
#'
#' Measure how strongly a sequence is reactivated over time.
#'
#' @param spikes Spike matrix
#' @param sequence Sequence object or cell indices
#' @param window_size Window size for computing reactivation
#'
#' @return Vector of reactivation strengths over time
#' @export
compute_reactivation_strength <- function(spikes, sequence, window_size = 10) {
  if (is.list(sequence)) {
    cells <- sequence$cells
    template_order <- rank(sequence$relative_times)
  } else {
    cells <- sequence
    template_order <- seq_along(cells)
  }

  n_time <- ncol(spikes)
  strengths <- rep(NA, n_time)

  for (t in seq(window_size, n_time)) {
    window <- spikes[cells, (t - window_size + 1):t, drop = FALSE]

    # Get activation order in window
    first_spikes <- sapply(seq_along(cells), function(i) {
      spk <- which(window[i, ] > 0)
      if (length(spk) > 0) spk[1] else NA
    })

    valid <- !is.na(first_spikes)
    if (sum(valid) < 3) {
      strengths[t] <- 0
      next
    }

    observed_order <- rank(first_spikes[valid])
    strengths[t] <- cor(template_order[valid], observed_order, method = "spearman")
  }

  strengths
}

#' Plot detected sequences
#'
#' @param sequences Neural sequences object
#' @param spikes Optional spike matrix for raster
#' @param n_sequences Number of sequences to plot
#'
#' @return ggplot object
#' @export
plot_sequences <- function(sequences, spikes = NULL, n_sequences = 5) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  n_plot <- min(n_sequences, length(sequences$sequences))

  if (n_plot == 0) {
    message("No sequences to plot")
    return(invisible(NULL))
  }

  # Sort by score
  scores <- sapply(sequences$sequences, `[[`, "score")
  top_idx <- order(scores, decreasing = TRUE)[1:n_plot]

  plots <- lapply(top_idx, function(i) {
    seq_i <- sequences$sequences[[i]]

    df <- data.frame(
      cell = seq_i$cells,
      time = if (!is.null(seq_i$times)) seq_i$times else seq_i$relative_times,
      order = seq_along(seq_i$cells)
    )

    ggplot2::ggplot(df, ggplot2::aes(x = time, y = reorder(factor(cell), order))) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_line(ggplot2::aes(group = 1), alpha = 0.5) +
      ggplot2::labs(
        title = sprintf("Sequence %d (score: %.2f)", i, seq_i$score),
        x = "Time", y = "Cell"
      ) +
      ggplot2::theme_minimal()
  })

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(plots, ncol = 1)
  } else {
    plots[[1]]
  }
}
