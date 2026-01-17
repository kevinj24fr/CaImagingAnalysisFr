#' Stimulus and Behavior Alignment Functions
#'
#' Functions for aligning calcium imaging data with stimulus presentations
#' and behavioral events.
#'
#' @name stimulus_alignment
NULL

#' Align Traces to Events
#'
#' Extract peri-event time histograms (PETHs) from calcium traces aligned to
#' stimulus or behavioral events.
#'
#' @param traces Data frame or matrix of calcium traces (time x cells)
#' @param events Vector of event times (in frames) or data frame with event info
#' @param pre_frames Number of frames before event
#' @param post_frames Number of frames after event
#' @param baseline_frames Frames to use for baseline (relative to event, negative)
#' @param normalize Method for normalization ("none", "zscore", "baseline", "percent")
#' @param frame_rate Frame rate in Hz (for time conversion)
#' @param verbose Show progress messages
#' @return List with aligned traces and statistics
#'
#' @details
#' This function extracts trial-aligned snippets of calcium activity centered on
#' events of interest (e.g., stimulus onset, reward delivery, movement initiation).
#'
#' @examples
#' \dontrun
#' # Align to stimulus onsets
#' stim_times <- c(100, 300, 500, 700)  # Frame indices
#' aligned <- align_to_events(traces, stim_times, pre_frames = 30, post_frames = 60)
#'
#' # Access trial-averaged response
#' mean_response <- aligned$mean_traces
#'
#' # Plot
#' plot_event_aligned(aligned)
#' }
#'
#' @export
align_to_events <- function(traces,
                            events,
                            pre_frames = 30,
                            post_frames = 60,
                            baseline_frames = c(-30, -10),
                            normalize = c("zscore", "baseline", "percent", "none"),
                            frame_rate = 30,
                            verbose = TRUE) {

  normalize <- match.arg(normalize)

  # Convert traces to matrix
  if (is.data.frame(traces)) {
    trace_mat <- as.matrix(traces[, sapply(traces, is.numeric)])
    cell_names <- colnames(trace_mat)
  } else {
    trace_mat <- as.matrix(traces)
    cell_names <- colnames(trace_mat)
    if (is.null(cell_names)) {
      cell_names <- paste0("Cell_", seq_len(ncol(trace_mat)))
    }
  }

  n_time <- nrow(trace_mat)
  n_cells <- ncol(trace_mat)

  # Process events
  if (is.data.frame(events)) {
    if ("frame" %in% names(events)) {
      event_frames <- events$frame
    } else if ("time" %in% names(events)) {
      event_frames <- round(events$time * frame_rate)
    } else {
      event_frames <- events[[1]]
    }
    event_info <- events
  } else {
    event_frames <- as.integer(events)
    event_info <- data.frame(frame = event_frames, trial = seq_along(event_frames))
  }

  # Filter valid events (enough pre/post frames)
  valid_events <- event_frames > pre_frames & event_frames <= (n_time - post_frames)
  event_frames <- event_frames[valid_events]
  n_events <- length(event_frames)

  if (n_events == 0) {
    stop("No valid events found within trace bounds")
  }

  if (verbose) {
    message("Aligning ", n_cells, " cells to ", n_events, " events")
    message("Window: -", pre_frames, " to +", post_frames, " frames")
  }

  # Time vector for aligned data
  time_vec <- seq(-pre_frames, post_frames)
  n_window <- length(time_vec)
  time_seconds <- time_vec / frame_rate

  # Initialize output arrays
  # aligned_traces: [time x cells x trials]
  aligned_traces <- array(NA, dim = c(n_window, n_cells, n_events))

  # Extract aligned traces
  for (trial in seq_len(n_events)) {
    event_frame <- event_frames[trial]
    frames <- (event_frame - pre_frames):(event_frame + post_frames)
    aligned_traces[, , trial] <- trace_mat[frames, ]
  }

  # Normalize if requested
  if (normalize != "none") {
    if (verbose) message("Applying ", normalize, " normalization")

    baseline_idx <- which(time_vec >= baseline_frames[1] & time_vec <= baseline_frames[2])

    for (cell in seq_len(n_cells)) {
      for (trial in seq_len(n_events)) {
        trace <- aligned_traces[, cell, trial]
        baseline <- trace[baseline_idx]
        baseline_mean <- mean(baseline, na.rm = TRUE)
        baseline_sd <- sd(baseline, na.rm = TRUE)

        aligned_traces[, cell, trial] <- switch(normalize,
          "zscore" = (trace - baseline_mean) / (baseline_sd + 1e-10),
          "baseline" = trace - baseline_mean,
          "percent" = 100 * (trace - baseline_mean) / (baseline_mean + 1e-10),
          trace
        )
      }
    }
  }

  # Compute statistics
  mean_traces <- apply(aligned_traces, c(1, 2), mean, na.rm = TRUE)
  sem_traces <- apply(aligned_traces, c(1, 2), function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  std_traces <- apply(aligned_traces, c(1, 2), sd, na.rm = TRUE)

  colnames(mean_traces) <- cell_names
  colnames(sem_traces) <- cell_names
  colnames(std_traces) <- cell_names

  result <- list(
    aligned_traces = aligned_traces,
    mean_traces = mean_traces,
    sem_traces = sem_traces,
    std_traces = std_traces,
    time_frames = time_vec,
    time_seconds = time_seconds,
    event_frames = event_frames,
    event_info = event_info[valid_events, , drop = FALSE],
    n_events = n_events,
    n_cells = n_cells,
    cell_names = cell_names,
    parameters = list(
      pre_frames = pre_frames,
      post_frames = post_frames,
      baseline_frames = baseline_frames,
      normalize = normalize,
      frame_rate = frame_rate
    )
  )

  class(result) <- c("event_aligned", "list")
  return(result)
}

#' Compute Event-Triggered Average
#'
#' Compute the average response aligned to events for each cell.
#'
#' @param traces Data frame or matrix of calcium traces
#' @param events Event times (frames)
#' @param window_frames Window size around event (single value or c(pre, post))
#' @param frame_rate Frame rate in Hz
#' @return Matrix of average responses (window x cells)
#'
#' @export
event_triggered_average <- function(traces, events, window_frames = 30, frame_rate = 30) {
  if (length(window_frames) == 1) {
    pre <- post <- window_frames
  } else {
    pre <- window_frames[1]
    post <- window_frames[2]
  }

  result <- align_to_events(traces, events, pre_frames = pre, post_frames = post,
                            normalize = "none", frame_rate = frame_rate, verbose = FALSE)
  return(result$mean_traces)
}

#' Detect Stimulus/Event Times from TTL
#'
#' Detect event times from a TTL (digital) signal.
#'
#' @param ttl_signal Digital signal vector (0s and 1s, or analog)
#' @param threshold Threshold for detecting onset (for analog signals)
#' @param edge Detection edge ("rising", "falling", "both")
#' @param min_interval Minimum interval between events (frames)
#' @param frame_rate Frame rate for time conversion
#' @return Data frame with event times
#'
#' @export
detect_events_ttl <- function(ttl_signal,
                              threshold = 0.5,
                              edge = c("rising", "falling", "both"),
                              min_interval = 1,
                              frame_rate = 30) {

  edge <- match.arg(edge)

  # Binarize signal
  binary <- as.numeric(ttl_signal > threshold)

  # Compute differences
  diff_signal <- diff(c(0, binary))

  # Find edges
  if (edge == "rising") {
    event_frames <- which(diff_signal == 1)
  } else if (edge == "falling") {
    event_frames <- which(diff_signal == -1)
  } else {
    event_frames <- which(diff_signal != 0)
  }

  if (length(event_frames) == 0) {
    return(data.frame(frame = integer(0), time = numeric(0), type = character(0)))
  }

  # Filter by minimum interval
  if (min_interval > 1 && length(event_frames) > 1) {
    keep <- c(TRUE, diff(event_frames) >= min_interval)
    event_frames <- event_frames[keep]
  }

  # Determine event type
  if (edge == "both") {
    event_types <- ifelse(diff_signal[event_frames] == 1, "onset", "offset")
  } else if (edge == "rising") {
    event_types <- rep("onset", length(event_frames))
  } else {
    event_types <- rep("offset", length(event_frames))
  }

  data.frame(
    frame = event_frames,
    time = event_frames / frame_rate,
    type = event_types
  )
}

#' Compute Response Metrics
#'
#' Compute response metrics for event-aligned data.
#'
#' @param aligned Result from align_to_events()
#' @param response_window Frames to consider as response window (relative to event)
#' @param baseline_window Frames to consider as baseline window (relative to event)
#' @param metrics Which metrics to compute
#' @return Data frame with metrics per cell
#'
#' @export
compute_response_metrics <- function(aligned,
                                     response_window = c(0, 30),
                                     baseline_window = c(-30, -5),
                                     metrics = c("peak", "mean", "auc", "latency", "reliability")) {

  time_vec <- aligned$time_frames
  n_cells <- aligned$n_cells
  n_trials <- aligned$n_events

  # Find window indices
  response_idx <- which(time_vec >= response_window[1] & time_vec <= response_window[2])
  baseline_idx <- which(time_vec >= baseline_window[1] & time_vec <= baseline_window[2])

  # Initialize results
  results <- data.frame(
    cell = aligned$cell_names
  )

  # Compute metrics for each cell
  for (cell in seq_len(n_cells)) {
    # Get trial responses
    trial_responses <- aligned$aligned_traces[, cell, ]  # time x trials

    # Baseline-subtracted responses
    baseline_means <- colMeans(trial_responses[baseline_idx, , drop = FALSE], na.rm = TRUE)
    subtracted <- sweep(trial_responses, 2, baseline_means, "-")

    # Mean across trials
    mean_response <- rowMeans(subtracted, na.rm = TRUE)
    response_period <- mean_response[response_idx]

    if ("peak" %in% metrics) {
      results$peak_response[cell] <- max(response_period, na.rm = TRUE)
      results$peak_time[cell] <- time_vec[response_idx[which.max(response_period)]]
    }

    if ("mean" %in% metrics) {
      results$mean_response[cell] <- mean(response_period, na.rm = TRUE)
    }

    if ("auc" %in% metrics) {
      # Area under curve (positive only)
      results$auc[cell] <- sum(pmax(response_period, 0), na.rm = TRUE)
    }

    if ("latency" %in% metrics) {
      # Time to half-max response
      half_max <- max(response_period, na.rm = TRUE) / 2
      above_half <- which(response_period >= half_max)
      if (length(above_half) > 0) {
        results$latency[cell] <- time_vec[response_idx[above_half[1]]]
      } else {
        results$latency[cell] <- NA
      }
    }

    if ("reliability" %in% metrics) {
      # Correlation between odd and even trials
      odd_trials <- seq(1, n_trials, by = 2)
      even_trials <- seq(2, n_trials, by = 2)

      if (length(odd_trials) > 1 && length(even_trials) > 1) {
        odd_mean <- rowMeans(subtracted[response_idx, odd_trials, drop = FALSE], na.rm = TRUE)
        even_mean <- rowMeans(subtracted[response_idx, even_trials, drop = FALSE], na.rm = TRUE)
        results$reliability[cell] <- cor(odd_mean, even_mean, use = "complete.obs")
      } else {
        results$reliability[cell] <- NA
      }
    }
  }

  return(results)
}

#' Statistical Test for Event Responses
#'
#' Test whether cells show significant responses to events.
#'
#' @param aligned Result from align_to_events()
#' @param response_window Response window (frames relative to event)
#' @param baseline_window Baseline window (frames relative to event)
#' @param test Statistical test ("ttest", "wilcox", "permutation")
#' @param n_permutations Number of permutations for permutation test
#' @param alpha Significance level
#' @return Data frame with test results
#'
#' @export
test_event_responses <- function(aligned,
                                 response_window = c(5, 30),
                                 baseline_window = c(-30, -5),
                                 test = c("ttest", "wilcox", "permutation"),
                                 n_permutations = 1000,
                                 alpha = 0.05) {

  test <- match.arg(test)
  time_vec <- aligned$time_frames
  n_cells <- aligned$n_cells
  n_trials <- aligned$n_events

  response_idx <- which(time_vec >= response_window[1] & time_vec <= response_window[2])
  baseline_idx <- which(time_vec >= baseline_window[1] & time_vec <= baseline_window[2])

  results <- data.frame(
    cell = aligned$cell_names,
    baseline_mean = NA,
    response_mean = NA,
    statistic = NA,
    p_value = NA,
    significant = NA
  )

  for (cell in seq_len(n_cells)) {
    trial_responses <- aligned$aligned_traces[, cell, ]

    # Compute baseline and response values per trial
    baseline_vals <- colMeans(trial_responses[baseline_idx, , drop = FALSE], na.rm = TRUE)
    response_vals <- colMeans(trial_responses[response_idx, , drop = FALSE], na.rm = TRUE)

    results$baseline_mean[cell] <- mean(baseline_vals, na.rm = TRUE)
    results$response_mean[cell] <- mean(response_vals, na.rm = TRUE)

    # Statistical test
    if (test == "ttest") {
      test_result <- t.test(response_vals, baseline_vals, paired = TRUE)
      results$statistic[cell] <- test_result$statistic
      results$p_value[cell] <- test_result$p.value

    } else if (test == "wilcox") {
      test_result <- wilcox.test(response_vals, baseline_vals, paired = TRUE)
      results$statistic[cell] <- test_result$statistic
      results$p_value[cell] <- test_result$p.value

    } else if (test == "permutation") {
      # Permutation test
      observed_diff <- mean(response_vals - baseline_vals, na.rm = TRUE)
      combined <- c(baseline_vals, response_vals)
      n_obs <- length(baseline_vals)

      perm_diffs <- numeric(n_permutations)
      for (p in seq_len(n_permutations)) {
        shuffled <- sample(combined)
        perm_diffs[p] <- mean(shuffled[(n_obs + 1):(2 * n_obs)] - shuffled[1:n_obs])
      }

      results$statistic[cell] <- observed_diff
      results$p_value[cell] <- mean(abs(perm_diffs) >= abs(observed_diff))
    }

    results$significant[cell] <- results$p_value[cell] < alpha
  }

  # Multiple comparison correction
  results$p_adjusted <- p.adjust(results$p_value, method = "fdr")
  results$significant_corrected <- results$p_adjusted < alpha

  return(results)
}

#' Plot Event-Aligned Data
#'
#' Create visualizations of event-aligned calcium activity.
#'
#' @param aligned Result from align_to_events()
#' @param cells Which cells to plot (indices or names)
#' @param type Plot type ("heatmap", "traces", "mean_sem", "raster")
#' @param show_baseline Show baseline period
#' @param ... Additional arguments passed to plotting functions
#' @return ggplot object
#'
#' @export
plot_event_aligned <- function(aligned,
                               cells = NULL,
                               type = c("heatmap", "traces", "mean_sem", "raster"),
                               show_baseline = TRUE,
                               ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }

  type <- match.arg(type)

  # Select cells
  if (is.null(cells)) {
    cells <- 1:min(10, aligned$n_cells)
  } else if (is.character(cells)) {
    cells <- match(cells, aligned$cell_names)
  }

  time_vec <- aligned$time_seconds

  if (type == "mean_sem") {
    # Mean +/- SEM for selected cells
    df_list <- lapply(cells, function(c) {
      data.frame(
        time = time_vec,
        mean = aligned$mean_traces[, c],
        sem = aligned$sem_traces[, c],
        cell = aligned$cell_names[c]
      )
    })
    df <- do.call(rbind, df_list)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = mean, color = cell, fill = cell)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = mean - sem, ymax = mean + sem),
                           alpha = 0.3, color = NA) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::labs(
        title = "Event-Aligned Responses",
        x = "Time from event (s)",
        y = "Response",
        color = "Cell",
        fill = "Cell"
      ) +
      ggplot2::theme_minimal()

  } else if (type == "heatmap") {
    # Heatmap of mean responses
    df <- expand.grid(
      time = time_vec,
      cell = factor(aligned$cell_names[cells], levels = aligned$cell_names[cells])
    )
    df$value <- as.vector(aligned$mean_traces[, cells])

    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = cell, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::geom_vline(xintercept = 0, color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(
        title = "Event-Aligned Heatmap",
        x = "Time from event (s)",
        y = "Cell",
        fill = "Response"
      ) +
      ggplot2::theme_minimal()

  } else if (type == "raster") {
    # Raster plot showing all trials
    cell_idx <- cells[1]  # Just first cell for raster

    df <- expand.grid(
      time = time_vec,
      trial = 1:aligned$n_events
    )
    df$value <- as.vector(aligned$aligned_traces[, cell_idx, ])

    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = trial, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::geom_vline(xintercept = 0, color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(
        title = paste("Trial Raster:", aligned$cell_names[cell_idx]),
        x = "Time from event (s)",
        y = "Trial",
        fill = "Response"
      ) +
      ggplot2::theme_minimal()

  } else if (type == "traces") {
    # Individual trial traces
    cell_idx <- cells[1]

    df <- data.frame(
      time = rep(time_vec, aligned$n_events),
      trial = rep(1:aligned$n_events, each = length(time_vec)),
      value = as.vector(aligned$aligned_traces[, cell_idx, ])
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = value, group = trial)) +
      ggplot2::geom_line(alpha = 0.3, color = "gray50") +
      ggplot2::geom_line(
        data = data.frame(time = time_vec, value = aligned$mean_traces[, cell_idx]),
        ggplot2::aes(x = time, y = value, group = 1),
        color = "red", linewidth = 1.2
      ) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title = paste("Trial Traces:", aligned$cell_names[cell_idx]),
        x = "Time from event (s)",
        y = "Response"
      ) +
      ggplot2::theme_minimal()
  }

  return(p)
}

#' Synchronize Behavior Data
#'
#' Synchronize behavioral data with imaging data using timestamps.
#'
#' @param behavior_data Data frame with behavior events and timestamps
#' @param imaging_timestamps Vector of imaging frame timestamps
#' @param behavior_time_column Name of time column in behavior_data
#' @param method Synchronization method ("nearest", "interpolate")
#' @return Data frame with behavior events matched to imaging frames
#'
#' @export
synchronize_behavior <- function(behavior_data,
                                 imaging_timestamps,
                                 behavior_time_column = "time",
                                 method = c("nearest", "interpolate")) {

  method <- match.arg(method)

  behavior_times <- behavior_data[[behavior_time_column]]
  n_frames <- length(imaging_timestamps)

  # Find matching frames for each behavior event
  if (method == "nearest") {
    matched_frames <- sapply(behavior_times, function(t) {
      which.min(abs(imaging_timestamps - t))
    })
    time_offsets <- behavior_times - imaging_timestamps[matched_frames]

  } else if (method == "interpolate") {
    # Linear interpolation to get fractional frame numbers
    matched_frames <- approx(imaging_timestamps, 1:n_frames, xout = behavior_times)$y
    time_offsets <- rep(0, length(behavior_times))
  }

  # Add frame information to behavior data
  result <- behavior_data
  result$imaging_frame <- matched_frames
  result$time_offset <- time_offsets

  # Flag events outside imaging range
  result$valid <- matched_frames >= 1 & matched_frames <= n_frames

  return(result)
}

#' Create PSTH (Peri-Stimulus Time Histogram)
#'
#' Create traditional PSTH from spike/event data.
#'
#' @param spike_times List of spike time vectors (one per trial or cell)
#' @param events Event times to align to
#' @param bin_size Bin size in seconds
#' @param pre_time Time before event (seconds)
#' @param post_time Time after event (seconds)
#' @return Data frame with PSTH
#'
#' @export
create_psth <- function(spike_times,
                        events,
                        bin_size = 0.05,
                        pre_time = 0.5,
                        post_time = 1.0) {

  # Create time bins
  bins <- seq(-pre_time, post_time, by = bin_size)
  n_bins <- length(bins) - 1
  bin_centers <- bins[-length(bins)] + bin_size / 2

  # Count spikes in each bin for each event
  n_events <- length(events)

  if (is.list(spike_times)) {
    # Multiple cells/trials
    n_units <- length(spike_times)
    counts <- matrix(0, nrow = n_bins, ncol = n_units)

    for (u in seq_len(n_units)) {
      spikes <- spike_times[[u]]
      for (e in events) {
        aligned_spikes <- spikes - e
        hist_result <- hist(aligned_spikes, breaks = bins, plot = FALSE)
        counts[, u] <- counts[, u] + hist_result$counts
      }
    }
    counts <- counts / n_events  # Average per event

  } else {
    # Single spike train
    counts <- numeric(n_bins)
    for (e in events) {
      aligned_spikes <- spike_times - e
      hist_result <- hist(aligned_spikes[aligned_spikes >= -pre_time & aligned_spikes <= post_time],
                          breaks = bins, plot = FALSE)
      counts <- counts + hist_result$counts
    }
    counts <- counts / n_events
  }

  # Convert to rate
  rate <- counts / bin_size

  data.frame(
    time = bin_centers,
    count = if (is.matrix(counts)) rowMeans(counts) else counts,
    rate = if (is.matrix(rate)) rowMeans(rate) else rate
  )
}
