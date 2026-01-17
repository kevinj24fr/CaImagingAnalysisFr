#' Pharmacological Treatment Analysis
#'
#' Tools for analyzing calcium imaging data from pharmacological experiments.
#' Designed for in vitro/ex vivo preparations comparing treatment conditions.
#'
#' @name pharmacological_analysis
#' @keywords internal
NULL

#' Detect calcium transients
#'
#' Identify individual calcium transients/events and extract their properties.
#'
#' @param traces Matrix of traces (cells x time) or single trace vector
#' @param frame_rate Acquisition frame rate in Hz
#' @param threshold_sd Threshold in standard deviations above baseline
#' @param min_duration Minimum event duration in seconds
#' @param baseline_method Method for baseline estimation: "rolling", "percentile"
#' @param baseline_window Window size for rolling baseline (seconds)
#'
#' @return List with transient times, amplitudes, durations, and kinetics
#' @export
#'
#' @examples
#' \dontrun{
#' events <- detect_transients(traces, frame_rate = 10, threshold_sd = 2.5)
#' summary(events)
#' }
detect_transients <- function(traces, frame_rate = 10, threshold_sd = 2.5,
                               min_duration = 0.3, baseline_method = "rolling",
                               baseline_window = 30) {
  if (is.vector(traces)) {
    traces <- matrix(traces, nrow = 1)
  }

  n_cells <- nrow(traces)
  n_frames <- ncol(traces)
  time_vec <- seq(0, (n_frames - 1) / frame_rate, length.out = n_frames)

  min_frames <- ceiling(min_duration * frame_rate)
  window_frames <- ceiling(baseline_window * frame_rate)

  all_events <- list()

  for (cell in seq_len(n_cells)) {
    trace <- traces[cell, ]

    # Estimate baseline
    if (baseline_method == "rolling") {
      baseline <- .rolling_percentile(trace, window_frames, 0.1)
    } else {
      baseline <- quantile(trace, 0.1)
    }

    # Compute dF/F
    dff <- (trace - baseline) / (baseline + 1e-10)

    # Threshold
    threshold <- threshold_sd * sd(dff)

    # Find events
    above <- dff > threshold
    events <- .find_events(above, dff, trace, time_vec, min_frames, frame_rate)

    if (!is.null(events) && nrow(events) > 0) {
      events$cell_id <- cell
      all_events[[cell]] <- events
    }
  }

  events_df <- do.call(rbind, all_events)

  if (is.null(events_df) || nrow(events_df) == 0) {
    events_df <- data.frame(
      cell_id = integer(),
      start_time = numeric(),
      end_time = numeric(),
      peak_time = numeric(),
      duration = numeric(),
      amplitude = numeric(),
      peak_dff = numeric(),
      rise_time = numeric(),
      decay_time = numeric(),
      auc = numeric()
    )
  }

  structure(
    list(
      events = events_df,
      n_events = nrow(events_df),
      n_cells = n_cells,
      frame_rate = frame_rate,
      threshold_sd = threshold_sd,
      per_cell = table(events_df$cell_id)
    ),
    class = "calcium_transients"
  )
}

.rolling_percentile <- function(x, window, p = 0.1) {
  n <- length(x)
  result <- rep(NA, n)
  half_w <- window %/% 2

  for (i in seq_len(n)) {
    start <- max(1, i - half_w)
    end <- min(n, i + half_w)
    result[i] <- quantile(x[start:end], p)
  }

  result
}

.find_events <- function(above, dff, trace, time_vec, min_frames, frame_rate) {
  events <- list()

  # Find contiguous regions above threshold
  starts <- which(diff(c(FALSE, above)) == 1)
  ends <- which(diff(c(above, FALSE)) == -1)

  if (length(starts) == 0) {
    return(data.frame(
      start_time = numeric(),
      end_time = numeric(),
      peak_time = numeric(),
      duration = numeric(),
      amplitude = numeric(),
      peak_dff = numeric(),
      rise_time = numeric(),
      decay_time = numeric(),
      auc = numeric()
    ))
  }

  for (i in seq_along(starts)) {
    if ((ends[i] - starts[i] + 1) < min_frames) next

    event_idx <- starts[i]:ends[i]
    event_dff <- dff[event_idx]
    event_trace <- trace[event_idx]

    peak_idx <- which.max(event_dff)
    peak_frame <- event_idx[peak_idx]

    # Kinetics
    rise_frames <- peak_idx
    decay_frames <- length(event_idx) - peak_idx

    events[[length(events) + 1]] <- data.frame(
      start_time = time_vec[starts[i]],
      end_time = time_vec[ends[i]],
      peak_time = time_vec[peak_frame],
      duration = (ends[i] - starts[i] + 1) / frame_rate,
      amplitude = max(event_trace) - min(event_trace[1:peak_idx]),
      peak_dff = max(event_dff),
      rise_time = rise_frames / frame_rate,
      decay_time = decay_frames / frame_rate,
      auc = sum(event_dff) / frame_rate
    )
  }

  do.call(rbind, events)
}

#' Summarize transient properties by condition
#'
#' Compare calcium transient properties across experimental conditions.
#'
#' @param transients_list Named list of transient objects (one per condition)
#' @param metrics Metrics to compare
#'
#' @return Summary data frame with statistics per condition
#' @export
summarize_transients_by_condition <- function(transients_list,
                                               metrics = c("frequency", "amplitude",
                                                          "duration", "rise_time")) {
  conditions <- names(transients_list)

  summaries <- lapply(conditions, function(cond) {
    tr <- transients_list[[cond]]
    events <- tr$events
    n_cells <- tr$n_cells
    total_time <- max(events$end_time, na.rm = TRUE)

    # Per-cell frequency
    cell_counts <- table(factor(events$cell_id, levels = 1:n_cells))
    frequencies <- as.numeric(cell_counts) / total_time * 60  # events per minute

    data.frame(
      condition = cond,
      n_cells = n_cells,
      n_events = tr$n_events,
      frequency_mean = mean(frequencies),
      frequency_sd = sd(frequencies),
      amplitude_mean = mean(events$amplitude, na.rm = TRUE),
      amplitude_sd = sd(events$amplitude, na.rm = TRUE),
      duration_mean = mean(events$duration, na.rm = TRUE),
      duration_sd = sd(events$duration, na.rm = TRUE),
      rise_time_mean = mean(events$rise_time, na.rm = TRUE),
      rise_time_sd = sd(events$rise_time, na.rm = TRUE),
      decay_time_mean = mean(events$decay_time, na.rm = TRUE),
      decay_time_sd = sd(events$decay_time, na.rm = TRUE),
      peak_dff_mean = mean(events$peak_dff, na.rm = TRUE),
      peak_dff_sd = sd(events$peak_dff, na.rm = TRUE)
    )
  })

  do.call(rbind, summaries)
}

#' Compare treatment conditions
#'
#' Statistical comparison of calcium activity between conditions.
#'
#' @param traces_list Named list of trace matrices (one per condition)
#' @param frame_rate Frame rate in Hz
#' @param test Statistical test: "wilcox", "t.test", "anova", "mixed"
#' @param paired Logical, are conditions paired (same cells)?
#' @param metrics Metrics to compare
#'
#' @return List with test results for each metric
#' @export
#'
#' @examples
#' \dontrun{
#' results <- compare_conditions(
#'   list(control = ctrl_traces, drug = drug_traces),
#'   frame_rate = 10,
#'   test = "wilcox"
#' )
#' }
compare_conditions <- function(traces_list, frame_rate = 10,
                                test = "wilcox", paired = FALSE,
                                metrics = c("mean_activity", "event_rate",
                                           "amplitude", "synchrony")) {
  conditions <- names(traces_list)
  n_conditions <- length(conditions)

  # Extract metrics for each condition
  all_metrics <- list()

  for (cond in conditions) {
    traces <- traces_list[[cond]]
    n_cells <- nrow(traces)

    # Detect transients
    transients <- detect_transients(traces, frame_rate = frame_rate)

    # Compute metrics per cell
    cell_metrics <- data.frame(
      condition = cond,
      cell_id = seq_len(n_cells),
      mean_activity = rowMeans(traces),
      sd_activity = apply(traces, 1, sd),
      max_activity = apply(traces, 1, max)
    )

    # Event rates
    total_time <- ncol(traces) / frame_rate / 60  # minutes
    event_counts <- table(factor(transients$events$cell_id, levels = 1:n_cells))
    cell_metrics$event_rate <- as.numeric(event_counts) / total_time

    # Mean amplitude per cell
    amp_by_cell <- tapply(transients$events$amplitude,
                          factor(transients$events$cell_id, levels = 1:n_cells),
                          mean, na.rm = TRUE)
    cell_metrics$amplitude <- as.numeric(amp_by_cell)
    cell_metrics$amplitude[is.na(cell_metrics$amplitude)] <- 0

    all_metrics[[cond]] <- cell_metrics
  }

  metrics_df <- do.call(rbind, all_metrics)

  # Add synchrony if requested
  if ("synchrony" %in% metrics) {
    for (cond in conditions) {
      traces <- traces_list[[cond]]
      sync <- compute_synchrony_index(traces)
      metrics_df$synchrony[metrics_df$condition == cond] <- sync$mean_synchrony
    }
  }

  # Perform statistical tests
  results <- list()

  for (metric in intersect(metrics, names(metrics_df))) {
    if (metric %in% c("condition", "cell_id")) next

    if (n_conditions == 2) {
      # Two-group comparison
      group1 <- metrics_df[[metric]][metrics_df$condition == conditions[1]]
      group2 <- metrics_df[[metric]][metrics_df$condition == conditions[2]]

      group1 <- group1[!is.na(group1)]
      group2 <- group2[!is.na(group2)]

      if (test == "wilcox") {
        test_result <- wilcox.test(group1, group2, paired = paired)
      } else if (test == "t.test") {
        test_result <- t.test(group1, group2, paired = paired)
      }

      results[[metric]] <- list(
        test = test,
        statistic = test_result$statistic,
        p.value = test_result$p.value,
        mean_diff = mean(group2, na.rm = TRUE) - mean(group1, na.rm = TRUE),
        effect_size = (mean(group2, na.rm = TRUE) - mean(group1, na.rm = TRUE)) /
          sqrt((var(group1, na.rm = TRUE) + var(group2, na.rm = TRUE)) / 2),
        group_means = c(mean(group1, na.rm = TRUE), mean(group2, na.rm = TRUE)),
        group_sds = c(sd(group1, na.rm = TRUE), sd(group2, na.rm = TRUE))
      )
      names(results[[metric]]$group_means) <- conditions
      names(results[[metric]]$group_sds) <- conditions

    } else {
      # Multi-group comparison
      if (test == "anova" || test == "kruskal") {
        formula_str <- paste(metric, "~ condition")
        if (test == "anova") {
          test_result <- summary(aov(as.formula(formula_str), data = metrics_df))
          results[[metric]] <- list(
            test = "anova",
            F_value = test_result[[1]]["condition", "F value"],
            p.value = test_result[[1]]["condition", "Pr(>F)"]
          )
        } else {
          test_result <- kruskal.test(as.formula(formula_str), data = metrics_df)
          results[[metric]] <- list(
            test = "kruskal",
            statistic = test_result$statistic,
            p.value = test_result$p.value
          )
        }
      }
    }
  }

  structure(
    list(
      results = results,
      metrics = metrics_df,
      conditions = conditions,
      test = test,
      paired = paired
    ),
    class = "condition_comparison"
  )
}

#' Compute synchrony index
#'
#' Measure synchronization of calcium activity across cells.
#'
#' @param traces Matrix of traces (cells x time)
#' @param method Synchrony method: "correlation", "phase", "event_coincidence"
#' @param window_size Window for sliding synchrony (NULL for global)
#'
#' @return List with synchrony measures
#' @export
compute_synchrony_index <- function(traces, method = "correlation",
                                     window_size = NULL) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  if (method == "correlation") {
    # Mean pairwise correlation
    cor_mat <- cor(t(traces))
    upper_tri <- cor_mat[upper.tri(cor_mat)]
    mean_sync <- mean(upper_tri, na.rm = TRUE)
    sd_sync <- sd(upper_tri, na.rm = TRUE)

  } else if (method == "phase") {
    # Phase synchrony using Hilbert transform
    if (!requireNamespace("signal", quietly = TRUE)) {
      warning("signal package not available, using correlation method")
      return(compute_synchrony_index(traces, method = "correlation"))
    }

    phases <- t(apply(traces, 1, function(x) {
      analytic <- signal::hilbert(x - mean(x))
      Arg(analytic)
    }))

    # Mean phase coherence
    phase_diff <- matrix(0, n_cells, n_cells)
    for (i in 1:(n_cells - 1)) {
      for (j in (i + 1):n_cells) {
        phase_diff[i, j] <- abs(mean(exp(1i * (phases[i, ] - phases[j, ]))))
      }
    }
    upper_tri <- phase_diff[upper.tri(phase_diff)]
    mean_sync <- mean(upper_tri)
    sd_sync <- sd(upper_tri)

  } else if (method == "event_coincidence") {
    # Event-based coincidence
    threshold <- apply(traces, 1, function(x) mean(x) + 2 * sd(x))
    events <- traces > threshold

    # Fraction of time multiple cells active
    n_active <- colSums(events)
    mean_sync <- mean(n_active > 1) / mean(n_active > 0)
    sd_sync <- NA
  }

  # Sliding window synchrony if requested
  if (!is.null(window_size)) {
    n_windows <- n_time - window_size + 1
    sync_over_time <- sapply(seq_len(n_windows), function(t) {
      window_traces <- traces[, t:(t + window_size - 1)]
      compute_synchrony_index(window_traces, method = method)$mean_synchrony
    })
  } else {
    sync_over_time <- NULL
  }

  list(
    mean_synchrony = mean_sync,
    sd_synchrony = sd_sync,
    method = method,
    synchrony_time_course = sync_over_time
  )
}

#' Detect network bursts
#'
#' Identify population-level burst events.
#'
#' @param traces Matrix of traces (cells x time)
#' @param frame_rate Frame rate in Hz
#' @param min_participation Minimum fraction of cells participating
#' @param threshold_sd Threshold for cell activation
#'
#' @return List with burst times and properties
#' @export
detect_network_bursts <- function(traces, frame_rate = 10,
                                   min_participation = 0.2,
                                   threshold_sd = 2) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Detect when each cell is active
  thresholds <- apply(traces, 1, function(x) mean(x) + threshold_sd * sd(x))
  active <- traces > thresholds

  # Compute participation over time
  participation <- colSums(active) / n_cells

  # Find bursts
  above <- participation > min_participation

  if (!any(above)) {
    return(structure(
      list(
        bursts = data.frame(),
        n_bursts = 0,
        participation = participation
      ),
      class = "network_bursts"
    ))
  }

  starts <- which(diff(c(FALSE, above)) == 1)
  ends <- which(diff(c(above, FALSE)) == -1)

  bursts <- data.frame(
    start_frame = starts,
    end_frame = ends,
    start_time = starts / frame_rate,
    end_time = ends / frame_rate,
    duration = (ends - starts + 1) / frame_rate,
    peak_participation = sapply(seq_along(starts), function(i) {
      max(participation[starts[i]:ends[i]])
    }),
    mean_participation = sapply(seq_along(starts), function(i) {
      mean(participation[starts[i]:ends[i]])
    })
  )

  # Compute inter-burst intervals
  if (nrow(bursts) > 1) {
    bursts$ibi <- c(NA, diff(bursts$start_time))
  } else {
    bursts$ibi <- NA
  }

  structure(
    list(
      bursts = bursts,
      n_bursts = nrow(bursts),
      participation = participation,
      frame_rate = frame_rate,
      burst_rate = nrow(bursts) / (n_time / frame_rate / 60)  # per minute
    ),
    class = "network_bursts"
  )
}

#' Dose-response analysis
#'
#' Fit dose-response curves for pharmacological data.
#'
#' @param doses Numeric vector of drug concentrations
#' @param responses Numeric vector or matrix of responses
#' @param model Model type: "hill", "logistic", "linear"
#'
#' @return Fitted dose-response model
#' @export
#'
#' @examples
#' \dontrun{
#' doses <- c(0, 0.1, 1, 10, 100)  # uM
#' responses <- c(1.0, 0.95, 0.7, 0.3, 0.1)  # normalized activity
#' fit <- fit_dose_response(doses, responses, model = "hill")
#' }
fit_dose_response <- function(doses, responses, model = "hill") {
  if (is.matrix(responses)) {
    responses <- rowMeans(responses)
  }

  df <- data.frame(dose = doses, response = responses)
  df <- df[order(df$dose), ]

  if (model == "hill") {
    # Hill equation: response = bottom + (top - bottom) / (1 + (IC50/dose)^hill)
    start_vals <- list(
      top = max(responses),
      bottom = min(responses),
      ic50 = median(doses[doses > 0]),
      hill = 1
    )

    fit <- tryCatch({
      nls(response ~ bottom + (top - bottom) / (1 + (ic50/dose)^hill),
          data = df[df$dose > 0, ],
          start = start_vals,
          control = list(maxiter = 100))
    }, error = function(e) {
      # Fall back to simpler model
      nls(response ~ bottom + (top - bottom) / (1 + (ic50/dose)),
          data = df[df$dose > 0, ],
          start = start_vals[1:3],
          control = list(maxiter = 100))
    })

    coefs <- coef(fit)

    result <- list(
      model = fit,
      type = "hill",
      ic50 = coefs["ic50"],
      ec50 = coefs["ic50"],  # Alias
      hill_coefficient = if ("hill" %in% names(coefs)) coefs["hill"] else 1,
      top = coefs["top"],
      bottom = coefs["bottom"],
      r_squared = 1 - sum(residuals(fit)^2) / sum((responses - mean(responses))^2)
    )

  } else if (model == "logistic") {
    # 4-parameter logistic
    fit <- nls(response ~ d + (a - d) / (1 + (dose/c)^b),
               data = df[df$dose > 0, ],
               start = list(a = max(responses), b = 1,
                          c = median(doses[doses > 0]), d = min(responses)))

    coefs <- coef(fit)
    result <- list(
      model = fit,
      type = "logistic",
      ic50 = coefs["c"],
      r_squared = 1 - sum(residuals(fit)^2) / sum((responses - mean(responses))^2)
    )

  } else if (model == "linear") {
    fit <- lm(response ~ log10(dose + min(doses[doses > 0])/10), data = df)
    result <- list(
      model = fit,
      type = "linear",
      slope = coef(fit)[2],
      r_squared = summary(fit)$r.squared
    )
  }

  result$data <- df
  structure(result, class = "dose_response")
}

#' Plot dose-response curve
#'
#' @param dr Dose-response object from fit_dose_response
#' @param log_scale Use log scale for dose axis
#'
#' @return ggplot object
#' @export
plot_dose_response <- function(dr, log_scale = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  df <- dr$data

  # Generate prediction curve
  if (log_scale) {
    dose_seq <- 10^seq(log10(min(df$dose[df$dose > 0]) / 10),
                       log10(max(df$dose) * 2), length.out = 100)
  } else {
    dose_seq <- seq(0, max(df$dose) * 1.1, length.out = 100)
  }

  pred_df <- data.frame(dose = dose_seq)
  pred_df$response <- predict(dr$model, newdata = pred_df)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = dose, y = response)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_line(data = pred_df, color = "blue") +
    ggplot2::labs(
      x = "Dose",
      y = "Response",
      title = sprintf("Dose-Response (IC50 = %.3g)", dr$ic50)
    ) +
    ggplot2::theme_minimal()

  if (log_scale) {
    p <- p + ggplot2::scale_x_log10()
  }

  # Add IC50 line
  if (!is.null(dr$ic50)) {
    mid_response <- (dr$top + dr$bottom) / 2
    p <- p +
      ggplot2::geom_hline(yintercept = mid_response, linetype = "dashed", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = dr$ic50, linetype = "dashed", alpha = 0.5)
  }

  p
}

#' Generate pharmacological analysis report
#'
#' Create summary report comparing treatment conditions.
#'
#' @param traces_list Named list of trace matrices
#' @param frame_rate Frame rate in Hz
#' @param conditions_info Optional data frame with condition metadata
#'
#' @return List with comprehensive analysis results
#' @export
pharmacological_report <- function(traces_list, frame_rate = 10,
                                   conditions_info = NULL) {
  conditions <- names(traces_list)

  # Transient analysis
  transients <- lapply(traces_list, function(tr) {
    detect_transients(tr, frame_rate = frame_rate)
  })
  names(transients) <- conditions

  transient_summary <- summarize_transients_by_condition(transients)

  # Network bursts
  bursts <- lapply(traces_list, function(tr) {
    detect_network_bursts(tr, frame_rate = frame_rate)
  })
  names(bursts) <- conditions

  burst_summary <- do.call(rbind, lapply(conditions, function(cond) {
    b <- bursts[[cond]]
    data.frame(
      condition = cond,
      n_bursts = b$n_bursts,
      burst_rate = b$burst_rate,
      mean_duration = mean(b$bursts$duration, na.rm = TRUE),
      mean_participation = mean(b$bursts$mean_participation, na.rm = TRUE)
    )
  }))

  # Synchrony
  synchrony <- lapply(traces_list, function(tr) {
    compute_synchrony_index(tr, method = "correlation")
  })
  names(synchrony) <- conditions

  sync_summary <- data.frame(
    condition = conditions,
    mean_synchrony = sapply(synchrony, `[[`, "mean_synchrony"),
    sd_synchrony = sapply(synchrony, `[[`, "sd_synchrony")
  )

  # Statistical comparisons
  if (length(conditions) >= 2) {
    comparisons <- compare_conditions(traces_list, frame_rate = frame_rate)
  } else {
    comparisons <- NULL
  }

  list(
    transient_summary = transient_summary,
    burst_summary = burst_summary,
    synchrony_summary = sync_summary,
    comparisons = comparisons,
    transients = transients,
    bursts = bursts,
    synchrony = synchrony
  )
}

#' Print method for condition comparison
#'
#' @param x Condition comparison object
#' @param ... Ignored
#'
#' @export
print.condition_comparison <- function(x, ...) {
  cat("Condition Comparison Results\n")
  cat("============================\n")
  cat("Conditions:", paste(x$conditions, collapse = " vs "), "\n")
  cat("Test:", x$test, ifelse(x$paired, "(paired)", "(unpaired)"), "\n\n")

  for (metric in names(x$results)) {
    res <- x$results[[metric]]
    cat(sprintf("%s:\n", metric))
    cat(sprintf("  p-value: %.4g %s\n", res$p.value,
                ifelse(res$p.value < 0.05, "*", "")))
    if (!is.null(res$effect_size)) {
      cat(sprintf("  Effect size (Cohen's d): %.2f\n", res$effect_size))
    }
    cat(sprintf("  Group means: %s\n",
                paste(sprintf("%s=%.3f", names(res$group_means), res$group_means),
                      collapse = ", ")))
    cat("\n")
  }

  invisible(x)
}
