#' Advanced Pharmacological Analysis for Cell Cultures
#'
#' Specialized analysis tools for in vitro calcium imaging with pharmacological
#' treatments. Includes responder analysis, wave propagation, baseline changes,
#' and temporal response characterization.
#'
#' @name advanced_pharmacology
#' @keywords internal
NULL

# ============================================================================
# Responder Analysis
# ============================================================================

#' Classify cells as responders or non-responders
#'
#' Identify which cells respond to treatment and characterize response types.
#'
#' @param traces_baseline Matrix of traces during baseline period (cells x time)
#' @param traces_treatment Matrix of traces during treatment (cells x time)
#' @param method Classification method: "threshold", "statistical", "ml"
#' @param threshold_sd SD threshold for response detection
#' @param alpha Significance level for statistical method
#'
#' @return List with responder classification and statistics
#' @export
#'
#' @examples
#' \dontrun{
#' resp <- classify_responders(baseline_traces, drug_traces)
#' table(resp$classification)  # responder counts
#' }
classify_responders <- function(traces_baseline, traces_treatment,
                                 method = "threshold", threshold_sd = 2,
                                 alpha = 0.05) {
  n_cells <- nrow(traces_baseline)

  # Compute baseline statistics per cell
  baseline_mean <- rowMeans(traces_baseline)
  baseline_sd <- apply(traces_baseline, 1, sd)

  # Compute treatment statistics per cell
  treatment_mean <- rowMeans(traces_treatment)
  treatment_max <- apply(traces_treatment, 1, max)
  treatment_min <- apply(traces_treatment, 1, min)

  if (method == "threshold") {
    # Threshold-based classification
    upper_threshold <- baseline_mean + threshold_sd * baseline_sd
    lower_threshold <- baseline_mean - threshold_sd * baseline_sd

    # Classify based on deviation from baseline
    response_magnitude <- treatment_mean - baseline_mean
    response_normalized <- response_magnitude / (baseline_sd + 1e-10)

    classification <- ifelse(
      treatment_max > upper_threshold, "activated",
      ifelse(treatment_min < lower_threshold, "inhibited", "non_responder")
    )

  } else if (method == "statistical") {
    # Statistical test per cell
    p_values <- sapply(seq_len(n_cells), function(i) {
      tryCatch({
        t.test(traces_treatment[i, ], traces_baseline[i, ])$p.value
      }, error = function(e) 1)
    })

    # Direction of change
    response_magnitude <- treatment_mean - baseline_mean
    response_normalized <- response_magnitude / (baseline_sd + 1e-10)

    classification <- ifelse(
      p_values < alpha & response_magnitude > 0, "activated",
      ifelse(p_values < alpha & response_magnitude < 0, "inhibited", "non_responder")
    )

    classification[p.adjust(p_values, method = "fdr") >= alpha] <- "non_responder"

  } else if (method == "ml") {
    # ML-based classification using features
    features <- extract_response_features(traces_baseline, traces_treatment)
    classification <- .ml_classify_responders(features)
    response_magnitude <- treatment_mean - baseline_mean
    response_normalized <- response_magnitude / (baseline_sd + 1e-10)
  }

  # Compute response metrics
  results <- data.frame(
    cell_id = seq_len(n_cells),
    classification = classification,
    baseline_mean = baseline_mean,
    baseline_sd = baseline_sd,
    treatment_mean = treatment_mean,
    response_magnitude = response_magnitude,
    response_normalized = response_normalized,
    fold_change = treatment_mean / (baseline_mean + 1e-10)
  )

  # Summary statistics
  responder_summary <- table(classification)
  responder_fraction <- sum(classification != "non_responder") / n_cells

  structure(
    list(
      results = results,
      classification = classification,
      summary = responder_summary,
      responder_fraction = responder_fraction,
      activated_fraction = mean(classification == "activated"),
      inhibited_fraction = mean(classification == "inhibited"),
      method = method,
      threshold_sd = threshold_sd
    ),
    class = "responder_analysis"
  )
}

#' Extract response features for classification
#'
#' @param traces_baseline Baseline traces
#' @param traces_treatment Treatment traces
#'
#' @return Feature matrix
#' @keywords internal
extract_response_features <- function(traces_baseline, traces_treatment) {
  n_cells <- nrow(traces_baseline)

  data.frame(
    # Baseline features
    baseline_mean = rowMeans(traces_baseline),
    baseline_sd = apply(traces_baseline, 1, sd),
    baseline_cv = apply(traces_baseline, 1, function(x) sd(x)/mean(x)),

    # Treatment features
    treatment_mean = rowMeans(traces_treatment),
    treatment_sd = apply(traces_treatment, 1, sd),
    treatment_max = apply(traces_treatment, 1, max),
    treatment_min = apply(traces_treatment, 1, min),

    # Change features
    mean_change = rowMeans(traces_treatment) - rowMeans(traces_baseline),
    max_change = apply(traces_treatment, 1, max) - rowMeans(traces_baseline),
    sd_change = apply(traces_treatment, 1, sd) - apply(traces_baseline, 1, sd),

    # Temporal features
    time_to_peak = apply(traces_treatment, 1, which.max),
    time_to_min = apply(traces_treatment, 1, which.min)
  )
}

.ml_classify_responders <- function(features) {
  # Simple clustering-based classification
  features_scaled <- scale(features)
  features_scaled[is.na(features_scaled)] <- 0

  km <- kmeans(features_scaled, centers = 3, nstart = 25)

  # Identify clusters by mean change
  cluster_means <- tapply(features$mean_change, km$cluster, mean)
  activated_cluster <- which.max(cluster_means)
  inhibited_cluster <- which.min(cluster_means)
  neutral_cluster <- setdiff(1:3, c(activated_cluster, inhibited_cluster))

  classification <- ifelse(
    km$cluster == activated_cluster, "activated",
    ifelse(km$cluster == inhibited_cluster, "inhibited", "non_responder")
  )

  classification
}

#' Subtype cells by response pattern
#'
#' Cluster cells into subtypes based on their response dynamics.
#'
#' @param traces Matrix of traces (cells x time)
#' @param n_subtypes Number of subtypes to identify
#' @param method Clustering method: "kmeans", "hierarchical", "dtw"
#' @param normalize Normalize traces before clustering
#'
#' @return List with subtype assignments and characteristics
#' @export
subtype_by_response <- function(traces, n_subtypes = 4, method = "kmeans",
                                 normalize = TRUE) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  if (normalize) {
    # Z-score normalize each cell
    traces_norm <- t(scale(t(traces)))
    traces_norm[is.na(traces_norm)] <- 0
  } else {
    traces_norm <- traces
  }

  if (method == "kmeans") {
    km <- kmeans(traces_norm, centers = n_subtypes, nstart = 25)
    assignments <- km$cluster
    centers <- km$centers

  } else if (method == "hierarchical") {
    dist_mat <- dist(traces_norm)
    hc <- hclust(dist_mat, method = "ward.D2")
    assignments <- cutree(hc, k = n_subtypes)
    centers <- t(sapply(1:n_subtypes, function(k) {
      colMeans(traces_norm[assignments == k, , drop = FALSE])
    }))

  } else if (method == "dtw") {
    if (!requireNamespace("dtw", quietly = TRUE)) {
      warning("dtw package not available, using kmeans")
      return(subtype_by_response(traces, n_subtypes, "kmeans", normalize))
    }

    # DTW distance matrix
    dist_mat <- matrix(0, n_cells, n_cells)
    for (i in 1:(n_cells - 1)) {
      for (j in (i + 1):n_cells) {
        d <- dtw::dtw(traces_norm[i, ], traces_norm[j, ])$distance
        dist_mat[i, j] <- d
        dist_mat[j, i] <- d
      }
    }

    hc <- hclust(as.dist(dist_mat), method = "ward.D2")
    assignments <- cutree(hc, k = n_subtypes)
    centers <- t(sapply(1:n_subtypes, function(k) {
      colMeans(traces_norm[assignments == k, , drop = FALSE])
    }))
  }

  # Characterize each subtype
  subtype_stats <- lapply(1:n_subtypes, function(k) {
    members <- which(assignments == k)
    subtype_traces <- traces[members, , drop = FALSE]

    list(
      n_cells = length(members),
      fraction = length(members) / n_cells,
      mean_response = mean(rowMeans(subtype_traces)),
      peak_time = which.max(colMeans(subtype_traces)),
      mean_trace = colMeans(subtype_traces)
    )
  })

  # Label subtypes by characteristic
  peak_times <- sapply(subtype_stats, `[[`, "peak_time")
  mean_responses <- sapply(subtype_stats, `[[`, "mean_response")

  labels <- paste0("subtype_", 1:n_subtypes)
  # Could add more sophisticated labeling based on characteristics

  structure(
    list(
      assignments = assignments,
      centers = centers,
      subtype_stats = subtype_stats,
      labels = labels,
      n_subtypes = n_subtypes,
      method = method
    ),
    class = "response_subtypes"
  )
}

# ============================================================================
# Baseline and Calcium Load Analysis
# ============================================================================

#' Analyze baseline calcium levels
#'
#' Compare resting calcium levels between conditions.
#'
#' @param traces Matrix of traces (cells x time)
#' @param method Baseline estimation: "mean", "median", "percentile"
#' @param percentile Percentile for baseline estimation
#' @param frame_rate Frame rate (for time conversion)
#'
#' @return List with baseline statistics
#' @export
analyze_baseline <- function(traces, method = "percentile", percentile = 10,
                              frame_rate = 10) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  baseline <- switch(method,
    mean = rowMeans(traces),
    median = apply(traces, 1, median),
    percentile = apply(traces, 1, quantile, probs = percentile/100)
  )

  # Baseline stability (CV of rolling baseline)
  window_frames <- min(100, n_time %/% 4)
  stability <- apply(traces, 1, function(x) {
    rolling_bl <- sapply(seq(1, n_time - window_frames, by = window_frames %/% 2),
                         function(t) mean(x[t:(t + window_frames - 1)]))
    sd(rolling_bl) / mean(rolling_bl)
  })

  # Baseline drift (slope of linear fit)
  time_vec <- seq_len(n_time) / frame_rate
  drift <- apply(traces, 1, function(x) {
    coef(lm(x ~ time_vec))[2]
  })

  list(
    baseline = baseline,
    mean_baseline = mean(baseline),
    sd_baseline = sd(baseline),
    stability = stability,
    drift = drift,
    method = method
  )
}

#' Compare baseline between conditions
#'
#' @param baseline1 Baseline values for condition 1
#' @param baseline2 Baseline values for condition 2
#' @param paired Are the measurements paired?
#'
#' @return Statistical comparison results
#' @export
compare_baseline <- function(baseline1, baseline2, paired = FALSE) {
  test_result <- wilcox.test(baseline1, baseline2, paired = paired)

  diff <- baseline2 - baseline1
  effect_size <- mean(diff) / sd(diff)

  list(
    condition1_mean = mean(baseline1),
    condition2_mean = mean(baseline2),
    difference = mean(diff),
    percent_change = 100 * mean(diff) / mean(baseline1),
    p_value = test_result$p.value,
    effect_size = effect_size
  )
}

#' Compute calcium load (AUC)
#'
#' Total calcium exposure over time.
#'
#' @param traces Matrix of traces (cells x time)
#' @param baseline Baseline values (or NULL to compute)
#' @param frame_rate Frame rate for time conversion
#' @param above_baseline Only count activity above baseline
#'
#' @return Calcium load per cell
#' @export
compute_calcium_load <- function(traces, baseline = NULL, frame_rate = 10,
                                  above_baseline = TRUE) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  if (is.null(baseline)) {
    baseline <- apply(traces, 1, quantile, probs = 0.1)
  }

  # Compute AUC
  dt <- 1 / frame_rate

  if (above_baseline) {
    # Only integrate positive deviations from baseline
    calcium_load <- sapply(seq_len(n_cells), function(i) {
      deviation <- traces[i, ] - baseline[i]
      sum(pmax(deviation, 0)) * dt
    })
  } else {
    # Total integral
    calcium_load <- rowSums(traces) * dt
  }

  list(
    calcium_load = calcium_load,
    mean_load = mean(calcium_load),
    sd_load = sd(calcium_load),
    load_per_minute = calcium_load / (n_time / frame_rate / 60)
  )
}

# ============================================================================
# Calcium Wave Propagation
# ============================================================================

#' Detect calcium waves
#'
#' Identify propagating waves of calcium activity.
#'
#' @param traces Matrix of traces (cells x time)
#' @param positions Matrix of cell positions (cells x 2: x, y coordinates)
#' @param frame_rate Frame rate in Hz
#' @param threshold_sd Threshold for wave detection
#' @param min_cells Minimum cells participating in wave
#' @param max_delay Maximum delay between cells (seconds)
#'
#' @return List with detected waves and propagation properties
#' @export
#'
#' @examples
#' \dontrun{
#' waves <- detect_calcium_waves(traces, cell_positions, frame_rate = 10)
#' }
detect_calcium_waves <- function(traces, positions, frame_rate = 10,
                                  threshold_sd = 2, min_cells = 5,
                                  max_delay = 2) {

  n_cells <- nrow(traces)
  n_time <- ncol(traces)
  max_delay_frames <- max_delay * frame_rate

  # Detect activation times for each cell
  thresholds <- apply(traces, 1, function(x) mean(x) + threshold_sd * sd(x))

  activation_times <- lapply(seq_len(n_cells), function(i) {
    above <- traces[i, ] > thresholds[i]
    starts <- which(diff(c(FALSE, above)) == 1)
    starts / frame_rate
  })

  # Find coordinated activations (potential waves)
  waves <- list()
  used_activations <- lapply(1:n_cells, function(i) logical(length(activation_times[[i]])))

  for (t in seq(0, n_time / frame_rate, by = 0.5)) {
    # Find cells activating around this time
    participating <- sapply(seq_len(n_cells), function(i) {
      times <- activation_times[[i]]
      any(times >= t & times < t + max_delay & !used_activations[[i]])
    })

    if (sum(participating) >= min_cells) {
      # Get exact activation times
      wave_times <- sapply(which(participating), function(i) {
        times <- activation_times[[i]]
        valid <- times >= t & times < t + max_delay & !used_activations[[i]]
        if (any(valid)) min(times[valid]) else NA
      })

      wave_cells <- which(participating)[!is.na(wave_times)]
      wave_times <- wave_times[!is.na(wave_times)]

      if (length(wave_cells) >= min_cells) {
        # Analyze propagation
        wave_positions <- positions[wave_cells, , drop = FALSE]
        wave_analysis <- .analyze_wave_propagation(wave_positions, wave_times)

        waves[[length(waves) + 1]] <- list(
          cells = wave_cells,
          activation_times = wave_times,
          start_time = min(wave_times),
          duration = max(wave_times) - min(wave_times),
          n_cells = length(wave_cells),
          speed = wave_analysis$speed,
          direction = wave_analysis$direction,
          origin = wave_analysis$origin,
          r_squared = wave_analysis$r_squared
        )

        # Mark activations as used
        for (i in seq_along(wave_cells)) {
          cell <- wave_cells[i]
          time <- wave_times[i]
          idx <- which.min(abs(activation_times[[cell]] - time))
          used_activations[[cell]][idx] <- TRUE
        }
      }
    }
  }

  structure(
    list(
      waves = waves,
      n_waves = length(waves),
      mean_speed = if (length(waves) > 0) mean(sapply(waves, `[[`, "speed"), na.rm = TRUE) else NA,
      mean_size = if (length(waves) > 0) mean(sapply(waves, `[[`, "n_cells")) else NA,
      wave_rate = length(waves) / (n_time / frame_rate / 60)  # waves per minute
    ),
    class = "calcium_waves"
  )
}

.analyze_wave_propagation <- function(positions, times) {
  # Fit linear model: time = a + b*x + c*y
  if (nrow(positions) < 4) {
    return(list(speed = NA, direction = NA, origin = NA, r_squared = NA))
  }

  df <- data.frame(
    time = times,
    x = positions[, 1],
    y = positions[, 2]
  )

  fit <- lm(time ~ x + y, data = df)
  coefs <- coef(fit)

  # Direction and speed from gradient
  dx <- coefs["x"]
  dy <- coefs["y"]

  speed <- 1 / sqrt(dx^2 + dy^2)  # units per second
  direction <- atan2(dy, dx) * 180 / pi  # degrees

  # Origin: where time = min(times)
  origin_cell <- which.min(times)
  origin <- positions[origin_cell, ]

  list(
    speed = speed,
    direction = direction,
    origin = origin,
    r_squared = summary(fit)$r.squared
  )
}

#' Identify wave initiation sites
#'
#' Find cells that frequently initiate calcium waves.
#'
#' @param waves Waves object from detect_calcium_waves
#' @param n_cells Total number of cells
#' @param positions Cell positions matrix
#'
#' @return Data frame with initiation site statistics
#' @export
identify_initiation_sites <- function(waves, n_cells, positions = NULL) {
  if (waves$n_waves == 0) {
    return(data.frame(
      cell_id = integer(),
      initiation_count = integer(),
      initiation_fraction = numeric()
    ))
  }

  # Count initiations per cell
  initiator_cells <- sapply(waves$waves, function(w) {
    w$cells[which.min(w$activation_times)]
  })

  initiation_counts <- table(factor(initiator_cells, levels = 1:n_cells))

  df <- data.frame(
    cell_id = 1:n_cells,
    initiation_count = as.numeric(initiation_counts),
    initiation_fraction = as.numeric(initiation_counts) / waves$n_waves
  )

  if (!is.null(positions)) {
    df$x <- positions[, 1]
    df$y <- positions[, 2]
  }

  # Sort by initiation count
  df <- df[order(df$initiation_count, decreasing = TRUE), ]

  df
}

#' Compute wave participation
#'
#' @param waves Waves object
#' @param n_cells Total number of cells
#'
#' @return Participation statistics per cell
#' @export
compute_wave_participation <- function(waves, n_cells) {
  if (waves$n_waves == 0) {
    return(list(
      participation_count = rep(0, n_cells),
      participation_fraction = rep(0, n_cells)
    ))
  }

  participation <- sapply(1:n_cells, function(i) {
    sum(sapply(waves$waves, function(w) i %in% w$cells))
  })

  list(
    participation_count = participation,
    participation_fraction = participation / waves$n_waves,
    mean_participation = mean(participation),
    highly_participating = which(participation > mean(participation) + sd(participation))
  )
}

# ============================================================================
# Temporal Drug Response
# ============================================================================

#' Analyze temporal response to drug
#'
#' Characterize how response evolves over time after drug application.
#'
#' @param traces Matrix of traces (cells x time)
#' @param drug_onset Frame number when drug was applied
#' @param frame_rate Frame rate in Hz
#' @param bin_size Time bin size in seconds for averaging
#'
#' @return List with temporal response characteristics
#' @export
analyze_temporal_response <- function(traces, drug_onset, frame_rate = 10,
                                       bin_size = 10) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Time relative to drug onset
  time_vec <- (seq_len(n_time) - drug_onset) / frame_rate

  # Bin the data
  bin_frames <- bin_size * frame_rate
  n_bins <- floor(n_time / bin_frames)

  binned_response <- matrix(NA, n_cells, n_bins)
  bin_times <- numeric(n_bins)

  for (b in seq_len(n_bins)) {
    start <- (b - 1) * bin_frames + 1
    end <- b * bin_frames
    binned_response[, b] <- rowMeans(traces[, start:end, drop = FALSE])
    bin_times[b] <- mean(time_vec[start:end])
  }

  # Population mean response
  mean_response <- colMeans(binned_response)
  sem_response <- apply(binned_response, 2, function(x) sd(x) / sqrt(length(x)))

  # Find key time points
  baseline_bins <- which(bin_times < 0)
  post_drug_bins <- which(bin_times >= 0)

  baseline_level <- if (length(baseline_bins) > 0) mean(mean_response[baseline_bins]) else mean_response[1]

  # Time to peak
  post_drug_response <- mean_response[post_drug_bins]
  peak_idx <- which.max(abs(post_drug_response - baseline_level))
  time_to_peak <- bin_times[post_drug_bins[peak_idx]]
  peak_response <- post_drug_response[peak_idx]

  # Time to half-max
  half_max <- baseline_level + (peak_response - baseline_level) / 2
  half_max_idx <- which(abs(post_drug_response - half_max) == min(abs(post_drug_response - half_max)))[1]
  time_to_half_max <- bin_times[post_drug_bins[half_max_idx]]

  # Sustained vs transient (ratio of late to peak response)
  late_bins <- tail(post_drug_bins, 3)
  late_response <- mean(mean_response[late_bins])
  sustained_ratio <- (late_response - baseline_level) / (peak_response - baseline_level + 1e-10)

  list(
    bin_times = bin_times,
    mean_response = mean_response,
    sem_response = sem_response,
    binned_response = binned_response,
    baseline_level = baseline_level,
    peak_response = peak_response,
    time_to_peak = time_to_peak,
    time_to_half_max = time_to_half_max,
    sustained_ratio = sustained_ratio,
    response_type = ifelse(sustained_ratio > 0.5, "sustained", "transient")
  )
}

#' Fit exponential response model
#'
#' Fit rise and decay kinetics of drug response.
#'
#' @param time_vec Time vector
#' @param response Response values
#' @param model_type Model type: "rise", "decay", "rise_decay"
#'
#' @return Fitted model parameters
#' @export
fit_response_kinetics <- function(time_vec, response, model_type = "rise_decay") {
  # Normalize response
  baseline <- mean(response[time_vec < 0])
  response_norm <- response - baseline

  # Fit based on model type
  if (model_type == "rise") {
    # Exponential rise: y = A * (1 - exp(-t/tau))
    positive_times <- time_vec >= 0
    fit <- tryCatch({
      nls(response_norm[positive_times] ~ A * (1 - exp(-time_vec[positive_times] / tau)),
          start = list(A = max(response_norm), tau = 1))
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      return(list(
        type = "rise",
        amplitude = coef(fit)["A"],
        tau_rise = coef(fit)["tau"],
        model = fit
      ))
    }

  } else if (model_type == "decay") {
    # Find peak and fit decay
    peak_idx <- which.max(response_norm)
    decay_data <- response_norm[peak_idx:length(response_norm)]
    decay_time <- time_vec[peak_idx:length(time_vec)] - time_vec[peak_idx]

    fit <- tryCatch({
      nls(decay_data ~ A * exp(-decay_time / tau) + B,
          start = list(A = decay_data[1], tau = 10, B = tail(decay_data, 1)))
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      return(list(
        type = "decay",
        amplitude = coef(fit)["A"],
        tau_decay = coef(fit)["tau"],
        plateau = coef(fit)["B"],
        model = fit
      ))
    }

  } else if (model_type == "rise_decay") {
    # Double exponential: y = A * (exp(-t/tau_d) - exp(-t/tau_r))
    positive_times <- time_vec >= 0
    positive_response <- response_norm[positive_times]
    positive_time <- time_vec[positive_times]

    fit <- tryCatch({
      nls(positive_response ~ A * (exp(-positive_time / tau_d) - exp(-positive_time / tau_r)),
          start = list(A = max(positive_response) * 2, tau_r = 0.5, tau_d = 10),
          control = list(maxiter = 100))
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      return(list(
        type = "rise_decay",
        amplitude = coef(fit)["A"],
        tau_rise = coef(fit)["tau_r"],
        tau_decay = coef(fit)["tau_d"],
        model = fit
      ))
    }
  }

  # Return NA if fitting failed
  list(type = model_type, amplitude = NA, tau_rise = NA, tau_decay = NA, model = NULL)
}

# ============================================================================
# Spatial Correlation
# ============================================================================

#' Compute distance-dependent correlation
#'
#' Analyze how correlation depends on cell-cell distance.
#'
#' @param traces Matrix of traces (cells x time)
#' @param positions Matrix of cell positions (cells x 2)
#' @param n_bins Number of distance bins
#'
#' @return List with distance-correlation relationship
#' @export
compute_spatial_correlation <- function(traces, positions, n_bins = 10) {
  n_cells <- nrow(traces)

  # Compute pairwise distances
  dist_mat <- as.matrix(dist(positions))

  # Compute pairwise correlations
  cor_mat <- cor(t(traces))

  # Bin by distance
  distances <- dist_mat[upper.tri(dist_mat)]
  correlations <- cor_mat[upper.tri(cor_mat)]

  # Remove self-correlations and NAs
  valid <- !is.na(correlations)
  distances <- distances[valid]
  correlations <- correlations[valid]

  # Create distance bins
  breaks <- quantile(distances, probs = seq(0, 1, length.out = n_bins + 1))
  bins <- cut(distances, breaks = breaks, labels = FALSE, include.lowest = TRUE)

  # Compute mean correlation per bin
  mean_dist <- tapply(distances, bins, mean)
  mean_corr <- tapply(correlations, bins, mean)
  sem_corr <- tapply(correlations, bins, function(x) sd(x) / sqrt(length(x)))

  # Fit exponential decay
  fit <- tryCatch({
    nls(mean_corr ~ a * exp(-mean_dist / lambda) + b,
        start = list(a = mean_corr[1], lambda = max(mean_dist) / 2, b = tail(mean_corr, 1)))
  }, error = function(e) NULL)

  space_constant <- if (!is.null(fit)) coef(fit)["lambda"] else NA

  list(
    distance_bins = mean_dist,
    mean_correlation = mean_corr,
    sem_correlation = sem_corr,
    space_constant = space_constant,
    fit = fit,
    all_distances = distances,
    all_correlations = correlations
  )
}

#' Test spatial correlation change with treatment
#'
#' @param spatial_corr1 Spatial correlation for condition 1
#' @param spatial_corr2 Spatial correlation for condition 2
#'
#' @return Statistical comparison
#' @export
compare_spatial_correlation <- function(spatial_corr1, spatial_corr2) {
  # Compare space constants
  if (!is.na(spatial_corr1$space_constant) && !is.na(spatial_corr2$space_constant)) {
    space_constant_change <- spatial_corr2$space_constant - spatial_corr1$space_constant
  } else {
    space_constant_change <- NA
  }

  # Compare overall correlation levels
  test <- wilcox.test(spatial_corr1$all_correlations, spatial_corr2$all_correlations)

  list(
    space_constant1 = spatial_corr1$space_constant,
    space_constant2 = spatial_corr2$space_constant,
    space_constant_change = space_constant_change,
    mean_corr1 = mean(spatial_corr1$all_correlations),
    mean_corr2 = mean(spatial_corr2$all_correlations),
    p_value = test$p.value
  )
}

# ============================================================================
# Recovery Analysis
# ============================================================================

#' Analyze recovery after treatment
#'
#' Characterize return to baseline after drug washout.
#'
#' @param traces Matrix of traces (cells x time)
#' @param baseline_period Frames for baseline period
#' @param treatment_period Frames for treatment period
#' @param recovery_period Frames for recovery period
#' @param frame_rate Frame rate in Hz
#'
#' @return List with recovery statistics
#' @export
analyze_recovery <- function(traces, baseline_period, treatment_period,
                              recovery_period, frame_rate = 10) {
  n_cells <- nrow(traces)

  # Extract periods
  baseline_traces <- traces[, baseline_period, drop = FALSE]
  treatment_traces <- traces[, treatment_period, drop = FALSE]
  recovery_traces <- traces[, recovery_period, drop = FALSE]

  # Compute means
  baseline_mean <- rowMeans(baseline_traces)
  treatment_mean <- rowMeans(treatment_traces)

  # Recovery over time
  recovery_time <- (seq_along(recovery_period) - 1) / frame_rate

  # Compute recovery fraction at each time point
  recovery_fraction <- t(sapply(seq_len(n_cells), function(i) {
    max_deviation <- treatment_mean[i] - baseline_mean[i]
    if (abs(max_deviation) < 1e-10) return(rep(1, length(recovery_period)))

    current_deviation <- recovery_traces[i, ] - baseline_mean[i]
    1 - current_deviation / max_deviation
  }))

  # Mean recovery across cells
  mean_recovery <- colMeans(recovery_fraction, na.rm = TRUE)

  # Time to X% recovery
  recovery_times <- sapply(c(0.5, 0.75, 0.9), function(threshold) {
    idx <- which(mean_recovery >= threshold)[1]
    if (is.na(idx)) return(NA)
    recovery_time[idx]
  })
  names(recovery_times) <- c("t50", "t75", "t90")

  # Final recovery level
  final_recovery <- mean(tail(mean_recovery, 10))

  # Per-cell recovery classification
  cell_final_recovery <- rowMeans(recovery_fraction[, tail(seq_len(ncol(recovery_fraction)), 10), drop = FALSE])
  recovery_classification <- ifelse(cell_final_recovery > 0.9, "full",
                                    ifelse(cell_final_recovery > 0.5, "partial", "none"))

  list(
    recovery_time = recovery_time,
    mean_recovery = mean_recovery,
    recovery_fraction = recovery_fraction,
    t50 = recovery_times["t50"],
    t75 = recovery_times["t75"],
    t90 = recovery_times["t90"],
    final_recovery = final_recovery,
    cell_recovery = cell_final_recovery,
    recovery_classification = recovery_classification,
    full_recovery_fraction = mean(recovery_classification == "full"),
    partial_recovery_fraction = mean(recovery_classification == "partial")
  )
}

# ============================================================================
# Excitability Metrics
# ============================================================================

#' Compute excitability metrics
#'
#' Measure threshold, sensitivity, and refractory properties.
#'
#' @param traces Matrix of traces (cells x time)
#' @param frame_rate Frame rate in Hz
#'
#' @return List with excitability measures
#' @export
compute_excitability <- function(traces, frame_rate = 10) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Detect events
  transients <- detect_transients(traces, frame_rate = frame_rate)
  events <- transients$events

  # Excitability metrics per cell
  excitability <- lapply(seq_len(n_cells), function(cell) {
    cell_events <- events[events$cell_id == cell, ]

    if (nrow(cell_events) < 2) {
      return(list(
        event_rate = nrow(cell_events) / (n_time / frame_rate / 60),
        threshold = NA,
        refractory_period = NA,
        mean_amplitude = if (nrow(cell_events) > 0) mean(cell_events$amplitude) else NA
      ))
    }

    # Estimate threshold (mean baseline + fraction of mean amplitude)
    trace <- traces[cell, ]
    baseline <- quantile(trace, 0.1)
    threshold <- baseline + 0.5 * mean(cell_events$amplitude)

    # Refractory period (minimum inter-event interval)
    ieis <- diff(cell_events$peak_time)
    refractory <- if (length(ieis) > 0) min(ieis) else NA

    list(
      event_rate = nrow(cell_events) / (n_time / frame_rate / 60),
      threshold = threshold,
      refractory_period = refractory,
      mean_amplitude = mean(cell_events$amplitude),
      mean_iei = mean(ieis),
      cv_iei = sd(ieis) / mean(ieis)
    )
  })

  # Summarize
  event_rates <- sapply(excitability, `[[`, "event_rate")
  thresholds <- sapply(excitability, `[[`, "threshold")
  refractory_periods <- sapply(excitability, `[[`, "refractory_period")

  list(
    cell_excitability = excitability,
    mean_event_rate = mean(event_rates, na.rm = TRUE),
    mean_threshold = mean(thresholds, na.rm = TRUE),
    mean_refractory = mean(refractory_periods, na.rm = TRUE),
    event_rates = event_rates,
    thresholds = thresholds,
    refractory_periods = refractory_periods
  )
}

# ============================================================================
# Transient Shape Classification
# ============================================================================

#' Classify transient shapes
#'
#' Categorize calcium transients by their waveform shape.
#'
#' @param traces Matrix of traces (cells x time)
#' @param frame_rate Frame rate in Hz
#' @param n_classes Number of shape classes
#'
#' @return List with shape classifications
#' @export
classify_transient_shapes <- function(traces, frame_rate = 10, n_classes = 4) {
  # Detect transients
  transients <- detect_transients(traces, frame_rate = frame_rate)
  events <- transients$events

  if (nrow(events) == 0) {
    return(list(classifications = character(), shape_stats = NULL))
  }

  # Extract waveforms for each event
  n_cells <- nrow(traces)
  window_before <- round(0.5 * frame_rate)  # 0.5 sec before

  window_after <- round(2 * frame_rate)  # 2 sec after

  waveforms <- list()
  event_info <- list()

  for (i in seq_len(nrow(events))) {
    cell <- events$cell_id[i]
    peak_frame <- round(events$peak_time[i] * frame_rate)

    start_frame <- max(1, peak_frame - window_before)
    end_frame <- min(ncol(traces), peak_frame + window_after)

    waveform <- traces[cell, start_frame:end_frame]

    # Normalize waveform
    waveform_norm <- (waveform - min(waveform)) / (max(waveform) - min(waveform) + 1e-10)

    # Resample to fixed length
    target_length <- window_before + window_after + 1
    if (length(waveform_norm) != target_length) {
      waveform_norm <- approx(seq_along(waveform_norm), waveform_norm,
                               n = target_length)$y
    }

    waveforms[[i]] <- waveform_norm
    event_info[[i]] <- events[i, ]
  }

  # Stack waveforms
  waveform_matrix <- do.call(rbind, waveforms)

  # Extract shape features
  shape_features <- t(apply(waveform_matrix, 1, function(w) {
    peak_idx <- which.max(w)
    rise_phase <- w[1:peak_idx]
    decay_phase <- w[peak_idx:length(w)]

    # Features
    c(
      rise_slope = if (length(rise_phase) > 1) max(diff(rise_phase)) else 0,
      decay_slope = if (length(decay_phase) > 1) min(diff(decay_phase)) else 0,
      rise_time = peak_idx / length(w),
      decay_time = (length(w) - peak_idx) / length(w),
      symmetry = abs(length(rise_phase) - length(decay_phase)) / length(w),
      plateau = mean(w[peak_idx:min(peak_idx + 5, length(w))]),
      n_peaks = sum(diff(sign(diff(w))) == -2),
      late_activity = mean(tail(w, 10))
    )
  }))

  # Cluster shapes
  km <- kmeans(scale(shape_features), centers = n_classes, nstart = 25)

  # Characterize each cluster
  cluster_waveforms <- lapply(1:n_classes, function(k) {
    colMeans(waveform_matrix[km$cluster == k, , drop = FALSE])
  })

  # Label clusters based on characteristics
  cluster_labels <- sapply(1:n_classes, function(k) {
    w <- cluster_waveforms[[k]]
    n_peaks <- sum(diff(sign(diff(w))) == -2)
    late_activity <- mean(tail(w, 10))
    decay_speed <- if (which.max(w) < length(w)) {
      (max(w) - tail(w, 1)) / (length(w) - which.max(w))
    } else 0

    if (n_peaks > 1) {
      "oscillatory"
    } else if (late_activity > 0.5) {
      "sustained"
    } else if (decay_speed > 0.05) {
      "fast_decay"
    } else {
      "slow_decay"
    }
  })

  # Add classifications to events
  events$shape_class <- km$cluster
  events$shape_label <- cluster_labels[km$cluster]

  list(
    events = events,
    cluster_assignments = km$cluster,
    cluster_labels = cluster_labels,
    cluster_waveforms = cluster_waveforms,
    shape_features = shape_features,
    shape_summary = table(cluster_labels),
    waveform_matrix = waveform_matrix
  )
}

#' Plot transient shapes
#'
#' @param shape_analysis Shape analysis result
#'
#' @return ggplot object
#' @export
plot_transient_shapes <- function(shape_analysis) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  n_classes <- length(shape_analysis$cluster_waveforms)
  n_time <- length(shape_analysis$cluster_waveforms[[1]])

  df <- do.call(rbind, lapply(1:n_classes, function(k) {
    data.frame(
      time = seq_len(n_time),
      amplitude = shape_analysis$cluster_waveforms[[k]],
      class = shape_analysis$cluster_labels[k],
      class_num = k
    )
  }))

  ggplot2::ggplot(df, ggplot2::aes(x = time, y = amplitude, color = class)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(~ class, scales = "free_y") +
    ggplot2::labs(
      title = "Transient Shape Classes",
      x = "Time (frames)",
      y = "Normalized Amplitude"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}

# ============================================================================
# Comprehensive Report
# ============================================================================

#' Generate comprehensive GBM analysis report
#'
#' @param baseline_traces Baseline period traces
#' @param treatment_traces Treatment period traces
#' @param recovery_traces Recovery period traces (optional)
#' @param positions Cell positions (optional)
#' @param frame_rate Frame rate in Hz
#'
#' @return Comprehensive analysis results
#' @export
gbm_analysis_report <- function(baseline_traces, treatment_traces,
                                 recovery_traces = NULL, positions = NULL,
                                 frame_rate = 10) {
  n_cells <- nrow(baseline_traces)

  message("Analyzing responders...")
  responders <- classify_responders(baseline_traces, treatment_traces)

  message("Analyzing baseline...")
  baseline_analysis <- analyze_baseline(baseline_traces, frame_rate = frame_rate)

  message("Computing calcium load...")
  baseline_load <- compute_calcium_load(baseline_traces, frame_rate = frame_rate)
  treatment_load <- compute_calcium_load(treatment_traces, frame_rate = frame_rate)

  message("Analyzing excitability...")
  baseline_excitability <- compute_excitability(baseline_traces, frame_rate)
  treatment_excitability <- compute_excitability(treatment_traces, frame_rate)

  message("Classifying transient shapes...")
  baseline_shapes <- classify_transient_shapes(baseline_traces, frame_rate)
  treatment_shapes <- classify_transient_shapes(treatment_traces, frame_rate)

  # Spatial analysis if positions provided
  if (!is.null(positions)) {
    message("Analyzing spatial correlations...")
    baseline_spatial <- compute_spatial_correlation(baseline_traces, positions)
    treatment_spatial <- compute_spatial_correlation(treatment_traces, positions)

    message("Detecting calcium waves...")
    baseline_waves <- detect_calcium_waves(baseline_traces, positions, frame_rate)
    treatment_waves <- detect_calcium_waves(treatment_traces, positions, frame_rate)
  } else {
    baseline_spatial <- treatment_spatial <- NULL
    baseline_waves <- treatment_waves <- NULL
  }

  # Recovery analysis if provided
  if (!is.null(recovery_traces)) {
    message("Analyzing recovery...")
    recovery <- analyze_recovery(
      cbind(baseline_traces, treatment_traces, recovery_traces),
      baseline_period = 1:ncol(baseline_traces),
      treatment_period = (ncol(baseline_traces) + 1):(ncol(baseline_traces) + ncol(treatment_traces)),
      recovery_period = (ncol(baseline_traces) + ncol(treatment_traces) + 1):
        (ncol(baseline_traces) + ncol(treatment_traces) + ncol(recovery_traces)),
      frame_rate = frame_rate
    )
  } else {
    recovery <- NULL
  }

  # Compile summary
  summary_stats <- data.frame(
    metric = c(
      "n_cells",
      "responder_fraction",
      "activated_fraction",
      "inhibited_fraction",
      "baseline_mean",
      "treatment_mean",
      "fold_change",
      "baseline_load",
      "treatment_load",
      "load_ratio",
      "baseline_event_rate",
      "treatment_event_rate"
    ),
    value = c(
      n_cells,
      responders$responder_fraction,
      responders$activated_fraction,
      responders$inhibited_fraction,
      mean(baseline_analysis$baseline),
      mean(rowMeans(treatment_traces)),
      mean(rowMeans(treatment_traces)) / mean(baseline_analysis$baseline),
      baseline_load$mean_load,
      treatment_load$mean_load,
      treatment_load$mean_load / baseline_load$mean_load,
      baseline_excitability$mean_event_rate,
      treatment_excitability$mean_event_rate
    )
  )

  structure(
    list(
      summary = summary_stats,
      responders = responders,
      baseline = baseline_analysis,
      calcium_load = list(baseline = baseline_load, treatment = treatment_load),
      excitability = list(baseline = baseline_excitability, treatment = treatment_excitability),
      transient_shapes = list(baseline = baseline_shapes, treatment = treatment_shapes),
      spatial = list(baseline = baseline_spatial, treatment = treatment_spatial),
      waves = list(baseline = baseline_waves, treatment = treatment_waves),
      recovery = recovery
    ),
    class = "gbm_analysis_report"
  )
}

#' Print GBM analysis report
#'
#' @param x GBM analysis report
#' @param ... Ignored
#'
#' @export
print.gbm_analysis_report <- function(x, ...) {
  cat("GBM Calcium Imaging Analysis Report\n")
  cat("====================================\n\n")

  cat("Cell Summary:\n")
  cat(sprintf("  Total cells: %d\n", x$summary$value[x$summary$metric == "n_cells"]))
  cat(sprintf("  Responders: %.1f%%\n", 100 * x$summary$value[x$summary$metric == "responder_fraction"]))
  cat(sprintf("    - Activated: %.1f%%\n", 100 * x$summary$value[x$summary$metric == "activated_fraction"]))
  cat(sprintf("    - Inhibited: %.1f%%\n", 100 * x$summary$value[x$summary$metric == "inhibited_fraction"]))

  cat("\nCalcium Dynamics:\n")
  cat(sprintf("  Fold change in mean: %.2f\n", x$summary$value[x$summary$metric == "fold_change"]))
  cat(sprintf("  Calcium load ratio: %.2f\n", x$summary$value[x$summary$metric == "load_ratio"]))

  cat("\nEvent Rates (per minute):\n")
  cat(sprintf("  Baseline: %.2f\n", x$summary$value[x$summary$metric == "baseline_event_rate"]))
  cat(sprintf("  Treatment: %.2f\n", x$summary$value[x$summary$metric == "treatment_event_rate"]))

  if (!is.null(x$waves$baseline) && !is.null(x$waves$treatment)) {
    cat("\nWave Analysis:\n")
    cat(sprintf("  Baseline waves: %d (%.2f/min)\n",
                x$waves$baseline$n_waves, x$waves$baseline$wave_rate))
    cat(sprintf("  Treatment waves: %d (%.2f/min)\n",
                x$waves$treatment$n_waves, x$waves$treatment$wave_rate))
  }

  if (!is.null(x$recovery)) {
    cat("\nRecovery:\n")
    cat(sprintf("  Time to 50%% recovery: %.1f sec\n", x$recovery$t50))
    cat(sprintf("  Final recovery: %.1f%%\n", 100 * x$recovery$final_recovery))
  }

  invisible(x)
}
