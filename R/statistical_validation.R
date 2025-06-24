#' Calculate Confidence Intervals for Spike Detection
#' 
#' Uses bootstrap resampling to calculate confidence intervals for spike detection results.
#' 
#' @param trace Numeric vector of fluorescence values
#' @param method Spike inference method
#' @param confidence_level Confidence level (default: 0.95)
#' @param n_bootstrap Number of bootstrap samples
#' @param verbose Whether to show progress messages
#' @return List with confidence intervals and statistics
#' @export
spike_confidence_intervals <- function(trace,
                                      method = "oasis",
                                      confidence_level = NULL,
                                      n_bootstrap = NULL,
                                      verbose = TRUE) {
  
  config <- get_config()
  
  if (is.null(confidence_level)) {
    confidence_level <- config$statistics$confidence_level
  }
  
  if (is.null(n_bootstrap)) {
    n_bootstrap <- config$statistics$bootstrap_samples
  }
  
  validate_trace(trace)
  
  if (verbose) {
    message("Calculating confidence intervals for spike detection...")
  }
  
  # Original spike detection
  original_spikes <- tryCatch({
    infer_spikes(trace, method = method, fallback = TRUE, verbose = FALSE)
  }, error = function(e) {
    warning("Original spike detection failed: ", e$message)
    return(data.frame(fit = rep(0, length(trace)), spike = rep(0, length(trace))))
  })
  
  # Bootstrap samples
  bootstrap_results <- list()
  spike_counts <- numeric(n_bootstrap)
  spike_rates <- numeric(n_bootstrap)
  
  for (i in seq_len(n_bootstrap)) {
    if (verbose && i %% 100 == 0) {
      message("  Bootstrap sample ", i, "/", n_bootstrap)
    }
    
    # Resample with replacement
    indices <- sample(seq_along(trace), replace = TRUE)
    bootstrap_trace <- trace[indices]
    
    # Detect spikes in bootstrap sample
    tryCatch({
      bootstrap_spikes <- infer_spikes(bootstrap_trace, method = method, fallback = TRUE, verbose = FALSE)
      
      spike_counts[i] <- sum(bootstrap_spikes$spike > 0, na.rm = TRUE)
      spike_rates[i] <- spike_counts[i] / length(bootstrap_trace)
      
      bootstrap_results[[i]] <- list(
        spike_count = spike_counts[i],
        spike_rate = spike_rates[i],
        trace_length = length(bootstrap_trace)
      )
    }, error = function(e) {
      spike_counts[i] <<- 0
      spike_rates[i] <<- 0
    })
  }
  
  # Calculate confidence intervals
  alpha <- 1 - confidence_level
  lower_quantile <- alpha / 2
  upper_quantile <- 1 - alpha / 2
  
  ci_spike_count <- quantile(spike_counts, c(lower_quantile, upper_quantile), na.rm = TRUE)
  ci_spike_rate <- quantile(spike_rates, c(lower_quantile, upper_quantile), na.rm = TRUE)
  
  # Original statistics
  original_spike_count <- sum(original_spikes$spike > 0, na.rm = TRUE)
  original_spike_rate <- original_spike_count / length(trace)
  
  if (verbose) {
    message("Confidence intervals calculated successfully")
  }
  
  list(
    original = list(
      spike_count = original_spike_count,
      spike_rate = original_spike_rate
    ),
    confidence_intervals = list(
      spike_count = ci_spike_count,
      spike_rate = ci_spike_rate
    ),
    bootstrap_samples = n_bootstrap,
    confidence_level = confidence_level,
    bootstrap_results = bootstrap_results
  )
}

#' Calculate Quality Metrics for Calcium Correction
#' 
#' @param raw_df Raw data frame
#' @param corrected_df Corrected data frame
#' @param verbose Whether to show progress messages
#' @return List of quality metrics
#' @export
calculate_correction_quality <- function(raw_df, corrected_df, verbose = TRUE) {
  
  config <- get_config()
  
  if (verbose) {
    message("Calculating correction quality metrics...")
  }
  
  # Get cell columns
  cell_cols <- names(raw_df)[grepl(config$cell_pattern, names(raw_df))]
  
  quality_metrics <- list()
  
  for (col in cell_cols) {
    raw_trace <- raw_df[[col]]
    corrected_trace <- corrected_df[[col]]
    
    # Signal-to-noise ratio improvement
    raw_snr <- mean(raw_trace, na.rm = TRUE) / sd(raw_trace, na.rm = TRUE)
    corrected_snr <- mean(corrected_trace, na.rm = TRUE) / sd(corrected_trace, na.rm = TRUE)
    snr_improvement <- (corrected_snr - raw_snr) / raw_snr * 100
    
    # Baseline stability
    baseline <- stats::filter(corrected_trace, rep(1/50, 50), sides = 2)
    baseline[is.na(baseline)] <- corrected_trace[is.na(baseline)]
    baseline_stability <- sd(baseline, na.rm = TRUE)
    
    # Outlier detection
    outliers <- detect_outliers(corrected_trace)
    outlier_fraction <- sum(outliers) / length(corrected_trace)
    
    # Dynamic range
    dynamic_range <- max(corrected_trace, na.rm = TRUE) - min(corrected_trace, na.rm = TRUE)
    
    quality_metrics[[col]] <- list(
      snr_improvement = snr_improvement,
      baseline_stability = baseline_stability,
      outlier_fraction = outlier_fraction,
      dynamic_range = dynamic_range,
      raw_snr = raw_snr,
      corrected_snr = corrected_snr
    )
  }
  
  # Overall quality score
  overall_scores <- sapply(quality_metrics, function(x) {
    # Combine metrics into overall score (0-1, higher is better)
    snr_score <- min(x$snr_improvement / 50, 1)  # Cap at 50% improvement
    stability_score <- max(0, 1 - x$baseline_stability / 0.1)  # Lower stability is better
    outlier_score <- max(0, 1 - x$outlier_fraction / 0.05)  # Lower outliers is better
    
    mean(c(snr_score, stability_score, outlier_score))
  })
  
  overall_quality <- mean(overall_scores, na.rm = TRUE)
  
  if (verbose) {
    message("Quality metrics calculated successfully")
  }
  
  list(
    cell_metrics = quality_metrics,
    overall_quality = overall_quality,
    quality_threshold = config$statistics$quality_threshold,
    passes_threshold = overall_quality >= config$statistics$quality_threshold
  )
}

#' Detect outliers in a trace
#' 
#' @param trace Numeric vector
#' @param method Outlier detection method
#' @param threshold Outlier threshold
#' @return Logical vector indicating outliers
#' @keywords internal
detect_outliers <- function(trace, method = NULL, threshold = NULL) {
  config <- get_config()
  
  if (is.null(method)) {
    method <- config$statistics$outlier_method
  }
  
  if (is.null(threshold)) {
    threshold <- config$statistics$outlier_threshold
  }
  
  if (method == "iqr") {
    q1 <- quantile(trace, 0.25, na.rm = TRUE)
    q3 <- quantile(trace, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower_bound <- q1 - threshold * iqr
    upper_bound <- q3 + threshold * iqr
    outliers <- trace < lower_bound | trace > upper_bound
    
  } else if (method == "zscore") {
    mean_val <- mean(trace, na.rm = TRUE)
    sd_val <- sd(trace, na.rm = TRUE)
    z_scores <- abs((trace - mean_val) / sd_val)
    outliers <- z_scores > threshold
    
  } else if (method == "mad") {
    median_val <- median(trace, na.rm = TRUE)
    mad_val <- median(abs(trace - median_val), na.rm = TRUE)
    modified_z_scores <- abs((trace - median_val) / mad_val)
    outliers <- modified_z_scores > threshold
  }
  
  outliers[is.na(outliers)] <- FALSE
  return(outliers)
}

#' Statistical validation of spike detection results
#' 
#' @param trace Numeric vector
#' @param spikes Spike detection results
#' @param method Detection method used
#' @param verbose Whether to show progress messages
#' @return List of validation results
#' @export
validate_spike_detection <- function(trace, spikes, method = "unknown", verbose = TRUE) {
  
  if (verbose) {
    message("Validating spike detection results...")
  }
  
  # Basic statistics
  n_spikes <- sum(spikes$spike > 0, na.rm = TRUE)
  spike_rate <- n_spikes / length(trace)
  
  # Spike interval analysis
  spike_times <- which(spikes$spike > 0)
  intervals <- diff(spike_times)
  
  # Minimum interval check
  config <- get_config()
  min_interval <- config$spike_detection$min_spike_interval
  
  if (length(intervals) > 0) {
    violations <- sum(intervals < min_interval)
    interval_violation_rate <- violations / length(intervals)
  } else {
    interval_violation_rate <- 0
  }
  
  # Amplitude analysis
  spike_amplitudes <- trace[spikes$spike > 0]
  min_amplitude <- config$spike_detection$amplitude_threshold
  
  if (length(spike_amplitudes) > 0) {
    weak_spikes <- sum(spike_amplitudes < min_amplitude)
    weak_spike_rate <- weak_spikes / length(spike_amplitudes)
  } else {
    weak_spike_rate <- 0
  }
  
  # Quality score calculation
  quality_score <- calculate_spike_quality_score(
    spike_rate, interval_violation_rate, weak_spike_rate
  )
  
  # Confidence intervals
  ci_results <- tryCatch({
    spike_confidence_intervals(trace, method = method, verbose = FALSE)
  }, error = function(e) {
    list(confidence_intervals = list(spike_rate = c(NA, NA)))
  })
  
  if (verbose) {
    message("Spike detection validation completed")
  }
  
  list(
    basic_stats = list(
      n_spikes = n_spikes,
      spike_rate = spike_rate,
      trace_length = length(trace)
    ),
    interval_analysis = list(
      mean_interval = ifelse(length(intervals) > 0, mean(intervals, na.rm = TRUE), NA),
      min_interval = ifelse(length(intervals) > 0, min(intervals, na.rm = TRUE), NA),
      violations = ifelse(length(intervals) > 0, violations, 0),
      violation_rate = interval_violation_rate
    ),
    amplitude_analysis = list(
      mean_amplitude = ifelse(length(spike_amplitudes) > 0, mean(spike_amplitudes, na.rm = TRUE), NA),
      weak_spikes = ifelse(length(spike_amplitudes) > 0, weak_spikes, 0),
      weak_spike_rate = weak_spike_rate
    ),
    quality_score = quality_score,
    confidence_intervals = ci_results$confidence_intervals,
    method = method,
    passes_validation = quality_score >= config$statistics$quality_threshold
  )
}

#' Calculate spike detection quality score
#' 
#' @param spike_rate Spike rate
#' @param interval_violation_rate Rate of interval violations
#' @param weak_spike_rate Rate of weak spikes
#' @return Quality score (0-1, higher is better)
#' @keywords internal
calculate_spike_quality_score <- function(spike_rate, interval_violation_rate, weak_spike_rate) {
  
  # Spike rate score (optimal around 2%)
  optimal_rate <- 0.02
  rate_score <- max(0, 1 - abs(spike_rate - optimal_rate) / optimal_rate)
  
  # Interval violation score (lower is better)
  interval_score <- max(0, 1 - interval_violation_rate)
  
  # Weak spike score (lower is better)
  weak_spike_score <- max(0, 1 - weak_spike_rate)
  
  # Combine scores (weighted average)
  overall_score <- 0.4 * rate_score + 0.3 * interval_score + 0.3 * weak_spike_score
  
  return(overall_score)
}

#' Comprehensive statistical analysis
#' 
#' @param raw_df Raw data frame
#' @param corrected_df Corrected data frame
#' @param spike_results List of spike detection results
#' @param verbose Whether to show progress messages
#' @return Comprehensive statistical analysis
#' @export
comprehensive_statistical_analysis <- function(raw_df, corrected_df, spike_results, verbose = TRUE) {
  
  if (verbose) {
    message("Performing comprehensive statistical analysis...")
  }
  
  # Correction quality
  correction_quality <- calculate_correction_quality(raw_df, corrected_df, verbose = FALSE)
  
  # Spike detection validation
  config <- get_config()
  cell_cols <- names(corrected_df)[grepl(config$cell_pattern, names(corrected_df))]
  
  spike_validations <- list()
  for (cell in cell_cols) {
    if (cell %in% names(spike_results)) {
      trace <- corrected_df[[cell]]
      spikes <- spike_results[[cell]]$spikes
      method <- spike_results[[cell]]$method
      
      spike_validations[[cell]] <- validate_spike_detection(
        trace, spikes, method, verbose = FALSE
      )
    }
  }
  
  # Overall assessment
  overall_correction_quality <- correction_quality$overall_quality
  overall_spike_quality <- mean(
    sapply(spike_validations, function(x) x$quality_score), 
    na.rm = TRUE
  )
  
  overall_quality <- (overall_correction_quality + overall_spike_quality) / 2
  
  if (verbose) {
    message("Statistical analysis completed")
  }
  
  list(
    correction_quality = correction_quality,
    spike_validations = spike_validations,
    overall_quality = overall_quality,
    passes_quality_threshold = overall_quality >= config$statistics$quality_threshold,
    recommendations = generate_recommendations(correction_quality, spike_validations)
  )
}

#' Generate recommendations based on analysis
#' 
#' @param correction_quality Correction quality results
#' @param spike_validations Spike validation results
#' @return List of recommendations
#' @keywords internal
generate_recommendations <- function(correction_quality, spike_validations) {
  
  recommendations <- list()
  
  # Correction recommendations
  if (correction_quality$overall_quality < 0.6) {
    recommendations$correction <- "Consider adjusting correction parameters or using different method"
  }
  
  # Spike detection recommendations
  spike_issues <- sapply(spike_validations, function(x) {
    issues <- c()
    if (x$interval_analysis$violation_rate > 0.1) {
      issues <- c(issues, "High interval violations")
    }
    if (x$amplitude_analysis$weak_spike_rate > 0.3) {
      issues <- c(issues, "Many weak spikes")
    }
    if (x$basic_stats$spike_rate < 0.001 || x$basic_stats$spike_rate > 0.1) {
      issues <- c(issues, "Unusual spike rate")
    }
    issues
  })
  
  if (any(sapply(spike_issues, length) > 0)) {
    recommendations$spike_detection <- "Consider adjusting spike detection parameters or method"
  }
  
  if (length(recommendations) == 0) {
    recommendations$general <- "Analysis quality is good, no major issues detected"
  }
  
  return(recommendations)
} 