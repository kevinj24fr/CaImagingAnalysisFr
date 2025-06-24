#' Optimize Calcium Correction Parameters
#' 
#' Uses cross-validation to find optimal parameters for calcium correction.
#' 
#' @param raw_df Input data frame
#' @param target_metric Optimization metric: "snr", "baseline_stability", "spike_detection"
#' @param span_range Range of span values to test
#' @param cv_folds Number of cross-validation folds
#' @param verbose Whether to show progress messages
#' @return List of optimization results
#' @export
optimize_correction_parameters <- function(raw_df,
                                          target_metric = c("snr", "baseline_stability", "spike_detection"),
                                          span_range = NULL,
                                          cv_folds = NULL,
                                          verbose = TRUE) {
  
  target_metric <- match.arg(target_metric)
  config <- get_config()
  
  if (is.null(span_range)) {
    span_range <- config$optimization$span_range
  }
  
  if (is.null(cv_folds)) {
    cv_folds <- config$optimization$cv_folds
  }
  
  # Validate input
  validate_data_frame(raw_df)
  
  if (verbose) {
    message("Optimizing correction parameters for metric: ", target_metric)
  }
  
  # Get cell columns
  cell_cols <- names(raw_df)[grepl(config$cell_pattern, names(raw_df))]
  
  # Create cross-validation folds
  n_samples <- nrow(raw_df)
  fold_size <- floor(n_samples / cv_folds)
  fold_indices <- lapply(seq_len(cv_folds), function(i) {
    start_idx <- (i - 1) * fold_size + 1
    end_idx <- ifelse(i == cv_folds, n_samples, i * fold_size)
    start_idx:end_idx
  })
  
  results <- list()
  
  for (span in span_range) {
    if (verbose) {
      message("  Testing span = ", span)
    }
    
    fold_scores <- numeric(cv_folds)
    
    for (fold in seq_len(cv_folds)) {
      # Split data
      test_indices <- fold_indices[[fold]]
      train_indices <- setdiff(seq_len(n_samples), test_indices)
      
      train_data <- raw_df[train_indices, ]
      test_data <- raw_df[test_indices, ]
      
      # Apply correction
      tryCatch({
        corrected_train <- calcium_correction(train_data, span = span, verbose = FALSE)
        corrected_test <- calcium_correction(test_data, span = span, verbose = FALSE)
        
        # Calculate metric
        fold_scores[fold] <- calculate_metric(corrected_test, target_metric)
      }, error = function(e) {
        fold_scores[fold] <<- NA
      })
    }
    
    # Average score across folds
    mean_score <- mean(fold_scores, na.rm = TRUE)
    sd_score <- sd(fold_scores, na.rm = TRUE)
    
    results[[as.character(span)]] <- list(
      span = span,
      mean_score = mean_score,
      sd_score = sd_score,
      fold_scores = fold_scores
    )
  }
  
  # Find optimal parameters
  scores <- sapply(results, function(x) x$mean_score)
  optimal_span <- span_range[which.max(scores)]
  
  if (verbose) {
    message("Optimal span: ", optimal_span, " (score: ", round(max(scores, na.rm = TRUE), 4), ")")
  }
  
  list(
    optimal_span = optimal_span,
    optimal_score = max(scores, na.rm = TRUE),
    all_results = results,
    target_metric = target_metric,
    span_range = span_range
  )
}

#' Calculate optimization metric
#' 
#' @param corrected_df Corrected data frame
#' @param metric Metric to calculate
#' @return Metric value
#' @keywords internal
calculate_metric <- function(corrected_df, metric) {
  config <- get_config()
  cell_cols <- names(corrected_df)[grepl(config$cell_pattern, names(corrected_df))]
  
  if (metric == "snr") {
    # Signal-to-noise ratio
    snr_values <- sapply(cell_cols, function(col) {
      trace <- corrected_df[[col]]
      signal <- mean(trace, na.rm = TRUE)
      noise <- sd(trace, na.rm = TRUE)
      signal / noise
    })
    return(mean(snr_values, na.rm = TRUE))
    
  } else if (metric == "baseline_stability") {
    # Baseline stability (lower variance is better, so we return negative)
    stability_values <- sapply(cell_cols, function(col) {
      trace <- corrected_df[[col]]
      # Calculate baseline using moving average
      baseline <- stats::filter(trace, rep(1/50, 50), sides = 2)
      baseline[is.na(baseline)] <- trace[is.na(baseline)]
      -sd(baseline, na.rm = TRUE)  # Negative because lower is better
    })
    return(mean(stability_values, na.rm = TRUE))
    
  } else if (metric == "spike_detection") {
    # Spike detection quality (using fallback method)
    spike_scores <- sapply(cell_cols, function(col) {
      trace <- corrected_df[[col]]
      tryCatch({
        spikes <- infer_spikes(trace, method = "oasis", fallback = TRUE, verbose = FALSE)
        # Calculate spike detection score based on reasonable spike rate
        spike_rate <- sum(spikes$spike > 0, na.rm = TRUE) / length(trace)
        # Penalize too high or too low spike rates
        if (spike_rate < 0.001 || spike_rate > 0.1) {
          return(0)
        } else {
          return(1 - abs(spike_rate - 0.02) / 0.02)  # Optimal around 2%
        }
      }, error = function(e) {
        return(0)
      })
    })
    return(mean(spike_scores, na.rm = TRUE))
  }
  
  return(NA)
}

#' Optimize Spike Detection Parameters
#' 
#' @param corrected_df Corrected data frame
#' @param method Spike detection method
#' @param threshold_range Range of threshold multipliers to test
#' @param cv_folds Number of cross-validation folds
#' @param verbose Whether to show progress messages
#' @return List of optimization results
#' @export
optimize_spike_parameters <- function(corrected_df,
                                     method = "oasis",
                                     threshold_range = NULL,
                                     cv_folds = NULL,
                                     verbose = TRUE) {
  
  config <- get_config()
  
  if (is.null(threshold_range)) {
    threshold_range <- config$optimization$threshold_range
  }
  
  if (is.null(cv_folds)) {
    cv_folds <- config$optimization$cv_folds
  }
  
  if (verbose) {
    message("Optimizing spike detection parameters for method: ", method)
  }
  
  # Get cell columns
  cell_cols <- names(corrected_df)[grepl(config$cell_pattern, names(corrected_df))]
  
  # Create cross-validation folds
  n_samples <- nrow(corrected_df)
  fold_size <- floor(n_samples / cv_folds)
  fold_indices <- lapply(seq_len(cv_folds), function(i) {
    start_idx <- (i - 1) * fold_size + 1
    end_idx <- ifelse(i == cv_folds, n_samples, i * fold_size)
    start_idx:end_idx
  })
  
  results <- list()
  
  for (threshold in threshold_range) {
    if (verbose) {
      message("  Testing threshold = ", threshold)
    }
    
    fold_scores <- numeric(cv_folds)
    
    for (fold in seq_len(cv_folds)) {
      # Split data
      test_indices <- fold_indices[[fold]]
      test_data <- corrected_df[test_indices, ]
      
      # Calculate spike detection score
      tryCatch({
        cell_scores <- sapply(cell_cols, function(col) {
          trace <- test_data[[col]]
          spikes <- infer_spikes(trace, method = method, fallback = TRUE, verbose = FALSE)
          
          # Calculate quality score
          spike_rate <- sum(spikes$spike > 0, na.rm = TRUE) / length(trace)
          
          # Penalize extreme spike rates
          if (spike_rate < 0.001 || spike_rate > 0.1) {
            return(0)
          } else {
            return(1 - abs(spike_rate - 0.02) / 0.02)
          }
        })
        
        fold_scores[fold] <- mean(cell_scores, na.rm = TRUE)
      }, error = function(e) {
        fold_scores[fold] <<- NA
      })
    }
    
    # Average score across folds
    mean_score <- mean(fold_scores, na.rm = TRUE)
    sd_score <- sd(fold_scores, na.rm = TRUE)
    
    results[[as.character(threshold)]] <- list(
      threshold = threshold,
      mean_score = mean_score,
      sd_score = sd_score,
      fold_scores = fold_scores
    )
  }
  
  # Find optimal parameters
  scores <- sapply(results, function(x) x$mean_score)
  optimal_threshold <- threshold_range[which.max(scores)]
  
  if (verbose) {
    message("Optimal threshold: ", optimal_threshold, " (score: ", round(max(scores, na.rm = TRUE), 4), ")")
  }
  
  list(
    optimal_threshold = optimal_threshold,
    optimal_score = max(scores, na.rm = TRUE),
    all_results = results,
    method = method,
    threshold_range = threshold_range
  )
}

#' Auto-optimize all parameters
#' 
#' @param raw_df Input data frame
#' @param optimize_correction Whether to optimize correction parameters
#' @param optimize_spikes Whether to optimize spike detection parameters
#' @param verbose Whether to show progress messages
#' @return List of optimization results
#' @export
auto_optimize_parameters <- function(raw_df,
                                    optimize_correction = TRUE,
                                    optimize_spikes = TRUE,
                                    verbose = TRUE) {
  
  if (verbose) {
    message("Starting automatic parameter optimization...")
  }
  
  results <- list()
  
  # Optimize correction parameters
  if (optimize_correction) {
    if (verbose) message("Optimizing correction parameters...")
    
    correction_results <- list()
    for (metric in c("snr", "baseline_stability")) {
      if (verbose) message("  Testing metric: ", metric)
      correction_results[[metric]] <- optimize_correction_parameters(
        raw_df, target_metric = metric, verbose = FALSE
      )
    }
    
    # Choose best metric based on overall performance
    best_metric <- names(correction_results)[which.max(
      sapply(correction_results, function(x) x$optimal_score)
    )]
    
    results$correction <- correction_results[[best_metric]]
    results$correction$best_metric <- best_metric
  }
  
  # Apply optimal correction if available
  if (optimize_correction && !is.null(results$correction)) {
    optimal_span <- results$correction$optimal_span
    if (verbose) message("Applying optimal correction (span = ", optimal_span, ")...")
    
    corrected_df <- calcium_correction(raw_df, span = optimal_span, verbose = FALSE)
  } else {
    corrected_df <- calcium_correction(raw_df, verbose = FALSE)
  }
  
  # Optimize spike detection parameters
  if (optimize_spikes) {
    if (verbose) message("Optimizing spike detection parameters...")
    
    spike_results <- list()
    for (method in c("oasis", "caiman")) {
      if (verbose) message("  Testing method: ", method)
      tryCatch({
        spike_results[[method]] <- optimize_spike_parameters(
          corrected_df, method = method, verbose = FALSE
        )
      }, error = function(e) {
        if (verbose) message("    Method ", method, " failed: ", e$message)
      })
    }
    
    # Choose best method
    if (length(spike_results) > 0) {
      best_method <- names(spike_results)[which.max(
        sapply(spike_results, function(x) x$optimal_score)
      )]
      results$spikes <- spike_results[[best_method]]
      results$spikes$best_method <- best_method
    }
  }
  
  if (verbose) {
    message("Parameter optimization completed!")
  }
  
  return(results)
} 