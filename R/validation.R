#' Validate input data frame
#' 
#' @param raw_df Data frame to validate
#' @param require_cells Whether to require cell columns
#' @param require_background Whether to require background columns
#' @return TRUE if valid, throws error otherwise
#' @keywords internal
validate_data_frame <- function(raw_df, require_cells = TRUE, require_background = TRUE) {
  config <- get_config()
  
  if (!is.data.frame(raw_df)) {
    stop("raw_df must be a data frame")
  }
  
  if (nrow(raw_df) == 0) {
    stop("raw_df cannot be empty")
  }
  
  if (ncol(raw_df) == 0) {
    stop("raw_df must have at least one column")
  }
  
  # Check for required column patterns
  cell_cols <- names(raw_df)[grepl(config$cell_pattern, names(raw_df))]
  bg_cols <- names(raw_df)[grepl(config$background_pattern, names(raw_df))]
  
  if (require_cells && length(cell_cols) == 0) {
    stop("No cell columns (", config$cell_pattern, ") found in data")
  }
  
  if (require_background && length(bg_cols) == 0) {
    stop("No background columns (", config$background_pattern, ") found in data")
  }
  
  # Check for numeric columns
  all_cols <- c(cell_cols, bg_cols)
  non_numeric <- sapply(raw_df[all_cols], function(x) !is.numeric(x))
  
  if (any(non_numeric)) {
    stop("All cell and background columns must be numeric")
  }
  
  # Check for missing values
  missing_vals <- sapply(raw_df[all_cols], function(x) any(is.na(x)))
  if (any(missing_vals)) {
    warning("Missing values detected in columns: ", 
            paste(names(missing_vals)[missing_vals], collapse = ", "))
  }
  
  # Quality control checks
  validate_quality_control(raw_df, cell_cols, bg_cols)
  
  TRUE
}

#' Validate quality control parameters
#' 
#' @param raw_df Input data frame
#' @param cell_cols Cell column names
#' @param bg_cols Background column names
#' @return TRUE if passes quality control
#' @keywords internal
validate_quality_control <- function(raw_df, cell_cols, bg_cols) {
  config <- get_config()
  qc <- config$quality_control
  
  # Check trace length
  if (nrow(raw_df) < qc$min_trace_length) {
    warning("Trace length (", nrow(raw_df), ") is below minimum recommended (", 
            qc$min_trace_length, ")")
  }
  
  # Check missing data
  for (col in c(cell_cols, bg_cols)) {
    missing_fraction <- sum(is.na(raw_df[[col]])) / nrow(raw_df)
    if (missing_fraction > qc$max_missing_fraction) {
      warning("Column ", col, " has ", round(missing_fraction * 100, 1), 
              "% missing data (above ", qc$max_missing_fraction * 100, "%)")
    }
    
    # Check for contiguous missing data
    if (missing_fraction > 0) {
      rle_missing <- rle(is.na(raw_df[[col]]))
      max_contiguous <- max(rle_missing$lengths[rle_missing$values])
      if (max_contiguous > qc$max_contiguous_missing) {
        warning("Column ", col, " has ", max_contiguous, 
                " contiguous missing values (above ", qc$max_contiguous_missing, ")")
      }
    }
  }
  
  # Check signal-to-noise ratio for cell traces
  for (col in cell_cols) {
    trace <- raw_df[[col]]
    if (all(is.na(trace))) next
    
    signal <- mean(trace, na.rm = TRUE)
    noise <- sd(trace, na.rm = TRUE)
    snr <- signal / noise
    
    if (snr < qc$min_signal_to_noise) {
      warning("Column ", col, " has low signal-to-noise ratio: ", 
              round(snr, 2), " (below ", qc$min_signal_to_noise, ")")
    }
  }
  
  TRUE
}

#' Validate numeric parameters
#' 
#' @param value Parameter value
#' @param min_val Minimum allowed value
#' @param max_val Maximum allowed value
#' @param param_name Parameter name for error messages
#' @return TRUE if valid, throws error otherwise
#' @keywords internal
validate_numeric_param <- function(value, min_val, max_val, param_name) {
  if (!is.numeric(value) || length(value) != 1) {
    stop(param_name, " must be a single numeric value")
  }
  
  if (value < min_val || value > max_val) {
    stop(param_name, " must be between ", min_val, " and ", max_val)
  }
  
  TRUE
}

#' Validate trace vector
#' 
#' @param trace Numeric vector to validate
#' @param param_name Parameter name for error messages
#' @return TRUE if valid, throws error otherwise
#' @keywords internal
validate_trace <- function(trace, param_name = "trace") {
  if (!is.numeric(trace)) {
    stop(param_name, " must be a numeric vector")
  }
  
  if (length(trace) == 0) {
    stop(param_name, " cannot be empty")
  }
  
  if (any(is.infinite(trace))) {
    stop(param_name, " contains infinite values")
  }
  
  # Check for reasonable fluorescence values (typically positive)
  if (any(trace < -1000, na.rm = TRUE)) {
    warning(param_name, " contains unusually low values (< -1000)")
  }
  
  if (any(trace > 10000, na.rm = TRUE)) {
    warning(param_name, " contains unusually high values (> 10000)")
  }
  
  TRUE
}

#' Validate temporal consistency
#' 
#' @param raw_df Input data frame
#' @return TRUE if temporally consistent
#' @keywords internal
validate_temporal_consistency <- function(raw_df) {
  config <- get_config()
  
  # Check for reasonable sampling rate (no extreme jumps)
  cell_cols <- names(raw_df)[grepl(config$cell_pattern, names(raw_df))]
  
  for (col in cell_cols) {
    trace <- raw_df[[col]]
    if (length(trace) < 2) next
    
    # Calculate first differences
    diff_trace <- diff(trace)
    
    # Check for extreme jumps (outliers)
    qc <- config$quality_control
    outlier_threshold <- config$statistics$outlier_threshold
    
    if (config$statistics$outlier_method == "iqr") {
      q1 <- quantile(diff_trace, 0.25, na.rm = TRUE)
      q3 <- quantile(diff_trace, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      outliers <- abs(diff_trace) > (q3 + outlier_threshold * iqr)
    } else if (config$statistics$outlier_method == "zscore") {
      mean_diff <- mean(diff_trace, na.rm = TRUE)
      sd_diff <- sd(diff_trace, na.rm = TRUE)
      outliers <- abs((diff_trace - mean_diff) / sd_diff) > outlier_threshold
    }
    
    if (sum(outliers, na.rm = TRUE) > 0) {
      warning("Column ", col, " contains ", sum(outliers, na.rm = TRUE), 
              " temporal outliers (extreme jumps)")
    }
  }
  
  TRUE
}

#' Check memory usage
#' 
#' @param data_size Estimated data size in MB
#' @return TRUE if within memory limits
#' @keywords internal
check_memory_usage <- function(data_size) {
  config <- get_config()
  
  if (data_size > config$memory_limit_mb) {
    warning("Estimated data size (", round(data_size, 1), " MB) exceeds memory limit (", 
            config$memory_limit_mb, " MB). Consider processing in chunks.")
    return(FALSE)
  }
  
  TRUE
} 