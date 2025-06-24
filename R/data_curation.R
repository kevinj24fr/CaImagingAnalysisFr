#' Data Curation Functions
#'
#' Data curation and validation utilities for calcium imaging data.
#'
#' @name data_curation
NULL

#' Data Validation and Quality Check
#'
#' Validate calcium imaging data and check for common issues.
#'
#' @param data Input data (matrix, data frame, or file path)
#' @param data_type Type of data ("calcium_traces", "spike_times", "metadata")
#' @param validation_rules List of validation rules to apply
#' @param output_format Output format for validation report ("summary", "detailed", "json")
#' @param ... Additional arguments
#' @return Validation results and recommendations
#' @export
validate_calcium_data <- function(data, data_type = "calcium_traces", validation_rules = NULL, output_format = "summary", ...) {
  # Load data if file path provided
  if (is.character(data) && file.exists(data)) {
    if (grepl("\\.csv$", data)) {
      data <- read.csv(data, stringsAsFactors = FALSE)
    } else if (grepl("\\.rds$", data)) {
      data <- readRDS(data)
    } else {
      stop("Unsupported file format. Please provide CSV or RDS file.")
    }
  }
  
  # Default validation rules
  if (is.null(validation_rules)) {
    validation_rules <- list(
      check_missing_values = TRUE,
      check_outliers = TRUE,
      check_data_types = TRUE,
      check_dimensions = TRUE,
      check_value_ranges = TRUE
    )
  }
  
  validation_results <- list()
  issues <- list()
  recommendations <- list()
  
  # Check data type and apply appropriate validations
  if (data_type == "calcium_traces") {
    # Validate calcium trace data
    if (validation_rules$check_dimensions) {
      if (!is.matrix(data) && !is.data.frame(data)) {
        issues$dimensions <- "Data should be a matrix or data frame"
        recommendations$dimensions <- "Convert data to matrix format"
      } else {
        validation_results$dimensions <- paste("Data dimensions:", nrow(data), "x", ncol(data))
      }
    }
    
    if (validation_rules$check_missing_values) {
      missing_count <- sum(is.na(data))
      missing_percent <- (missing_count / length(data)) * 100
      
      if (missing_count > 0) {
        issues$missing_values <- paste("Found", missing_count, "missing values (", round(missing_percent, 2), "%)")
        recommendations$missing_values <- "Consider imputation or removal of missing values"
      } else {
        validation_results$missing_values <- "No missing values found"
      }
    }
    
    if (validation_rules$check_outliers) {
      # Simple outlier detection using IQR method
      outliers <- list()
      for (i in 1:ncol(data)) {
        col_data <- data[, i]
        q1 <- quantile(col_data, 0.25, na.rm = TRUE)
        q3 <- quantile(col_data, 0.75, na.rm = TRUE)
        iqr <- q3 - q1
        lower_bound <- q1 - 1.5 * iqr
        upper_bound <- q3 + 1.5 * iqr
        
        outlier_count <- sum(col_data < lower_bound | col_data > upper_bound, na.rm = TRUE)
        if (outlier_count > 0) {
          outliers[[paste0("Column_", i)]] <- outlier_count
        }
      }
      
      if (length(outliers) > 0) {
        issues$outliers <- paste("Found outliers in", length(outliers), "columns")
        recommendations$outliers <- "Review and handle outliers appropriately"
      } else {
        validation_results$outliers <- "No significant outliers detected"
      }
    }
    
    if (validation_rules$check_value_ranges) {
      # Check for reasonable calcium signal ranges
      min_val <- min(data, na.rm = TRUE)
      max_val <- max(data, na.rm = TRUE)
      
      if (min_val < -10 || max_val > 10) {
        issues$value_ranges <- paste("Values outside typical range: [", min_val, ",", max_val, "]")
        recommendations$value_ranges <- "Check if data is properly normalized (Delta F/F)"
      } else {
        validation_results$value_ranges <- paste("Value range: [", round(min_val, 3), ",", round(max_val, 3), "]")
      }
    }
    
  } else if (data_type == "metadata") {
    # Validate metadata
    required_fields <- c("experiment_id", "date", "animal_id", "session_id")
    missing_fields <- setdiff(required_fields, names(data))
    
    if (length(missing_fields) > 0) {
      issues$missing_fields <- paste("Missing required fields:", paste(missing_fields, collapse = ", "))
      recommendations$missing_fields <- "Add missing metadata fields"
    } else {
      validation_results$required_fields <- "All required metadata fields present"
    }
  }
  
  # Compute quality_score (simple: 1 if no issues, 0 if any issues)
  quality_score <- if (length(issues) == 0) 1 else 0

  # Compute basic_stats (mean, sd, min, max for each column if data is matrix/data.frame)
  basic_stats <- NULL
  if (is.matrix(data) || is.data.frame(data)) {
    basic_stats <- data.frame(
      column = colnames(data),
      mean = apply(data, 2, function(x) mean(x, na.rm = TRUE)),
      sd = apply(data, 2, function(x) sd(x, na.rm = TRUE)),
      min = apply(data, 2, function(x) min(x, na.rm = TRUE)),
      max = apply(data, 2, function(x) max(x, na.rm = TRUE))
    )
  } else {
    basic_stats <- data.frame()
  }

  # Create summary
  summary <- list(
    data_type = data_type,
    validation_passed = length(issues) == 0,
    n_issues = length(issues),
    issues = issues,
    recommendations = recommendations,
    validation_results = validation_results,
    quality_score = quality_score,
    basic_stats = basic_stats
  )
  
  # Format output
  if (output_format == "summary") {
    return(summary)
  } else if (output_format == "detailed") {
    return(list(
      summary = summary,
      raw_data = data,
      validation_rules = validation_rules
    ))
  } else if (output_format == "json") {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("jsonlite is required for JSON output")
    }
    return(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE))
  }
  
  return(summary)
}

#' Metadata Management System
#'
#' Create and manage metadata for calcium imaging experiments.
#'
#' @param experiment_info List containing experiment information
#' @param metadata_template Template for metadata structure (default: NULL)
#' @param auto_generate Whether to auto-generate missing fields (default: TRUE)
#' @param ... Additional arguments
#' @return Structured metadata object
#' @export
create_experiment_metadata <- function(experiment_info, metadata_template = NULL, auto_generate = TRUE, ...) {
  # Default metadata template
  if (is.null(metadata_template)) {
    metadata_template <- list(
      experiment_id = NA_character_,
      date = NA_character_,
      time = NA_character_,
      animal_id = NA_character_,
      session_id = NA_character_,
      experimenter = NA_character_,
      protocol = NA_character_,
      imaging_parameters = list(
        frame_rate = NA_real_,
        resolution = NA_character_,
        magnification = NA_real_,
        excitation_wavelength = NA_real_,
        emission_wavelength = NA_real_
      ),
      data_processing = list(
        preprocessing_steps = character(),
        spike_detection_method = NA_character_,
        correction_method = NA_character_
      ),
      notes = NA_character_
    )
  }
  
  # Create metadata object
  metadata <- metadata_template
  
  # Fill in provided information
  for (field in names(experiment_info)) {
    if (field %in% names(metadata)) {
      metadata[[field]] <- experiment_info[[field]]
    }
  }
  
  # Auto-generate missing fields
  if (auto_generate) {
    if (is.na(metadata$date)) {
      metadata$date <- as.character(Sys.Date())
    }
    
    if (is.na(metadata$time)) {
      metadata$time <- format(Sys.time(), "%H:%M:%S")
    }
    
    if (is.na(metadata$experiment_id) && !is.na(metadata$animal_id)) {
      metadata$experiment_id <- paste0(metadata$animal_id, "_", metadata$date)
    }
    
    if (is.na(metadata$session_id)) {
      metadata$session_id <- paste0("session_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    }
  }
  
  # Add metadata class
  class(metadata) <- "calcium_imaging_metadata"
  
  return(metadata)
}

#' Automated Data Preprocessing Pipeline
#'
#' Automatically preprocess calcium imaging data based on data characteristics.
#'
#' @param data Input calcium imaging data
#' @param preprocessing_steps Vector of preprocessing steps to apply
#' @param parameters List of preprocessing parameters
#' @param quality_threshold Quality threshold for automatic decisions (default: 0.8)
#' @param ... Additional arguments
#' @return Preprocessed data and preprocessing report
#' @export
automated_preprocessing <- function(data, preprocessing_steps = NULL, parameters = NULL, quality_threshold = 0.8, ...) {
  # Default preprocessing steps
  if (is.null(preprocessing_steps)) {
    preprocessing_steps <- c("detect_outliers", "remove_trend", "normalize", "smooth")
  }
  
  # Default parameters
  if (is.null(parameters)) {
    parameters <- list(
      outlier_method = "iqr",
      trend_method = "linear",
      normalization_method = "zscore",
      smoothing_method = "gaussian",
      smoothing_window = 5
    )
  }
  
  preprocessing_report <- list(
    steps_applied = character(),
    parameters_used = parameters,
    quality_metrics = list(),
    warnings = character()
  )
  
  processed_data <- data
  
  # Apply preprocessing steps
  for (step in preprocessing_steps) {
    if (step == "detect_outliers") {
      result <- detect_and_handle_outliers(processed_data, method = parameters$outlier_method)
      processed_data <- result$data
      preprocessing_report$steps_applied <- c(preprocessing_report$steps_applied, "outlier_detection")
      preprocessing_report$quality_metrics$outliers_removed <- result$outliers_removed
      
    } else if (step == "remove_trend") {
      result <- remove_trend(processed_data, method = parameters$trend_method)
      processed_data <- result$data
      preprocessing_report$steps_applied <- c(preprocessing_report$steps_applied, "trend_removal")
      preprocessing_report$quality_metrics$trend_removed <- result$trend_removed
      
    } else if (step == "normalize") {
      result <- normalize_data(processed_data, method = parameters$normalization_method)
      processed_data <- result$data
      preprocessing_report$steps_applied <- c(preprocessing_report$steps_applied, "normalization")
      preprocessing_report$quality_metrics$normalization_method <- result$method
      
    } else if (step == "smooth") {
      result <- smooth_data(processed_data, method = parameters$smoothing_method, window = parameters$smoothing_window)
      processed_data <- result$data
      preprocessing_report$steps_applied <- c(preprocessing_report$steps_applied, "smoothing")
      preprocessing_report$quality_metrics$smoothing_applied <- TRUE
    }
  }
  
  # Calculate overall quality score
  quality_score <- calculate_preprocessing_quality(processed_data, data)
  preprocessing_report$quality_metrics$overall_quality <- quality_score
  
  if (quality_score < quality_threshold) {
    preprocessing_report$warnings <- c(preprocessing_report$warnings, 
                                     "Quality score below threshold - review preprocessing steps")
  }
  
  return(list(
    processed_data = processed_data,
    preprocessing_report = preprocessing_report,
    original_data = data
  ))
}

#' Data Format Standardization
#'
#' Standardize calcium imaging data to a common format.
#'
#' @param data Input data in various formats
#' @param target_format Target format ("matrix", "data.frame", "list")
#' @param metadata Metadata to include (optional)
#' @param ... Additional arguments
#' @return Standardized data object
#' @export
standardize_data_format <- function(data, target_format = "matrix", metadata = NULL, ...) {
  # Detect input format
  input_format <- class(data)[1]
  
  # Convert to target format
  if (target_format == "matrix") {
    if (input_format == "data.frame") {
      standardized_data <- as.matrix(data)
    } else if (input_format == "list") {
      # Assume list contains traces
      standardized_data <- do.call(cbind, data)
    } else if (input_format == "matrix") {
      standardized_data <- data
    } else {
      stop("Unsupported input format for matrix conversion")
    }
    
    # Add column names if missing
    if (is.null(colnames(standardized_data))) {
      colnames(standardized_data) <- paste0("Cell_", 1:ncol(standardized_data))
    }
    
  } else if (target_format == "data.frame") {
    if (input_format == "matrix") {
      standardized_data <- as.data.frame(data)
    } else if (input_format == "list") {
      standardized_data <- as.data.frame(do.call(cbind, data))
    } else if (input_format == "data.frame") {
      standardized_data <- data
    } else {
      stop("Unsupported input format for data.frame conversion")
    }
    
  } else if (target_format == "list") {
    if (input_format == "matrix" || input_format == "data.frame") {
      standardized_data <- as.list(as.data.frame(data))
    } else if (input_format == "list") {
      standardized_data <- data
    } else {
      stop("Unsupported input format for list conversion")
    }
  }
  
  # Add metadata if provided
  if (!is.null(metadata)) {
    attr(standardized_data, "metadata") <- metadata
  }
  
  # Add format information
  attr(standardized_data, "format_info") <- list(
    original_format = input_format,
    target_format = target_format,
    conversion_date = Sys.time()
  )
  
  return(standardized_data)
}

#' Data Export and Archiving
#'
#' Export processed data and metadata in various formats.
#'
#' @param data Processed data object
#' @param metadata Metadata object
#' @param output_dir Output directory (default: "exports")
#' @param formats Vector of export formats ("csv", "rds", "json", "mat")
#' @param include_metadata Whether to include metadata in export (default: TRUE)
#' @param compression Whether to compress files (default: TRUE)
#' @param ... Additional arguments
#' @return List of exported file paths
#' @export
export_processed_data <- function(data, metadata = NULL, output_dir = "exports", formats = c("csv", "rds"), include_metadata = TRUE, compression = TRUE, ...) {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Generate timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  exported_files <- list()
  
  # Export in requested formats
  for (format in formats) {
    if (format == "csv") {
      filename <- file.path(output_dir, paste0("calcium_data_", timestamp, ".csv"))
      write.csv(data, filename, row.names = FALSE)
      exported_files$csv <- filename
      
    } else if (format == "rds") {
      filename <- file.path(output_dir, paste0("calcium_data_", timestamp, ".rds"))
      if (include_metadata && !is.null(metadata)) {
        save_data <- list(data = data, metadata = metadata)
      } else {
        save_data <- data
      }
      saveRDS(save_data, filename, compress = compression)
      exported_files$rds <- filename
      
    } else if (format == "json") {
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        warning("jsonlite not available, skipping JSON export")
        next
      }
      filename <- file.path(output_dir, paste0("calcium_data_", timestamp, ".json"))
      export_data <- list(data = as.data.frame(data))
      if (include_metadata && !is.null(metadata)) {
        export_data$metadata <- metadata
      }
      jsonlite::write_json(export_data, filename, pretty = TRUE)
      exported_files$json <- filename
      
    } else if (format == "mat") {
      if (!requireNamespace("R.matlab", quietly = TRUE)) {
        warning("R.matlab not available, skipping MAT export")
        next
      }
      filename <- file.path(output_dir, paste0("calcium_data_", timestamp, ".mat"))
      export_data <- list(calcium_data = data)
      if (include_metadata && !is.null(metadata)) {
        export_data$metadata <- metadata
      }
      R.matlab::writeMat(filename, ... = export_data)
      exported_files$mat <- filename
    }
  }
  
  # Create export summary
  summary_file <- file.path(output_dir, paste0("export_summary_", timestamp, ".txt"))
  summary_content <- c(
    "=== Calcium Imaging Data Export Summary ===",
    "",
    paste("Export Date:", Sys.time()),
    paste("Output Directory:", output_dir),
    "",
    "Exported Files:",
    paste(names(exported_files), unlist(exported_files), sep = ": ")
  )
  writeLines(summary_content, summary_file)
  exported_files$summary <- summary_file
  
  message("Data exported successfully to: ", output_dir)
  return(exported_files)
}

#' Helper Functions for Preprocessing
#' @keywords internal

detect_and_handle_outliers <- function(data, method = "iqr") {
  outliers_removed <- 0
  
  if (method == "iqr") {
    for (i in 1:ncol(data)) {
      col_data <- data[, i]
      q1 <- quantile(col_data, 0.25, na.rm = TRUE)
      q3 <- quantile(col_data, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      lower_bound <- q1 - 1.5 * iqr
      upper_bound <- q3 + 1.5 * iqr
      
      outlier_indices <- which(col_data < lower_bound | col_data > upper_bound)
      if (length(outlier_indices) > 0) {
        data[outlier_indices, i] <- NA
        outliers_removed <- outliers_removed + length(outlier_indices)
      }
    }
  }
  
  return(list(data = data, outliers_removed = outliers_removed))
}

remove_trend <- function(data, method = "linear") {
  trend_removed <- FALSE
  
  if (method == "linear") {
    for (i in 1:ncol(data)) {
      col_data <- data[, i]
      time_points <- 1:length(col_data)
      
      # Fit linear trend
      trend_model <- lm(col_data ~ time_points, na.action = na.exclude)
      trend <- predict(trend_model, newdata = data.frame(time_points = time_points))
      
      # Remove trend
      data[, i] <- col_data - trend
      trend_removed <- TRUE
    }
  }
  
  return(list(data = data, trend_removed = trend_removed))
}

normalize_data <- function(data, method = "zscore") {
  if (method == "zscore") {
    for (i in 1:ncol(data)) {
      col_data <- data[, i]
      mean_val <- mean(col_data, na.rm = TRUE)
      sd_val <- sd(col_data, na.rm = TRUE)
      data[, i] <- (col_data - mean_val) / sd_val
    }
  } else if (method == "minmax") {
    for (i in 1:ncol(data)) {
      col_data <- data[, i]
      min_val <- min(col_data, na.rm = TRUE)
      max_val <- max(col_data, na.rm = TRUE)
      data[, i] <- (col_data - min_val) / (max_val - min_val)
    }
  }
  
  return(list(data = data, method = method))
}

smooth_data <- function(data, method = "gaussian", window = 5) {
  if (method == "gaussian") {
    # Simple moving average as placeholder
    for (i in 1:ncol(data)) {
      col_data <- data[, i]
      smoothed <- stats::filter(col_data, rep(1/window, window), sides = 2)
      data[, i] <- smoothed
    }
  }
  
  return(list(data = data, method = method))
}

calculate_preprocessing_quality <- function(processed_data, original_data) {
  # Simple quality metric based on signal-to-noise ratio improvement
  original_snr <- mean(apply(original_data, 2, function(x) mean(x, na.rm = TRUE) / sd(x, na.rm = TRUE)))
  processed_snr <- mean(apply(processed_data, 2, function(x) mean(x, na.rm = TRUE) / sd(x, na.rm = TRUE)))
  
  quality_score <- min(1.0, processed_snr / original_snr)
  return(quality_score)
}

#' Unified Data Curation Interface
#'
#' Curate and validate calcium imaging data with metadata.
#'
#' @param data Input data (matrix, data frame, or file path)
#' @param metadata Optional metadata information
#' @param validation_rules List of validation rules to apply
#' @param auto_preprocess Whether to apply automated preprocessing
#' @param ... Additional arguments
#' @return Curated data with metadata and validation results
#' @export
curate_data <- function(data, metadata = NULL, validation_rules = NULL, auto_preprocess = TRUE, ...) {
  message("Curating calcium imaging data")
  
  # Validate the data
  validation_results <- validate_calcium_data(data, validation_rules = validation_rules, ...)
  
  # Create metadata if not provided
  if (is.null(metadata)) {
    metadata <- create_experiment_metadata(list(
      experiment_id = paste0("exp_", format(Sys.time(), "%Y%m%d_%H%M%S")),
      date = format(Sys.Date(), "%Y-%m-%d"),
      time = format(Sys.time(), "%H:%M:%S")
    ), auto_generate = TRUE)
  }
  
  # Apply automated preprocessing if requested
  if (auto_preprocess && validation_results$validation_passed) {
    data <- automated_preprocessing(data, ...)
  }
  
  # Standardize data format
  data <- standardize_data_format(data, ...)
  
  # Detect and handle outliers
  data <- detect_and_handle_outliers(data, ...)
  
  return(list(
    data = data,
    metadata = metadata,
    validation_results = validation_results,
    curation_timestamp = Sys.time()
  ))
} 