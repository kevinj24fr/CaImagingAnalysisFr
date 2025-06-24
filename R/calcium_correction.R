#' Calcium Trace Correction
#' 
#' Performs background subtraction, normalization, and baseline correction
#' on calcium imaging traces. This is the main correction function that
#' consolidates both modern and legacy approaches.
#' 
#' @param raw_df A data frame containing calcium traces with columns named
#'   `Cell_*` for cell traces and optionally `BG_*` for background traces
#' @param method Correction method: "modern" (default) or "legacy"
#' @param span Smoothing span for baseline estimation (0.01-1.0). Default 0.45
#' @param normalize Whether to normalize traces by their mean. Default TRUE
#' @param verbose Whether to show progress messages. Default TRUE
#' 
#' @return A data frame with corrected traces and a Time column
#' 
#' @examples
#' # Basic usage with both cell and background columns
#' raw <- generate_synthetic_data(3, 500)
#' corrected <- calcium_correction(raw)
#' 
#' # Usage with only cell columns (no background correction)
#' cell_traces <- raw[, grep("^Cell_", names(raw))]
#' corrected <- calcium_correction(cell_traces)
#' 
#' # Custom parameters
#' corrected <- calcium_correction(raw, span = 0.3, normalize = FALSE)
#' 
#' # Legacy method
#' corrected <- calcium_correction(raw, method = "legacy")
#' 
#' @export
calcium_correction <- function(raw_df, 
                              method = c("modern", "legacy"),
                              span = NULL,
                              normalize = TRUE,
                              verbose = TRUE) {
  
  # Validate method
  method <- match.arg(method)
  
  # Get default span from config if not provided
  if (is.null(span)) {
    span <- get_config()$default_span
  }
  
  # Check if background columns are present
  config <- get_config()
  bg_cols <- names(raw_df)[grepl(config$background_pattern, names(raw_df))]
  has_background <- length(bg_cols) > 0
  
  # Validate inputs - make background optional
  validate_data_frame(raw_df, require_cells = TRUE, require_background = FALSE)
  validate_numeric_param(span, get_config()$min_span, get_config()$max_span, "span")
  
  if (verbose) {
    if (has_background) {
      message("Processing ", ncol(raw_df), " traces using ", method, " method...")
    } else {
      message("Processing ", ncol(raw_df), " cell traces (no background correction)...")
    }
  }
  
  # Apply appropriate correction method
  if (method == "legacy") {
    result <- calcium_correction_legacy(raw_df, span, verbose)
  } else {
    result <- calcium_correction_modern(raw_df, span, normalize, verbose)
  }
  
  if (verbose) {
    message("Correction completed successfully")
  }
  
  return(result)
}

#' Modern calcium correction
#' 
#' @param raw_df Input data frame
#' @param span Smoothing span
#' @param normalize Whether to normalize
#' @param verbose Progress messages
#' @return Corrected data frame
#' @keywords internal
calcium_correction_modern <- function(raw_df, span, normalize, verbose) {
  config <- get_config()
  
  # Extract background and cell columns
  bg_cols <- names(raw_df)[grepl(config$background_pattern, names(raw_df))]
  cell_cols <- names(raw_df)[grepl(config$cell_pattern, names(raw_df))]
  
  # If no cell columns found, assume all columns are cell columns
  if (length(cell_cols) == 0) {
    cell_cols <- names(raw_df)
    if (verbose) message("  No cell columns found, treating all columns as cell traces...")
  }
  
  # Background correction (only if background columns exist)
  if (length(bg_cols) > 0) {
    if (verbose) message("  Performing background correction...")
    bg <- raw_df[, bg_cols, drop = FALSE]
    bg_avg <- rowMeans(bg, na.rm = TRUE)
    corrected <- raw_df[, cell_cols, drop = FALSE] - bg_avg
  } else {
    if (verbose) message("  No background columns found, skipping background correction...")
    corrected <- raw_df[, cell_cols, drop = FALSE]
  }
  
  # Normalization (optional)
  if (normalize) {
    if (verbose) message("  Normalizing traces...")
    norm_means <- colMeans(corrected, na.rm = TRUE)
    corrected <- sweep(corrected, 2, norm_means, "/")
  }
  
  # Baseline correction using loess
  if (verbose) message("  Performing baseline correction...")
  time <- seq_len(nrow(corrected))
  baseline <- apply(corrected, 2, function(x) {
    tryCatch({
      # Check if we have enough non-NA values
      if (sum(!is.na(x)) < 3) {
        warning("Insufficient non-NA values for baseline correction, using zeros")
        return(rep(0, length(x)))
      }
      stats::predict(stats::loess(x ~ time, span = span, na.action = na.exclude))
    }, error = function(e) {
      warning("Loess fitting failed for a trace, using linear trend instead")
      tryCatch({
        stats::predict(stats::lm(x ~ time, na.action = na.exclude))
      }, error = function(e2) {
        warning("Linear fitting also failed, using zeros")
        rep(0, length(x))
      })
    })
  })
  
  # Apply baseline correction
  norm <- corrected - baseline
  norm$Time <- time
  
  return(norm)
}

#' Legacy calcium correction
#' 
#' @param raw_df Input data frame
#' @param span Smoothing span
#' @param verbose Progress messages
#' @return Corrected data frame
#' @keywords internal
calcium_correction_legacy <- function(raw_df, span, verbose) {
  config <- get_config()
  
  # Extract background and cell columns
  bg_cols <- names(raw_df)[grepl(config$background_pattern, names(raw_df))]
  cell_cols <- names(raw_df)[grepl(config$cell_pattern, names(raw_df))]
  
  # If no cell columns found, assume all columns are cell columns
  if (length(cell_cols) == 0) {
    cell_cols <- names(raw_df)
    if (verbose) message("  No cell columns found, treating all columns as cell traces...")
  }
  
  if (verbose) message("  Performing legacy background correction...")
  
  # Background correction (legacy method) - only if background columns exist
  if (length(bg_cols) > 0) {
    background <- raw_df[, bg_cols, drop = FALSE]
    background$average <- rowMeans(background, na.rm = TRUE)
    bgcorrectedtraces <- raw_df[, cell_cols, drop = FALSE] - background$average
  } else {
    if (verbose) message("  No background columns found, skipping background correction...")
    bgcorrectedtraces <- raw_df[, cell_cols, drop = FALSE]
  }
  
  # Normalization (legacy method)
  if (verbose) message("  Performing legacy normalization...")
  normalizationvalues <- colMeans(bgcorrectedtraces, na.rm = TRUE)
  normalizedtraces <- sweep(bgcorrectedtraces, 2, normalizationvalues, "/")
  
  # Baseline correction (legacy method)
  if (verbose) message("  Performing legacy baseline correction...")
  loessmatrix <- data.frame(matrix(NA, nrow = nrow(bgcorrectedtraces), 
                                   ncol = ncol(bgcorrectedtraces)))
  Time <- seq_len(nrow(bgcorrectedtraces))
  
  for (i in seq_len(ncol(bgcorrectedtraces))) {
    Testintermediate <- as.matrix(bgcorrectedtraces[, i])
    loessintermediate <- stats::loess(Testintermediate ~ Time, span = span, 
                                     na.action = na.exclude)
    smoothed <- stats::predict(loessintermediate)
    loessmatrix[, i] <- smoothed - 2  # Legacy offset
  }
  
  correctedmatrix <- bgcorrectedtraces - loessmatrix
  correctedmatrix$Time <- Time
  
  return(correctedmatrix)
}
