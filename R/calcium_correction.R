#' Calcium Trace Correction
#' 
#' Performs background subtraction, normalization, and baseline correction
#' on calcium imaging traces. This is the main correction function that
#' consolidates both modern and legacy approaches.
#' 
#' @param raw_df A data frame containing calcium traces with columns named
#'   `Cell_*` for cell traces and `BG_*` for background traces
#' @param method Correction method: "modern" (default) or "legacy"
#' @param span Smoothing span for baseline estimation (0.01-1.0). Default 0.45
#' @param normalize Whether to normalize traces by their mean. Default TRUE
#' @param verbose Whether to show progress messages. Default TRUE
#' 
#' @return A data frame with corrected traces and a Time column
#' 
#' @examples
#' # Basic usage
#' raw <- generate_synthetic_data(3, 500)
#' corrected <- calcium_correction(raw)
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
  
  # Validate inputs
  validate_data_frame(raw_df)
  validate_numeric_param(span, get_config()$min_span, get_config()$max_span, "span")
  
  if (verbose) {
    message("Processing ", ncol(raw_df), " traces using ", method, " method...")
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
  
  # Background correction
  if (verbose) message("  Performing background correction...")
  bg <- raw_df[, bg_cols, drop = FALSE]
  bg_avg <- rowMeans(bg, na.rm = TRUE)
  corrected <- raw_df[, cell_cols, drop = FALSE] - bg_avg
  
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
  
  if (verbose) message("  Performing legacy background correction...")
  
  # Background correction (legacy method)
  background <- raw_df[, bg_cols, drop = FALSE]
  background$average <- rowMeans(background, na.rm = TRUE)
  bgcorrectedtraces <- raw_df[, cell_cols, drop = FALSE] - background$average
  
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
