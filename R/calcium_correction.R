#' Baseline correction for calcium imaging traces
#'
#' @description Applies background subtraction, normalization, and a LOESS
#'   baseline correction to calcium imaging data.
#'
#' @param rawtrace Data frame containing cell columns named `Cell_*` and
#'   background columns named `BG_*`.
#' @param n_bg Number of background columns to average. Defaults to 3.
#' @param span Smoothing parameter passed to `loess` for drift correction.
#'
#' @return Data frame of corrected cell traces with an added `Time` column.
#' @examples
#' dummy <- data.frame(
#'   Cell_1 = 1:10,
#'   BG_1 = rnorm(10),
#'   BG_2 = rnorm(10),
#'   BG_3 = rnorm(10)
#' )
#' calcium_correction(dummy)
#' @export
calcium_correction <- function(rawtrace, n_bg = 3, span = 0.45) {
  bg_cols <- grep('^BG', names(rawtrace))
  cell_cols <- grep('^Cell', names(rawtrace))

  bg_avg <- rowMeans(rawtrace[, bg_cols, drop = FALSE][, seq_len(n_bg)])
  bg_corrected <- rawtrace[, cell_cols, drop = FALSE] - bg_avg

  normalization <- colMeans(bg_corrected)
  normalized <- sweep(bg_corrected, 2, normalization, '/')

  time <- seq_len(nrow(bg_corrected))
  loess_mat <- vapply(
    seq_len(ncol(bg_corrected)),
    function(i) {
      smoothed <- predict(loess(bg_corrected[, i] ~ time, span = span))
      smoothed - 2
    },
    numeric(nrow(bg_corrected))
  )

  corrected <- bg_corrected - loess_mat
  corrected <- as.data.frame(corrected)
  names(corrected) <- names(rawtrace)[cell_cols]
  corrected$Time <- time
  corrected
}

# Backwards compatibility
calciumcorrection <- calcium_correction
