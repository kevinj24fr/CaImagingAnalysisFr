#' Background and baseline correction
#'
#' @param raw_df Data frame with Cell_* and BG_* columns
#' @param span Smoothing span for baseline estimation
#' @return Data frame with corrected traces and Time column
#' @export
calcium_correction <- function(raw_df, span = 0.45) {
  bg <- raw_df[, grepl("^BG_", names(raw_df))]
  bg_avg <- rowMeans(bg)

  cells <- raw_df[, grepl("^Cell_", names(raw_df))]
  corrected <- cells - bg_avg

  time <- seq_len(nrow(corrected))
  baseline <- apply(corrected, 2, function(x) {
    stats::predict(stats::loess(x ~ time, span = span))
  })

  norm <- corrected - baseline
  norm$Time <- time
  norm
}
