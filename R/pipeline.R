#' Generate synthetic calcium imaging dataset
#'
#' @param n_cells Number of cells to simulate
#' @param n_time  Number of time points
#' @param spike_prob Probability of a spike at any frame
#' @return data.frame with cell and background traces
#' @export
generate_synthetic_data <- function(n_cells = 5, n_time = 1000, spike_prob = 0.02) {
  set.seed(1)
  bg <- replicate(3, cumsum(rnorm(n_time, sd = 0.001)) + rnorm(n_time, sd = 0.02))
  bg <- as.data.frame(bg)
  names(bg) <- paste0("BG_", seq_len(ncol(bg)))

  cells <- replicate(n_cells, {
    trace <- cumsum(rnorm(n_time, sd = 0.001)) + rnorm(n_time, sd = 0.05)
    spikes <- rbinom(n_time, 1, spike_prob) == 1
    for (t in which(spikes)) {
      amp <- runif(1, 0.5, 2)
      decay <- amp * exp(-seq_len(50) / 10)
      end <- min(n_time, t + length(decay) - 1)
      trace[t:end] <- trace[t:end] + decay[seq_len(end - t + 1)]
    }
    trace
  })
  cells <- as.data.frame(cells)
  names(cells) <- paste0("Cell_", seq_len(ncol(cells)))

  cbind(cells, bg)
}

#' Background and baseline correction
#'
#' @param raw_df Data frame containing Cell_* and BG_* columns
#' @param span Smoothing span for baseline estimation
#' @return data frame with corrected traces and Time column
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

#' Spike inference using OASIS via reticulate
#'
#' @param trace Numeric vector of fluorescence values
#' @return data frame with deconvolved activity and estimated spikes
#' @export
infer_spikes <- function(trace) {
  if (!reticulate::py_module_available("oasis")) {
    tryCatch({
      reticulate::py_install("oasis")
    }, error = function(e) {
      warning("oasis package not available; spikes will be zero")
      return(data.frame(fit = rep(0, length(trace)), spike = rep(0, length(trace))))
    })
  }
  if (!reticulate::py_module_available("oasis")) {
    return(data.frame(fit = rep(0, length(trace)), spike = rep(0, length(trace))))
  }
  oasis <- reticulate::import("oasis.functions")
  result <- oasis$deconvolve(trace)
  data.frame(fit = result[[1]], spike = result[[2]])
}

#' Quick visualization of corrected traces with detected spikes
#'
#' @param corrected_df Output of calcium_correction
#' @param cell Column name of a single cell
#' @return ggplot object
#' @export
plot_cell_trace <- function(corrected_df, cell) {
  dat <- data.frame(Time = corrected_df$Time,
                    Signal = corrected_df[[cell]])
  spikes <- infer_spikes(dat$Signal)
  dat$Deconvolved <- spikes$fit
  dat$Spike <- spikes$spike > 0

  ggplot2::ggplot(dat, ggplot2::aes(Time, Signal)) +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::geom_line(ggplot2::aes(y = Deconvolved), color = "orange", alpha = 0.7) +
    ggplot2::geom_point(data = dat[dat$Spike, ], ggplot2::aes(y = Signal), color = "red", size = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = cell, y = "dF/F")
}
