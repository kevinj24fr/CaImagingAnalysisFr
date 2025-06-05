#' Generate synthetic calcium imaging dataset
#'
#' @param n_cells Number of cells to simulate
#' @param n_time Number of time points
#' @param spike_prob Probability of a spike at any frame
#' @return Data frame with cell and background traces
#' @examples
#' df <- generate_synthetic_data(3, 500)
#' head(df)
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
