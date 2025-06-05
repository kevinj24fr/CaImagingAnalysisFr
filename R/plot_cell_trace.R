#' Plot corrected trace with detected spikes
#'
#' @param corrected_df Output of [calcium_correction()]
#' @param cell Name of a cell column
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
