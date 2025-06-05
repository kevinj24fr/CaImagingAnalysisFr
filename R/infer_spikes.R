#' Spike inference via OASIS
#'
#' This uses the Python package `oasis` to deconvolve calcium traces. If the
#' package is missing, it attempts installation and returns zeros on failure.
#'
#' @param trace Numeric vector of fluorescence values
#' @return Data frame with deconvolved activity and estimated spikes
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
