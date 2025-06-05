#' Spike inference with multiple backends
#'
#' Provides access to additional algorithms such as CaImAn and Suite2p
#' via `reticulate`.
#'
#' @param trace Numeric vector of fluorescence values
#' @param method One of "oasis", "caiman", or "suite2p"
#' @return Data frame with deconvolved activity and estimated spikes
#' @export
infer_spikes_advanced <- function(trace, method = c("oasis", "caiman", "suite2p")) {
  method <- match.arg(method)
  if (method == "oasis") {
    return(infer_spikes(trace))
  }

  if (method == "caiman") {
    if (!reticulate::py_module_available("caiman")) {
      tryCatch(reticulate::py_install("caiman"), error = function(e) NULL)
    }
    if (!reticulate::py_module_available("caiman")) {
      warning("caiman not available; falling back to oasis")
      return(infer_spikes(trace))
    }
    cm <- reticulate::import("caiman")
    # Simple wrapper using OASIS from CaImAn
    res <- cm$deconvolution$cd_oasis(trace)
    return(data.frame(fit = res[[1]], spike = res[[2]]))
  }

  if (method == "suite2p") {
    if (!reticulate::py_module_available("suite2p")) {
      tryCatch(reticulate::py_install("suite2p"), error = function(e) NULL)
    }
    if (!reticulate::py_module_available("suite2p")) {
      warning("suite2p not available; falling back to oasis")
      return(infer_spikes(trace))
    }
    s2p <- reticulate::import("suite2p")
    res <- s2p$deconv$deconvolve(trace)
    return(data.frame(fit = res[["C_dec" ]], spike = res[["spks" ]]))
  }
}
