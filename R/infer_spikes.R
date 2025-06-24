#' Spike Inference via OASIS, CaImAn, Suite2p, or Deep Learning
#' 
#' Performs spike inference on calcium traces using the OASIS, CaImAn, Suite2p, or deep learning (CASCADE) algorithm.
#' 
#' @param trace Numeric vector of fluorescence values
#' @param method Inference method: "oasis" (default), "caiman", "suite2p", or "deep"
#' @param fallback Whether to use fallback method if primary fails (default: TRUE)
#' @param model_path Path to pretrained deep learning model (for method = "deep")
#' @param verbose Whether to show progress messages (default: TRUE)
#' 
#' @return Data frame with deconvolved activity and estimated spikes
#' 
#' @examples
#' # Basic usage
#' raw <- generate_synthetic_data(1, 500)
#' corrected <- calcium_correction(raw)
#' spikes <- infer_spikes(corrected$Cell_1)
#' 
#' # Try deep learning method (requires model)
#' # spikes <- infer_spikes(corrected$Cell_1, method = "deep", model_path = "cascade_model.h5")
#' 
#' @export
infer_spikes <- function(trace, 
                        method = c("oasis", "caiman", "suite2p", "deep"),
                        fallback = TRUE,
                        model_path = NULL,
                        verbose = TRUE) {
  
  # Validate inputs
  validate_trace(trace)
  method <- match.arg(method)
  
  if (verbose) {
    message("Performing spike inference using ", method, " method...")
  }
  
  # Try primary method
  result <- tryCatch({
    if (method == "deep") {
      infer_spikes_deep(trace, model_path = model_path, verbose = verbose)
    } else {
      infer_spikes_primary(trace, method, verbose)
    }
  }, error = function(e) {
    if (fallback && method != "deep") {
      if (verbose) {
        warning("Primary method failed, using fallback: ", e$message)
      }
      infer_spikes_fallback(trace, verbose)
    } else {
      stop("Spike inference failed: ", e$message)
    }
  })
  
  if (verbose) {
    message("Spike inference completed successfully")
  }
  
  return(result)
}

#' Primary spike inference methods
#' 
#' @param trace Input trace
#' @param method Method to use
#' @param verbose Progress messages
#' @return Spike inference results
#' @keywords internal
infer_spikes_primary <- function(trace, method, verbose) {
  
  if (method == "oasis") {
    return(infer_spikes_oasis(trace, verbose))
  } else if (method == "caiman") {
    return(infer_spikes_caiman(trace, verbose))
  } else if (method == "suite2p") {
    return(infer_spikes_suite2p(trace, verbose))
  }
}

#' OASIS spike inference
#' 
#' @param trace Input trace
#' @param verbose Progress messages
#' @return OASIS results
#' @keywords internal
infer_spikes_oasis <- function(trace, verbose) {
  if (verbose) message("  Using OASIS algorithm...")
  
  # Check Python dependencies
  deps <- manage_python_dependencies(packages = "oasis", install_missing = TRUE, verbose = FALSE)
  
  if (!deps$oasis$available) {
    stop("OASIS package is not available and could not be installed")
  }
  
  # Import and run OASIS
  oasis <- reticulate::import("oasis.functions")
  result <- oasis$deconvolve(trace)
  
  # Format results
  data.frame(
    fit = as.numeric(result[[1]]),
    spike = as.numeric(result[[2]]),
    stringsAsFactors = FALSE
  )
}

#' CaImAn spike inference
#' 
#' @param trace Input trace
#' @param verbose Progress messages
#' @return CaImAn results
#' @keywords internal
infer_spikes_caiman <- function(trace, verbose) {
  if (verbose) message("  Using CaImAn algorithm...")
  
  # Check Python dependencies
  deps <- manage_python_dependencies(packages = "caiman", install_missing = TRUE, verbose = FALSE)
  
  if (!deps$caiman$available) {
    stop("CaImAn package is not available and could not be installed")
  }
  
  # Import and run CaImAn
  cm <- reticulate::import("caiman")
  res <- cm$deconvolution$cd_oasis(trace)
  
  # Format results
  data.frame(
    fit = as.numeric(res[[1]]),
    spike = as.numeric(res[[2]]),
    stringsAsFactors = FALSE
  )
}

#' Suite2p spike inference
#' 
#' @param trace Input trace
#' @param verbose Progress messages
#' @return Suite2p results
#' @keywords internal
infer_spikes_suite2p <- function(trace, verbose) {
  if (verbose) message("  Using Suite2p algorithm...")
  
  # Check Python dependencies
  deps <- manage_python_dependencies(packages = "suite2p", install_missing = TRUE, verbose = FALSE)
  
  if (!deps$suite2p$available) {
    stop("Suite2p package is not available and could not be installed")
  }
  
  # Import and run Suite2p
  s2p <- reticulate::import("suite2p")
  res <- s2p$deconv$deconvolve(trace)
  
  # Format results
  data.frame(
    fit = as.numeric(res[["C_dec"]]),
    spike = as.numeric(res[["spks"]]),
    stringsAsFactors = FALSE
  )
}

#' Fallback spike inference
#' 
#' Simple threshold-based spike detection as fallback
#' 
#' @param trace Input trace
#' @param verbose Progress messages
#' @return Fallback results
#' @keywords internal
infer_spikes_fallback <- function(trace, verbose) {
  if (verbose) message("  Using fallback threshold method...")
  
  config <- get_config()
  
  # Simple threshold-based detection
  threshold <- config$default_threshold_multiplier * sd(trace, na.rm = TRUE)
  spikes <- trace > threshold
  
  # Simple smoothing for fit
  fit <- stats::filter(trace, rep(1/5, 5), sides = 2)
  fit[is.na(fit)] <- trace[is.na(fit)]
  
  data.frame(
    fit = as.numeric(fit),
    spike = as.numeric(spikes),
    stringsAsFactors = FALSE
  )
}
