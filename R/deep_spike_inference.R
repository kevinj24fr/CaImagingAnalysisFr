#' Deep Learning-Based Spike Inference (CASCADE)
#'
#' Performs spike inference using a deep neural network model (CASCADE) via Python.
#' Requires the 'cascade' Python package and a pretrained model.
#'
#' @param trace Numeric vector of fluorescence values
#' @param model_path Path to the pretrained CASCADE model (.h5 or .pt)
#' @param install_missing Whether to install missing Python dependencies
#' @param verbose Whether to show progress messages
#' @return Data frame with deconvolved activity and estimated spikes
#' @export
infer_spikes_deep <- function(trace, model_path = NULL, install_missing = TRUE, verbose = TRUE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for deep learning spike inference.")
  }
  
  # Check Python dependencies
  deps <- manage_python_dependencies(packages = c("cascade"), install_missing = install_missing, verbose = verbose)
  if (!deps$cascade$available) {
    stop("CASCADE Python package is not available and could not be installed.")
  }
  
  if (is.null(model_path)) {
    stop("A path to a pretrained CASCADE model must be provided.")
  }
  if (!file.exists(model_path)) {
    stop("Model file not found: ", model_path)
  }
  
  if (verbose) message("Loading CASCADE model from ", model_path)
  cascade <- reticulate::import("cascade")
  model <- cascade$load_model(model_path)
  
  # Prepare input
  input_trace <- as.numeric(trace)
  input_trace <- matrix(input_trace, nrow = 1)  # shape: (1, T)
  
  if (verbose) message("Running deep learning spike inference...")
  result <- model$predict(input_trace)
  
  # Assume result is a 2D array: (1, T, 2) or (1, T)
  if (length(dim(result)) == 3) {
    fit <- as.numeric(result[1, , 1])
    spike <- as.numeric(result[1, , 2])
  } else {
    fit <- as.numeric(result[1, ])
    spike <- fit  # If only one output, treat as spike
  }
  
  data.frame(
    fit = fit,
    spike = spike,
    stringsAsFactors = FALSE
  )
}

#' Check if deep learning spike inference is available
#' @return TRUE if available, FALSE otherwise
#' @keywords internal
is_deep_spike_inference_available <- function() {
  requireNamespace("reticulate", quietly = TRUE) &&
    reticulate::py_module_available("cascade")
} 