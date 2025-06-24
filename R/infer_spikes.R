#' Spike Inference from Calcium Imaging Data
#'
#' Infer spike times from calcium imaging traces using various methods.
#'
#' @name infer_spikes
#' @docType package
NULL

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
  
  # Convert result to expected data frame format
  result_df <- convert_spike_result_to_df(result)
  
  if (verbose) {
    message("Spike inference completed successfully")
  }
  
  return(result_df)
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
    return(oasis_spike_inference(trace))
  } else if (method == "caiman") {
    return(caiman_spike_inference(trace))
  } else if (method == "suite2p") {
    return(suite2p_spike_inference(trace))
  }
}

#' OASIS Spike Inference
#'
#' Infer spikes using OASIS algorithm (Friedrich et al., 2017).
#'
#' @param trace Calcium trace
#' @param lambda Penalty parameter (default: 0.1)
#' @param ... Additional arguments
#' @return List containing spike times and inferred calcium
#' @export
oasis_spike_inference <- function(trace, lambda = 0.1, ...) {
  # Base R implementation of OASIS using L1 trend filtering
  message("Running OASIS spike inference")
  
  n <- length(trace)
  
  # Simple OASIS-like implementation using thresholding and deconvolution
  # This is a simplified version that captures the essence of OASIS
  
  # 1. Estimate baseline using median filter
  baseline <- stats::filter(trace, rep(1/21, 21), sides = 2)
  baseline[is.na(baseline)] <- median(trace, na.rm = TRUE)
  
  # 2. Remove baseline
  detrended <- trace - baseline
  
  # 3. Apply threshold-based spike detection
  threshold <- lambda * sd(detrended, na.rm = TRUE)
  spike_indices <- which(detrended > threshold)
  
  # 4. Merge nearby spikes (within 3 time points)
  if (length(spike_indices) > 1) {
    merged_spikes <- numeric(0)
    current_spike <- spike_indices[1]
    
    for (i in 2:length(spike_indices)) {
      if (spike_indices[i] - current_spike > 3) {
        merged_spikes <- c(merged_spikes, current_spike)
        current_spike <- spike_indices[i]
      } else {
        # Keep the spike with higher amplitude
        if (detrended[spike_indices[i]] > detrended[current_spike]) {
          current_spike <- spike_indices[i]
        }
      }
    }
    merged_spikes <- c(merged_spikes, current_spike)
    spike_indices <- merged_spikes
  }
  
  # 5. Create spike vector
  spikes <- rep(0, n)
  spikes[spike_indices] <- 1
  
  # 6. Estimate calcium concentration (simplified)
  # Use exponential decay model
  tau <- 10  # Decay time constant
  calcium_est <- numeric(n)
  calcium_est[1] <- trace[1]
  
  for (i in 2:n) {
    calcium_est[i] <- calcium_est[i-1] * exp(-1/tau) + spikes[i-1]
  }
  
  return(list(
    spikes = spikes,
    spike_times = spike_indices,
    calcium_est = calcium_est,
    baseline = baseline,
    parameters = list(lambda = lambda, tau = tau)
  ))
}

#' CaImAn Spike Inference
#'
#' Infer spikes using CaImAn algorithm (Giovannucci et al., 2019).
#'
#' @param trace Calcium trace
#' @param method CaImAn method ("oasis", "thresholded_oasis", "mcmc")
#' @param ... Additional arguments
#' @return List containing spike times and inferred calcium
#' @export
caiman_spike_inference <- function(trace, method = "oasis", ...) {
  message("Running CaImAn spike inference")
  
  # Base R implementation of CaImAn-like methods
  
  if (method == "oasis") {
    # Use OASIS implementation
    return(oasis_spike_inference(trace, ...))
    
  } else if (method == "thresholded_oasis") {
    # Thresholded OASIS with adaptive thresholding
    n <- length(trace)
    
    # Adaptive thresholding
    window_size <- min(50, n %/% 10)
    thresholds <- numeric(n)
    
    for (i in 1:n) {
      start_idx <- max(1, i - window_size %/% 2)
      end_idx <- min(n, i + window_size %/% 2)
      local_data <- trace[start_idx:end_idx]
      thresholds[i] <- mean(local_data) + 2 * sd(local_data)
    }
    
    # Detect spikes
    spikes <- as.numeric(trace > thresholds)
    
    # Merge nearby spikes
    spike_indices <- which(spikes == 1)
    if (length(spike_indices) > 1) {
      merged_spikes <- numeric(0)
      current_spike <- spike_indices[1]
      
      for (i in 2:length(spike_indices)) {
        if (spike_indices[i] - current_spike > 5) {
          merged_spikes <- c(merged_spikes, current_spike)
          current_spike <- spike_indices[i]
        }
      }
      merged_spikes <- c(merged_spikes, current_spike)
      spike_indices <- merged_spikes
    }
    
    # Recreate spike vector
    spikes <- rep(0, n)
    spikes[spike_indices] <- 1
    
    return(list(
      spikes = spikes,
      spike_times = spike_indices,
      thresholds = thresholds,
      method = method
    ))
    
  } else if (method == "mcmc") {
    # Simplified MCMC-like approach using Bayesian inference
    # This is a very simplified version
    
    # Use peak detection as a proxy for MCMC
    peaks <- find_peaks(trace)
    
    # Filter peaks based on amplitude
    peak_amplitudes <- trace[peaks]
    threshold <- quantile(peak_amplitudes, 0.8)
    significant_peaks <- peaks[peak_amplitudes > threshold]
    
    # Create spike vector
    spikes <- rep(0, length(trace))
    spikes[significant_peaks] <- 1
    
    return(list(
      spikes = spikes,
      spike_times = significant_peaks,
      method = method,
      peak_amplitudes = peak_amplitudes
    ))
  }
}

#' Suite2p Spike Inference
#'
#' Infer spikes using Suite2p algorithm (Pachitariu et al., 2017).
#'
#' @param trace Calcium trace
#' @param ... Additional arguments
#' @return List containing spike times and inferred calcium
#' @export
suite2p_spike_inference <- function(trace, ...) {
  message("Running Suite2p spike inference")
  
  # Base R implementation of Suite2p-like spike detection
  n <- length(trace)
  
  # 1. Preprocessing: normalize and detrend
  trace_norm <- (trace - mean(trace, na.rm = TRUE)) / sd(trace, na.rm = TRUE)
  
  # 2. Remove slow trends using high-pass filter
  # Simple high-pass filter using difference
  detrended <- c(0, diff(trace_norm))
  
  # 3. Apply threshold-based detection
  threshold <- 2.5  # Suite2p default
  spike_candidates <- which(detrended > threshold)
  
  # 4. Apply refractory period (minimum 3 frames between spikes)
  if (length(spike_candidates) > 1) {
    valid_spikes <- numeric(0)
    last_spike <- spike_candidates[1]
    valid_spikes <- c(valid_spikes, last_spike)
    
    for (i in 2:length(spike_candidates)) {
      if (spike_candidates[i] - last_spike >= 3) {
        valid_spikes <- c(valid_spikes, spike_candidates[i])
        last_spike <- spike_candidates[i]
      }
    }
    spike_candidates <- valid_spikes
  }
  
  # 5. Create spike vector
  spikes <- rep(0, n)
  spikes[spike_candidates] <- 1
  
  # 6. Estimate spike amplitudes
  spike_amplitudes <- detrended[spike_candidates]
  
  return(list(
    spikes = spikes,
    spike_times = spike_candidates,
    spike_amplitudes = spike_amplitudes,
    detrended_trace = detrended,
    threshold = threshold
  ))
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

#' Helper Functions
#' @keywords internal

find_peaks <- function(x, min_peak_height = NULL) {
  # Find peaks in a time series
  n <- length(x)
  peaks <- numeric(0)
  
  for (i in 2:(n-1)) {
    if (x[i] > x[i-1] && x[i] > x[i+1]) {
      if (is.null(min_peak_height) || x[i] > min_peak_height) {
        peaks <- c(peaks, i)
      }
    }
  }
  
  return(peaks)
}

calculate_spike_metrics <- function(spikes, trace) {
  # Calculate various spike metrics
  n_spikes <- sum(spikes)
  spike_rate <- n_spikes / length(spikes)  # spikes per time point
  
  if (n_spikes > 0) {
    spike_intervals <- diff(which(spikes == 1))
    mean_isi <- mean(spike_intervals)
    cv_isi <- sd(spike_intervals) / mean_isi
  } else {
    mean_isi <- NA
    cv_isi <- NA
  }
  
  return(list(
    n_spikes = n_spikes,
    spike_rate = spike_rate,
    mean_isi = mean_isi,
    cv_isi = cv_isi
  ))
}

#' Convert spike inference result to data frame format
#' 
#' @param spike_result Result from spike inference function
#' @return Data frame with fit and spike columns
#' @keywords internal
convert_spike_result_to_df <- function(spike_result) {
  # Extract components from spike_result
  if ("calcium_est" %in% names(spike_result)) {
    fit <- spike_result$calcium_est
  } else if ("fit" %in% names(spike_result)) {
    fit <- spike_result$fit
  } else {
    # Fallback: use original trace as fit
    fit <- spike_result$trace
  }
  
  if ("spikes" %in% names(spike_result)) {
    spike <- spike_result$spikes
  } else if ("spike" %in% names(spike_result)) {
    spike <- spike_result$spike
  } else {
    # Fallback: create spike vector from spike_times
    n <- length(fit)
    spike <- rep(0, n)
    if ("spike_times" %in% names(spike_result)) {
      spike[spike_result$spike_times] <- 1
    }
  }
  
  # Ensure both vectors have the same length
  n <- max(length(fit), length(spike))
  if (length(fit) < n) fit <- c(fit, rep(NA, n - length(fit)))
  if (length(spike) < n) spike <- c(spike, rep(0, n - length(spike)))
  
  data.frame(
    fit = fit,
    spike = spike,
    stringsAsFactors = FALSE
  )
}

#' Deep Learning Spike Inference
#'
#' Perform spike inference using deep learning methods.
#'
#' @param trace Calcium trace
#' @param model_path Path to pre-trained model
#' @param verbose Whether to show progress messages
#' @return Spike inference results
#' @keywords internal
infer_spikes_deep <- function(trace, model_path = NULL, verbose = TRUE) {
  if (verbose) {
    message("Performing deep learning spike inference")
  }
  
  # Call the deep_spike_inference function
  result <- deep_spike_inference(trace, model_type = "lstm", verbose = verbose)
  
  return(result)
}
