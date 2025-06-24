#' Bayesian and Probabilistic Modeling for Calcium Imaging
#' 
#' This module provides Bayesian inference methods for calcium imaging data,
#' including Bayesian spike inference, parameter estimation, and uncertainty quantification.
#' 
#' @name bayesian_modeling
#' @docType package
NULL

#' Bayesian Spike Inference using PyMC
#' 
#' Performs Bayesian spike inference using PyMC for uncertainty quantification.
#' 
#' @param calcium_trace Numeric vector of calcium fluorescence values
#' @param sampling_rate Numeric, sampling rate in Hz (default: 30)
#' @param n_samples Integer, number of MCMC samples (default: 1000)
#' @param n_chains Integer, number of MCMC chains (default: 4)
#' @param prior_spike_rate Numeric, prior spike rate (default: 0.1)
#' @param prior_decay_rate Numeric, prior decay rate (default: 0.95)
#' @param prior_noise_sd Numeric, prior noise standard deviation (default: 0.1)
#' @param random_seed Integer, random seed for reproducibility (default: 42)
#' @param ... Additional arguments passed to PyMC
#' 
#' @return List containing:
#'   - spikes: Posterior mean spike probabilities
#'   - spike_samples: MCMC samples of spike indicators
#'   - parameters: Posterior parameter estimates
#'   - diagnostics: MCMC diagnostics
#'   - model_info: Model information
#' 
#' @examples
#' \dontrun{
#' # Generate synthetic data
#' data <- generate_synthetic_data(n_cells = 1, n_frames = 1000)
#' 
#' # Bayesian spike inference
#' result <- bayesian_spike_inference(data$calcium_traces[, 1])
#' 
#' # Plot results
#' plot_bayesian_results(result)
#' }
#' 
#' @export
bayesian_spike_inference <- function(calcium_trace,
                                   sampling_rate = 30,
                                   n_samples = 1000,
                                   n_chains = 4,
                                   prior_spike_rate = 0.1,
                                   prior_decay_rate = 0.95,
                                   prior_noise_sd = 0.1,
                                   random_seed = 42,
                                   ...) {
  
  # Input validation
  validate_calcium_trace(calcium_trace)
  validate_positive_numeric(sampling_rate, "sampling_rate")
  validate_positive_integer(n_samples, "n_samples")
  validate_positive_integer(n_chains, "n_chains")
  validate_probability(prior_spike_rate, "prior_spike_rate")
  validate_probability(prior_decay_rate, "prior_decay_rate")
  validate_positive_numeric(prior_noise_sd, "prior_noise_sd")
  validate_positive_integer(random_seed, "random_seed")
  
  tryCatch({
    # Check Python dependencies
    manage_python_dependencies(c("pymc", "numpy", "arviz"))
    
    # Prepare data
    n_frames <- length(calcium_trace)
    calcium_trace <- as.numeric(calcium_trace)
    
    # Call Python function
    result <- reticulate::py_call(
      reticulate::py_eval("bayesian_spike_inference_python"),
      calcium_trace = calcium_trace,
      sampling_rate = sampling_rate,
      n_samples = n_samples,
      n_chains = n_chains,
      prior_spike_rate = prior_spike_rate,
      prior_decay_rate = prior_decay_rate,
      prior_noise_sd = prior_noise_sd,
      random_seed = random_seed,
      ...
    )
    
    # Convert to R format
    result$spikes <- as.numeric(result$spikes)
    result$spike_samples <- as.matrix(result$spike_samples)
    result$parameters <- as.list(result$parameters)
    result$diagnostics <- as.list(result$diagnostics)
    
    # Add metadata
    result$model_info <- list(
      method = "bayesian_spike_inference",
      n_frames = n_frames,
      sampling_rate = sampling_rate,
      n_samples = n_samples,
      n_chains = n_chains,
      timestamp = Sys.time()
    )
    
    class(result) <- "bayesian_spike_result"
    return(result)
    
  }, error = function(e) {
    stop("Bayesian spike inference failed: ", e$message)
  })
}

#' Bayesian Parameter Estimation
#' 
#' Estimates model parameters using Bayesian inference with uncertainty quantification.
#' 
#' @param calcium_trace Numeric vector of calcium fluorescence values
#' @param model_type Character, model type ("exponential", "gamma", "custom")
#' @param priors List of prior distributions
#' @param n_samples Integer, number of MCMC samples (default: 2000)
#' @param n_chains Integer, number of MCMC chains (default: 4)
#' @param random_seed Integer, random seed (default: 42)
#' @param ... Additional arguments
#' 
#' @return List containing posterior parameter estimates and diagnostics
#' 
#' @export
bayesian_parameter_estimation <- function(calcium_trace,
                                        model_type = "exponential",
                                        priors = NULL,
                                        n_samples = 2000,
                                        n_chains = 4,
                                        random_seed = 42,
                                        ...) {
  
  # Input validation
  validate_calcium_trace(calcium_trace)
  validate_character(model_type, "model_type")
  validate_positive_integer(n_samples, "n_samples")
  validate_positive_integer(n_chains, "n_chains")
  validate_positive_integer(random_seed, "random_seed")
  
  # Set default priors if not provided
  if (is.null(priors)) {
    priors <- list(
      decay_rate = list(dist = "beta", alpha = 9.5, beta = 0.5),
      noise_sd = list(dist = "half_normal", sigma = 0.1),
      baseline = list(dist = "normal", mu = 0, sigma = 1)
    )
  }
  
  tryCatch({
    # Check Python dependencies
    manage_python_dependencies(c("pymc", "numpy", "arviz"))
    
    # Call Python function
    result <- reticulate::py_call(
      reticulate::py_eval("bayesian_parameter_estimation_python"),
      calcium_trace = as.numeric(calcium_trace),
      model_type = model_type,
      priors = priors,
      n_samples = n_samples,
      n_chains = n_chains,
      random_seed = random_seed,
      ...
    )
    
    # Convert to R format
    result$posterior_samples <- as.matrix(result$posterior_samples)
    result$summary <- as.data.frame(result$summary)
    result$diagnostics <- as.list(result$diagnostics)
    
    class(result) <- "bayesian_parameter_result"
    return(result)
    
  }, error = function(e) {
    stop("Bayesian parameter estimation failed: ", e$message)
  })
}

#' Uncertainty Quantification
#' 
#' Quantifies uncertainty in spike inference and parameter estimates.
#' 
#' @param bayesian_result Result from bayesian_spike_inference or bayesian_parameter_estimation
#' @param confidence_level Numeric, confidence level (default: 0.95)
#' @param method Character, uncertainty method ("credible_interval", "bootstrap", "ensemble")
#' @param ... Additional arguments
#' 
#' @return List containing uncertainty measures
#' 
#' @export
uncertainty_quantification <- function(bayesian_result,
                                     confidence_level = 0.95,
                                     method = "credible_interval",
                                     ...) {
  
  # Input validation
  if (!inherits(bayesian_result, c("bayesian_spike_result", "bayesian_parameter_result"))) {
    stop("bayesian_result must be from bayesian_spike_inference or bayesian_parameter_estimation")
  }
  validate_probability(confidence_level, "confidence_level")
  validate_character(method, "method")
  
  tryCatch({
    # Check Python dependencies
    manage_python_dependencies(c("numpy", "scipy"))
    
    # Call Python function
    result <- reticulate::py_call(
      reticulate::py_eval("uncertainty_quantification_python"),
      bayesian_result = bayesian_result,
      confidence_level = confidence_level,
      method = method,
      ...
    )
    
    # Convert to R format
    result$credible_intervals <- as.data.frame(result$credible_intervals)
    result$uncertainty_metrics <- as.list(result$uncertainty_metrics)
    
    class(result) <- "uncertainty_result"
    return(result)
    
  }, error = function(e) {
    stop("Uncertainty quantification failed: ", e$message)
  })
}

#' Probabilistic Spike Detection
#' 
#' Detects spikes using probabilistic methods with uncertainty estimates.
#' 
#' @param calcium_trace Numeric vector of calcium fluorescence values
#' @param method Character, method ("threshold", "changepoint", "hmm")
#' @param threshold Numeric, detection threshold (default: NULL, auto-detect)
#' @param min_spike_interval Integer, minimum frames between spikes (default: 5)
#' @param uncertainty_estimation Logical, estimate uncertainty (default: TRUE)
#' @param ... Additional arguments
#' 
#' @return List containing spike detection results and uncertainty
#' 
#' @export
probabilistic_spike_detection <- function(calcium_trace,
                                        method = "threshold",
                                        threshold = NULL,
                                        min_spike_interval = 5,
                                        uncertainty_estimation = TRUE,
                                        ...) {
  
  # Input validation
  validate_calcium_trace(calcium_trace)
  validate_character(method, "method")
  validate_positive_integer(min_spike_interval, "min_spike_interval")
  validate_logical(uncertainty_estimation, "uncertainty_estimation")
  
  if (!is.null(threshold)) {
    validate_positive_numeric(threshold, "threshold")
  }
  
  tryCatch({
    # Check Python dependencies
    manage_python_dependencies(c("numpy", "scipy", "scikit_learn"))
    
    # Call Python function
    result <- reticulate::py_call(
      reticulate::py_eval("probabilistic_spike_detection_python"),
      calcium_trace = as.numeric(calcium_trace),
      method = method,
      threshold = threshold,
      min_spike_interval = min_spike_interval,
      uncertainty_estimation = uncertainty_estimation,
      ...
    )
    
    # Convert to R format
    result$spike_times <- as.integer(result$spike_times)
    result$spike_probabilities <- as.numeric(result$spike_probabilities)
    if (!is.null(result$uncertainty)) {
      result$uncertainty <- as.numeric(result$uncertainty)
    }
    
    class(result) <- "probabilistic_spike_result"
    return(result)
    
  }, error = function(e) {
    stop("Probabilistic spike detection failed: ", e$message)
  })
}

#' Hierarchical Bayesian Modeling
#' 
#' Fits hierarchical Bayesian models for multi-cell calcium imaging data.
#' 
#' @param calcium_traces Matrix of calcium traces (cells x time)
#' @param group_structure List defining group structure (optional)
#' @param model_type Character, model type ("hierarchical", "multilevel")
#' @param n_samples Integer, number of MCMC samples (default: 2000)
#' @param n_chains Integer, number of MCMC chains (default: 4)
#' @param random_seed Integer, random seed (default: 42)
#' @param ... Additional arguments
#' 
#' @return List containing hierarchical model results
#' 
#' @export
hierarchical_bayesian_modeling <- function(calcium_traces,
                                         group_structure = NULL,
                                         model_type = "hierarchical",
                                         n_samples = 2000,
                                         n_chains = 4,
                                         random_seed = 42,
                                         ...) {
  
  # Input validation
  validate_calcium_traces_matrix(calcium_traces)
  validate_character(model_type, "model_type")
  validate_positive_integer(n_samples, "n_samples")
  validate_positive_integer(n_chains, "n_chains")
  validate_positive_integer(random_seed, "random_seed")
  
  tryCatch({
    # Check Python dependencies
    manage_python_dependencies(c("pymc", "numpy", "arviz"))
    
    # Call Python function
    result <- reticulate::py_call(
      reticulate::py_eval("hierarchical_bayesian_modeling_python"),
      calcium_traces = as.matrix(calcium_traces),
      group_structure = group_structure,
      model_type = model_type,
      n_samples = n_samples,
      n_chains = n_chains,
      random_seed = random_seed,
      ...
    )
    
    # Convert to R format
    result$group_effects <- as.data.frame(result$group_effects)
    result$individual_effects <- as.data.frame(result$individual_effects)
    result$diagnostics <- as.list(result$diagnostics)
    
    class(result) <- "hierarchical_bayesian_result"
    return(result)
    
  }, error = function(e) {
    stop("Hierarchical Bayesian modeling failed: ", e$message)
  })
}

#' Plot Bayesian Results
#' 
#' Creates comprehensive plots for Bayesian analysis results.
#' 
#' @param bayesian_result Result from Bayesian analysis functions
#' @param plot_type Character, type of plot ("trace", "posterior", "uncertainty", "all")
#' @param confidence_level Numeric, confidence level for intervals (default: 0.95)
#' @param ... Additional plotting arguments
#' 
#' @return ggplot object or list of plots
#' 
#' @export
plot_bayesian_results <- function(bayesian_result,
                                plot_type = "all",
                                confidence_level = 0.95,
                                ...) {
  
  # Input validation
  if (!inherits(bayesian_result, c("bayesian_spike_result", "bayesian_parameter_result", 
                                  "uncertainty_result", "probabilistic_spike_result", 
                                  "hierarchical_bayesian_result"))) {
    stop("Invalid bayesian_result object")
  }
  validate_character(plot_type, "plot_type")
  validate_probability(confidence_level, "confidence_level")
  
  tryCatch({
    # Check Python dependencies
    manage_python_dependencies(c("matplotlib", "seaborn"))
    
    # Call Python function
    plots <- reticulate::py_call(
      reticulate::py_eval("plot_bayesian_results_python"),
      bayesian_result = bayesian_result,
      plot_type = plot_type,
      confidence_level = confidence_level,
      ...
    )
    
    return(plots)
    
  }, error = function(e) {
    stop("Plotting Bayesian results failed: ", e$message)
  })
}

#' Bayesian Model Comparison
#' 
#' Compares different Bayesian models using information criteria.
#' 
#' @param model_results List of Bayesian model results
#' @param comparison_metrics Character vector, metrics to use ("waic", "loo", "dic")
#' @param ... Additional arguments
#' 
#' @return Data frame with model comparison results
#' 
#' @export
bayesian_model_comparison <- function(model_results,
                                    comparison_metrics = c("waic", "loo"),
                                    ...) {
  
  # Input validation
  if (!is.list(model_results) || length(model_results) < 2) {
    stop("model_results must be a list with at least 2 models")
  }
  validate_character_vector(comparison_metrics, "comparison_metrics")
  
  tryCatch({
    # Check Python dependencies
    manage_python_dependencies(c("arviz", "numpy"))
    
    # Call Python function
    result <- reticulate::py_call(
      reticulate::py_eval("bayesian_model_comparison_python"),
      model_results = model_results,
      comparison_metrics = comparison_metrics,
      ...
    )
    
    # Convert to R format
    result$comparison_table <- as.data.frame(result$comparison_table)
    result$model_weights <- as.numeric(result$model_weights)
    
    return(result)
    
  }, error = function(e) {
    stop("Bayesian model comparison failed: ", e$message)
  })
}

#' Bayesian Predictive Checks
#' 
#' Performs posterior predictive checks for model validation.
#' 
#' @param bayesian_result Result from Bayesian analysis
#' @param n_predictions Integer, number of predictions to generate (default: 100)
#' @param test_statistics Character vector, test statistics to compute
#' @param ... Additional arguments
#' 
#' @return List containing predictive check results
#' 
#' @export
bayesian_predictive_checks <- function(bayesian_result,
                                     n_predictions = 100,
                                     test_statistics = c("mean", "std", "min", "max"),
                                     ...) {
  
  # Input validation
  if (!inherits(bayesian_result, c("bayesian_spike_result", "bayesian_parameter_result"))) {
    stop("Invalid bayesian_result object")
  }
  validate_positive_integer(n_predictions, "n_predictions")
  validate_character_vector(test_statistics, "test_statistics")
  
  tryCatch({
    # Check Python dependencies
    manage_python_dependencies(c("pymc", "numpy", "arviz"))
    
    # Call Python function
    result <- reticulate::py_call(
      reticulate::py_eval("bayesian_predictive_checks_python"),
      bayesian_result = bayesian_result,
      n_predictions = n_predictions,
      test_statistics = test_statistics,
      ...
    )
    
    # Convert to R format
    result$predictions <- as.matrix(result$predictions)
    result$test_statistics <- as.data.frame(result$test_statistics)
    result$p_values <- as.numeric(result$p_values)
    
    class(result) <- "predictive_check_result"
    return(result)
    
  }, error = function(e) {
    stop("Bayesian predictive checks failed: ", e$message)
  })
}

# Helper functions for validation
validate_calcium_trace <- function(trace) {
  if (!is.numeric(trace) || length(trace) < 10) {
    stop("calcium_trace must be a numeric vector with at least 10 elements")
  }
  if (any(is.na(trace)) || any(is.infinite(trace))) {
    stop("calcium_trace contains NA or infinite values")
  }
}

validate_calcium_traces_matrix <- function(traces) {
  if (!is.matrix(traces) || nrow(traces) < 1 || ncol(traces) < 10) {
    stop("calcium_traces must be a matrix with at least 1 row and 10 columns")
  }
  if (any(is.na(traces)) || any(is.infinite(traces))) {
    stop("calcium_traces contains NA or infinite values")
  }
}

validate_probability <- function(x, name) {
  if (!is.numeric(x) || x <= 0 || x >= 1) {
    stop(name, " must be a probability between 0 and 1")
  }
}

validate_positive_numeric <- function(x, name) {
  if (!is.numeric(x) || x <= 0) {
    stop(name, " must be a positive numeric value")
  }
}

validate_positive_integer <- function(x, name) {
  if (!is.numeric(x) || x <= 0 || x != round(x)) {
    stop(name, " must be a positive integer")
  }
}

validate_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1) {
    stop(name, " must be a single character string")
  }
}

validate_character_vector <- function(x, name) {
  if (!is.character(x) || length(x) < 1) {
    stop(name, " must be a character vector with at least one element")
  }
}

validate_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1) {
    stop(name, " must be a single logical value")
  }
} 