#' Bayesian Modeling for Calcium Imaging Data
#'
#' Perform Bayesian analysis and modeling of calcium imaging data.
#'
#' @name bayesian_modeling
#' @docType package
NULL

#' Validate Probability
#'
#' Validate that a value is a probability between 0 and 1.
#'
#' @param value Value to check
#' @param name Name of the parameter (for error message)
#' @return NULL if valid, throws error if invalid
#' @keywords internal
validate_probability <- function(value, name = "value") {
  if (!is.numeric(value) || length(value) != 1 || value < 0 || value > 1) {
    stop(sprintf("%s must be a probability between 0 and 1", name))
  }
}

#' Validate Calcium Trace
#'
#' Validate input calcium trace for Bayesian analysis.
#'
#' @param trace Calcium trace to validate
#' @return NULL if valid, throws error if invalid
#' @export
validate_calcium_trace <- function(trace) {
  if (is.null(trace)) {
    stop("calcium_trace must be a numeric vector")
  }
  if (!is.numeric(trace)) {
    stop("calcium_trace must be a numeric vector")
  }
  if (length(trace) < 10) {
    stop("calcium_trace must be a numeric vector with at least 10 elements")
  }
  if (any(is.na(trace) | is.infinite(trace))) {
    stop("calcium_trace contains NA or infinite values")
  }
}

#' Validate Calcium Traces Matrix
#'
#' Validate input matrix of calcium traces for hierarchical Bayesian analysis.
#'
#' @param traces Matrix of calcium traces
#' @return NULL if valid, throws error if invalid
#' @export
validate_calcium_traces_matrix <- function(traces) {
  if (is.null(traces) || !is.matrix(traces)) {
    stop("calcium_traces must be a matrix")
  }
  if (nrow(traces) < 1 || ncol(traces) < 10) {
    stop("calcium_traces must be a matrix with at least 1 row and 10 columns")
  }
  if (any(is.na(traces) | is.infinite(traces))) {
    stop("calcium_traces contains NA or infinite values")
  }
}

#' Validate Positive Numeric
#'
#' Validate that a value is a positive numeric value.
#'
#' @param value Value to check
#' @param name Name of the parameter (for error message)
#' @return NULL if valid, throws error if invalid
#' @keywords internal
validate_positive_numeric <- function(value, name = "value") {
  if (!is.numeric(value) || length(value) != 1 || value <= 0) {
    stop(sprintf("%s must be a positive numeric value", name))
  }
}

#' Validate Positive Integer
#'
#' Validate that a value is a positive integer.
#'
#' @param value Value to check
#' @param name Name of the parameter (for error message)
#' @return NULL if valid, throws error if invalid
#' @keywords internal
validate_positive_integer <- function(value, name = "value") {
  if (!is.numeric(value) || length(value) != 1 || value <= 0 || value != as.integer(value)) {
    stop(sprintf("%s must be a positive integer", name))
  }
}

#' Validate Character
#'
#' Validate that a value is a single character string.
#'
#' @param value Value to check
#' @param name Name of the parameter (for error message)
#' @return NULL if valid, throws error if invalid
#' @keywords internal
validate_character <- function(value, name = "value") {
  if (!is.character(value) || length(value) != 1) {
    stop(sprintf("%s must be a single character string", name))
  }
}

#' Validate Character Vector
#'
#' Validate that a value is a character vector with at least one element.
#'
#' @param value Value to check
#' @param name Name of the parameter (for error message)
#' @return NULL if valid, throws error if invalid
#' @keywords internal
validate_character_vector <- function(value, name = "value") {
  if (!is.character(value)) {
    stop(sprintf("%s must be a character vector", name))
  }
  if (length(value) == 0) {
    stop(sprintf("%s must be a character vector with at least one element", name))
  }
}

#' Validate Logical
#'
#' Validate that a value is a single logical value.
#'
#' @param value Value to check
#' @param name Name of the parameter (for error message)
#' @return NULL if valid, throws error if invalid
#' @keywords internal
validate_logical <- function(value, name = "value") {
  if (!is.logical(value) || length(value) != 1) {
    stop(sprintf("%s must be a single logical value", name))
  }
}

#' Bayesian Spike Inference
#'
#' Infer spikes using Bayesian methods with MCMC sampling.
#'
#' @param trace Calcium trace
#' @param model_type Model type ("poisson", "gaussian", "hierarchical")
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param sampling_rate Optional sampling rate (must be positive)
#' @param prior_spike_rate Optional prior spike rate (must be between 0 and 1)
#' @param ... Additional arguments
#' @return List containing MCMC samples and inference results
#' @export
bayesian_spike_inference <- function(trace, model_type = "poisson", 
                                   n_samples = 1000, burnin = 100, 
                                   sampling_rate = NULL, prior_spike_rate = NULL, ...) {
  message("Running Bayesian spike inference")
  
  # Validate inputs
  validate_calcium_trace(trace)
  
  if (!is.character(model_type) || length(model_type) != 1) {
    stop("model_type must be a single character string")
  }
  
  if (!is.numeric(n_samples) || n_samples <= 0) {
    stop("n_samples must be a positive integer")
  }
  
  if (!is.numeric(burnin) || burnin < 0) {
    stop("burnin must be a non-negative integer")
  }
  
  if (!is.null(sampling_rate) && (!is.numeric(sampling_rate) || sampling_rate <= 0)) {
    stop("sampling_rate must be a positive numeric value")
  }
  
  if (!is.null(prior_spike_rate)) {
    validate_probability(prior_spike_rate, "prior_spike_rate")
  }
  
  # Base R implementation of Bayesian spike inference
  # This uses MCMC methods to sample from posterior distributions
  
  if (model_type == "poisson") {
    return(poisson_spike_model(trace, n_samples, burnin, ...))
  } else if (model_type == "gaussian") {
    return(gaussian_spike_model(trace, n_samples, burnin, ...))
  } else if (model_type == "hierarchical") {
    return(hierarchical_spike_model(trace, n_samples, burnin, ...))
  } else {
    stop("Unknown model type: ", model_type)
  }
}

#' Poisson Spike Model
#'
#' Bayesian spike inference using Poisson model.
#'
#' @param trace Calcium trace
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return MCMC results
#' @keywords internal
poisson_spike_model <- function(trace, n_samples, burnin, ...) {
  n <- length(trace)
  
  # Initialize parameters
  lambda <- mean(trace, na.rm = TRUE)  # Baseline firing rate
  spikes <- rep(0, n)  # Spike indicators
  tau <- 10  # Decay time constant
  
  # Storage for MCMC samples
  lambda_samples <- numeric(n_samples)
  spike_samples <- matrix(0, n_samples, n)
  tau_samples <- numeric(n_samples)
  
  # MCMC sampling
  for (iter in 1:(burnin + n_samples)) {
    # Sample lambda (baseline rate)
    lambda <- rgamma(1, shape = 1 + sum(spikes), rate = 1 + n)
    
    # Sample tau (decay time constant)
    tau <- rgamma(1, shape = 1 + n/2, rate = 1 + sum((trace - lambda)^2)/2)
    
    # Sample spikes using Metropolis-Hastings
    for (i in 1:n) {
      # Propose new spike state
      spike_prop <- 1 - spikes[i]  # Flip current state
      
      # Calculate acceptance ratio
      if (spike_prop == 1) {
        # Proposing to add spike
        log_ratio <- log(lambda) - log(1 - lambda) + 
                    (trace[i] - lambda)^2 / (2 * tau) - 
                    (trace[i] - lambda - 1)^2 / (2 * tau)
      } else {
        # Proposing to remove spike
        log_ratio <- log(1 - lambda) - log(lambda) + 
                    (trace[i] - lambda - 1)^2 / (2 * tau) - 
                    (trace[i] - lambda)^2 / (2 * tau)
      }
      
      # Accept or reject
      if (log(runif(1)) < log_ratio) {
        spikes[i] <- spike_prop
      }
    }
    
    # Store samples after burn-in
    if (iter > burnin) {
      sample_idx <- iter - burnin
      lambda_samples[sample_idx] <- lambda
      spike_samples[sample_idx, ] <- spikes
      tau_samples[sample_idx] <- tau
    }
  }
  
  # Compute posterior summaries
  posterior_spikes <- colMeans(spike_samples)
  posterior_lambda <- mean(lambda_samples)
  posterior_tau <- mean(tau_samples)
  
  return(list(
    lambda_samples = lambda_samples,
    spike_samples = spike_samples,
    tau_samples = tau_samples,
    posterior_spikes = posterior_spikes,
    posterior_lambda = posterior_lambda,
    posterior_tau = posterior_tau,
    model_type = "poisson",
    parameters = list(n_samples = n_samples, burnin = burnin)
  ))
}

#' Gaussian Spike Model
#'
#' Bayesian spike inference using Gaussian model.
#'
#' @param trace Calcium trace
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return MCMC results
#' @keywords internal
gaussian_spike_model <- function(trace, n_samples, burnin, ...) {
  n <- length(trace)
  
  # Initialize parameters
  mu <- mean(trace, na.rm = TRUE)  # Baseline calcium
  sigma <- sd(trace, na.rm = TRUE)  # Noise standard deviation
  spikes <- rep(0, n)  # Spike indicators
  spike_amplitude <- 1.0  # Spike amplitude
  
  # Storage for MCMC samples
  mu_samples <- numeric(n_samples)
  sigma_samples <- numeric(n_samples)
  spike_samples <- matrix(0, n_samples, n)
  amplitude_samples <- numeric(n_samples)
  
  # MCMC sampling
  for (iter in 1:(burnin + n_samples)) {
    # Sample mu (baseline calcium)
    mu <- rnorm(1, mean = mean(trace - spike_amplitude * spikes), 
                sd = sigma / sqrt(n))
    
    # Sample sigma (noise standard deviation)
    ss <- sum((trace - mu - spike_amplitude * spikes)^2)
    sigma <- sqrt(1 / rgamma(1, shape = n/2, rate = ss/2))
    
    # Sample spike amplitude
    spike_sum <- sum(spikes)
    if (spike_sum > 0) {
      spike_mean <- sum(spikes * (trace - mu)) / spike_sum
      spike_var <- sigma^2 / spike_sum
      spike_amplitude <- rnorm(1, mean = spike_mean, sd = sqrt(spike_var))
    }
    
    # Sample spikes using Metropolis-Hastings
    for (i in 1:n) {
      # Propose new spike state
      spike_prop <- 1 - spikes[i]  # Flip current state
      
      # Calculate acceptance ratio
      if (spike_prop == 1) {
        # Proposing to add spike
        log_ratio <- log(0.1) - log(0.9) + 
                    (trace[i] - mu)^2 / (2 * sigma^2) - 
                    (trace[i] - mu - spike_amplitude)^2 / (2 * sigma^2)
      } else {
        # Proposing to remove spike
        log_ratio <- log(0.9) - log(0.1) + 
                    (trace[i] - mu - spike_amplitude)^2 / (2 * sigma^2) - 
                    (trace[i] - mu)^2 / (2 * sigma^2)
      }
      
      # Accept or reject
      if (log(runif(1)) < log_ratio) {
        spikes[i] <- spike_prop
      }
    }
    
    # Store samples after burn-in
    if (iter > burnin) {
      sample_idx <- iter - burnin
      mu_samples[sample_idx] <- mu
      sigma_samples[sample_idx] <- sigma
      spike_samples[sample_idx, ] <- spikes
      amplitude_samples[sample_idx] <- spike_amplitude
    }
  }
  
  # Compute posterior summaries
  posterior_spikes <- colMeans(spike_samples)
  posterior_mu <- mean(mu_samples)
  posterior_sigma <- mean(sigma_samples)
  posterior_amplitude <- mean(amplitude_samples)
  
  return(list(
    mu_samples = mu_samples,
    sigma_samples = sigma_samples,
    spike_samples = spike_samples,
    amplitude_samples = amplitude_samples,
    posterior_spikes = posterior_spikes,
    posterior_mu = posterior_mu,
    posterior_sigma = posterior_sigma,
    posterior_amplitude = posterior_amplitude,
    model_type = "gaussian",
    parameters = list(n_samples = n_samples, burnin = burnin)
  ))
}

#' Hierarchical Spike Model
#'
#' Bayesian spike inference using hierarchical model.
#'
#' @param trace Calcium trace
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return MCMC results
#' @keywords internal
hierarchical_spike_model <- function(trace, n_samples, burnin, ...) {
  n <- length(trace)
  
  # Initialize parameters
  mu <- mean(trace, na.rm = TRUE)
  sigma <- sd(trace, na.rm = TRUE)
  spikes <- rep(0, n)
  spike_amplitude <- 1.0
  
  # Hyperparameters
  mu_prior_mean <- 0
  mu_prior_sd <- 10
  sigma_prior_shape <- 1
  sigma_prior_rate <- 1
  
  # Storage for MCMC samples
  mu_samples <- numeric(n_samples)
  sigma_samples <- numeric(n_samples)
  spike_samples <- matrix(0, n_samples, n)
  amplitude_samples <- numeric(n_samples)
  
  # MCMC sampling
  for (iter in 1:(burnin + n_samples)) {
    # Sample mu with hierarchical prior
    mu_precision <- 1/mu_prior_sd^2 + n/sigma^2
    mu_mean <- (mu_prior_mean/mu_prior_sd^2 + sum(trace - spike_amplitude * spikes)/sigma^2) / mu_precision
    mu <- rnorm(1, mean = mu_mean, sd = sqrt(1/mu_precision))
    
    # Sample sigma with hierarchical prior
    ss <- sum((trace - mu - spike_amplitude * spikes)^2)
    sigma <- sqrt(1 / rgamma(1, shape = sigma_prior_shape + n/2, 
                            rate = sigma_prior_rate + ss/2))
    
    # Sample spike amplitude
    spike_sum <- sum(spikes)
    if (spike_sum > 0) {
      spike_mean <- sum(spikes * (trace - mu)) / spike_sum
      spike_var <- sigma^2 / spike_sum
      spike_amplitude <- rnorm(1, mean = spike_mean, sd = sqrt(spike_var))
    }
    
    # Sample spikes using Metropolis-Hastings
    for (i in 1:n) {
      # Propose new spike state
      spike_prop <- 1 - spikes[i]
      
      # Calculate acceptance ratio with hierarchical prior
      if (spike_prop == 1) {
        log_ratio <- log(0.1) - log(0.9) + 
                    (trace[i] - mu)^2 / (2 * sigma^2) - 
                    (trace[i] - mu - spike_amplitude)^2 / (2 * sigma^2)
      } else {
        log_ratio <- log(0.9) - log(0.1) + 
                    (trace[i] - mu - spike_amplitude)^2 / (2 * sigma^2) - 
                    (trace[i] - mu)^2 / (2 * sigma^2)
      }
      
      # Accept or reject
      if (log(runif(1)) < log_ratio) {
        spikes[i] <- spike_prop
      }
    }
    
    # Store samples after burn-in
    if (iter > burnin) {
      sample_idx <- iter - burnin
      mu_samples[sample_idx] <- mu
      sigma_samples[sample_idx] <- sigma
      spike_samples[sample_idx, ] <- spikes
      amplitude_samples[sample_idx] <- spike_amplitude
    }
  }
  
  # Compute posterior summaries
  posterior_spikes <- colMeans(spike_samples)
  posterior_mu <- mean(mu_samples)
  posterior_sigma <- mean(sigma_samples)
  posterior_amplitude <- mean(amplitude_samples)
  
  return(list(
    mu_samples = mu_samples,
    sigma_samples = sigma_samples,
    spike_samples = spike_samples,
    amplitude_samples = amplitude_samples,
    posterior_spikes = posterior_spikes,
    posterior_mu = posterior_mu,
    posterior_sigma = posterior_sigma,
    posterior_amplitude = posterior_amplitude,
    model_type = "hierarchical",
    parameters = list(n_samples = n_samples, burnin = burnin)
  ))
}

#' Bayesian Model Comparison
#'
#' Compare different Bayesian models using model selection criteria.
#'
#' @param model_results List of model results to compare
#' @param comparison_metrics Optional character vector of metrics
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return Model comparison results
#' @export
bayesian_model_comparison <- function(model_results, comparison_metrics = NULL, 
                                    n_samples = 1000, burnin = 100, ...) {
  message("Comparing Bayesian models")
  
  # Validate model_results argument first to match test expectations
  if (!is.list(model_results) || length(model_results) < 2) {
    stop("model_results must be a list with at least 2 models")
  }
  
  if (!is.null(comparison_metrics) && !is.character(comparison_metrics)) {
    stop("comparison_metrics must be a character vector")
  }
  
  if (!is.numeric(n_samples) || n_samples <= 0 || n_samples != as.integer(n_samples)) {
    stop("n_samples must be a positive integer")
  }
  
  if (!is.numeric(burnin) || burnin < 0 || burnin != as.integer(burnin)) {
    stop("burnin must be a non-negative integer")
  }
  
  # Extract model names and validate model results
  model_names <- names(model_results)
  if (is.null(model_names)) {
    model_names <- paste0("model", seq_along(model_results))
  }
  
  # Validate that each model result has required components
  for (i in seq_along(model_results)) {
    if (!is.list(model_results[[i]]) || !("spikes" %in% names(model_results[[i]]))) {
      stop("Each model result must be a list with a 'spikes' component")
    }
  }
  
  # Compute model selection criteria
  comparison <- data.frame(
    model = model_names,
    log_likelihood = sapply(model_results, function(x) {
      if ("log_likelihood" %in% names(x)) x$log_likelihood else NA
    }),
    aic = sapply(model_results, function(x) {
      if ("aic" %in% names(x)) x$aic else NA
    }),
    bic = sapply(model_results, function(x) {
      if ("bic" %in% names(x)) x$bic else NA
    }),
    stringsAsFactors = FALSE
  )
  
  # Add DIC (Deviance Information Criterion) if available
  comparison$dic <- sapply(model_results, function(x) {
    if ("dic" %in% names(x)) x$dic else NA
  })
  
  # Determine best model based on available criteria
  best_model <- NULL
  if (!all(is.na(comparison$dic))) {
    best_model <- comparison$model[which.min(comparison$dic)]
  } else if (!all(is.na(comparison$aic))) {
    best_model <- comparison$model[which.min(comparison$aic)]
  } else if (!all(is.na(comparison$bic))) {
    best_model <- comparison$model[which.min(comparison$bic)]
  }
  
  return(list(
    models = model_results,
    comparison = comparison,
    best_model = best_model
  ))
}

#' Bayesian Parameter Estimation
#'
#' Estimate parameters using Bayesian methods.
#'
#' @param data Input data
#' @param model Model specification
#' @param priors Prior distributions
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param model_type Optional model type (must be character)
#' @param ... Additional arguments
#' @return Parameter estimation results
#' @export
bayesian_parameter_estimation <- function(data, model = "gaussian", 
                                        priors = NULL, n_samples = 1000, 
                                        burnin = 100, model_type = NULL, ...) {
  message("Running Bayesian parameter estimation")
  
  # Validate inputs
  validate_calcium_trace(data)
  
  if (!is.character(model) || length(model) != 1) {
    stop("model must be a single character string")
  }
  
  if (!is.null(model_type) && (!is.character(model_type) || length(model_type) != 1)) {
    stop("model_type must be a single character string")
  }
  
  if (!is.numeric(n_samples) || n_samples <= 0) {
    stop("n_samples must be a positive integer")
  }
  
  if (!is.numeric(burnin) || burnin < 0) {
    stop("burnin must be a non-negative integer")
  }
  
  # Base R implementation of Bayesian parameter estimation
  # This uses MCMC methods for parameter inference
  
  if (model == "gaussian") {
    return(gaussian_parameter_estimation(data, priors, n_samples, burnin, ...))
  } else if (model == "poisson") {
    return(poisson_parameter_estimation(data, priors, n_samples, burnin, ...))
  } else {
    stop("Unknown model: ", model)
  }
}

#' Gaussian Parameter Estimation
#'
#' Estimate Gaussian model parameters using MCMC.
#'
#' @param data Input data
#' @param priors Prior specifications
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return Parameter estimation results
#' @keywords internal
gaussian_parameter_estimation <- function(data, priors, n_samples, burnin, ...) {
  n <- length(data)
  
  # Set default priors if not provided
  if (is.null(priors)) {
    priors <- list(
      mu_mean = 0,
      mu_sd = 10,
      sigma_shape = 1,
      sigma_rate = 1
    )
  }
  
  # Initialize parameters
  mu <- mean(data, na.rm = TRUE)
  sigma <- sd(data, na.rm = TRUE)
  
  # Storage for MCMC samples
  mu_samples <- numeric(n_samples)
  sigma_samples <- numeric(n_samples)
  
  # MCMC sampling
  for (iter in 1:(burnin + n_samples)) {
    # Sample mu
    mu_precision <- 1/priors$mu_sd^2 + n/sigma^2
    mu_mean <- (priors$mu_mean/priors$mu_sd^2 + sum(data)/sigma^2) / mu_precision
    mu <- rnorm(1, mean = mu_mean, sd = sqrt(1/mu_precision))
    
    # Sample sigma
    ss <- sum((data - mu)^2)
    sigma <- sqrt(1 / rgamma(1, shape = priors$sigma_shape + n/2, 
                            rate = priors$sigma_rate + ss/2))
    
    # Store samples after burn-in
    if (iter > burnin) {
      sample_idx <- iter - burnin
      mu_samples[sample_idx] <- mu
      sigma_samples[sample_idx] <- sigma
    }
  }
  
  # Compute posterior summaries
  posterior_summary <- list(
    mu = list(
      mean = mean(mu_samples),
      sd = sd(mu_samples),
      quantiles = quantile(mu_samples, c(0.025, 0.25, 0.5, 0.75, 0.975))
    ),
    sigma = list(
      mean = mean(sigma_samples),
      sd = sd(sigma_samples),
      quantiles = quantile(sigma_samples, c(0.025, 0.25, 0.5, 0.75, 0.975))
    )
  )
  
  return(list(
    mu_samples = mu_samples,
    sigma_samples = sigma_samples,
    posterior_summary = posterior_summary,
    model = "gaussian",
    parameters = list(n_samples = n_samples, burnin = burnin)
  ))
}

#' Poisson Parameter Estimation
#'
#' Estimate Poisson model parameters using MCMC.
#'
#' @param data Input data
#' @param priors Prior specifications
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return Parameter estimation results
#' @keywords internal
poisson_parameter_estimation <- function(data, priors, n_samples, burnin, ...) {
  n <- length(data)
  
  # Set default priors if not provided
  if (is.null(priors)) {
    priors <- list(
      lambda_shape = 1,
      lambda_rate = 1
    )
  }
  
  # Initialize parameters
  lambda <- mean(data, na.rm = TRUE)
  
  # Storage for MCMC samples
  lambda_samples <- numeric(n_samples)
  
  # MCMC sampling
  for (iter in 1:(burnin + n_samples)) {
    # Sample lambda
    lambda <- rgamma(1, shape = priors$lambda_shape + sum(data), 
                    rate = priors$lambda_rate + n)
    
    # Store samples after burn-in
    if (iter > burnin) {
      sample_idx <- iter - burnin
      lambda_samples[sample_idx] <- lambda
    }
  }
  
  # Compute posterior summaries
  posterior_summary <- list(
    lambda = list(
      mean = mean(lambda_samples),
      sd = sd(lambda_samples),
      quantiles = quantile(lambda_samples, c(0.025, 0.25, 0.5, 0.75, 0.975))
    )
  )
  
  return(list(
    lambda_samples = lambda_samples,
    posterior_summary = posterior_summary,
    model = "poisson",
    parameters = list(n_samples = n_samples, burnin = burnin)
  ))
}

#' Bayesian Hypothesis Testing
#'
#' Perform Bayesian hypothesis testing.
#'
#' @param data1 First dataset
#' @param data2 Second dataset
#' @param test_type Type of test ("mean", "variance", "proportion")
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return Hypothesis testing results
#' @export
bayesian_hypothesis_testing <- function(data1, data2, test_type = "mean", 
                                      n_samples = 1000, burnin = 100, ...) {
  message("Running Bayesian hypothesis testing")
  
  if (test_type == "mean") {
    return(bayesian_mean_test(data1, data2, n_samples, burnin, ...))
  } else if (test_type == "variance") {
    return(bayesian_variance_test(data1, data2, n_samples, burnin, ...))
  } else if (test_type == "proportion") {
    return(bayesian_proportion_test(data1, data2, n_samples, burnin, ...))
  } else {
    stop("Unknown test type: ", test_type)
  }
}

#' Bayesian Mean Test
#'
#' Test for difference in means using Bayesian methods.
#'
#' @param data1 First dataset
#' @param data2 Second dataset
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return Mean test results
#' @keywords internal
bayesian_mean_test <- function(data1, data2, n_samples, burnin, ...) {
  # Fit models to both datasets
  fit1 <- gaussian_parameter_estimation(data1, NULL, n_samples, burnin)
  fit2 <- gaussian_parameter_estimation(data2, NULL, n_samples, burnin)
  
  # Compute difference in means
  mean_diff <- fit1$mu_samples - fit2$mu_samples
  
  # Compute posterior probability of difference
  prob_positive <- mean(mean_diff > 0)
  prob_negative <- mean(mean_diff < 0)
  
  # Compute credible interval
  credible_interval <- quantile(mean_diff, c(0.025, 0.975))
  
  return(list(
    mean_diff_samples = mean_diff,
    posterior_mean = mean(mean_diff),
    posterior_sd = sd(mean_diff),
    prob_positive = prob_positive,
    prob_negative = prob_negative,
    credible_interval = credible_interval,
    test_type = "mean"
  ))
}

#' Bayesian Variance Test
#'
#' Test for difference in variances using Bayesian methods.
#'
#' @param data1 First dataset
#' @param data2 Second dataset
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return Variance test results
#' @keywords internal
bayesian_variance_test <- function(data1, data2, n_samples, burnin, ...) {
  # Fit models to both datasets
  fit1 <- gaussian_parameter_estimation(data1, NULL, n_samples, burnin)
  fit2 <- gaussian_parameter_estimation(data2, NULL, n_samples, burnin)
  
  # Compute ratio of variances
  var_ratio <- (fit1$sigma_samples^2) / (fit2$sigma_samples^2)
  
  # Compute posterior probability
  prob_ratio_gt_1 <- mean(var_ratio > 1)
  prob_ratio_lt_1 <- mean(var_ratio < 1)
  
  # Compute credible interval
  credible_interval <- quantile(var_ratio, c(0.025, 0.975))
  
  return(list(
    var_ratio_samples = var_ratio,
    posterior_mean = mean(var_ratio),
    posterior_sd = sd(var_ratio),
    prob_ratio_gt_1 = prob_ratio_gt_1,
    prob_ratio_lt_1 = prob_ratio_lt_1,
    credible_interval = credible_interval,
    test_type = "variance"
  ))
}

#' Bayesian Proportion Test
#'
#' Test for difference in proportions using Bayesian methods.
#'
#' @param data1 First dataset (binary)
#' @param data2 Second dataset (binary)
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return Proportion test results
#' @keywords internal
bayesian_proportion_test <- function(data1, data2, n_samples, burnin, ...) {
  # Fit models to both datasets
  fit1 <- poisson_parameter_estimation(data1, NULL, n_samples, burnin)
  fit2 <- poisson_parameter_estimation(data2, NULL, n_samples, burnin)
  
  # Compute difference in proportions
  prop_diff <- fit1$lambda_samples - fit2$lambda_samples
  
  # Compute posterior probability
  prob_positive <- mean(prop_diff > 0)
  prob_negative <- mean(prop_diff < 0)
  
  # Compute credible interval
  credible_interval <- quantile(prop_diff, c(0.025, 0.975))
  
  return(list(
    prop_diff_samples = prop_diff,
    posterior_mean = mean(prop_diff),
    posterior_sd = sd(prop_diff),
    prob_positive = prob_positive,
    prob_negative = prob_negative,
    credible_interval = credible_interval,
    test_type = "proportion"
  ))
}

#' Helper Functions
#' @keywords internal

compute_log_likelihood <- function(trace, model_result) {
  # Compute log-likelihood for model comparison
  if (model_result$model_type == "poisson") {
    lambda <- model_result$posterior_lambda
    return(sum(dpois(trace, lambda, log = TRUE)))
  } else if (model_result$model_type == "gaussian") {
    mu <- model_result$posterior_mu
    sigma <- model_result$posterior_sigma
    return(sum(dnorm(trace, mu, sigma, log = TRUE)))
  } else {
    return(0)  # Placeholder
  }
}

compute_aic <- function(trace, model_result) {
  # Compute AIC
  log_lik <- compute_log_likelihood(trace, model_result)
  n_params <- length(model_result$posterior_spikes) + 2  # Simplified
  return(2 * n_params - 2 * log_lik)
}

compute_bic <- function(trace, model_result) {
  # Compute BIC
  log_lik <- compute_log_likelihood(trace, model_result)
  n_params <- length(model_result$posterior_spikes) + 2  # Simplified
  n_obs <- length(trace)
  return(log(n_obs) * n_params - 2 * log_lik)
}

compute_dic <- function(trace, model_result) {
  # Compute DIC (Deviance Information Criterion)
  # This is a simplified version
  log_lik <- compute_log_likelihood(trace, model_result)
  return(-2 * log_lik)  # Simplified DIC
}

#' Uncertainty Quantification
#'
#' Quantify uncertainty in Bayesian spike inference results.
#'
#' @param bayesian_result Result from bayesian_spike_inference
#' @param confidence_level Confidence level for intervals
#' @param method Method for uncertainty quantification
#' @param ... Additional arguments
#' @return Uncertainty quantification results
#' @export
uncertainty_quantification <- function(bayesian_result, confidence_level = 0.95, 
                                     method = "credible_interval", ...) {
  if (is.null(bayesian_result) || !is.list(bayesian_result)) {
    stop("bayesian_result must be from bayesian_spike_inference")
  }
  
  if (!is.numeric(confidence_level) || confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be a probability between 0 and 1")
  }
  
  if (!is.character(method) || length(method) != 1) {
    stop("method must be a single character string")
  }
  
  message("Quantifying uncertainty in Bayesian results")
  
  # Extract samples
  if ("spike_samples" %in% names(bayesian_result)) {
    spike_samples <- bayesian_result$spike_samples
  } else {
    stop("Invalid bayesian_result object")
  }
  
  # Compute credible intervals
  alpha <- 1 - confidence_level
  lower_quantile <- alpha / 2
  upper_quantile <- 1 - alpha / 2
  
  credible_intervals <- t(apply(spike_samples, 2, quantile, 
                               probs = c(lower_quantile, upper_quantile)))
  
  # Compute posterior means and standard deviations
  posterior_means <- colMeans(spike_samples)
  posterior_sds <- apply(spike_samples, 2, sd)
  
  return(list(
    credible_intervals = credible_intervals,
    posterior_means = posterior_means,
    posterior_sds = posterior_sds,
    confidence_level = confidence_level,
    method = method
  ))
}

#' Probabilistic Spike Detection
#'
#' Detect spikes using probabilistic methods.
#'
#' @param trace Calcium trace
#' @param method Detection method
#' @param min_spike_interval Minimum interval between spikes
#' @param uncertainty_estimation Whether to estimate uncertainty
#' @param ... Additional arguments
#' @return Probabilistic spike detection results
#' @export
probabilistic_spike_detection <- function(trace, method = "threshold", 
                                        min_spike_interval = 5, 
                                        uncertainty_estimation = TRUE, ...) {
  validate_calcium_trace(trace)
  
  if (!is.character(method) || length(method) != 1) {
    stop("method must be a single character string")
  }
  
  if (!is.numeric(min_spike_interval) || min_spike_interval <= 0) {
    stop("min_spike_interval must be a positive integer")
  }
  
  if (!is.logical(uncertainty_estimation) || length(uncertainty_estimation) != 1) {
    stop("uncertainty_estimation must be a single logical value")
  }
  
  message("Running probabilistic spike detection")
  
  # Simple threshold-based detection
  threshold <- mean(trace) + 2 * sd(trace)
  spike_probs <- ifelse(trace > threshold, 0.8, 0.1)
  
  # Apply minimum interval constraint
  if (min_spike_interval > 1) {
    for (i in 1:length(spike_probs)) {
      if (spike_probs[i] > 0.5) {
        # Suppress nearby spikes
        start_idx <- max(1, i - min_spike_interval)
        end_idx <- min(length(spike_probs), i + min_spike_interval)
        spike_probs[start_idx:end_idx] <- pmin(spike_probs[start_idx:end_idx], 0.3)
        spike_probs[i] <- 0.8  # Keep current spike
      }
    }
  }
  
  # Estimate uncertainty if requested
  uncertainty <- NULL
  if (uncertainty_estimation) {
    uncertainty <- list(
      standard_error = sd(spike_probs) / sqrt(length(spike_probs)),
      confidence_interval = quantile(spike_probs, c(0.025, 0.975))
    )
  }
  
  return(list(
    spike_probabilities = spike_probs,
    threshold = threshold,
    method = method,
    uncertainty = uncertainty
  ))
}

#' Hierarchical Bayesian Modeling
#'
#' Perform hierarchical Bayesian modeling on multiple calcium traces.
#'
#' @param calcium_traces Matrix of calcium traces (rows = cells, cols = time)
#' @param model_type Type of hierarchical model
#' @param n_samples Number of MCMC samples
#' @param burnin Number of burn-in samples
#' @param ... Additional arguments
#' @return Hierarchical modeling results
#' @export
hierarchical_bayesian_modeling <- function(calcium_traces, model_type = "gaussian", 
                                         n_samples = 1000, burnin = 100, ...) {
  validate_calcium_traces_matrix(calcium_traces)
  
  if (!is.character(model_type) || length(model_type) != 1) {
    stop("model_type must be a single character string")
  }
  
  if (!is.numeric(n_samples) || n_samples <= 0) {
    stop("n_samples must be a positive integer")
  }
  
  message("Running hierarchical Bayesian modeling")
  
  n_cells <- nrow(calcium_traces)
  n_timepoints <- ncol(calcium_traces)
  
  # Fit individual models
  individual_results <- list()
  for (i in 1:n_cells) {
    individual_results[[i]] <- bayesian_spike_inference(
      calcium_traces[i, ], model_type, n_samples, burnin, ...
    )
  }
  
  # Compute population-level statistics
  population_spikes <- matrix(0, n_cells, n_timepoints)
  for (i in 1:n_cells) {
    population_spikes[i, ] <- individual_results[[i]]$posterior_spikes
  }
  
  # Population mean and variance
  population_mean <- colMeans(population_spikes)
  population_variance <- apply(population_spikes, 2, var)
  
  return(list(
    individual_results = individual_results,
    population_spikes = population_spikes,
    population_mean = population_mean,
    population_variance = population_variance,
    model_type = model_type,
    n_cells = n_cells,
    n_timepoints = n_timepoints
  ))
}

#' Plot Bayesian Results
#'
#' Create plots for Bayesian analysis results.
#'
#' @param bayesian_result Result from Bayesian analysis
#' @param plot_type Type of plot to create
#' @param confidence_level Confidence level for intervals
#' @param ... Additional arguments
#' @return Plot object
#' @export
plot_bayesian_results <- function(bayesian_result, plot_type = "spikes", 
                                confidence_level = 0.95, ...) {
  if (is.null(bayesian_result) || !is.list(bayesian_result)) {
    stop("Invalid bayesian_result object")
  }
  
  if (!is.character(plot_type) || length(plot_type) != 1) {
    stop("plot_type must be a single character string")
  }
  
  if (!is.numeric(confidence_level) || confidence_level <= 0 || confidence_level >= 1) {
    stop("confidence_level must be a probability between 0 and 1")
  }
  
  message("Creating Bayesian results plot")
  
  if (plot_type == "spikes") {
    # Plot spike probabilities
    if ("posterior_spikes" %in% names(bayesian_result)) {
      plot(bayesian_result$posterior_spikes, type = "l", 
           xlab = "Time", ylab = "Spike Probability",
           main = "Posterior Spike Probabilities")
    }
  } else if (plot_type == "parameters") {
    # Plot parameter distributions
    if ("mu_samples" %in% names(bayesian_result)) {
      hist(bayesian_result$mu_samples, main = "Posterior Distribution of μ",
           xlab = "μ", ylab = "Frequency")
    }
  } else if (plot_type == "trace") {
    # Plot MCMC trace
    if ("lambda_samples" %in% names(bayesian_result)) {
      plot(bayesian_result$lambda_samples, type = "l",
           xlab = "MCMC Iteration", ylab = "λ",
           main = "MCMC Trace for λ")
    }
  }
  
  return(invisible(NULL))
}

#' Bayesian Predictive Checks
#'
#' Perform Bayesian predictive checks for model validation.
#'
#' @param bayesian_result Result from Bayesian analysis
#' @param n_predictions Number of predictions to generate
#' @param test_statistics Statistics to compute
#' @param ... Additional arguments
#' @return Predictive check results
#' @export
bayesian_predictive_checks <- function(bayesian_result, n_predictions = 100, 
                                     test_statistics = c("mean", "variance"), ...) {
  if (is.null(bayesian_result) || !is.list(bayesian_result)) {
    stop("Invalid bayesian_result object")
  }
  
  if (!is.numeric(n_predictions) || n_predictions <= 0) {
    stop("n_predictions must be a positive integer")
  }
  
  if (!is.character(test_statistics)) {
    stop("test_statistics must be a character vector")
  }
  
  message("Running Bayesian predictive checks")
  
  # Generate predictions based on posterior samples
  predictions <- matrix(0, n_predictions, length(bayesian_result$posterior_spikes))
  
  if ("mu_samples" %in% names(bayesian_result) && "sigma_samples" %in% names(bayesian_result)) {
    # Gaussian model predictions
    for (i in 1:n_predictions) {
      mu_idx <- sample(length(bayesian_result$mu_samples), 1)
      sigma_idx <- sample(length(bayesian_result$sigma_samples), 1)
      
      predictions[i, ] <- rnorm(length(bayesian_result$posterior_spikes),
                               mean = bayesian_result$mu_samples[mu_idx],
                               sd = bayesian_result$sigma_samples[sigma_idx])
    }
  } else if ("lambda_samples" %in% names(bayesian_result)) {
    # Poisson model predictions
    for (i in 1:n_predictions) {
      lambda_idx <- sample(length(bayesian_result$lambda_samples), 1)
      predictions[i, ] <- rpois(length(bayesian_result$posterior_spikes),
                               lambda = bayesian_result$lambda_samples[lambda_idx])
    }
  }
  
  # Compute test statistics
  observed_stats <- list()
  predicted_stats <- list()
  
  for (stat in test_statistics) {
    if (stat == "mean") {
      observed_stats[[stat]] <- mean(bayesian_result$posterior_spikes)
      predicted_stats[[stat]] <- apply(predictions, 1, mean)
    } else if (stat == "variance") {
      observed_stats[[stat]] <- var(bayesian_result$posterior_spikes)
      predicted_stats[[stat]] <- apply(predictions, 1, var)
    }
  }
  
  # Compute p-values
  p_values <- list()
  for (stat in test_statistics) {
    p_values[[stat]] <- mean(predicted_stats[[stat]] >= observed_stats[[stat]])
  }
  
  return(list(
    predictions = predictions,
    observed_statistics = observed_stats,
    predicted_statistics = predicted_stats,
    p_values = p_values,
    test_statistics = test_statistics
  ))
}

#' Posterior Summary Statistics
#'
#' Compute summary statistics for posterior samples.
#'
#' @param samples Matrix or vector of posterior samples
#' @param quantiles Quantiles to compute (default: c(0.025, 0.5, 0.975))
#' @param ... Additional arguments
#' @return List with summary statistics
#' @export
posterior_summary <- function(samples, quantiles = c(0.025, 0.5, 0.975), ...) {
  # Flatten all values if matrix
  if (is.matrix(samples) || is.data.frame(samples)) {
    samples <- as.numeric(samples)
  }
  
  # Remove NA and infinite values
  samples <- samples[is.finite(samples)]
  
  if (length(samples) == 0) {
    return(list(
      mean = NA,
      sd = NA,
      lower = NA,
      upper = NA,
      median = NA
    ))
  }
  
  q_values <- quantile(samples, probs = quantiles, na.rm = TRUE)
  
  return(list(
    mean = mean(samples, na.rm = TRUE),
    sd = sd(samples, na.rm = TRUE),
    lower = q_values[1],
    upper = q_values[length(q_values)],
    median = q_values[2]
  ))
} 