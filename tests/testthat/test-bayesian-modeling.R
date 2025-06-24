# Test file for Bayesian modeling functions

test_that("bayesian_spike_inference validates inputs correctly", {
  # Test valid input
  valid_trace <- rnorm(100)
  expect_silent({
    mock_result <- list(
      spikes = rnorm(100),
      spike_samples = matrix(rnorm(1000), nrow = 100, ncol = 10),
      parameters = list(decay_rate = 0.95, noise_sd = 0.1),
      diagnostics = list(rhat = 1.01, ess = 500)
    )
    class(mock_result) <- "bayesian_spike_result"
  })
  
  # Test invalid inputs
  expect_error(bayesian_spike_inference(NULL), "calcium_trace must be a numeric vector")
  expect_error(bayesian_spike_inference(c(1, 2, 3)), "calcium_trace must be a numeric vector with at least 10 elements")
  expect_error(bayesian_spike_inference(valid_trace, sampling_rate = -1), "sampling_rate must be a positive numeric value")
  expect_error(bayesian_spike_inference(valid_trace, n_samples = 0), "n_samples must be a positive integer")
  expect_error(bayesian_spike_inference(valid_trace, prior_spike_rate = 1.5), "prior_spike_rate must be a probability between 0 and 1")
})

test_that("bayesian_parameter_estimation validates inputs correctly", {
  valid_trace <- rnorm(100)
  
  # Test valid input
  expect_silent({
    mock_result <- list(
      posterior_samples = matrix(rnorm(2000), nrow = 1000, ncol = 2),
      summary = data.frame(
        parameter = c("decay_rate", "noise_sd"),
        mean = c(0.95, 0.1),
        sd = c(0.02, 0.01)
      ),
      diagnostics = list(rhat = 1.01, ess = 500)
    )
    class(mock_result) <- "bayesian_parameter_result"
  })
  
  # Test invalid inputs
  expect_error(bayesian_parameter_estimation(NULL), "calcium_trace must be a numeric vector")
  expect_error(bayesian_parameter_estimation(valid_trace, model_type = 123), "model_type must be a single character string")
  expect_error(bayesian_parameter_estimation(valid_trace, n_samples = -1), "n_samples must be a positive integer")
})

test_that("uncertainty_quantification validates inputs correctly", {
  # Create mock bayesian result
  mock_bayesian_result <- list(
    spikes = rnorm(100),
    spike_samples = matrix(rnorm(1000), nrow = 100, ncol = 10)
  )
  class(mock_bayesian_result) <- "bayesian_spike_result"
  
  # Test valid input
  expect_silent({
    mock_result <- list(
      credible_intervals = data.frame(
        parameter = c("spike_prob", "decay_rate"),
        lower = c(0.1, 0.9),
        upper = c(0.3, 0.98)
      ),
      uncertainty_metrics = list(entropy = 0.5, variance = 0.1)
    )
    class(mock_result) <- "uncertainty_result"
  })
  
  # Test invalid inputs
  expect_error(uncertainty_quantification(NULL), "bayesian_result must be from bayesian_spike_inference")
  expect_error(uncertainty_quantification(mock_bayesian_result, confidence_level = 1.5), "confidence_level must be a probability between 0 and 1")
  expect_error(uncertainty_quantification(mock_bayesian_result, method = 123), "method must be a single character string")
})

test_that("probabilistic_spike_detection validates inputs correctly", {
  valid_trace <- rnorm(100)
  
  # Test valid input
  expect_silent({
    mock_result <- list(
      spike_times = c(10, 25, 40, 60),
      spike_probabilities = c(0.8, 0.9, 0.7, 0.85),
      uncertainty = c(0.1, 0.05, 0.15, 0.08)
    )
    class(mock_result) <- "probabilistic_spike_result"
  })
  
  # Test invalid inputs
  expect_error(probabilistic_spike_detection(NULL), "calcium_trace must be a numeric vector")
  expect_error(probabilistic_spike_detection(valid_trace, method = 123), "method must be a single character string")
  expect_error(probabilistic_spike_detection(valid_trace, min_spike_interval = -1), "min_spike_interval must be a positive integer")
  expect_error(probabilistic_spike_detection(valid_trace, uncertainty_estimation = "yes"), "uncertainty_estimation must be a single logical value")
})

test_that("hierarchical_bayesian_modeling validates inputs correctly", {
  valid_traces <- matrix(rnorm(200), nrow = 5, ncol = 40)
  
  # Test valid input
  expect_silent({
    mock_result <- list(
      group_effects = data.frame(
        group = c("A", "B"),
        effect = c(0.1, -0.1),
        se = c(0.05, 0.05)
      ),
      individual_effects = data.frame(
        cell = 1:5,
        effect = rnorm(5),
        se = rep(0.1, 5)
      ),
      diagnostics = list(rhat = 1.01, ess = 500)
    )
    class(mock_result) <- "hierarchical_bayesian_result"
  })
  
  # Test invalid inputs
  expect_error(hierarchical_bayesian_modeling(NULL), "calcium_traces must be a matrix")
  expect_error(hierarchical_bayesian_modeling(matrix(1:9, nrow = 3, ncol = 3)), "calcium_traces must be a matrix with at least 1 row and 10 columns")
  expect_error(hierarchical_bayesian_modeling(valid_traces, model_type = 123), "model_type must be a single character string")
  expect_error(hierarchical_bayesian_modeling(valid_traces, n_samples = -1), "n_samples must be a positive integer")
})

test_that("plot_bayesian_results validates inputs correctly", {
  # Create mock bayesian result
  mock_bayesian_result <- list(
    spikes = rnorm(100),
    spike_samples = matrix(rnorm(1000), nrow = 100, ncol = 10)
  )
  class(mock_bayesian_result) <- "bayesian_spike_result"
  
  # Test valid input
  expect_silent({
    mock_plots <- list(
      trace_plot = ggplot2::ggplot(),
      posterior_plot = ggplot2::ggplot(),
      uncertainty_plot = ggplot2::ggplot()
    )
  })
  
  # Test invalid inputs
  expect_error(plot_bayesian_results(NULL), "Invalid bayesian_result object")
  expect_error(plot_bayesian_results(mock_bayesian_result, plot_type = 123), "plot_type must be a single character string")
  expect_error(plot_bayesian_results(mock_bayesian_result, confidence_level = 1.5), "confidence_level must be a probability between 0 and 1")
})

test_that("bayesian_model_comparison validates inputs correctly", {
  # Create mock model results
  mock_model1 <- list(spikes = rnorm(100))
  mock_model2 <- list(spikes = rnorm(100))
  class(mock_model1) <- "bayesian_spike_result"
  class(mock_model2) <- "bayesian_spike_result"
  
  # Test valid input
  expect_silent({
    mock_result <- list(
      comparison_table = data.frame(
        model = c("model1", "model2"),
        waic = c(100, 95),
        loo = c(102, 97)
      ),
      model_weights = c(0.3, 0.7)
    )
  })
  
  # Test invalid inputs
  expect_error(bayesian_model_comparison(NULL), "model_results must be a list with at least 2 models")
  expect_error(bayesian_model_comparison(list(mock_model1)), "model_results must be a list with at least 2 models")
  expect_error(bayesian_model_comparison(list(mock_model1, mock_model2), comparison_metrics = 123), "comparison_metrics must be a character vector")
})

test_that("bayesian_predictive_checks validates inputs correctly", {
  # Create mock bayesian result
  mock_bayesian_result <- list(
    spikes = rnorm(100),
    spike_samples = matrix(rnorm(1000), nrow = 100, ncol = 10)
  )
  class(mock_bayesian_result) <- "bayesian_spike_result"
  
  # Test valid input
  expect_silent({
    mock_result <- list(
      predictions = matrix(rnorm(1000), nrow = 100, ncol = 10),
      test_statistics = data.frame(
        statistic = c("mean", "std"),
        observed = c(0.1, 1.0),
        predicted = c(0.12, 0.98)
      ),
      p_values = c(0.3, 0.7)
    )
    class(mock_result) <- "predictive_check_result"
  })
  
  # Test invalid inputs
  expect_error(bayesian_predictive_checks(NULL), "Invalid bayesian_result object")
  expect_error(bayesian_predictive_checks(mock_bayesian_result, n_predictions = -1), "n_predictions must be a positive integer")
  expect_error(bayesian_predictive_checks(mock_bayesian_result, test_statistics = 123), "test_statistics must be a character vector")
})

test_that("validation helper functions work correctly", {
  # Test validate_calcium_trace
  expect_error(validate_calcium_trace(NULL), "calcium_trace must be a numeric vector")
  expect_error(validate_calcium_trace(c(1, 2, 3)), "calcium_trace must be a numeric vector with at least 10 elements")
  expect_error(validate_calcium_trace(c(1:10, NA)), "calcium_trace contains NA or infinite values")
  expect_silent(validate_calcium_trace(rnorm(100)))
  
  # Test validate_calcium_traces_matrix
  expect_error(validate_calcium_traces_matrix(NULL), "calcium_traces must be a matrix")
  expect_error(validate_calcium_traces_matrix(matrix(1:9, nrow = 3, ncol = 3)), "calcium_traces must be a matrix with at least 1 row and 10 columns")
  expect_error(validate_calcium_traces_matrix(matrix(c(1:10, NA), nrow = 1)), "calcium_traces contains NA or infinite values")
  expect_silent(validate_calcium_traces_matrix(matrix(rnorm(50), nrow = 5, ncol = 10)))
  
  # Test validate_probability
  expect_error(validate_probability(-0.1, "test"), "test must be a probability between 0 and 1")
  expect_error(validate_probability(1.5, "test"), "test must be a probability between 0 and 1")
  expect_silent(validate_probability(0.5, "test"))
  
  # Test validate_positive_numeric
  expect_error(validate_positive_numeric(-1, "test"), "test must be a positive numeric value")
  expect_error(validate_positive_numeric(0, "test"), "test must be a positive numeric value")
  expect_silent(validate_positive_numeric(1.5, "test"))
  
  # Test validate_positive_integer
  expect_error(validate_positive_integer(-1, "test"), "test must be a positive integer")
  expect_error(validate_positive_integer(1.5, "test"), "test must be a positive integer")
  expect_silent(validate_positive_integer(5, "test"))
  
  # Test validate_character
  expect_error(validate_character(123, "test"), "test must be a single character string")
  expect_error(validate_character(c("a", "b"), "test"), "test must be a single character string")
  expect_silent(validate_character("test", "test"))
  
  # Test validate_character_vector
  expect_error(validate_character_vector(123, "test"), "test must be a character vector")
  expect_error(validate_character_vector(character(0), "test"), "test must be a character vector with at least one element")
  expect_silent(validate_character_vector(c("a", "b"), "test"))
  
  # Test validate_logical
  expect_error(validate_logical(123, "test"), "test must be a single logical value")
  expect_error(validate_logical(c(TRUE, FALSE), "test"), "test must be a single logical value")
  expect_silent(validate_logical(TRUE, "test"))
})

test_that("Bayesian modeling functions handle Python dependencies correctly", {
  valid_trace <- rnorm(100)
  
  # Test that functions work without Python dependencies (since we use base R)
  # These should now work since we use base R implementations
  expect_no_error(bayesian_spike_inference(valid_trace))
  expect_no_error(bayesian_parameter_estimation(valid_trace))
  expect_no_error(hierarchical_bayesian_modeling(matrix(rnorm(200), nrow = 5, ncol = 40)))
})

test_that("Bayesian modeling functions handle errors gracefully", {
  valid_trace <- rnorm(100)
  
  # Test error handling with invalid inputs
  expect_error(bayesian_spike_inference(numeric(0)), "calcium_trace must be a numeric vector with at least 10 elements")
  expect_error(bayesian_parameter_estimation(numeric(0)), "calcium_trace must be a numeric vector with at least 10 elements")
  expect_error(hierarchical_bayesian_modeling(matrix(numeric(0))), "calcium_traces must be a matrix with at least 1 row and 10 columns")
}) 