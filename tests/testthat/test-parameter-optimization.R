test_that("Parameter optimization functions work", {
  # Generate test data
  test_data <- generate_synthetic_data(n_cells = 2, n_time = 100, verbose = FALSE)
  
  # Test correction parameter optimization
  opt_results <- optimize_correction_parameters(
    test_data, 
    target_metric = "snr", 
    span_range = c(0.3, 0.5), 
    cv_folds = 2, 
    verbose = FALSE
  )
  
  expect_true(is.list(opt_results))
  expect_true("optimal_span" %in% names(opt_results))
  expect_true("optimal_score" %in% names(opt_results))
  expect_true(is.numeric(opt_results$optimal_span))
  expect_true(is.numeric(opt_results$optimal_score))
  
  # Test spike parameter optimization
  corrected_data <- calcium_correction(test_data, verbose = FALSE)
  
  spike_opt_results <- optimize_spike_parameters(
    corrected_data,
    method = "oasis",
    threshold_range = c(2.0, 2.5),
    cv_folds = 2,
    verbose = FALSE
  )
  
  expect_true(is.list(spike_opt_results))
  expect_true("optimal_threshold" %in% names(spike_opt_results))
  expect_true("optimal_score" %in% names(spike_opt_results))
})

test_that("Auto-optimization works", {
  test_data <- generate_synthetic_data(n_cells = 1, n_time = 50, verbose = FALSE)
  
  auto_results <- auto_optimize_parameters(
    test_data,
    optimize_correction = TRUE,
    optimize_spikes = FALSE,  # Skip spikes to avoid Python dependencies
    verbose = FALSE
  )
  
  expect_true(is.list(auto_results))
  expect_true("correction" %in% names(auto_results))
})

test_that("Optimization handles edge cases", {
  # Test with very small dataset
  small_data <- generate_synthetic_data(n_cells = 1, n_time = 10, verbose = FALSE)
  
  expect_warning(
    optimize_correction_parameters(small_data, cv_folds = 5, verbose = FALSE)
  )
}) 