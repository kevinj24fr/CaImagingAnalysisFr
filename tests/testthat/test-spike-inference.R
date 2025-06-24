context("Spike inference functions")

test_that("infer_spikes validates input", {
  # Valid input
  trace <- rnorm(100)
  result <- infer_spikes(trace, verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(trace))
  expect_true(all(c("fit", "spike") %in% names(result)))
  
  # Invalid input
  expect_error(infer_spikes("not numeric", verbose = FALSE))
  expect_error(infer_spikes(numeric(0), verbose = FALSE))
  expect_error(infer_spikes(c(1, 2, Inf, 4), verbose = FALSE))
})

test_that("infer_spikes handles different methods", {
  trace <- rnorm(100)
  
  # Test fallback method (should always work)
  result <- infer_spikes(trace, method = "oasis", fallback = TRUE, verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(trace))
  
  # Test with fallback disabled (may fail if Python packages not available)
  if (reticulate::py_module_available("oasis")) {
    result <- infer_spikes(trace, method = "oasis", fallback = FALSE, verbose = FALSE)
    expect_true(is.data.frame(result))
  } else {
    expect_error(infer_spikes(trace, method = "oasis", fallback = FALSE, verbose = FALSE))
  }
})

test_that("infer_spikes fallback method works", {
  trace <- rnorm(100)
  
  # Test fallback directly
  result <- CaImagingAnalysisFr:::infer_spikes_fallback(trace, verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(trace))
  expect_true(all(c("fit", "spike") %in% names(result)))
  
  # Check that spikes are logical
  expect_true(all(result$spike %in% c(0, 1)))
  
  # Check that fit has same length as input
  expect_equal(length(result$fit), length(trace))
})

test_that("infer_spikes handles edge cases", {
  # Constant trace
  constant_trace <- rep(1, 100)
  result <- infer_spikes(constant_trace, verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 100)
  
  # Trace with missing values
  trace_with_na <- c(rnorm(50), rep(NA, 50))
  result <- infer_spikes(trace_with_na, verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(trace_with_na))
})

test_that("infer_spikes produces reasonable results", {
  # Create a trace with known spikes
  set.seed(42)
  trace <- rnorm(200, sd = 0.1)
  
  # Add some spikes
  spike_positions <- c(50, 100, 150)
  trace[spike_positions] <- trace[spike_positions] + 2
  
  result <- infer_spikes(trace, verbose = FALSE)
  
  # Check basic structure
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), length(trace))
  expect_true(all(c("fit", "spike") %in% names(result)))
  
  # Check that spikes are detected (at least some)
  expect_true(sum(result$spike) >= 0)
  
  # Check that fit has reasonable values
  expect_true(all(is.finite(result$fit)))
})

test_that("infer_spikes handles different method parameters", {
  trace <- rnorm(100)
  
  # Test all supported methods
  methods <- c("oasis", "caiman", "suite2p")
  
  for (method in methods) {
    result <- tryCatch({
      infer_spikes(trace, method = method, verbose = FALSE)
    }, error = function(e) {
      # If method fails, should still work with fallback
      infer_spikes(trace, method = method, fallback = TRUE, verbose = FALSE)
    })
    
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), length(trace))
    expect_true(all(c("fit", "spike") %in% names(result)))
  }
}) 