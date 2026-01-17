test_that("neuropil_correct with coefficient method", {
  set.seed(123)
  n_cells <- 5
  n_time <- 100

  cell_traces <- matrix(rnorm(n_cells * n_time), nrow = n_cells)
  neuropil_traces <- matrix(rnorm(n_cells * n_time, sd = 0.5), nrow = n_cells)

  # Add neuropil contamination
  contaminated <- cell_traces + 0.7 * neuropil_traces

  result <- neuropil_correct(
    contaminated,
    neuropil_traces,
    method = "coefficient",
    coefficient = 0.7
  )

  expect_true(is.list(result))
  expect_true("corrected" %in% names(result))
  expect_equal(dim(result$corrected), dim(contaminated))
})

test_that("neuropil_correct with regression method", {
  set.seed(456)
  cell_traces <- matrix(rnorm(3 * 50), nrow = 3)
  neuropil_traces <- matrix(rnorm(3 * 50), nrow = 3)

  result <- neuropil_correct(
    cell_traces,
    neuropil_traces,
    method = "regression"
  )

  expect_true(is.list(result))
  expect_true("coefficients" %in% names(result))
  expect_equal(length(result$coefficients), 3)
})

test_that("estimate_neuropil_coefficient returns valid values", {
  set.seed(789)
  cell_trace <- rnorm(100)
  neuropil_trace <- rnorm(100)

  coef <- estimate_neuropil_coefficient(cell_trace, neuropil_trace)

  expect_true(is.list(coef))
  expect_true("regression" %in% names(coef))
  expect_true(is.numeric(coef$regression))
})

test_that("neuropil_correct handles single cell", {
  cell_trace <- matrix(rnorm(100), nrow = 1)
  neuropil_trace <- matrix(rnorm(100), nrow = 1)

  result <- neuropil_correct(cell_trace, neuropil_trace, method = "coefficient")

  expect_equal(nrow(result$corrected), 1)
})
