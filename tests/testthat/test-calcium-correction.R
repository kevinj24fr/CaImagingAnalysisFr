context("Calcium correction functions")

test_that("calcium_correction works with valid input", {
  # Create test data
  test_data <- data.frame(
    Cell_1 = 1:10 + rnorm(10, sd = 0.1),
    Cell_2 = 2:11 + rnorm(10, sd = 0.1),
    BG_1 = rep(1, 10) + rnorm(10, sd = 0.05),
    BG_2 = rep(2, 10) + rnorm(10, sd = 0.05)
  )
  
  # Test modern method
  result_modern <- calcium_correction(test_data, method = "modern", verbose = FALSE)
  expect_true(is.data.frame(result_modern))
  expect_equal(nrow(result_modern), nrow(test_data))
  expect_true("Time" %in% names(result_modern))
  expect_true("Cell_1" %in% names(result_modern))
  expect_true("Cell_2" %in% names(result_modern))
  
  # Test legacy method
  result_legacy <- calcium_correction(test_data, method = "legacy", verbose = FALSE)
  expect_true(is.data.frame(result_legacy))
  expect_equal(nrow(result_legacy), nrow(test_data))
  expect_true("Time" %in% names(result_legacy))
})

test_that("calcium_correction handles edge cases", {
  # Test with missing values
  test_data_missing <- data.frame(
    Cell_1 = c(1:9, NA),
    BG_1 = 1:10
  )
  
  expect_warning(result <- calcium_correction(test_data_missing, verbose = FALSE))
  expect_true(is.data.frame(result))
  
  # Test with single cell
  test_data_single <- data.frame(
    Cell_1 = 1:10,
    BG_1 = 1:10
  )
  
  result <- calcium_correction(test_data_single, verbose = FALSE)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 10)
})

test_that("calcium_correction validates parameters", {
  test_data <- data.frame(
    Cell_1 = 1:10,
    BG_1 = 1:10
  )
  
  # Invalid span
  expect_error(calcium_correction(test_data, span = -1, verbose = FALSE))
  expect_error(calcium_correction(test_data, span = 2, verbose = FALSE))
  
  # Invalid method
  expect_error(calcium_correction(test_data, method = "invalid", verbose = FALSE))
})

test_that("calcium_correction produces different results for different methods", {
  test_data <- data.frame(
    Cell_1 = 1:20 + rnorm(20, sd = 0.1),
    BG_1 = rep(1, 20) + rnorm(20, sd = 0.05)
  )
  
  result_modern <- calcium_correction(test_data, method = "modern", verbose = FALSE)
  result_legacy <- calcium_correction(test_data, method = "legacy", verbose = FALSE)
  
  # Results should be different due to different normalization approaches
  expect_false(identical(result_modern$Cell_1, result_legacy$Cell_1))
})

test_that("calcium_correction handles normalization parameter", {
  test_data <- data.frame(
    Cell_1 = 1:10,
    BG_1 = 1:10
  )
  
  result_normalized <- calcium_correction(test_data, normalize = TRUE, verbose = FALSE)
  result_not_normalized <- calcium_correction(test_data, normalize = FALSE, verbose = FALSE)
  
  # Results should be different when normalization is toggled
  expect_false(identical(result_normalized$Cell_1, result_not_normalized$Cell_1))
}) 