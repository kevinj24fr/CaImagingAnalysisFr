context("Input validation functions")

test_that("validate_data_frame works correctly", {
  # Valid data frame
  valid_df <- data.frame(
    Cell_1 = 1:10,
    Cell_2 = 2:11,
    BG_1 = rep(1, 10),
    BG_2 = rep(2, 10)
  )
  expect_true(validate_data_frame(valid_df))
  
  # Invalid: not a data frame
  expect_error(validate_data_frame(matrix(1:10, ncol = 2)))
  
  # Invalid: empty data frame
  expect_error(validate_data_frame(data.frame()))
  
  # Invalid: no columns
  expect_error(validate_data_frame(data.frame(1:10)[, FALSE]))
  
  # Invalid: no cell columns
  expect_error(validate_data_frame(data.frame(BG_1 = 1:10)))
  
  # Invalid: no background columns
  expect_error(validate_data_frame(data.frame(Cell_1 = 1:10)))
  
  # Invalid: non-numeric columns
  invalid_df <- data.frame(
    Cell_1 = letters[1:10],
    BG_1 = 1:10
  )
  expect_error(validate_data_frame(invalid_df))
  
  # Warning for missing values
  missing_df <- data.frame(
    Cell_1 = c(1:9, NA),
    BG_1 = 1:10
  )
  expect_warning(validate_data_frame(missing_df))
})

test_that("validate_numeric_param works correctly", {
  # Valid parameters
  expect_true(validate_numeric_param(0.5, 0, 1, "test_param"))
  expect_true(validate_numeric_param(10, 1, 100, "test_param"))
  
  # Invalid: not numeric
  expect_error(validate_numeric_param("0.5", 0, 1, "test_param"))
  
  # Invalid: not single value
  expect_error(validate_numeric_param(c(0.5, 0.6), 0, 1, "test_param"))
  
  # Invalid: below minimum
  expect_error(validate_numeric_param(-1, 0, 1, "test_param"))
  
  # Invalid: above maximum
  expect_error(validate_numeric_param(2, 0, 1, "test_param"))
})

test_that("validate_trace works correctly", {
  # Valid traces
  expect_true(validate_trace(1:10))
  expect_true(validate_trace(rnorm(100)))
  
  # Invalid: not numeric
  expect_error(validate_trace(letters[1:10]))
  
  # Invalid: empty
  expect_error(validate_trace(numeric(0)))
  
  # Invalid: infinite values
  expect_error(validate_trace(c(1, 2, Inf, 4)))
}) 