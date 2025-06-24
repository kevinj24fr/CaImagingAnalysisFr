context("Synthetic data generation")

test_that("generate_synthetic_data produces correct structure", {
  # Test with default parameters
  df <- generate_synthetic_data(verbose = FALSE)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 1000)  # default n_time
  expect_equal(ncol(df), 8)     # 5 cells + 3 background
  
  # Check column names
  expect_true(all(grepl("^Cell_", names(df)[1:5])))
  expect_true(all(grepl("^BG_", names(df)[6:8])))
})

test_that("generate_synthetic_data works with custom parameters", {
  # Test with custom parameters
  df <- generate_synthetic_data(n_cells = 3, n_time = 100, spike_prob = 0.01, verbose = FALSE)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 100)
  expect_equal(ncol(df), 6)  # 3 cells + 3 background
  
  # Check column names
  expect_true(all(grepl("^Cell_", names(df)[1:3])))
  expect_true(all(grepl("^BG_", names(df)[4:6])))
})

test_that("generate_synthetic_data validates parameters", {
  # Invalid n_cells
  expect_error(generate_synthetic_data(n_cells = 0, verbose = FALSE))
  expect_error(generate_synthetic_data(n_cells = 1001, verbose = FALSE))
  
  # Invalid n_time
  expect_error(generate_synthetic_data(n_time = 5, verbose = FALSE))
  expect_error(generate_synthetic_data(n_time = 100001, verbose = FALSE))
  
  # Invalid spike_prob
  expect_error(generate_synthetic_data(spike_prob = 0, verbose = FALSE))
  expect_error(generate_synthetic_data(spike_prob = 0.5, verbose = FALSE))
})

test_that("generate_synthetic_data is reproducible", {
  # Test reproducibility with same seed
  df1 <- generate_synthetic_data(seed = 42, verbose = FALSE)
  df2 <- generate_synthetic_data(seed = 42, verbose = FALSE)
  expect_identical(df1, df2)
  
  # Test different seeds produce different results
  df3 <- generate_synthetic_data(seed = 43, verbose = FALSE)
  expect_false(identical(df1, df3))
})

test_that("generate_synthetic_data produces realistic data", {
  df <- generate_synthetic_data(n_cells = 2, n_time = 500, verbose = FALSE)
  
  # Check that data is numeric
  expect_true(all(sapply(df, is.numeric)))
  
  # Check that data has reasonable ranges
  for (col in names(df)) {
    expect_true(all(is.finite(df[[col]])))
    expect_true(sd(df[[col]]) > 0)  # Some variation
  }
  
  # Check that background has lower variance than cells
  bg_cols <- names(df)[grepl("^BG_", names(df))]
  cell_cols <- names(df)[grepl("^Cell_", names(df))]
  
  bg_var <- mean(sapply(df[bg_cols], var))
  cell_var <- mean(sapply(df[cell_cols], var))
  
  # Background should generally have lower variance
  expect_true(bg_var < cell_var)
})

test_that("generate_synthetic_data handles edge cases", {
  # Single cell
  df <- generate_synthetic_data(n_cells = 1, n_time = 50, verbose = FALSE)
  expect_equal(ncol(df), 4)  # 1 cell + 3 background
  
  # Minimum valid parameters
  df <- generate_synthetic_data(n_cells = 1, n_time = 10, spike_prob = 0.001, verbose = FALSE)
  expect_equal(nrow(df), 10)
  expect_equal(ncol(df), 4)
}) 