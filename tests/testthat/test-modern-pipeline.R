context('Modern pipeline functions')

test_that('Synthetic data generation returns expected format', {
  df <- generate_synthetic_data(n_cells = 3, n_time = 100)
  expect_true(is.data.frame(df))
  expect_true(all(grepl('Cell_', names(df)[1:3])))
  expect_true(all(grepl('BG_', names(df)[4:6])))
  expect_equal(nrow(df), 100)
})

test_that('calcium_correction returns normalized traces with Time', {
  df <- generate_synthetic_data(n_cells = 2, n_time = 50)
  corrected <- calcium_correction(df)
  expect_true('Time' %in% names(corrected))
  expect_equal(nrow(corrected), 50)
})
