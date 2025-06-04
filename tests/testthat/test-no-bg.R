source(file.path('..','..','R','pipeline.R'))

test_that('calcium_correction handles missing background', {
  df <- data.frame(Cell_1 = rnorm(10))
  corrected <- calcium_correction(df)
  expect_equal(nrow(corrected), 10)
  expect_true('Time' %in% names(corrected))
})
