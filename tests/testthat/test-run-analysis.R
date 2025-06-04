source(file.path('..','..','R','pipeline.R'))

test_that('run_analysis returns expected list components', {
  df <- generate_synthetic_data(n_cells = 2, n_time = 30)
  res <- run_analysis(df, penalty = 1, ar_order = 1)
  expect_true(is.list(res))
  expect_true(all(c('corrected','cell_plot','raster_plot') %in% names(res)))
})
