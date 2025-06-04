source(file.path('..','..','R', 'pipeline.R'))

test_that('calciumcorrection returns expected structure', {
  dummy <- data.frame(
    Cell_1 = 1:10,
    Cell_2 = seq(2,20,by=2),
    BG_1   = rep(1,10),
    BG_2   = rep(2,10),
    BG_3   = rep(3,10)
  )

  result <- suppressWarnings(calcium_correction(dummy))

  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(dummy))
  expect_equal(ncol(result), 3)
  expect_true('Time' %in% names(result))
})
