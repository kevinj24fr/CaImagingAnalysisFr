
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

test_that('calcium_correction handles multiple cells and noisy background', {
  set.seed(1)
  dummy <- data.frame(
    Cell_1 = rnorm(10),
    Cell_2 = rnorm(10),
    Cell_3 = rnorm(10),
    BG_1   = rnorm(10, 1, 0.1),
    BG_2   = rnorm(10, 1, 0.1),
    BG_3   = rnorm(10, 1, 0.1)
  )

  res <- suppressWarnings(calcium_correction(dummy))

  expect_equal(nrow(res), nrow(dummy))
  expect_equal(ncol(res), 4)
  expect_true(all(c('Cell_1','Cell_2','Cell_3','Time') %in% names(res)))
})
