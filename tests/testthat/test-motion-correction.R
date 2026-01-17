test_that("shift_image performs basic shifts", {
  img <- matrix(0, 10, 10)
  img[5, 5] <- 1

  # Shift by 1 pixel
  shifted <- shift_image(img, dy = 1, dx = 0)

  expect_equal(dim(shifted), dim(img))
  # The bright pixel should have moved
  expect_true(shifted[6, 5] > shifted[5, 5])
})

test_that("shift_image handles zero shift", {
  img <- matrix(runif(100), 10, 10)

  shifted <- shift_image(img, dy = 0, dx = 0)

  expect_equal(shifted, img, tolerance = 1e-10)
})

test_that("phase_correlation returns valid shifts", {
  # Create template
  template <- matrix(0, 20, 20)
  template[8:12, 8:12] <- 1

  # Create shifted image
  shifted <- matrix(0, 20, 20)
  shifted[10:14, 9:13] <- 1  # Shifted by (2, 1)

  result <- phase_correlation(shifted, template, max_shift = 5)

  expect_true(is.list(result))
  expect_true("dy" %in% names(result))
  expect_true("dx" %in% names(result))
  expect_true(abs(result$dy) <= 5)
  expect_true(abs(result$dx) <= 5)
})

test_that("motion_correct works on small movie", {
  # Create small movie with simulated motion
  n_frames <- 5
  movie <- array(0, dim = c(20, 20, n_frames))

  for (i in 1:n_frames) {
    movie[8:12, 8:12, i] <- 1
  }

  # Add some artificial shift to one frame
  movie[, , 3] <- shift_image(movie[, , 1], dy = 1, dx = 1)

  result <- motion_correct(movie, method = "template", max_shift = 3, verbose = FALSE)

  expect_true(is.list(result))
  expect_true("corrected" %in% names(result))
  expect_true("shifts" %in% names(result))
  expect_equal(dim(result$corrected), dim(movie))
})

test_that("assess_motion_correction returns quality metrics", {
  original <- array(rnorm(20 * 20 * 5), dim = c(20, 20, 5))
  corrected <- original  # No actual correction

  metrics <- assess_motion_correction(original, corrected)

  expect_true(is.list(metrics))
  expect_true("correlation_improvement" %in% names(metrics))
})
