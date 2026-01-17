test_that("read_calcium_imaging handles CSV files", {
  # Create temporary CSV file
 tmp_file <- tempfile(fileext = ".csv")
  test_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  write.csv(test_data, tmp_file, row.names = FALSE)

  result <- read_calcium_imaging(tmp_file, format = "csv", verbose = FALSE)

  expect_true(is.matrix(result) || is.data.frame(result))
  expect_equal(nrow(result), 10)

  unlink(tmp_file)
})

test_that("read_calcium_imaging auto-detects format", {
  tmp_file <- tempfile(fileext = ".csv")
  test_data <- matrix(rnorm(50), nrow = 5, ncol = 10)
  write.csv(test_data, tmp_file, row.names = FALSE)

  result <- read_calcium_imaging(tmp_file, format = "auto", verbose = FALSE)

  expect_true(!is.null(result))

  unlink(tmp_file)
})

test_that("extract_traces_from_rois works with simple ROIs", {
  # Create synthetic movie
  movie <- array(rnorm(50 * 50 * 20), dim = c(50, 50, 20))

  # Create simple circular ROIs
  rois <- list()
  for (i in 1:3) {
    mask <- matrix(0, 50, 50)
    center_x <- 10 + (i - 1) * 15
    center_y <- 25
    for (x in 1:50) {
      for (y in 1:50) {
        if ((x - center_x)^2 + (y - center_y)^2 < 25) {
          mask[x, y] <- 1
        }
      }
    }
    rois[[i]] <- list(mask = mask)
  }

  traces <- extract_traces_from_rois(movie, rois, method = "mean", verbose = FALSE)

  expect_true(is.matrix(traces))
  expect_equal(nrow(traces), 3)
  expect_equal(ncol(traces), 20)
})

test_that("extract_traces_from_rois handles coordinate-based ROIs", {
  movie <- array(rnorm(30 * 30 * 10), dim = c(30, 30, 10))

  # Coordinate-based ROIs
  rois <- list(
    matrix(c(5, 5, 5, 6, 6, 5, 6, 6), ncol = 2, byrow = TRUE),
    matrix(c(15, 15, 15, 16, 16, 15), ncol = 2, byrow = TRUE)
  )

  traces <- extract_traces_from_rois(movie, rois, method = "mean", verbose = FALSE)

  expect_true(is.matrix(traces))
  expect_equal(nrow(traces), 2)
})
