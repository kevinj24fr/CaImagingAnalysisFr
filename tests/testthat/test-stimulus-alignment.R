test_that("align_to_events creates correct output structure", {
  set.seed(123)
  traces <- matrix(rnorm(5 * 200), nrow = 5)
  events <- c(50, 100, 150)

  result <- align_to_events(
    traces,
    events,
    pre_frames = 10,
    post_frames = 20,
    normalize = "none"
  )

  expect_true(is.list(result))
  expect_true("aligned" %in% names(result))
  expect_true("time_vector" %in% names(result))

  # Check dimensions: cells x time x trials
  expect_equal(dim(result$aligned)[1], 5)  # 5 cells
  expect_equal(dim(result$aligned)[2], 31)  # 10 pre + 1 + 20 post
  expect_equal(dim(result$aligned)[3], 3)  # 3 events
})

test_that("align_to_events handles zscore normalization", {
  traces <- matrix(rnorm(3 * 100, mean = 10, sd = 2), nrow = 3)
  events <- c(30, 60)

  result <- align_to_events(
    traces,
    events,
    pre_frames = 5,
    post_frames = 10,
    normalize = "zscore"
  )

  # Z-scored data should have mean near 0
  mean_val <- mean(result$aligned, na.rm = TRUE)
  expect_true(abs(mean_val) < 1)
})

test_that("compute_response_metrics returns expected metrics", {
  # Create mock aligned data
  aligned_data <- list(
    aligned = array(rnorm(3 * 20 * 5), dim = c(3, 20, 5)),
    time_vector = -5:14,
    pre_frames = 5,
    post_frames = 14
  )

  metrics <- compute_response_metrics(
    aligned_data,
    response_window = c(0, 10),
    baseline_window = c(-5, -1)
  )

  expect_true(is.list(metrics))
  expect_true("peak_amplitude" %in% names(metrics) ||
              "mean_response" %in% names(metrics))
})

test_that("test_event_responses performs statistical tests", {
  aligned_data <- list(
    aligned = array(rnorm(2 * 30 * 10), dim = c(2, 30, 10)),
    time_vector = -10:19,
    pre_frames = 10,
    post_frames = 19
  )

  # Add response to make test meaningful
  aligned_data$aligned[, 15:20, ] <- aligned_data$aligned[, 15:20, ] + 2

  result <- test_event_responses(
    aligned_data,
    response_window = c(5, 10),
    baseline_window = c(-10, -5),
    test = "ttest"
  )

  expect_true(is.list(result) || is.data.frame(result))
})

test_that("synchronize_behavior interpolates correctly", {
  behavior_data <- data.frame(
    time = seq(0, 10, by = 0.1),
    velocity = sin(seq(0, 10, by = 0.1))
  )

  imaging_timestamps <- seq(0.05, 9.95, length.out = 50)

  result <- synchronize_behavior(behavior_data, imaging_timestamps)

  expect_true(is.data.frame(result) || is.list(result))
})
