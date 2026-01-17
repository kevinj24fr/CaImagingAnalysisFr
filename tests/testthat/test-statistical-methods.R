test_that("pca_representation returns correct structure", {
  set.seed(123)
  traces <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)

  result <- pca_representation(traces, embedding_dim = 5)

  expect_true(is.list(result))
  expect_true("embeddings" %in% names(result))
  expect_true("explained_variance" %in% names(result))
  expect_true(ncol(result$embeddings) <= 5)
})

test_that("threshold_spike_detection finds spikes", {
  set.seed(456)
  # Create trace with clear spikes
  traces <- matrix(rnorm(3 * 100), nrow = 3)
  # Add spikes
  traces[1, c(20, 50, 80)] <- 5
  traces[2, c(30, 70)] <- 4

  result <- threshold_spike_detection(traces, threshold_sd = 2)

  expect_true(is.list(result))
  expect_true("spike_predictions" %in% names(result))
  expect_equal(dim(result$spike_predictions), dim(traces))

  # Should detect the added spikes
  expect_true(sum(result$spike_predictions[1, ]) >= 3)
})

test_that("adaptive_threshold_detection works with varying baseline", {
  # Create trace with drifting baseline
  trace <- sin(seq(0, 4 * pi, length.out = 200)) * 2 + rnorm(200, sd = 0.5)
  traces <- matrix(trace, nrow = 1)

  result <- adaptive_threshold_detection(traces, window_size = 30)

  expect_true(is.list(result))
  expect_equal(nrow(result$spike_predictions), 1)
})

test_that("pca_kmeans_clustering produces valid clusters", {
  set.seed(789)
  # Create traces with distinct patterns
  traces <- rbind(
    matrix(rnorm(5 * 50, mean = 0), nrow = 5),
    matrix(rnorm(5 * 50, mean = 3), nrow = 5)
  )

  result <- pca_kmeans_clustering(traces, n_clusters = 2)

  expect_true(is.list(result))
  expect_true("cluster_assignments" %in% names(result))
  expect_equal(length(result$cluster_assignments), 10)
  expect_true(all(result$cluster_assignments %in% 1:2))
})

test_that("wavelet_denoise reduces noise", {
  set.seed(101)
  clean_signal <- sin(seq(0, 4 * pi, length.out = 100))
  noisy_signal <- clean_signal + rnorm(100, sd = 0.5)
  traces <- matrix(noisy_signal, nrow = 1)

  result <- wavelet_denoise(traces, method = "smooth", smooth_window = 5)

  expect_true(is.list(result))
  expect_true("denoised_traces" %in% names(result))

  # Denoised should be smoother (lower variance of differences)
  original_diff_var <- var(diff(traces[1, ]))
  denoised_diff_var <- var(diff(result$denoised_traces[1, ]))
  expect_true(denoised_diff_var < original_diff_var)
})

test_that("benchmark_spike_methods compares methods", {
  traces <- matrix(rnorm(3 * 100), nrow = 3)

  result <- benchmark_spike_methods(
    traces,
    methods = c("threshold", "adaptive")
  )

  expect_true(is.data.frame(result))
  expect_true("method" %in% names(result))
  expect_true(nrow(result) >= 2)
})

test_that("statistical_spike_inference works with all methods", {
  trace <- rnorm(100)

  for (method in c("rolling", "convolution", "correlation")) {
    result <- statistical_spike_inference(trace, method = method)

    expect_true(is.list(result), info = paste("Failed for method:", method))
    expect_true("spikes" %in% names(result))
    expect_true("spike_times" %in% names(result))
  }
})

test_that("deprecated functions emit warnings", {
  traces <- matrix(rnorm(5 * 50), nrow = 5)

  expect_warning(
    self_supervised_learning(traces),
    "deprecated"
  )

  expect_warning(
    contrastive_clustering(traces),
    "deprecated"
  )
})
