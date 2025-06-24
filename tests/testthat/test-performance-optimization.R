test_that("Performance optimization functions work", {
  # Generate test data
  test_data <- generate_synthetic_data(n_cells = 2, n_time = 100, verbose = FALSE)
  
  # Test chunked processing
  process_func <- function(chunk) {
    calcium_correction(chunk, verbose = FALSE)
  }
  
  combine_func <- function(results) {
    do.call(rbind, results[!sapply(results, is.null)])
  }
  
  chunked_result <- process_in_chunks(
    test_data,
    chunk_size = 50,
    process_function = process_func,
    combine_function = combine_func,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(chunked_result))
  expect_equal(nrow(chunked_result), nrow(test_data))
  
  # Test chunked correction
  chunked_corrected <- calcium_correction_chunked(
    test_data,
    chunk_size = 50,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(chunked_corrected))
  expect_equal(nrow(chunked_corrected), nrow(test_data))
  
  # Test data structure optimization
  optimized_data <- optimize_data_structure(
    test_data,
    optimize_types = TRUE,
    remove_na = FALSE,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(optimized_data))
  expect_equal(nrow(optimized_data), nrow(test_data))
  expect_equal(ncol(optimized_data), ncol(test_data))
})

test_that("Memory monitoring works", {
  # Test memory monitoring
  mem_info <- monitor_memory_usage("test_operation", verbose = FALSE)
  
  expect_true(is.list(mem_info))
  expect_true("operation" %in% names(mem_info))
  expect_true("memory_mb" %in% names(mem_info))
})

test_that("Performance benchmarking works", {
  test_data <- generate_synthetic_data(n_cells = 1, n_time = 50, verbose = FALSE)
  
  # Define simple operations for testing
  operations <- list(
    correction = function(x) calcium_correction(x, verbose = FALSE)
  )
  
  benchmark_results <- benchmark_performance(
    test_data,
    operations = operations,
    n_repeats = 2,  # Small number for testing
    verbose = FALSE
  )
  
  expect_true(is.list(benchmark_results))
  expect_true("correction" %in% names(benchmark_results))
  expect_true("mean_time" %in% names(benchmark_results$correction))
  expect_true("mean_memory" %in% names(benchmark_results$correction))
})

test_that("Chunked spike inference works", {
  test_data <- generate_synthetic_data(n_cells = 2, n_time = 100, verbose = FALSE)
  corrected_data <- calcium_correction(test_data, verbose = FALSE)
  
  chunked_spikes <- infer_spikes_chunked(
    corrected_data,
    method = "oasis",
    chunk_size = 50,
    verbose = FALSE
  )
  
  expect_true(is.list(chunked_spikes))
  expect_true("Cell_1" %in% names(chunked_spikes))
  expect_true("Cell_2" %in% names(chunked_spikes))
  expect_true(is.data.frame(chunked_spikes$Cell_1))
  expect_true("fit" %in% names(chunked_spikes$Cell_1))
  expect_true("spike" %in% names(chunked_spikes$Cell_1))
}) 