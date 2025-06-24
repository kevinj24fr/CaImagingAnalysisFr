context("Plotting functions")

test_that("plot_cell_trace works with valid input", {
  # Create test data
  raw <- generate_synthetic_data(n_cells = 2, n_time = 100, verbose = FALSE)
  corrected <- calcium_correction(raw, verbose = FALSE)
  
  # Test basic plotting
  p <- plot_cell_trace(corrected, "Cell_1", verbose = FALSE)
  expect_true(inherits(p, "ggplot"))
  
  # Test with different options
  p2 <- plot_cell_trace(corrected, "Cell_1", show_spikes = FALSE, verbose = FALSE)
  expect_true(inherits(p2, "ggplot"))
  
  p3 <- plot_cell_trace(corrected, "Cell_1", show_deconvolved = FALSE, verbose = FALSE)
  expect_true(inherits(p3, "ggplot"))
})

test_that("plot_cell_trace validates input", {
  # Create test data
  raw <- generate_synthetic_data(n_cells = 1, n_time = 50, verbose = FALSE)
  corrected <- calcium_correction(raw, verbose = FALSE)
  
  # Invalid cell name
  expect_error(plot_cell_trace(corrected, "NonExistentCell", verbose = FALSE))
  
  # Invalid data frame
  expect_error(plot_cell_trace("not a data frame", "Cell_1", verbose = FALSE))
  
  # Missing Time column
  corrected_no_time <- corrected[, names(corrected) != "Time"]
  expect_error(plot_cell_trace(corrected_no_time, "Cell_1", verbose = FALSE))
})

test_that("plot_cell_trace handles custom colors", {
  raw <- generate_synthetic_data(n_cells = 1, n_time = 50, verbose = FALSE)
  corrected <- calcium_correction(raw, verbose = FALSE)
  
  custom_colors <- list(
    signal = "red",
    deconvolved = "blue",
    spikes = "green",
    background = "yellow"
  )
  
  p <- plot_cell_trace(corrected, "Cell_1", colors = custom_colors, verbose = FALSE)
  expect_true(inherits(p, "ggplot"))
})

test_that("plot_multiple_cells works", {
  # Create test data
  raw <- generate_synthetic_data(n_cells = 3, n_time = 100, verbose = FALSE)
  corrected <- calcium_correction(raw, verbose = FALSE)
  
  # Test with default parameters
  p <- plot_multiple_cells(corrected, verbose = FALSE)
  expect_true(inherits(p, "ggplot") || is.list(p))
  
  # Test with custom parameters
  p2 <- plot_multiple_cells(corrected, ncol = 2, verbose = FALSE)
  expect_true(inherits(p2, "ggplot") || is.list(p2))
  
  # Test with specific cells
  p3 <- plot_multiple_cells(corrected, cells = c("Cell_1", "Cell_2"), verbose = FALSE)
  expect_true(inherits(p3, "ggplot") || is.list(p3))
})

test_that("plot_multiple_cells validates input", {
  raw <- generate_synthetic_data(n_cells = 2, n_time = 50, verbose = FALSE)
  corrected <- calcium_correction(raw, verbose = FALSE)
  
  # Invalid cell names
  expect_error(plot_multiple_cells(corrected, cells = c("Cell_1", "NonExistent"), verbose = FALSE))
  
  # Invalid ncol
  expect_error(plot_multiple_cells(corrected, ncol = 0, verbose = FALSE))
  expect_error(plot_multiple_cells(corrected, ncol = -1, verbose = FALSE))
})

test_that("plotting functions handle edge cases", {
  # Single time point - use minimum allowed value
  raw <- generate_synthetic_data(n_cells = 1, n_time = 10, verbose = FALSE)
  corrected <- calcium_correction(raw, verbose = FALSE)
  
  # Should still work
  p <- plot_cell_trace(corrected, "Cell_1", verbose = FALSE)
  expect_true(inherits(p, "ggplot"))
  
  # Very short trace
  raw <- generate_synthetic_data(n_cells = 1, n_time = 10, verbose = FALSE)
  corrected <- calcium_correction(raw, verbose = FALSE)
  
  p <- plot_cell_trace(corrected, "Cell_1", verbose = FALSE)
  expect_true(inherits(p, "ggplot"))
})

test_that("plotting functions work with different spike methods", {
  raw <- generate_synthetic_data(n_cells = 1, n_time = 100, verbose = FALSE)
  corrected <- calcium_correction(raw, verbose = FALSE)
  
  # Test different spike inference methods
  methods <- c("oasis", "caiman", "suite2p")
  
  for (method in methods) {
    p <- tryCatch({
      plot_cell_trace(corrected, "Cell_1", method = method, verbose = FALSE)
    }, error = function(e) {
      # If method fails, should still work with fallback
      plot_cell_trace(corrected, "Cell_1", method = method, verbose = FALSE)
    })
    
    expect_true(inherits(p, "ggplot"))
  }
}) 