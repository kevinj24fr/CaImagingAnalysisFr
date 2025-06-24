#' Generate Synthetic Calcium Imaging Dataset
#' 
#' Creates realistic synthetic calcium imaging data with configurable
#' parameters for testing and demonstration purposes.
#' 
#' @param n_cells Number of cells to simulate (default: 5)
#' @param n_time Number of time points (default: 1000)
#' @param spike_prob Probability of a spike at any frame (default: 0.02)
#' @param seed Random seed for reproducibility (default: 1)
#' @param verbose Whether to show progress messages (default: TRUE)
#' 
#' @return Data frame with cell and background traces
#' 
#' @examples
#' # Basic usage
#' df <- generate_synthetic_data(3, 500)
#' head(df)
#' 
#' # Custom parameters
#' df <- generate_synthetic_data(n_cells = 10, n_time = 2000, spike_prob = 0.01)
#' 
#' @export
generate_synthetic_data <- function(n_cells = NULL, 
                                   n_time = NULL, 
                                   spike_prob = NULL,
                                   seed = 1,
                                   verbose = TRUE) {
  
  # Get defaults from config if not provided
  config <- get_config()
  if (is.null(n_cells)) n_cells <- config$default_n_cells
  if (is.null(n_time)) n_time <- config$default_n_time
  if (is.null(spike_prob)) spike_prob <- config$default_spike_prob
  
  # Validate inputs
  validate_numeric_param(n_cells, 1, 1000, "n_cells")
  validate_numeric_param(n_time, 10, 100000, "n_time")
  validate_numeric_param(spike_prob, config$min_spike_prob, config$max_spike_prob, "spike_prob")
  
  # Set seed for reproducibility
  set.seed(seed)
  
  if (verbose) {
    message("Generating synthetic data: ", n_cells, " cells, ", n_time, " time points")
  }
  
  # Generate background traces with slow drift
  if (verbose) message("  Generating background traces...")
  bg <- generate_background_traces(n_time, verbose)
  
  # Generate cell traces with spikes
  if (verbose) message("  Generating cell traces...")
  cells <- generate_cell_traces(n_cells, n_time, spike_prob, verbose)
  
  # Combine and return
  result <- cbind(cells, bg)
  
  if (verbose) {
    message("Synthetic data generation completed")
  }
  
  return(result)
}

#' Generate background traces
#' 
#' @param n_time Number of time points
#' @param verbose Progress messages
#' @return Background traces data frame
#' @keywords internal
generate_background_traces <- function(n_time, verbose) {
  # Vectorized background generation
  bg_matrix <- matrix(
    cumsum(rnorm(n_time * 3, sd = 0.001)) + rnorm(n_time * 3, sd = 0.02),
    ncol = 3
  )
  
  bg <- as.data.frame(bg_matrix)
  names(bg) <- paste0("BG_", seq_len(ncol(bg)))
  
  return(bg)
}

#' Generate cell traces with spikes
#' 
#' @param n_cells Number of cells
#' @param n_time Number of time points
#' @param spike_prob Spike probability
#' @param verbose Progress messages
#' @return Cell traces data frame
#' @keywords internal
generate_cell_traces <- function(n_cells, n_time, spike_prob, verbose) {
  # Pre-allocate matrix for efficiency
  cells_matrix <- matrix(0, nrow = n_time, ncol = n_cells)
  
  # Generate base traces with drift
  cells_matrix <- matrix(
    cumsum(rnorm(n_time * n_cells, sd = 0.001)) + rnorm(n_time * n_cells, sd = 0.05),
    ncol = n_cells
  )
  
  # Generate spike times for all cells at once
  spike_times <- lapply(seq_len(n_cells), function(i) {
    which(rbinom(n_time, 1, spike_prob) == 1)
  })
  
  # Apply spikes efficiently
  if (verbose) message("    Applying spikes...")
  
  for (cell_idx in seq_len(n_cells)) {
    spike_locs <- spike_times[[cell_idx]]
    
    if (length(spike_locs) > 0) {
      # Generate spike parameters
      amps <- runif(length(spike_locs), 0.5, 2)
      
      # Apply each spike
      for (spike_idx in seq_along(spike_locs)) {
        t <- spike_locs[spike_idx]
        amp <- amps[spike_idx]
        
        # Generate decay kernel
        decay_length <- 50
        decay <- amp * exp(-seq_len(decay_length) / 10)
        
        # Apply decay to trace
        end_idx <- min(n_time, t + decay_length - 1)
        actual_length <- end_idx - t + 1
        cells_matrix[t:end_idx, cell_idx] <- 
          cells_matrix[t:end_idx, cell_idx] + decay[seq_len(actual_length)]
      }
    }
  }
  
  # Convert to data frame
  cells <- as.data.frame(cells_matrix)
  names(cells) <- paste0("Cell_", seq_len(ncol(cells)))
  
  return(cells)
}
