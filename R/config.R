#' Package Configuration
#' 
#' Centralized configuration for the CaImagingAnalysisFr package.
#' 
#' @return List of configuration parameters
#' @export
get_config <- function() {
  list(
    # Default parameters
    default_span = 0.45,
    default_spike_prob = 0.02,
    default_n_cells = 5,
    default_n_time = 1000,
    default_threshold_multiplier = 2,
    
    # Python packages with version requirements
    python_packages = list(
      oasis = ">=1.0.0",
      caiman = ">=1.9.0", 
      suite2p = ">=0.12.0",
      numpy = ">=1.20.0",
      scipy = ">=1.7.0"
    ),
    supported_methods = c("oasis", "caiman", "suite2p"),
    
    # Column name patterns
    cell_pattern = "^Cell_",
    background_pattern = "^BG_",
    
    # Validation parameters
    min_span = 0.01,
    max_span = 1.0,
    min_spike_prob = 0.001,
    max_spike_prob = 0.1,
    
    # Performance settings
    chunk_size = 1000,
    progress_bar = TRUE,
    memory_limit_mb = 1024,  # 1GB memory limit
    
    # Parameter optimization settings
    optimization = list(
      cv_folds = 5,
      span_range = seq(0.1, 0.9, by = 0.1),
      threshold_range = seq(1.5, 3.0, by = 0.25),
      max_iterations = 100,
      tolerance = 1e-6
    ),
    
    # Statistical validation parameters
    statistics = list(
      confidence_level = 0.95,
      bootstrap_samples = 1000,
      quality_threshold = 0.8,
      outlier_method = "iqr",  # "iqr", "zscore", "mad"
      outlier_threshold = 1.5
    ),
    
    # Spike detection parameters
    spike_detection = list(
      min_spike_interval = 5,  # minimum frames between spikes
      decay_time_constant = 10,  # exponential decay parameter
      amplitude_threshold = 0.1,  # minimum spike amplitude
      smoothing_window = 5  # window for spike smoothing
    ),
    
    # Quality control parameters
    quality_control = list(
      min_signal_to_noise = 2.0,
      max_missing_fraction = 0.1,
      min_trace_length = 100,
      max_contiguous_missing = 20
    )
  )
} 