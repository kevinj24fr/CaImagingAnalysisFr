# CaImagingAnalysisFr

A comprehensive R package for calcium imaging analysis with advanced features for parameter optimization, statistical validation, and performance optimization.

## Features

### Core Analysis
- **Calcium Trace Correction**: Background subtraction, normalization, and baseline correction
- **Spike Inference**: Multiple algorithms (OASIS, CaImAn, Suite2p) with fallback mechanisms
- **Synthetic Data Generation**: Realistic calcium imaging data for testing and validation
- **Interactive Visualization**: Shiny app for real-time parameter adjustment and exploration

### Advanced Features

#### Parameter Optimization
- **Cross-validation based optimization** for correction parameters
- **Automatic parameter selection** based on multiple metrics (SNR, baseline stability, spike detection)
- **Multi-method comparison** to find optimal spike detection algorithms

```r
# Auto-optimize all parameters
optimization_results <- auto_optimize_parameters(raw_data)

# Optimize specific parameters
correction_opt <- optimize_correction_parameters(raw_data, target_metric = "snr")
spike_opt <- optimize_spike_parameters(corrected_data, method = "oasis")
```

#### Statistical Validation
- **Confidence intervals** for spike detection using bootstrap resampling
- **Quality metrics** for correction and spike detection
- **Comprehensive statistical analysis** with recommendations
- **Outlier detection** using multiple methods (IQR, Z-score, MAD)

```r
# Calculate confidence intervals
ci_results <- spike_confidence_intervals(trace, method = "oasis")

# Quality assessment
quality <- calculate_correction_quality(raw_data, corrected_data)
validation <- validate_spike_detection(trace, spikes, method = "oasis")

# Comprehensive analysis
analysis <- comprehensive_statistical_analysis(raw_data, corrected_data, spike_results)
```

#### Performance Optimization
- **Chunked processing** for large datasets
- **Memory management** and monitoring
- **Stream processing** for very large files
- **Performance benchmarking** tools

```r
# Process large datasets efficiently
chunked_result <- calcium_correction_chunked(large_data, chunk_size = 1000)
chunked_spikes <- infer_spikes_chunked(corrected_data, method = "oasis")

# Monitor memory usage
mem_info <- monitor_memory_usage("my_analysis")

# Benchmark performance
benchmark_results <- benchmark_performance(test_data)
```

#### Python Dependency Management
- **Automatic installation** of required Python packages
- **Version checking** and compatibility validation
- **Environment management** for production deployment

```r
# Install all dependencies
install_python_dependencies()

# Check specific packages
deps <- manage_python_dependencies(packages = c("oasis", "caiman"))

# Get environment info
python_info <- get_python_info()
```

## Installation

### Prerequisites
- R (>= 4.0.0)
- Python (>= 3.8) with pip
- Required R packages: `reticulate`, `ggplot2`, `shiny`, `plotly`, `targets`

### Install R Package
```r
# Install from GitHub
devtools::install_github("kevinj24fr/CaImagingAnalysisFr")

# Or install locally
devtools::install_local("path/to/CaImagingAnalysisFr")
```

### Install Python Dependencies
```r
library(CaImagingAnalysisFr)

# Install all required Python packages
install_python_dependencies()
```

## Quick Start

### Basic Analysis
```r
library(CaImagingAnalysisFr)

# Generate synthetic data
raw_data <- generate_synthetic_data(n_cells = 5, n_time = 1000)

# Apply correction
corrected_data <- calcium_correction(raw_data)

# Detect spikes
spikes <- infer_spikes(corrected_data$Cell_1)

# Visualize results
plot_cell_trace(corrected_data, "Cell_1")
```

### Advanced Analysis with Optimization
```r
# Auto-optimize parameters
optimization <- auto_optimize_parameters(raw_data)

# Apply optimized correction
corrected_opt <- calcium_correction(raw_data, span = optimization$correction$optimal_span)

# Statistical validation
analysis <- comprehensive_statistical_analysis(raw_data, corrected_opt, spike_results)

# Check quality
if (analysis$passes_quality_threshold) {
  message("Analysis quality is good!")
} else {
  message("Consider adjusting parameters: ", analysis$recommendations)
}
```

### Large Dataset Processing
```r
# Process large datasets efficiently
chunked_corrected <- calcium_correction_chunked(large_data, chunk_size = 5000)

# Monitor memory usage
monitor_memory_usage("large_analysis")

# Stream process very large files
stream_results <- stream_process_data(
  "large_file.csv",
  chunk_size = 1000,
  process_function = function(chunk) calcium_correction(chunk, verbose = FALSE)
)
```

### Interactive Analysis
```r
# Launch interactive viewer
launch_interactive_viewer(raw_data)
```

## Pipeline Workflow

### Using targets for Reproducible Analysis
```r
# Define pipeline
pipeline <- calcium_pipeline(
  n_cells = 10,
  n_time = 2000,
  correction_method = "modern",
  spike_method = "oasis"
)

# Run pipeline
library(targets)
tar_make()
```

## Configuration

The package uses a centralized configuration system:

```r
# Get current configuration
config <- get_config()

# Key configuration areas:
# - Default parameters for all functions
# - Python package requirements
# - Quality control thresholds
# - Performance settings
# - Statistical validation parameters
```

## Quality Control

The package includes comprehensive quality control:

- **Signal-to-noise ratio** assessment
- **Missing data** detection and handling
- **Temporal consistency** validation
- **Outlier detection** and flagging
- **Statistical validation** of results

## Performance Considerations

### Memory Management
- **Chunked processing** for datasets > 1GB
- **Memory monitoring** and warnings
- **Data structure optimization**
- **Stream processing** for very large files

### Optimization Tips
- Use `calcium_correction_chunked()` for large datasets
- Enable `optimize_data_structure()` for memory efficiency
- Monitor memory usage with `monitor_memory_usage()`
- Use `benchmark_performance()` to identify bottlenecks

## Statistical Rigor

The package provides tools for statistically rigorous analysis:

- **Bootstrap confidence intervals** for spike detection
- **Cross-validation** for parameter optimization
- **Quality metrics** with thresholds
- **Comprehensive validation** reports
- **Recommendations** for parameter adjustment

## Troubleshooting

### Python Issues
```r
# Check Python environment
python_info <- get_python_info()

# Reinstall dependencies
install_python_dependencies(verbose = TRUE)

# Check specific package
deps <- manage_python_dependencies(packages = "oasis")
```

### Memory Issues
```r
# Monitor memory usage
monitor_memory_usage("before_analysis")

# Use chunked processing
result <- calcium_correction_chunked(data, chunk_size = 500)

# Optimize data structure
optimized_data <- optimize_data_structure(data)
```

### Quality Issues
```r
# Check data quality
validate_data_frame(raw_data)

# Assess correction quality
quality <- calculate_correction_quality(raw_data, corrected_data)

# Get recommendations
analysis <- comprehensive_statistical_analysis(raw_data, corrected_data, spikes)
print(analysis$recommendations)
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Testing

Run the test suite:
```r
devtools::test()
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```
CaImagingAnalysisFr: A Comprehensive R Package for Calcium Imaging Analysis
Author, Year
```

## Support

For issues and questions:
- Check the troubleshooting section
- Review the vignettes
- Open an issue on GitHub
- Contact the maintainers

## Changelog

### Version 2.0.0 (Current)
- Added parameter optimization with cross-validation
- Implemented statistical validation and confidence intervals
- Added performance optimization and memory management
- Enhanced Python dependency management
- Added comprehensive quality control
- Improved error handling and fallback mechanisms
- Added performance benchmarking tools
- Enhanced documentation and examples

### Version 1.0.0
- Initial release with basic calcium imaging analysis
- Core correction and spike detection functions
- Basic visualization and interactive app
- Pipeline workflow with targets
