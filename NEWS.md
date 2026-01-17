# CaImagingAnalysisFr 0.5.0

## New Features - Classical Machine Learning

### Cell Type Classification
* `classify_cell_types()` - Classify neurons by activity patterns using ranger, xgboost, SVM, or glmnet
* `predict_cell_types()` - Predict cell types for new data
* Automatic feature extraction from traces (statistical, temporal, spectral)
* Feature importance ranking for model interpretation

### Neural Decoding
* `decode_behavior()` - Predict behavioral variables from neural activity
* Methods: ridge, lasso, elastic_net, ranger, xgboost
* `create_lagged_features()` - Include temporal context with time lags
* `cv_decode_behavior()` - Cross-validated decoder evaluation
* Returns R-squared, correlation, and decoder weights

### Event/Spike Classification
* `train_event_classifier()` - Train ML classifier on annotated events
* `detect_events_ml()` - Apply trained classifier to new data
* Window-based feature extraction around each time point
* Class balancing for imbalanced event data

### Ensemble Methods
* `ensemble_spike_detection()` - Combine multiple spike detection methods
* Voting ensemble (majority rule)
* Stacking ensemble with ranger or xgboost meta-learner

### Feature Extraction
* `extract_trace_features()` - Comprehensive feature extraction
* Statistical features: mean, SD, skewness, kurtosis, IQR, CV
* Temporal features: autocorrelation, derivative stats, peak statistics
* Spectral features: centroid, spread, entropy, frequency band power

## Dependencies
* New suggested packages: ranger, xgboost, e1071, glmnet

---

# CaImagingAnalysisFr 0.4.0

## New Features - R-Native Advantages

### S7 Class System
* Added type-safe S7 classes: `CalciumMovie`, `CalciumTraces`, `SpikeTrains`, `ROISet`, `ExperimentSession`
* Converter functions: `as_calcium_traces()`, `as_calcium_movie()`, `as_roi_set()`
* Runtime validation and method dispatch with S7

### data.table Integration
* High-performance trace operations: `traces_to_dt()`, `dt_to_traces()`
* Fast dF/F computation: `dt_compute_dff()` with rolling baseline support
* Grouped operations: `dt_zscore()`, `dt_detect_events()`, `dt_cell_summary()`
* Efficient correlation: `dt_pairwise_cor()`
* Trace manipulation: `dt_bin_traces()`, `dt_transform()`, `dt_smooth()`

### Parallel Processing (future/furrr)
* Backend configuration: `setup_parallel()`, `shutdown_parallel()`, `parallel_info()`
* Parallel analysis: `parallel_apply_traces()`, `parallel_motion_correct()`
* Parallel detection: `parallel_spike_detection()`, `parallel_extract_traces()`
* Batch processing: `parallel_batch_process()`, `parallel_cross_correlation()`

### Formula Interface
* Formula-based analysis: `analyze_traces(dff ~ time | cell_id, data)`
* Response modeling: `fit_response_model(amplitude ~ stimulus | cell_id, data)`
* Pipeline definition: `define_pipeline(result ~ motion_correct + extract_traces + compute_dff)`
* Spike detection: `detect_spikes_formula(spikes ~ trace + threshold(2.5))`

### Custom ggplot2 Visualization
* Calcium-specific geoms: `geom_calcium_heatmap()`, `geom_spike_train()`, `geom_aligned_trace()`
* Color scales: `scale_fill_calcium()`, `scale_fill_correlation()`
* Faceting: `facet_cell()` for multi-cell plots
* Convenience plots: `plot_traces()`, `plot_correlation_matrix()`, `plot_spike_raster()`, `plot_aligned_responses()`

### Arrow/DuckDB Backend
* Parquet I/O: `save_traces_parquet()`, `load_traces_parquet()`
* SQL queries: `duckdb_traces()`, `query_traces()` for large datasets
* Streaming: `stream_process_traces()` for memory-efficient processing
* Partitioned storage: `create_partitioned_dataset()` for multi-session data
* Summary statistics: `arrow_cell_summary()` on out-of-core data

### Rcpp Acceleration
* C++ implementations for hot paths with R fallbacks
* Motion correction: `shift_image_fast()`, `phase_correlation_fast()`, `motion_correct_fast()`
* Image processing: `gaussian_blur_fast()`, `median_filter_fast()`, `threshold_image_fast()`
* Segmentation: `connected_components_fast()`, `distance_transform_fast()`
* Trace extraction: `extract_trace_fast()`, `local_correlation_fast()`

## Improvements
* All new functions have pure R fallbacks when optional packages are not installed
* Package structure supports optional dependencies (S7, data.table, Rcpp, future, arrow, duckdb)
* Comprehensive NAMESPACE exports for all new functionality

## Dependencies
* New suggested packages: S7, data.table, Rcpp, RcppArmadillo, future, furrr, arrow, duckdb, DBI, dplyr, tidyr, lme4, viridisLite
* LinkingTo: Rcpp, RcppArmadillo (for compiled code)

---

# CaImagingAnalysisFr 0.3.0

## New Features

* Added TIFF stack reading and writing support (`read_tiff_stack()`, `write_tiff_stack()`)
* Added multi-format file I/O including HDF5, NWB, and MATLAB formats (`read_calcium_imaging()`)
* Added motion correction with phase correlation (`motion_correct()`, `phase_correlation()`)
* Added neuropil contamination correction (`neuropil_correct()`, `extract_neuropil_traces()`)
* Added cross-session cell registration with Hungarian algorithm matching (`register_cells()`)
* Added stimulus-aligned analysis tools (`align_to_events()`, `compute_response_metrics()`, `test_event_responses()`)
* Added PCA-based representation learning (`pca_representation()`)
* Added multiple spike detection methods (`threshold_spike_detection()`, `adaptive_threshold_detection()`)

## Improvements

* Fixed all stub functions in segmentation module with proper implementations:
  - `label_connected_components()` now uses BFS flood-fill algorithm

  - `distance_transform()` implements two-pass Euclidean distance transform
  - `find_local_maxima()` includes non-maximum suppression
  - `apply_convolution()` performs proper 2D convolution with zero-padding
  - All QC metric functions now return meaningful results

* Renamed misleading function names to accurately reflect implementations:
  - `self_supervised_learning()` -> `pca_representation()` (uses PCA, not deep learning)
  - `transformer_spike_inference()` -> `threshold_spike_detection()`
  - `contrastive_clustering()` -> `pca_kmeans_clustering()`
  - Old names retained as deprecated aliases for backwards compatibility

* Modernized pipeline infrastructure with proper variable scoping for targets
* Updated NAMESPACE with all new exports

## Breaking Changes
* Minimum R version increased to 4.1.0 for native pipe support
* Deprecated functions now emit warnings directing users to new names

---
# CaImagingAnalysisFr 0.2.0

## New Features

* Added Bayesian spike inference methods
* Added batch effect correction (ComBat, MNN)
* Added dynamic network analysis
* Added interactive visualization tools
* Added automated QC reporting

## Improvements

* Improved calcium correction algorithms
* Enhanced validation functions
* Better error handling throughout

---

# CaImagingAnalysisFr 0.1.0

* Initial release
* Basic calcium correction and spike inference
* Cell segmentation framework
* Network analysis tools
* Plotting utilities
