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
