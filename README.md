# CaImagingAnalysisFr

A comprehensive R package for calcium imaging analysis, from raw image processing to advanced pharmacological characterization.

[![R-CMD-check](https://github.com/kevinj24fr/CaImagingAnalysisFr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kevinj24fr/CaImagingAnalysisFr/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/kevinj24fr/CaImagingAnalysisFr/branch/main/graph/badge.svg)](https://codecov.io/gh/kevinj24fr/CaImagingAnalysisFr)

## Overview

CaImagingAnalysisFr provides a complete toolkit for analyzing calcium imaging data in R. The package covers the full analysis pipeline from raw image processing to statistical inference and is designed for both general neuroscience applications and specialized pharmacological studies (e.g., GBM cell cultures with drug treatments).

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("kevinj24fr/CaImagingAnalysisFr")
```

## Core Features

### Image Processing and Preprocessing
- TIFF/HDF5/NWB file I/O
- Motion correction with subpixel precision
- Neuropil contamination correction
- Background subtraction and dF/F computation

### Cell Detection and Tracking
- Multiple segmentation methods (threshold, watershed, k-means)
- ROI quality control and curation
- Cross-session cell registration using Hungarian algorithm matching

### Signal Analysis
- Spike inference (OASIS, Bayesian, threshold-based, ML ensemble)
- Transient detection with kinetic characterization
- Stimulus-aligned analysis (PETH, response metrics, statistical testing)

### Network Analysis
- Functional connectivity (correlation, partial correlation, Granger causality)
- Graph metrics and community detection
- Dynamic network analysis and stability metrics

## Advanced Analysis Modules

### Pharmacological Analysis
Tools for in vitro/ex vivo drug studies:
- Responder classification (activated/inhibited/non-responder)
- Dose-response curve fitting (Hill equation)
- Temporal drug response characterization
- Recovery analysis after washout

### Calcium Wave Analysis
- Wave detection with propagation speed and direction
- Initiation site identification
- Wave participation metrics per cell

### Information Theory
- Mutual information between cell pairs
- Transfer entropy for directed connectivity
- Active information storage and entropy rate

### Oscillation and Spectral Analysis
- Power spectral density (Welch, multitaper)
- Coherence and phase-locking analysis
- Frequency band power quantification

### State Space Analysis
- Population trajectory visualization (PCA, UMAP, t-SNE)
- Fixed point detection
- Trajectory comparison (Procrustes, DTW, Frechet distance)

### Neural Assemblies
- Assembly detection (PCA, ICA, NMF, correlation-based)
- Assembly activation tracking
- Sequence detection and replay analysis

### Population Dynamics (GPFA & dPCA)
Modern latent variable methods for trial-structured data:
- **GPFA**: Extract smooth low-dimensional trajectories with GP priors
- **dPCA**: Demix variance by task parameters (stimulus, time, decision)
- Cross-validation for dimensionality selection
- Decoding from demixed components

### Tuning Curve Analysis
Characterize neural selectivity:
- Orientation/direction tuning (von Mises fits, OSI/DSI)
- Contrast response functions (Naka-Rushton)
- Place field fitting (2D Gaussian)
- Spatial information in bits/spike

## R-Native Performance Features

### S7 Type System
Type-safe classes for calcium imaging data:
```r
traces <- CalciumTraces(data = trace_matrix, frame_rate = 10)
movie <- CalciumMovie(data = movie_array, frame_rate = 30)
```

### data.table Integration
High-performance operations for large datasets:
```r
dt <- traces_to_dt(traces)
dt_compute_dff(dt, baseline_frames = 1:100)
dt_detect_events(dt, threshold = 2.5)
```

### Parallel Processing
Multi-core analysis via future/furrr:
```r
setup_parallel(workers = 8)
results <- parallel_spike_detection(traces, method = "oasis")
```

### Out-of-Core Processing
Arrow/DuckDB backend for datasets exceeding RAM:
```r
save_traces_parquet(traces, "experiment.parquet")
results <- query_traces("SELECT * FROM traces WHERE cell_id < 100")
```

### Rcpp Acceleration
C++ implementations for computationally intensive operations with automatic R fallbacks.

## Quick Start

```r
library(CaImagingAnalysisFr)

# Load or generate data
data <- generate_synthetic_data(n_cells = 50, n_time = 3000)

# Basic analysis pipeline
corrected <- calcium_correction(data, method = "modern")
spikes <- infer_spikes(corrected, method = "oasis")
connectivity <- functional_connectivity(corrected)

# Pharmacological analysis (for drug treatment experiments)
responders <- classify_responders(baseline_traces, treatment_traces)
report <- gbm_analysis_report(
  baseline_traces = baseline,
  treatment_traces = treatment,
  positions = cell_positions,
  frame_rate = 10
)
```

## Documentation

- Function reference: `help(package = "CaImagingAnalysisFr")`
- Vignettes: `browseVignettes("CaImagingAnalysisFr")`
- Package website: https://kevinj24fr.github.io/CaImagingAnalysisFr/

## Requirements

- R >= 4.1.0
- Core dependencies: ggplot2, stats, MASS, igraph
- Optional packages for extended functionality:
  - S7, data.table, Rcpp (performance features)
  - future, furrr (parallelization)
  - arrow, duckdb (out-of-core processing)
  - ranger, xgboost, e1071, glmnet (machine learning)

## Citation

```bibtex
@software{CaImagingAnalysisFr,
  title = {CaImagingAnalysisFr: Comprehensive Calcium Imaging Analysis in R},
  author = {Kevin J and contributors},
  year = {2026},
  url = {https://github.com/kevinj24fr/CaImagingAnalysisFr}
}
```

## Contributing

Contributions are welcome. Please open an issue to discuss proposed changes or submit a pull request.

## License

MIT License. See [LICENSE](LICENSE) for details.
