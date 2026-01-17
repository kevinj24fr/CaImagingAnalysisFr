# CaImagingAnalysisFr

Comprehensive R package for calcium imaging analysis.

[![R-CMD-check](https://github.com/kevinj24fr/CaImagingAnalysisFr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kevinj24fr/CaImagingAnalysisFr/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/kevinj24fr/CaImagingAnalysisFr/branch/main/graph/badge.svg)](https://codecov.io/gh/kevinj24fr/CaImagingAnalysisFr)

## Overview

CaImagingAnalysisFr provides a complete toolkit for analyzing calcium imaging data in R. The package covers the full analysis pipeline from raw image processing to statistical inference, without requiring Python dependencies for core functionality.

## Features

**Image Processing**
- TIFF stack reading and writing
- Motion correction with subpixel precision
- Neuropil contamination correction

**Cell Detection and Tracking**
- Multiple segmentation methods (threshold, k-means, watershed)
- ROI quality control metrics
- Cross-session cell registration

**Signal Analysis**
- Background correction with multiple methods
- Spike inference (OASIS, Bayesian, threshold-based)
- Stimulus-aligned analysis (PETH, response metrics)

**Network Analysis**
- Functional connectivity estimation
- Graph metrics and community detection
- Dynamic network analysis

**Batch Processing**
- Batch effect correction (ComBat, MNN)
- Reproducible workflows via targets
- Automated QC reporting

## Installation

Install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("kevinj24fr/CaImagingAnalysisFr")
```

## Quick Start

```r
library(CaImagingAnalysisFr)

# Generate example data
data <- generate_synthetic_data(n_cells = 20, n_time = 1000)

# Correct calcium traces
corrected <- calcium_correction(data, method = "modern")

# Infer spikes
spikes <- infer_spikes(corrected, method = "oasis")

# Compute functional connectivity
connectivity <- functional_connectivity(corrected, method = "correlation")

# Generate QC report
qc <- generate_qc_report(data = corrected)
```

## Motion Correction

```r
# Read imaging data
movie <- read_tiff_stack("recording.tif")

# Correct motion
corrected_movie <- motion_correct(movie, method = "template")

# Assess correction quality
quality <- assess_motion_correction(movie, corrected_movie)
```
## Stimulus-Aligned Analysis

```r
# Align traces to stimulus events
aligned <- align_to_events(
  traces = cell_traces,
  events = stimulus_times,
  pre_frames = 30,
  post_frames = 60
)

# Compute response metrics
metrics <- compute_response_metrics(aligned, response_window = c(5, 30))

# Statistical testing
stats <- test_event_responses(aligned, test = "wilcox")
```

## Cross-Session Registration

```r
# Register cells across sessions
registration <- register_cells(
  sessions = list(session1_rois, session2_rois),
  method = "hungarian",
  max_distance = 10
)

# Extract matched traces
matched <- extract_registered_traces(registration, traces_list)
```

## Documentation

- Function reference: `help(package = "CaImagingAnalysisFr")`
- Vignettes: `browseVignettes("CaImagingAnalysisFr")`
- Package website: https://kevinj24fr.github.io/CaImagingAnalysisFr/

## Requirements

- R >= 4.1.0
- See DESCRIPTION for package dependencies

## Citation

If you use this package in your research, please cite:

```
@software{CaImagingAnalysisFr,
  title = {CaImagingAnalysisFr: Comprehensive Calcium Imaging Analysis in R},
  author = {Kevin J and contributors},
  year = {2024},
  url = {https://github.com/kevinj24fr/CaImagingAnalysisFr}
}
```

## Contributing

Contributions are welcome. Please open an issue to discuss proposed changes or submit a pull request.

## License

MIT License. See [LICENSE](LICENSE) for details.
