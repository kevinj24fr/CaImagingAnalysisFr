# CaImagingAnalysisFr

This package provides basic utilities for processing calcium imaging traces. Functions are organized as an R package for easier reuse and testing. Key features include:

- Background correction and normalization
- Spike inference via OASIS with optional CaImAn or Suite2p backends
- Interactive visualization using Shiny and Plotly
- Reproducible workflows using the `targets` package

## Installation

```
# From the package root
R CMD build .
R CMD INSTALL CaImagingAnalysisFr_0.1.0.tar.gz
```

To manage R dependencies, you can initialize `renv`:

```
renv::init()
```

Python modules are accessed via `reticulate`. You can create a dedicated virtual environment:

```
reticulate::virtualenv_create('caimaging-env')
reticulate::virtualenv_install('caimaging-env', c('oasis', 'caiman', 'suite2p'))
```

## Usage

Generate synthetic data and view it interactively:

```r
library(CaImagingAnalysisFr)
raw <- generate_synthetic_data(3, 500)
launch_interactive_viewer(raw)
```

## Continuous Integration

Tests in `tests/testthat` are executed automatically via GitHub Actions using the workflow found in `.github/workflows/test.yaml`.
