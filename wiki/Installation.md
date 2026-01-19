# Installation Guide

## Prerequisites

- **R version:** 4.1.0 or higher
- **Operating System:** Windows, macOS, or Linux

## Basic Installation

### From GitHub (Recommended)

```r
# Install remotes if you don't have it
install.packages("remotes")

# Install CaImagingAnalysisFr
remotes::install_github("kevinj24fr/CaImagingAnalysisFr")
```

### From Source

```bash
git clone https://github.com/kevinj24fr/CaImagingAnalysisFr.git
cd CaImagingAnalysisFr
R CMD INSTALL .
```

## Installing Optional Dependencies

The package has many optional dependencies that enable additional features. Install based on your needs:

### Performance Enhancement

For high-performance data processing and large datasets:

```r
install.packages(c(
  "S7",           # Type-safe object system
  "data.table",   # High-performance data frames
  "Rcpp",         # C++ acceleration
  "RcppArmadillo" # Linear algebra acceleration
))
```

### Parallel Computing

For parallel processing of large datasets:

```r
install.packages(c(
  "future",       # Parallel backend
  "furrr"         # Parallel purrr operations
))
```

### Big Data / Out-of-Core Processing

For datasets larger than RAM:

```r
install.packages(c(
  "arrow",        # Apache Arrow backend

  "duckdb",       # DuckDB SQL engine
  "DBI"           # Database interface
))
```

### File Format Support

```r
# HDF5 files
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rhdf5")

# TIFF stacks
install.packages("tiff")

# MATLAB files
install.packages("R.matlab")
```

### Machine Learning

```r
install.packages(c(
  "ranger",       # Random forests
  "xgboost",      # Gradient boosting
  "e1071",        # SVM and Naive Bayes
  "glmnet"        # Regularized regression
))
```

### Dimensionality Reduction

```r
install.packages(c(
  "Rtsne",        # t-SNE
  "uwot",         # UMAP
  "fastICA"       # Independent Component Analysis
))
```

### Signal Processing

```r
install.packages(c(
  "signal",       # Signal processing
  "waveslim",     # Wavelet analysis
  "multitaper",   # Spectral analysis
  "dtw"           # Dynamic time warping
))
```

### Interactive Visualization

```r
install.packages(c(
  "plotly",       # Interactive plots
  "shiny"         # Web applications
))
```

### Statistical Modeling

```r
install.packages(c(
  "lme4",         # Mixed-effects models
  "sva"           # Batch correction
))
```

## Python Integration (Optional)

Some advanced features use Python packages via reticulate. To enable these:

### 1. Install reticulate

```r
install.packages("reticulate")
```

### 2. Create a Python Environment

```r
library(reticulate)

# Create a dedicated conda environment
conda_create("caimaging", python = "3.9")
use_condaenv("caimaging")

# Or use virtualenv
virtualenv_create("caimaging")
use_virtualenv("caimaging")
```

### 3. Install Python Packages

```r
# Spike inference
py_install("oasis-deconv", envname = "caimaging")

# Motion correction and cell segmentation
py_install("caiman", envname = "caimaging")
py_install("suite2p", envname = "caimaging")

# Additional packages
py_install(c("numpy", "scipy", "scikit-learn"), envname = "caimaging")
```

## Verifying Installation

```r
library(CaImagingAnalysisFr)

# Check package version
packageVersion("CaImagingAnalysisFr")

# Run a simple test
data <- matrix(rnorm(1000), nrow = 10, ncol = 100)
experiment <- CaExperiment(data, frame_rate = 30)
print(experiment)
```

## Troubleshooting

### Common Issues

**Issue:** Installation fails with compilation errors
```r
# Try installing without compiled code first
remotes::install_github("kevinj24fr/CaImagingAnalysisFr",
                         INSTALL_opts = "--no-multiarch")
```

**Issue:** Missing system dependencies on Linux
```bash
# Ubuntu/Debian
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev

# CentOS/RHEL
sudo yum install libcurl-devel openssl-devel libxml2-devel
```

**Issue:** rhdf5 installation fails
```r
# Make sure BiocManager is up to date
BiocManager::install("rhdf5", force = TRUE)
```

### Getting Help

If you encounter issues:

1. Check the [FAQ](FAQ) page
2. Search existing [GitHub Issues](https://github.com/kevinj24fr/CaImagingAnalysisFr/issues)
3. Open a new issue with:
   - Your R version (`R.version`)
   - Package version (`packageVersion("CaImagingAnalysisFr")`)
   - Operating system
   - Full error message
