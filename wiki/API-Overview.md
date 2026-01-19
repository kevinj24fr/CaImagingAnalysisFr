# API Overview

A comprehensive overview of the main modules and functions in CaImagingAnalysisFr.

## Data Input/Output

### Reading Data

| Function | Description |
|----------|-------------|
| `read_tiff_stack()` | Load multi-page TIFF files |
| `read_hdf5_traces()` | Load traces from HDF5 files |
| `read_nwb()` | Load NWB (Neurodata Without Borders) files |
| `read_matlab_traces()` | Load MATLAB .mat files |
| `read_suite2p()` | Load Suite2P output |
| `read_caiman()` | Load CaImAn output |

### Writing Data

| Function | Description |
|----------|-------------|
| `write_hdf5()` | Export to HDF5 format |
| `write_nwb()` | Export to NWB format |
| `export_traces()` | Export traces to CSV |
| `ExportCommands()` | Export analysis pipeline as R script |

## Preprocessing

### Motion Correction

```r
experiment <- RunMotionCorrection(experiment,
  method = "template",    # "template", "piecewise", "optical_flow"
  max_shift = 20,
  template = "mean"
)
```

### Signal Correction

```r
# Neuropil correction
experiment <- RunCorrection(experiment,
  method = "neuropil",
  neuropil_factor = 0.7
)

# Background subtraction
experiment <- RunCorrection(experiment,
  method = "background",
  percentile = 8
)
```

### dF/F Computation

```r
experiment <- RunDFF(experiment,
  method = "rolling",     # "rolling", "baseline", "percentile"
  window = 300,           # frames
  percentile = 8
)
```

## Spike Inference

### Threshold-Based

```r
experiment <- RunSpikes(experiment,
  method = "threshold",
  threshold = 2.5,        # SD units
  min_distance = 3        # frames
)
```

### Deconvolution

```r
experiment <- RunSpikes(experiment,
  method = "oasis",       # Requires Python
  tau = 1.0,              # decay time constant
  smin = 0.5              # minimum spike size
)
```

### Statistical Methods

```r
experiment <- RunSpikes(experiment,
  method = "mlspike",
  prior = "sparse"
)
```

## Dimensionality Reduction

### Standard Methods

```r
# PCA
experiment <- RunPCA(experiment, n_components = 50)

# UMAP
experiment <- RunUMAP(experiment, n_neighbors = 15, min_dist = 0.1)

# t-SNE
experiment <- RunTSNE(experiment, perplexity = 30)

# ICA
experiment <- RunICA(experiment, n_components = 20)
```

### Advanced Methods

```r
# Gaussian Process Factor Analysis
experiment <- RunGPFA(experiment, n_latents = 10)

# Demixed PCA
experiment <- RunDPCA(experiment, conditions = condition_labels)

# Tensor Decomposition
experiment <- RunTensorDecomp(experiment,
  method = "cp",          # "cp" or "tucker"
  rank = 10
)
```

## Network Analysis

### Connectivity

```r
# Correlation-based
experiment <- RunConnectivity(experiment,
  method = "correlation",
  threshold = 0.3
)

# Partial correlation
experiment <- RunConnectivity(experiment,
  method = "partial",
  regularization = 0.01
)

# Transfer entropy
experiment <- RunConnectivity(experiment,
  method = "transfer_entropy",
  lag = 1
)
```

### Graph Analysis

```r
# Compute graph metrics
metrics <- graph_metrics(GetGraph(experiment, "connectivity"))

# Community detection
experiment <- RunCommunities(experiment, method = "louvain")

# Graph visualization
plot_network(experiment, layout = "fr")
```

### Causal Discovery

```r
experiment <- RunCausalDiscovery(experiment,
  method = "fci",         # "fci", "pc"
  alpha = 0.05
)
```

## Neural Assemblies

```r
# Detect assemblies
experiment <- RunAssemblies(experiment,
  method = "pca",         # "pca", "ica", "nmf"
  n_assemblies = 10,
  threshold = 2
)

# Get assembly patterns
patterns <- GetAssemblies(experiment)

# Get activation time series
activations <- GetAssemblyActivations(experiment)
```

## State Space Models

### Hidden Markov Models

```r
experiment <- RunHMM(experiment,
  n_states = 3,
  emission = "gaussian"
)

# Get state sequence
states <- GetStates(experiment, model = "hmm")
```

### Switching Linear Dynamical Systems

```r
experiment <- RunSLDS(experiment,
  n_states = 3,
  latent_dim = 5
)
```

## Advanced Analysis

### Topological Data Analysis

```r
experiment <- RunTDA(experiment,
  max_dimension = 2,
  max_scale = 1.0
)
```

### Recurrence Quantification

```r
experiment <- RunRQA(experiment,
  embedding_dim = 3,
  delay = 1,
  threshold = 0.1
)
```

### Changepoint Detection

```r
experiment <- RunChangepoint(experiment,
  method = "pelt",        # "pelt", "binseg", "bocpd"
  penalty = "bic"
)
```

### Optimal Transport

```r
# Compare distributions
distance <- compute_wasserstein(dist1, dist2)

# Manifold alignment
experiment <- RunOTCompare(experiment, reference = ref_experiment)
```

## Machine Learning

### Classification

```r
# Train classifier
model <- train_classifier(experiment,
  labels = cell_types,
  method = "rf",          # "rf", "xgboost", "svm"
  cv_folds = 5
)

# Predict
predictions <- predict(model, new_data)
```

### Clustering

```r
experiment <- RunClustering(experiment,
  method = "kmeans",      # "kmeans", "hierarchical", "dbscan"
  k = 5
)
```

### Conformal Prediction

```r
# Get prediction intervals with coverage guarantees
predictions <- conformal_predict(model, new_data, alpha = 0.1)
```

## Pharmacological Analysis

```r
# Dose-response curve
dr <- fit_dose_response(concentrations, responses,
  model = "4pl"           # 4-parameter logistic
)

# Get EC50
ec50 <- get_ec50(dr)

# Plot dose-response
plot_dose_response(dr)
```

## Visualization

### Trace Plots

```r
plot_traces(experiment, cells = 1:10, assay = "dff")
plot_traces_heatmap(experiment)
```

### Raster Plots

```r
plot_raster(experiment)
plot_raster(experiment, sort_by = "firing_rate")
```

### Network Plots

```r
plot_network(experiment, graph = "connectivity")
plot_adjacency(experiment, graph = "connectivity")
```

### Embedding Plots

```r
plot_embedding(experiment, reduction = "umap", color_by = "cell_type")
```

### Publication Themes

```r
# Apply Cell journal theme
p + theme_cell()

# Apply Nature journal theme
p + theme_nature()
```

## Performance Features

### Parallel Processing

```r
# Enable parallel backend
library(future)
plan(multisession, workers = 4)

# Functions automatically use parallel processing
experiment <- RunConnectivity(experiment, parallel = TRUE)
```

### Out-of-Core Processing

```r
# For datasets larger than RAM
experiment <- CaExperiment(traces, frame_rate = 30, backend = "arrow")

# Process in chunks
experiment <- RunDFF(experiment, chunk_size = 10000)
```

### data.table Operations

```r
# Fast data operations (automatic when data.table is installed)
dt_result <- fast_correlation(traces)
```

## Utility Functions

| Function | Description |
|----------|-------------|
| `get_config()` | Get package configuration |
| `validate_traces()` | Validate trace matrix |
| `downsample_traces()` | Reduce temporal resolution |
| `align_to_stimulus()` | Stimulus-triggered averaging |
| `compute_firing_rate()` | Get firing rates per cell |
| `bin_spikes()` | Bin spike times |
