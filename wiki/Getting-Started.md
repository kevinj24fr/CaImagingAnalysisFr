# Getting Started

This guide walks you through your first calcium imaging analysis with CaImagingAnalysisFr.

## Loading the Package

```r
library(CaImagingAnalysisFr)
```

## Loading Your Data

### From a TIFF Stack

```r
# Load a multi-page TIFF file
movie <- read_tiff_stack("path/to/recording.tif")

# Check dimensions (pixels x pixels x frames)
dim(movie)
```

### From Pre-extracted Traces

If you already have extracted calcium traces (e.g., from Suite2P or CaImAn):

```r
# Load from CSV (cells x timepoints)
traces <- read.csv("traces.csv", row.names = 1)
traces <- as.matrix(traces)

# Or create synthetic data for testing
traces <- matrix(rnorm(500 * 1000), nrow = 50, ncol = 1000)
```

### From HDF5/NWB Files

```r
# Requires rhdf5 package
traces <- read_hdf5_traces("recording.h5", dataset = "/traces")
```

## Creating a CaExperiment Object

The `CaExperiment` object is the central data structure that holds all your data and analysis results:

```r
# Create experiment with frame rate
experiment <- CaExperiment(
  traces = traces,
  frame_rate = 30  # Hz
)

# View summary
print(experiment)
```

## Basic Analysis Pipeline

### Step 1: Signal Correction

Remove background and neuropil contamination:

```r
experiment <- experiment |>
  RunCorrection(method = "neuropil", neuropil_factor = 0.7)
```

### Step 2: Compute dF/F

Convert raw fluorescence to normalized dF/F:

```r
experiment <- experiment |>
  RunDFF(method = "rolling", window = 300)  # 10s window at 30Hz
```

### Step 3: Spike Inference

Detect neural activity events:
```r
experiment <- experiment |>
  RunSpikes(method = "threshold", threshold = 2.5)  # 2.5 SD threshold
```

### Step 4: Dimensionality Reduction

Reduce dimensions for visualization and analysis:

```r
experiment <- experiment |>
  RunPCA(n_components = 20) |>
  RunUMAP(n_neighbors = 15)
```

### Step 5: Network Analysis

Compute functional connectivity:

```r
experiment <- experiment |>
  RunConnectivity(method = "correlation", threshold = 0.3)
```

## Accessing Results

### Get Processed Traces

```r
# Get dF/F traces
dff <- GetTraces(experiment, assay = "dff")

# Get raw traces
raw <- GetTraces(experiment, assay = "raw")
```

### Get Spike Data

```r
spikes <- GetSpikes(experiment)
```

### Get Embeddings

```r
pca_coords <- GetEmbedding(experiment, reduction = "pca")
umap_coords <- GetEmbedding(experiment, reduction = "umap")
```

### Get Network Graph

```r
graph <- GetGraph(experiment, graph_name = "connectivity")
```

## Basic Visualization

### Plot Traces

```r
# Plot first 5 cells
plot_traces(experiment, cells = 1:5)
```

### Plot Raster

```r
# Spike raster plot
plot_raster(experiment)
```

### Plot UMAP

```r
# Colored by some metadata
plot_embedding(experiment, reduction = "umap")
```

### Plot Network

```r
# Functional connectivity network
plot_network(experiment, graph = "connectivity")
```

## Complete Pipeline Example

Here's a complete analysis pipeline in one code block:

```r
library(CaImagingAnalysisFr)

# Load and analyze
experiment <- CaExperiment(traces, frame_rate = 30) |>
  # Preprocessing
  RunCorrection(method = "neuropil") |>
  RunDFF(method = "rolling") |>

  # Spike detection
  RunSpikes(method = "threshold") |>

  # Dimensionality reduction
  RunPCA(n_components = 20) |>
  RunUMAP() |>


  # Network analysis
  RunConnectivity(method = "correlation") |>

  # Neural assemblies
  RunAssemblies(method = "pca", n_assemblies = 5)

# Export pipeline for reproducibility
ExportCommands(experiment, "my_pipeline.R")
```

## Adding Metadata

You can add cell and experiment metadata:

```r
# Add cell-level metadata
experiment <- AddMetadata(experiment,
  cell_type = c(rep("excitatory", 40), rep("inhibitory", 10)),
  layer = sample(c("L2/3", "L4", "L5"), 50, replace = TRUE)
)

# Add experiment-level metadata
experiment@misc$condition <- "control"
experiment@misc$animal_id <- "mouse_001"
```

## Saving and Loading

```r
# Save experiment
saveRDS(experiment, "my_experiment.rds")

# Load experiment
experiment <- readRDS("my_experiment.rds")
```

## Next Steps

- Learn about [Core Concepts](Core-Concepts) to understand the data structure
- Explore the [API Overview](API-Overview) for all available functions
- See [Examples](Examples) for real-world analysis workflows
- Check the [FAQ](FAQ) for common questions
