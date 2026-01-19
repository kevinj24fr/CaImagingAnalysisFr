# Core Concepts

Understanding the fundamental concepts and data structures in CaImagingAnalysisFr.

## The CaExperiment Object

The `CaExperiment` object is inspired by Seurat's design philosophy - a single object that contains all your data and analysis results. This enables:

- **Unified Interface:** One object, consistent methods
- **Reproducibility:** Track all processing steps
- **Convenience:** Easy data access and manipulation

### Structure

```
CaExperiment
├── assays           # Different data representations
│   ├── raw          # Raw fluorescence traces
│   ├── corrected    # After neuropil correction
│   └── dff          # Normalized dF/F
├── spikes           # Inferred spike data
├── reductions       # Dimensionality reductions
│   ├── pca          # Principal components
│   └── umap         # UMAP coordinates
├── graphs           # Network graphs
│   └── connectivity # Functional connectivity
├── assemblies       # Neural assembly patterns
├── metadata         # Cell-level annotations
├── misc             # Experiment-level data
└── commands         # Processing history
```

### Creating a CaExperiment

```r
# Basic creation
experiment <- CaExperiment(traces, frame_rate = 30)

# With metadata
experiment <- CaExperiment(
  traces = traces,
  frame_rate = 30,
  metadata = data.frame(
    cell_type = cell_types,
    roi_area = areas
  )
)
```

## S7 Type-Safe Classes

For performance-critical applications, CaImagingAnalysisFr provides S7 classes with strict type validation:

### CalciumMovie

Represents raw imaging data:

```r
movie <- CalciumMovie(
  data = array_3d,      # height x width x frames
  frame_rate = 30,
  pixel_size = 0.5      # microns
)
```

### CalciumTraces

Represents extracted fluorescence traces:

```r
traces <- CalciumTraces(
  data = matrix_data,   # cells x timepoints
  frame_rate = 30,
  cell_ids = paste0("cell_", 1:n_cells)
)
```

### SpikeTrains

Represents inferred spike times:

```r
spikes <- SpikeTrains(
  spike_times = list_of_vectors,  # spike times per cell
  cell_ids = cell_ids,
  duration = total_duration
)
```

### ROISet

Represents regions of interest:

```r
rois <- ROISet(
  masks = list_of_masks,
  centroids = centroid_matrix,
  areas = area_vector
)
```

## Assays

Assays store different representations of your trace data:

| Assay | Description |
|-------|-------------|
| `raw` | Original fluorescence values |
| `corrected` | After background/neuropil correction |
| `dff` | Normalized dF/F values |
| `zscore` | Z-scored traces |
| `denoised` | After denoising (e.g., wavelet) |

### Working with Assays

```r
# Get data from specific assay
dff_data <- GetTraces(experiment, assay = "dff")

# Set default assay
DefaultAssay(experiment) <- "dff"

# List available assays
Assays(experiment)
```

## Reductions

Dimensionality reductions for visualization and downstream analysis:

| Reduction | Description |
|-----------|-------------|
| `pca` | Principal Component Analysis |
| `umap` | UMAP embedding |
| `tsne` | t-SNE embedding |
| `ica` | Independent Component Analysis |
| `gpfa` | Gaussian Process Factor Analysis |

### Working with Reductions

```r
# Run reduction
experiment <- RunPCA(experiment, n_components = 50)

# Get coordinates
pca_coords <- GetEmbedding(experiment, reduction = "pca")

# Get loadings (for PCA)
loadings <- GetLoadings(experiment, reduction = "pca")

# List reductions
Reductions(experiment)
```

## Graphs

Network representations of neural relationships:

| Graph | Description |
|-------|-------------|
| `connectivity` | Functional connectivity (correlation-based) |
| `granger` | Granger causality network |
| `transfer_entropy` | Information-theoretic connectivity |
| `partial_corr` | Partial correlation network |

### Working with Graphs

```r
# Compute connectivity
experiment <- RunConnectivity(experiment, method = "correlation")

# Get igraph object
g <- GetGraph(experiment, "connectivity")

# Compute graph metrics
metrics <- graph_metrics(g)
```

## Neural Assemblies

Co-activation patterns detected in your data:

```r
# Detect assemblies
experiment <- RunAssemblies(experiment, method = "pca", n_assemblies = 10)

# Get assembly patterns
patterns <- GetAssemblies(experiment)

# Get assembly activations over time
activations <- GetAssemblyActivations(experiment)
```

## Metadata

Cell-level annotations stored with your experiment:

```r
# Add metadata
experiment <- AddMetadata(experiment,
  cell_type = types,
  region = regions
)

# Access metadata
meta <- experiment@metadata

# Use in plots
plot_embedding(experiment, color_by = "cell_type")
```

## The Pipe Workflow

CaImagingAnalysisFr is designed for pipe-based workflows:

```r
experiment <- CaExperiment(traces, frame_rate = 30) |>
  RunCorrection() |>
  RunDFF() |>
  RunSpikes() |>
  RunPCA() |>
  RunConnectivity()
```

Each `Run*` function:
1. Takes a CaExperiment as input
2. Performs an analysis step
3. Stores results in the appropriate slot
4. Returns the modified CaExperiment
5. Logs the command for reproducibility

## Command History

All processing steps are logged:

```r
# View command history
Commands(experiment)

# Export as reproducible script
ExportCommands(experiment, "pipeline.R")
```

## Data Access Pattern Summary

| Data Type | Accessor Function |
|-----------|-------------------|
| Traces | `GetTraces(exp, assay = "...")` |
| Spikes | `GetSpikes(exp)` |
| Embeddings | `GetEmbedding(exp, reduction = "...")` |
| Graphs | `GetGraph(exp, graph_name = "...")` |
| Assemblies | `GetAssemblies(exp)` |
| Metadata | `exp@metadata` or `GetMetadata(exp)` |

## Best Practices

1. **Always specify frame_rate:** Required for temporal calculations
2. **Use pipes for clarity:** Chain operations for readable code
3. **Check intermediate results:** Use `print(experiment)` to see current state
4. **Export commands:** Save your pipeline for reproducibility
5. **Use appropriate assays:** Run analyses on dF/F, not raw traces
