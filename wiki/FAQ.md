# Frequently Asked Questions

## General Questions

### What is CaImagingAnalysisFr?

CaImagingAnalysisFr is an R package for calcium imaging analysis that provides a unified, Seurat-style workflow. It covers the entire analysis pipeline from raw data to publication figures, eliminating the need to switch between multiple tools.

### Why should I use this instead of Suite2P or CaImAn?

CaImagingAnalysisFr complements rather than replaces these tools:

- **Suite2P/CaImAn:** Excel at cell segmentation and motion correction
- **CaImagingAnalysisFr:** Excels at downstream analysis (network analysis, assemblies, state-space models, visualization)

You can import Suite2P or CaImAn output directly into CaImagingAnalysisFr for downstream analysis.

### What R version do I need?

R 4.1.0 or higher is required.

### Is this package on CRAN?

Not yet. Install from GitHub using `remotes::install_github()`.

## Data & Input

### What data formats are supported?

- **TIFF stacks:** Multi-page TIFF files
- **HDF5:** Including NWB format
- **MATLAB:** .mat files
- **CSV:** Standard comma-separated values
- **Suite2P output:** Direct import of F.npy, etc.
- **CaImAn output:** Direct import of estimates

### How do I load my data?

```r
# TIFF stack
traces <- read_tiff_stack("recording.tif")

# CSV file
traces <- as.matrix(read.csv("traces.csv", row.names = 1))

# HDF5
traces <- read_hdf5_traces("data.h5", dataset = "/traces")

# Suite2P
data <- read_suite2p("suite2p/plane0/")
```

### What should my data look like?

For `CaExperiment`, provide a matrix where:
- **Rows:** Cells
- **Columns:** Timepoints

```r
dim(traces)
# [1]   50 10000  # 50 cells, 10000 frames
```

### How do I handle multiple imaging planes?

Process each plane separately, then combine:

```r
experiments <- lapply(planes, function(p) {
  CaExperiment(p, frame_rate = 30) |>
    RunDFF(method = "rolling")
})

# Combine later for population analysis
combined <- merge_experiments(experiments)
```

## Analysis Questions

### What is dF/F and why do I need it?

dF/F (delta F over F) normalizes fluorescence changes relative to baseline, making signals comparable across cells with different absolute fluorescence levels.

```r
experiment <- RunDFF(experiment, method = "rolling", window = 300)
```

### Which spike inference method should I use?

| Method | Best For | Speed |
|--------|----------|-------|
| `threshold` | Quick analysis, high SNR data | Fast |
| `oasis` | Accurate spike times | Medium |
| `mlspike` | Sparse firing, single spikes | Slow |

Start with `threshold` for exploration, then use `oasis` for final analysis.

### How do I choose the number of PCA components?

```r
experiment <- RunPCA(experiment, n_components = 50)

# Check variance explained
var_explained <- GetVarianceExplained(experiment, reduction = "pca")
plot(cumsum(var_explained), type = "l")
abline(h = 0.9)  # Find elbow or 90% variance
```

### What connectivity method should I use?

| Method | Use When |
|--------|----------|
| `correlation` | Simple, interpretable |
| `partial` | Want to remove indirect connections |
| `transfer_entropy` | Looking for directional/causal relationships |
| `granger` | Time-lagged causal analysis |

### How do I detect neural assemblies?

```r
experiment <- RunAssemblies(experiment,
  method = "pca",        # or "ica", "nmf"
  n_assemblies = 10
)
```

Choose `n_assemblies` based on the number of PCA components with significant variance, or use cross-validation.

## Performance Questions

### My analysis is slow. How can I speed it up?

1. **Enable parallel processing:**
   ```r
   library(future)
   plan(multisession, workers = 4)
   ```

2. **Install performance packages:**
   ```r
   install.packages(c("data.table", "Rcpp"))
   ```

3. **Reduce data if exploring:**
   ```r
   # Downsample temporally
   traces_ds <- downsample_traces(traces, factor = 2)
   ```

### My dataset is larger than RAM. What can I do?

Use the Arrow backend for out-of-core processing:

```r
install.packages(c("arrow", "duckdb"))

experiment <- CaExperiment(traces, frame_rate = 30, backend = "arrow")
```

### How do I process multiple recordings efficiently?

Use parallel batch processing:

```r
library(future)
plan(multisession, workers = 4)

results <- future_map(files, analyze_recording)
```

## Visualization Questions

### How do I make publication-quality figures?

Use the built-in themes:

```r
p <- plot_network(experiment) + theme_cell()
ggsave("figure.pdf", p, width = 180, height = 150, units = "mm")
```

### Can I customize the plots?

Yes, all plots return ggplot2 objects:

```r
p <- plot_traces(experiment, cells = 1:5) +
  scale_color_viridis_d() +
  labs(title = "My Custom Title") +
  theme(axis.text = element_text(size = 12))
```

### How do I create multi-panel figures?

Use the patchwork package:

```r
library(patchwork)

p1 <- plot_traces(experiment)
p2 <- plot_raster(experiment)
p3 <- plot_network(experiment)

combined <- p1 / (p2 | p3)
```

## Troubleshooting

### Error: "could not find function"

The function might be from a suggested package. Install it:

```r
# Check which package provides the function
?function_name

# Install the required package
install.packages("package_name")
```

### Error: "Python module not found"

For Python-dependent features (OASIS, etc.):

```r
library(reticulate)
conda_create("caimaging")
py_install("oasis-deconv", envname = "caimaging")
use_condaenv("caimaging")
```

### Error: "object of type 'closure' is not subsettable"

You're probably trying to subset a function instead of its result:

```r
# Wrong
GetTraces[1:10, ]

# Right
GetTraces(experiment)[1:10, ]
```

### My connectivity matrix has NaN values

This usually happens with constant traces (no variance). Remove or filter these cells:

```r
# Find cells with variance
variances <- apply(GetTraces(experiment, "dff"), 1, var)
good_cells <- which(variances > 0.001)

# Subset experiment
experiment <- subset(experiment, cells = good_cells)
```

### HMM fitting doesn't converge

Try:
1. Reduce the number of states
2. Initialize with k-means
3. Increase iterations

```r
experiment <- RunHMM(experiment,
  n_states = 3,           # Fewer states
  init = "kmeans",        # Better initialization
  n_iter = 200            # More iterations
)
```

## Contributing

### How can I contribute?

1. Report bugs via [GitHub Issues](https://github.com/kevinj24fr/CaImagingAnalysisFr/issues)
2. Submit pull requests for fixes or features
3. Improve documentation
4. Share example analyses

### Where can I get help?

1. Check this FAQ
2. Read the [documentation](https://kevinj24fr.github.io/CaImagingAnalysisFr/)
3. Search [GitHub Issues](https://github.com/kevinj24fr/CaImagingAnalysisFr/issues)
4. Open a new issue with a reproducible example
