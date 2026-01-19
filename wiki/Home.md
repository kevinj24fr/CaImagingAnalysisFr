# CaImagingAnalysisFr Wiki

Welcome to the **CaImagingAnalysisFr** wiki - your comprehensive guide to calcium imaging analysis in R.

## What is CaImagingAnalysisFr?

CaImagingAnalysisFr is an R package providing a complete, unified workflow for calcium imaging analysis in neuroscience research. It solves the fragmentation problem where calcium imaging analysis traditionally requires multiple tools (Python, MATLAB, various R packages) with different data formats and interfaces.

**Core Philosophy:** One R package, one workflow - from raw pixels to publication-ready figures.

## Key Features

| Category | Capabilities |
|----------|-------------|
| **Signal Processing** | Motion correction, neuropil subtraction, spike inference, dF/F computation |
| **Network Analysis** | Functional connectivity, graph metrics, community detection, causal discovery |
| **Population Analysis** | GPFA, demixed PCA, neural assemblies, tensor decomposition |
| **Advanced Methods** | HMM/SLDS, optimal transport, topological data analysis, neural ODEs |
| **Machine Learning** | Classification, regression, clustering, conformal prediction |
| **Pharmacology** | Dose-response analysis, concentration effects, drug interactions |
| **Visualization** | Publication-ready plots with Cell/Nature journal themes |
| **Performance** | S7 classes, data.table, Rcpp acceleration, Arrow/DuckDB for big data |

## Quick Start

```r
library(CaImagingAnalysisFr)

# Create unified analysis object with Seurat-style workflow
experiment <- CaExperiment(traces, frame_rate = 30) |>
  RunCorrection(method = "neuropil") |>
  RunDFF(method = "rolling") |>
  RunSpikes(method = "threshold") |>
  RunPCA(n_components = 20) |>
  RunConnectivity(method = "correlation")

# Access results
GetTraces(experiment, assay = "dff")
GetSpikes(experiment)
```

## Getting Help

- **[Installation Guide](Installation)** - How to install the package
- **[Getting Started](Getting-Started)** - Your first analysis workflow
- **[Core Concepts](Core-Concepts)** - Understanding the CaExperiment object
- **[API Overview](API-Overview)** - Module-by-module function reference
- **[Examples](Examples)** - Real-world analysis examples
- **[FAQ](FAQ)** - Frequently asked questions

## Requirements

- R >= 4.1.0
- Core dependencies: ggplot2 (>= 3.4.0), igraph, MASS
- Optional: S7, data.table, Rcpp, arrow, duckdb for enhanced performance

## Citation

If you use CaImagingAnalysisFr in your research, please cite:

```
@software{caimaginganalysisfr,
  title = {CaImagingAnalysisFr: Unified Calcium Imaging Analysis in R},
  author = {Kevin J},
  year = {2024},
  url = {https://github.com/kevinj24fr/CaImagingAnalysisFr}
}
```

## License

This project is open source. See the repository for license details.
