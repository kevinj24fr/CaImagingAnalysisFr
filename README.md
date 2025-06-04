# CaImagingAnalysisFr

This repository contains a lightweight collection of functions for working with
calcium imaging traces.  The package provides utilities to generate synthetic
traces, perform baseline correction, infer spikes using the OASIS algorithm, and
visualise results.
The analysis functions are fully parameterised so you can control deconvolution settings and smoothing spans. An example R Markdown report (`analysis_report.Rmd`) demonstrates the workflow and produces publication-quality plots.


## Quick start

```R
# load package functions
source('R/pipeline.R')

# generate synthetic data
df <- generate_synthetic_data(n_cells = 3, n_time = 200)

# run the full analysis with spike inference (AR(2) model)
res <- run_analysis(df, penalty = 1, ar_order = 2)
res$cell_plot
res$raster_plot
```
