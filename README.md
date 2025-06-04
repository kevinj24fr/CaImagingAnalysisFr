# CaImagingAnalysisFr

This repository contains a lightweight collection of functions for working with
calcium imaging traces.  The package provides utilities to generate synthetic
traces, perform baseline correction, infer spikes using the OASIS algorithm, and
visualise results.

## Quick start

```R
# load package functions
source('R/pipeline.R')

# generate synthetic data
df <- generate_synthetic_data(n_cells = 3, n_time = 200)

# correct and plot
df_corrected <- calcium_correction(df)
plot_cell_trace(df_corrected, 'Cell_1')
```
