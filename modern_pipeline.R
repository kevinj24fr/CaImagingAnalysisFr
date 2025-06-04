# Modern Calcium Imaging Analysis Pipeline - Demo

# This demonstration script loads the package functions and runs the
# basic workflow on the bundled synthetic dataset.  It illustrates how
# to generate traces, perform baseline correction, spike inference, and
# visualize a single cell.

library(ggplot2)

# Load functions if running interactively without installation
if (!requireNamespace("CaImagingAnalysisFr", quietly = TRUE)) {
  source(file.path("R", "pipeline.R"))
}

# Use bundled synthetic data if available
if (file.exists("data/synthetic_data.csv")) {
  raw <- read.csv("data/synthetic_data.csv")
} else {
  raw <- generate_synthetic_data(n_cells = 5, n_time = 500)
}

# Run the full analysis (AR(2) model)
res <- run_analysis(raw, penalty = 1, ar_order = 2)

# Display plots
print(res$cell_plot)
print(res$raster_plot)
