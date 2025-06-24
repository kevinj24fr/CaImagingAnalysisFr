## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load---------------------------------------------------------------------
library(CaImagingAnalysisFr)

## ----data---------------------------------------------------------------------
# Generate synthetic data for demonstration
set.seed(123)
synthetic_data <- generate_synthetic_data(
  n_cells = 5,
  n_time = 300,
  spike_prob = 0.05
)

# Extract traces
traces <- as.matrix(synthetic_data[, grep("^Cell_", names(synthetic_data))])
print(paste("Working with", ncol(traces), "cells and", nrow(traces), "timepoints"))

## ----preprocessing------------------------------------------------------------
# Apply calcium correction
corrected_traces <- calcium_correction(synthetic_data, method = "modern")

# Compare raw vs corrected
par(mfrow = c(2, 1))
plot(as.numeric(traces[1, ]), type = "l", main = "Raw Trace", ylab = "Fluorescence")
plot(as.numeric(corrected_traces[1, ]), type = "l", main = "Corrected Trace", ylab = "Delta F/F")

