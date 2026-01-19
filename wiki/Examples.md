# Examples

Real-world examples demonstrating CaImagingAnalysisFr workflows.

## Example 1: Basic Calcium Imaging Analysis

A complete workflow from raw traces to basic characterization.

```r
library(CaImagingAnalysisFr)

# Load data
traces <- read_tiff_stack("recording.tif")
# Or if you have pre-extracted traces:
# traces <- as.matrix(read.csv("traces.csv", row.names = 1))

# Create experiment and run standard pipeline
experiment <- CaExperiment(traces, frame_rate = 30) |>
  # Preprocessing
  RunCorrection(method = "neuropil", neuropil_factor = 0.7) |>
  RunDFF(method = "rolling", window = 300) |>

  # Spike detection
  RunSpikes(method = "threshold", threshold = 2.5) |>

  # Dimensionality reduction
  RunPCA(n_components = 20) |>
  RunUMAP(n_neighbors = 15)

# Visualize results
plot_traces(experiment, cells = 1:5)
plot_raster(experiment)
plot_embedding(experiment, reduction = "umap")

# Get summary statistics
firing_rates <- compute_firing_rate(experiment)
summary(firing_rates)
```

## Example 2: Functional Connectivity Analysis

Analyze neural network structure and communities.

```r
library(CaImagingAnalysisFr)
library(igraph)

# Load and preprocess
experiment <- CaExperiment(traces, frame_rate = 30) |>
  RunCorrection(method = "neuropil") |>
  RunDFF(method = "rolling")

# Compute multiple connectivity measures
experiment <- experiment |>
  RunConnectivity(method = "correlation", threshold = 0.3) |>
  RunConnectivity(method = "partial", threshold = 0.2)

# Get the correlation-based graph
g <- GetGraph(experiment, "connectivity")

# Compute graph metrics
metrics <- graph_metrics(g)
print(metrics)

# Detect communities
experiment <- RunCommunities(experiment, method = "louvain")
communities <- GetCommunities(experiment)

# Visualize network
plot_network(experiment,
  graph = "connectivity",
  color_by = "community",
  layout = "fr"
)

# Plot adjacency matrix sorted by community
plot_adjacency(experiment,
  graph = "connectivity",
  sort_by = "community"
)

# Identify hub neurons
hub_scores <- hub_score(g)$vector
top_hubs <- order(hub_scores, decreasing = TRUE)[1:5]
cat("Top hub neurons:", top_hubs, "\n")
```

## Example 3: Neural Assembly Detection

Find co-activation patterns in population activity.

```r
library(CaImagingAnalysisFr)

# Load and preprocess
experiment <- CaExperiment(traces, frame_rate = 30) |>
  RunCorrection(method = "neuropil") |>
  RunDFF(method = "rolling") |>
  RunSpikes(method = "threshold")

# Detect assemblies using PCA
experiment <- RunAssemblies(experiment,
  method = "pca",
  n_assemblies = 10,
  threshold = 2.0
)

# Get assembly patterns (cells x assemblies)
patterns <- GetAssemblies(experiment)

# Get assembly activations over time
activations <- GetAssemblyActivations(experiment)

# Visualize assembly patterns
plot_assemblies(experiment)

# Plot assembly activation time series
plot_assembly_activations(experiment, assemblies = 1:3)

# Identify cells belonging to each assembly
for (i in 1:3) {
  member_cells <- which(abs(patterns[, i]) > 0.3)
  cat(sprintf("Assembly %d members: %s\n", i,
              paste(member_cells, collapse = ", ")))
}
```

## Example 4: Stimulus-Response Analysis

Analyze neural responses to experimental stimuli.

```r
library(CaImagingAnalysisFr)

# Load data with stimulus timing
experiment <- CaExperiment(traces, frame_rate = 30) |>
  RunCorrection(method = "neuropil") |>
  RunDFF(method = "rolling")

# Define stimulus times (in frames)
stim_times <- c(100, 400, 700, 1000, 1300)

# Compute stimulus-triggered averages
sta <- align_to_stimulus(experiment,
  stim_times = stim_times,
  pre_frames = 30,    # 1 second before
  post_frames = 90    # 3 seconds after
)

# Plot average response
plot_stimulus_response(sta)

# Identify responsive cells
responsive <- identify_responsive_cells(experiment,
  stim_times = stim_times,
  method = "ttest",
  alpha = 0.05
)

cat("Responsive cells:", sum(responsive), "out of", length(responsive), "\n")

# Compute tuning curves (if multiple stimulus types)
# stim_values <- c(0, 45, 90, 135, 180, 225, 270, 315)  # orientations
# tuning <- compute_tuning_curves(experiment, stim_times, stim_values)
# plot_tuning_curves(tuning, cells = 1:5)
```

## Example 5: State Space Analysis with HMM

Identify hidden neural states using Hidden Markov Models.

```r
library(CaImagingAnalysisFr)

# Load and preprocess
experiment <- CaExperiment(traces, frame_rate = 30) |>
  RunCorrection(method = "neuropil") |>
  RunDFF(method = "rolling") |>
  RunPCA(n_components = 10)

# Fit HMM on PCA-reduced data
experiment <- RunHMM(experiment,
  n_states = 4,
  emission = "gaussian",
  n_iter = 100
)

# Get state sequence
states <- GetStates(experiment, model = "hmm")

# Plot state transitions
plot_states(experiment, model = "hmm")

# Get transition matrix
trans_mat <- GetTransitionMatrix(experiment, model = "hmm")
print(round(trans_mat, 3))

# Compute state occupancy
state_occupancy <- table(states) / length(states)
print(state_occupancy)

# Visualize states in PCA space
plot_embedding(experiment,
  reduction = "pca",
  dims = c(1, 2),
  color_by = "state"
)
```

## Example 6: Cross-Session Cell Registration

Match cells across multiple recording sessions.

```r
library(CaImagingAnalysisFr)

# Load data from two sessions
session1 <- CaExperiment(traces1, frame_rate = 30)
session2 <- CaExperiment(traces2, frame_rate = 30)

# If you have ROI information
rois1 <- ROISet(masks = masks1, centroids = centroids1)
rois2 <- ROISet(masks = masks2, centroids = centroids2)

# Register cells across sessions
registration <- register_cells(
  session1, session2,
  rois1 = rois1,
  rois2 = rois2,
  method = "centroid",      # or "correlation", "shape"
  max_distance = 10         # microns
)

# Get matched cell pairs
matched_pairs <- registration$matches
cat("Matched", nrow(matched_pairs), "cells across sessions\n")

# Create aligned dataset
aligned <- align_sessions(session1, session2, registration)
```

## Example 7: Pharmacological Analysis

Analyze drug effects on neural activity.

```r
library(CaImagingAnalysisFr)

# Load baseline and drug condition data
baseline <- CaExperiment(traces_baseline, frame_rate = 30) |>
  RunDFF(method = "rolling") |>
  RunSpikes(method = "threshold")

drug <- CaExperiment(traces_drug, frame_rate = 30) |>
  RunDFF(method = "rolling") |>
  RunSpikes(method = "threshold")

# Compare firing rates
fr_baseline <- compute_firing_rate(baseline)
fr_drug <- compute_firing_rate(drug)

# Statistical comparison
wilcox.test(fr_baseline, fr_drug, paired = TRUE)

# Dose-response analysis (if multiple concentrations)
concentrations <- c(0, 0.1, 1, 10, 100)  # nM
responses <- c(1.0, 0.95, 0.7, 0.3, 0.1)  # normalized activity

# Fit dose-response curve
dr_fit <- fit_dose_response(concentrations, responses, model = "4pl")

# Get IC50
ic50 <- get_ic50(dr_fit)
cat("IC50:", ic50, "nM\n")

# Plot dose-response curve
plot_dose_response(dr_fit)
```

## Example 8: Publication-Ready Figures

Create figures suitable for journal submission.

```r
library(CaImagingAnalysisFr)
library(ggplot2)
library(patchwork)

# Prepare data
experiment <- CaExperiment(traces, frame_rate = 30) |>
  RunDFF(method = "rolling") |>
  RunSpikes(method = "threshold") |>
  RunPCA(n_components = 20) |>
  RunUMAP() |>
  RunConnectivity(method = "correlation")

# Create individual panels
p1 <- plot_traces(experiment, cells = 1:3) +
  theme_cell() +
  labs(title = "A")

p2 <- plot_raster(experiment) +
  theme_cell() +
  labs(title = "B")

p3 <- plot_embedding(experiment, reduction = "umap") +
  theme_cell() +
  labs(title = "C")

p4 <- plot_network(experiment, graph = "connectivity") +
  theme_cell() +
  labs(title = "D")

# Combine into multi-panel figure
figure <- (p1 | p2) / (p3 | p4)

# Save at publication resolution
ggsave("figure1.pdf", figure, width = 180, height = 150, units = "mm")
ggsave("figure1.png", figure, width = 180, height = 150, units = "mm", dpi = 300)
```

## Example 9: Batch Processing Multiple Recordings

Process many recordings with consistent parameters.

```r
library(CaImagingAnalysisFr)
library(future)

# Enable parallel processing
plan(multisession, workers = 4)

# List of recording files
files <- list.files("data/", pattern = "*.tif", full.names = TRUE)

# Define analysis function
analyze_recording <- function(file) {
  traces <- read_tiff_stack(file)

  experiment <- CaExperiment(traces, frame_rate = 30) |>
    RunCorrection(method = "neuropil") |>
    RunDFF(method = "rolling") |>
    RunSpikes(method = "threshold") |>
    RunConnectivity(method = "correlation")

  # Extract summary metrics
  list(
    file = basename(file),
    n_cells = ncol(GetTraces(experiment)),
    mean_firing_rate = mean(compute_firing_rate(experiment)),
    mean_connectivity = mean(GetAdjacency(experiment, "connectivity"))
  )
}

# Process all files
results <- future_map(files, analyze_recording)

# Combine results
summary_df <- bind_rows(results)
print(summary_df)
```

## Example 10: Integration with targets

Create a reproducible analysis pipeline with the targets package.

```r
# _targets.R file
library(targets)
library(CaImagingAnalysisFr)

tar_option_set(packages = c("CaImagingAnalysisFr"))

list(
  tar_target(raw_data, read_tiff_stack("recording.tif")),

  tar_target(experiment, {
    CaExperiment(raw_data, frame_rate = 30) |>
      RunCorrection(method = "neuropil") |>
      RunDFF(method = "rolling")
  }),

  tar_target(with_spikes, RunSpikes(experiment, method = "threshold")),

  tar_target(with_network, RunConnectivity(with_spikes, method = "correlation")),

  tar_target(summary_stats, {
    list(
      n_cells = ncol(GetTraces(with_network)),
      firing_rates = compute_firing_rate(with_network),
      graph_metrics = graph_metrics(GetGraph(with_network, "connectivity"))
    )
  }),

  tar_target(figure, {
    p <- plot_network(with_network)
    ggsave("figures/network.pdf", p)
    "figures/network.pdf"
  })
)
```

Run with:
```r
targets::tar_make()
targets::tar_visnetwork()  # Visualize pipeline
```
