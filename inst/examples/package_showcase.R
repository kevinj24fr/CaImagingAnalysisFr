# =============================================================================
# CaImagingAnalysisFr Package Showcase
# =============================================================================
#
# This script demonstrates the key features of CaImagingAnalysisFr using
# synthetic calcium imaging data. Run sections interactively to explore
# the package capabilities.
#
# =============================================================================

library(CaImagingAnalysisFr)

# -----------------------------------------------------------------------------
# 1. DATA GENERATION
# -----------------------------------------------------------------------------

cat("\n========== 1. GENERATING SYNTHETIC DATA ==========\n")

# Generate a larger dataset for demonstration
set.seed(42)
data <- generate_synthetic_data(
  n_cells = 30,
  n_time = 6000,    # 10 minutes at 10 Hz
  spike_prob = 0.015
)

# Extract traces matrix (cells x time)
cell_cols <- grep("^Cell_", names(data))
traces <- t(as.matrix(data[, cell_cols]))

cat(sprintf("Generated data: %d cells, %d timepoints\n", nrow(traces), ncol(traces)))
cat(sprintf("Recording duration: %.1f minutes at 10 Hz\n", ncol(traces) / 10 / 60))

# -----------------------------------------------------------------------------
# 2. BASIC PREPROCESSING
# -----------------------------------------------------------------------------

cat("\n========== 2. PREPROCESSING ==========\n")

# Apply calcium correction (dF/F)
corrected <- calcium_correction(data, method = "modern")
cat("Applied dF/F correction\n")

# Convert back to matrix format
if (is.data.frame(corrected)) {
  traces_dff <- t(as.matrix(corrected[, cell_cols]))
} else {
  traces_dff <- corrected
}

cat(sprintf("Trace range: [%.3f, %.3f]\n", min(traces_dff), max(traces_dff)))

# -----------------------------------------------------------------------------
# 3. SPIKE INFERENCE
# -----------------------------------------------------------------------------

cat("\n========== 3. SPIKE INFERENCE ==========\n")

# Detect spikes using threshold method
spikes <- threshold_spike_detection(traces_dff, threshold_sd = 2.0)
cat("Detected spikes using threshold method\n")

# Count events per cell
n_spikes <- spikes$n_events
cat(sprintf("Total events detected: %d\n", sum(n_spikes)))
cat(sprintf("Mean events per cell: %.1f\n", mean(n_spikes)))

# -----------------------------------------------------------------------------
# 4. TRANSIENT DETECTION
# -----------------------------------------------------------------------------

cat("\n========== 4. TRANSIENT DETECTION ==========\n")

# Detect calcium transients with kinetics
transients <- detect_transients(traces_dff, frame_rate = 10)

cat(sprintf("Detected %d transients\n", nrow(transients$events)))
if (nrow(transients$events) > 0) {
  cat(sprintf("Mean amplitude: %.3f\n", mean(transients$events$amplitude)))
  cat(sprintf("Mean duration: %.2f seconds\n", mean(transients$events$duration)))
}

# -----------------------------------------------------------------------------
# 5. NETWORK ANALYSIS
# -----------------------------------------------------------------------------

cat("\n========== 5. NETWORK ANALYSIS ==========\n")

# Compute functional connectivity
connectivity <- functional_connectivity(traces_dff, method = "correlation")
cat("Computed correlation-based connectivity\n")

# Threshold to get significant connections
threshold <- 0.3
n_connections <- sum(abs(connectivity) > threshold &
                     upper.tri(connectivity))
cat(sprintf("Connections above r=%.1f: %d\n", threshold, n_connections))

# Graph metrics
metrics <- graph_metrics(connectivity)
cat(sprintf("Network density: %.3f\n", metrics$density))
cat(sprintf("Mean clustering coefficient: %.3f\n", metrics$clustering_coefficient))

# Community detection
communities <- community_detection(connectivity, method = "louvain")
n_communities <- length(unique(communities$membership))
cat(sprintf("Detected %d communities\n", n_communities))

# -----------------------------------------------------------------------------
# 6. PHARMACOLOGICAL ANALYSIS (Simulated Drug Response)
# -----------------------------------------------------------------------------

cat("\n========== 6. PHARMACOLOGICAL ANALYSIS ==========\n")

# Split data into baseline and "treatment" phases
# (In real experiments, drug would be applied)
n_time <- ncol(traces_dff)
baseline_period <- 1:(n_time %/% 2)
treatment_period <- (n_time %/% 2 + 1):n_time

baseline_traces <- traces_dff[, baseline_period]
treatment_traces <- traces_dff[, treatment_period]

cat(sprintf("Baseline: %d frames, Treatment: %d frames\n",
            length(baseline_period), length(treatment_period)))

# Classify responders
responders <- classify_responders(
  baseline_traces,
  treatment_traces,
  method = "threshold",
  threshold_sd = 2
)

cat("\nResponder Classification:\n")
print(responders$summary)
cat(sprintf("Responder fraction: %.1f%%\n", responders$responder_fraction * 100))

# Analyze baseline calcium
baseline_analysis <- analyze_baseline(baseline_traces, frame_rate = 10)
cat(sprintf("\nMean baseline fluorescence: %.4f\n", baseline_analysis$mean_baseline))

# Compute calcium load
baseline_load <- compute_calcium_load(baseline_traces, frame_rate = 10)
treatment_load <- compute_calcium_load(treatment_traces, frame_rate = 10)
cat(sprintf("Calcium load ratio (treatment/baseline): %.2f\n",
            treatment_load$mean_load / baseline_load$mean_load))

# -----------------------------------------------------------------------------
# 7. TEMPORAL RESPONSE ANALYSIS
# -----------------------------------------------------------------------------

cat("\n========== 7. TEMPORAL RESPONSE ANALYSIS ==========\n")

# Analyze temporal response (using full trace with midpoint as "drug onset")
full_traces <- traces_dff
drug_onset <- n_time %/% 2

temporal <- analyze_temporal_response(
  full_traces,
  drug_onset = drug_onset,
  frame_rate = 10,
  bin_size = 30  # 30-second bins
)

cat(sprintf("Baseline level: %.4f\n", temporal$baseline_level))
cat(sprintf("Peak response: %.4f\n", temporal$peak_response))
cat(sprintf("Time to peak: %.1f seconds\n", temporal$time_to_peak))
cat(sprintf("Response type: %s\n", temporal$response_type))

# -----------------------------------------------------------------------------
# 8. EXCITABILITY METRICS
# -----------------------------------------------------------------------------

cat("\n========== 8. EXCITABILITY METRICS ==========\n")

# Compute excitability for baseline period
excitability <- compute_excitability(baseline_traces, frame_rate = 10)

cat(sprintf("Mean event rate: %.2f per minute\n", excitability$mean_event_rate))
cat(sprintf("Mean refractory period: %.2f seconds\n",
            excitability$mean_refractory))

# -----------------------------------------------------------------------------
# 9. INFORMATION THEORY
# -----------------------------------------------------------------------------

cat("\n========== 9. INFORMATION THEORY ==========\n")

# Compute mutual information between cell pairs (subset for speed)
subset_traces <- traces_dff[1:10, ]
mi_matrix <- mutual_information(subset_traces)

cat("Computed mutual information matrix (10x10 subset)\n
")
cat(sprintf("Mean MI: %.4f bits\n", mean(mi_matrix[upper.tri(mi_matrix)])))
cat(sprintf("Max MI: %.4f bits\n", max(mi_matrix[upper.tri(mi_matrix)])))

# -----------------------------------------------------------------------------
# 10. SPECTRAL ANALYSIS
# -----------------------------------------------------------------------------

cat("\n========== 10. SPECTRAL ANALYSIS ==========\n")

# Compute power spectrum for first cell
spectrum <- compute_power_spectrum(
  traces_dff[1, ],
  frame_rate = 10,
  method = "welch"
)

cat("Computed power spectrum\n")
if (!is.null(spectrum$peak_frequency)) {
  cat(sprintf("Dominant frequency: %.3f Hz\n", spectrum$peak_frequency))
}

# Compute band power
bands <- list(
  slow = c(0.01, 0.1),
  medium = c(0.1, 1),
  fast = c(1, 5)
)
band_power <- compute_band_power(traces_dff[1, ], frame_rate = 10, bands = bands)
cat("Band power distribution:\n")
for (band in names(band_power$band_power)) {
  cat(sprintf("  %s: %.1f%%\n", band, band_power$relative_power[[band]] * 100))
}

# -----------------------------------------------------------------------------
# 11. NEURAL ASSEMBLIES
# -----------------------------------------------------------------------------

cat("\n========== 11. NEURAL ASSEMBLIES ==========\n")

# Detect neural assemblies
assemblies <- detect_assemblies(traces_dff, method = "pca", n_assemblies = 5)

cat(sprintf("Detected %d assemblies\n", assemblies$n_assemblies))
for (i in seq_len(min(3, assemblies$n_assemblies))) {
  members <- which(abs(assemblies$weights[, i]) > 0.3)
  cat(sprintf("Assembly %d: %d core members\n", i, length(members)))
}

# -----------------------------------------------------------------------------
# 12. VARIABILITY ANALYSIS
# -----------------------------------------------------------------------------

cat("\n========== 12. VARIABILITY ANALYSIS ==========\n")

# Compute Fano factor
fano <- compute_fano_factor(traces_dff, window_size = 100)
cat(sprintf("Mean Fano factor: %.3f\n", mean(fano$fano_factors, na.rm = TRUE)))

# Compute coefficient of variation
cv_result <- compute_cv(traces_dff)
cat(sprintf("Mean CV: %.3f\n", mean(cv_result$cv, na.rm = TRUE)))

# Noise correlations (subset for speed)
noise_corr <- compute_noise_correlations(subset_traces)
cat(sprintf("Mean noise correlation: %.3f\n",
            mean(noise_corr[upper.tri(noise_corr)], na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 13. TRANSIENT SHAPE CLASSIFICATION
# -----------------------------------------------------------------------------

cat("\n========== 13. TRANSIENT SHAPE CLASSIFICATION ==========\n")

# Classify transient shapes
shapes <- classify_transient_shapes(traces_dff, frame_rate = 10, n_classes = 4)

if (!is.null(shapes$shape_summary)) {
  cat("Transient shape distribution:\n")
  print(shapes$shape_summary)
}

# -----------------------------------------------------------------------------
# 14. STATE SPACE ANALYSIS
# -----------------------------------------------------------------------------

cat("\n========== 14. STATE SPACE ANALYSIS ==========\n")

# Compute population trajectories
trajectories <- compute_trajectories(traces_dff, method = "pca", n_dims = 3)

cat(sprintf("Trajectory dimensions: %d\n", ncol(trajectories$coords)))
cat(sprintf("Variance explained (3 PCs): %.1f%%\n",
            sum(trajectories$var_explained[1:3]) * 100))

# Estimate intrinsic dimensionality
dim_estimate <- estimate_dimensionality(traces_dff)
cat(sprintf("Estimated intrinsic dimensionality: %d\n", dim_estimate$dimensionality))

# -----------------------------------------------------------------------------
# 15. SUMMARY STATISTICS
# -----------------------------------------------------------------------------

cat("\n========== 15. ANALYSIS SUMMARY ==========\n")

cat("\nDataset Overview:\n")
cat(sprintf("  Cells: %d\n", nrow(traces_dff)))
cat(sprintf("  Duration: %.1f minutes\n", ncol(traces_dff) / 10 / 60))
cat(sprintf("  Frame rate: 10 Hz\n"))

cat("\nActivity Metrics:\n")
cat(sprintf("  Total transients: %d\n", nrow(transients$events)))
cat(sprintf("  Mean rate: %.2f events/cell/minute\n",
            nrow(transients$events) / nrow(traces_dff) / (ncol(traces_dff) / 10 / 60)))

cat("\nNetwork Properties:\n")
cat(sprintf("  Communities: %d\n", n_communities))
cat(sprintf("  Clustering: %.3f\n", metrics$clustering_coefficient))

cat("\nPharmacology (simulated):\n")
cat(sprintf("  Responders: %.1f%%\n", responders$responder_fraction * 100))

cat("\n========== SHOWCASE COMPLETE ==========\n")

# -----------------------------------------------------------------------------
# OPTIONAL: VISUALIZATION (requires graphics device)
# -----------------------------------------------------------------------------

if (interactive()) {
  cat("\nGenerating plots...\n")

  # Plot 1: Sample traces
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

  # Raw traces
  time_sec <- seq_len(ncol(traces_dff)) / 10
  plot(time_sec[1:1000], traces_dff[1, 1:1000], type = "l",
       xlab = "Time (s)", ylab = "dF/F",
       main = "Cell 1 Trace (first 100s)", col = "darkblue")

  # Correlation matrix
  image(connectivity[1:15, 1:15],
        main = "Connectivity (15 cells)",
        col = hcl.colors(50, "RdBu", rev = TRUE))

  # Event rate histogram
  hist(n_spikes, breaks = 15, main = "Events per Cell",
       xlab = "Number of Events", col = "steelblue")

  # Trajectory (first 2 PCs)
  plot(trajectories$coords[, 1], trajectories$coords[, 2],
       type = "l", xlab = "PC1", ylab = "PC2",
       main = "Population Trajectory", col = "darkgreen")
  points(trajectories$coords[1, 1], trajectories$coords[1, 2],
         pch = 19, col = "red", cex = 2)

  par(mfrow = c(1, 1))
}

cat("\nRun this script interactively to see visualizations.\n")
cat("For more examples, see: vignette('getting-started')\n")
