# =============================================================================
# CaImagingAnalysisFr Package Showcase
# =============================================================================
#
# Complete demonstration of CaImagingAnalysisFr capabilities:
# - CaExperiment unified workflow (Seurat-style)
# - Advanced statistical methods (HMM, SLDS, TDA, Neural ODEs, etc.)
# - Publication-ready visualizations
#
# =============================================================================

library(CaImagingAnalysisFr)

# =============================================================================
# PART 1: THE CAEXPERIMENT WORKFLOW
# =============================================================================
#
# CaExperiment is the recommended way to use this package. It provides:
# - Unified data container for all analysis results
# - Automatic provenance tracking
# - Pipe-friendly interface
# - Intelligent subsetting
#
# =============================================================================

cat("\n")
cat("===========================================================================\n")
cat("  CaImagingAnalysisFr - Complete Package Showcase\n")
cat("===========================================================================\n\n")

# -----------------------------------------------------------------------------
# 1.1 Generate Synthetic Data
# -----------------------------------------------------------------------------

cat(">>> 1. GENERATING SYNTHETIC DATA\n")
cat("---------------------------------------------------------------------------\n")

set.seed(42)
raw_data <- generate_synthetic_data(

  n_cells = 50,
  n_time = 9000,    # 15 minutes at 10 Hz
  spike_prob = 0.02
)

# Extract traces matrix
cell_cols <- grep("^Cell_", names(raw_data))
traces <- t(as.matrix(raw_data[, cell_cols]))
frame_rate <- 10

cat(sprintf("   Cells: %d\n", nrow(traces)))
cat(sprintf("   Frames: %d (%.1f minutes at %d Hz)\n",
            ncol(traces), ncol(traces)/frame_rate/60, frame_rate))

# -----------------------------------------------------------------------------
# 1.2 Create CaExperiment and Run Full Pipeline
# -----------------------------------------------------------------------------

cat("\n>>> 2. CAEXPERIMENT UNIFIED WORKFLOW\n")
cat("---------------------------------------------------------------------------\n")
cat("   Creating CaExperiment and running analysis pipeline...\n\n")

# The recommended way: pipe-based workflow
# Note: RunCorrection is for neuropil subtraction (requires neuropil traces)
# For synthetic data, we skip it and go straight to dF/F
ca <- CaExperiment(traces, frame_rate = frame_rate) |>
  RunDFF(method = "rolling", window_size = 300) |>
  RunSpikes(method = "threshold", threshold_sd = 2.5) |>
  RunPCA(n_components = 15) |>
  RunConnectivity(method = "correlation") |>
  RunGraphMetrics() |>
  RunAssemblies(method = "pca", n_assemblies = 5) |>
  RunTransients(threshold_sd = 2.0)

cat("   Pipeline complete! Results stored in CaExperiment object.\n\n")

# Show what's in the object
print(ca)

# -----------------------------------------------------------------------------
# 1.3 Access Results
# -----------------------------------------------------------------------------

cat("\n>>> 3. ACCESSING RESULTS\n")
cat("---------------------------------------------------------------------------\n")

# Get processed traces
dff_traces <- GetTraces(ca, assay = "dff")
cat(sprintf("   dF/F traces: %d x %d\n", nrow(dff_traces), ncol(dff_traces)))

# Get spikes
spikes <- GetSpikes(ca)
total_spikes <- sum(spikes$spike_predictions)
cat(sprintf("   Total spikes detected: %d\n", total_spikes))

# Get PCA reduction
pca <- GetReduction(ca, "pca")
var_explained <- sum(pca$sdev[1:3]^2) / sum(pca$sdev^2) * 100
cat(sprintf("   PCA variance (3 PCs): %.1f%%\n", var_explained))

# Get connectivity graph
conn <- GetGraph(ca, "connectivity_correlation")
cat(sprintf("   Connectivity matrix: %d x %d\n", nrow(conn), ncol(conn)))

# Get assemblies
assemblies <- GetAssemblies(ca)
cat(sprintf("   Assemblies detected: %d\n", assemblies$n_assemblies))

# Get transients
transients <- GetTransients(ca)
cat(sprintf("   Transients detected: %d\n", nrow(transients$events)))

# -----------------------------------------------------------------------------
# 1.4 Command Logging / Provenance
# -----------------------------------------------------------------------------

cat("\n>>> 4. PROVENANCE TRACKING\n")
cat("---------------------------------------------------------------------------\n")

commands <- GetCommands(ca)
cat("   Analysis history:\n")
for (i in seq_along(commands)) {
  cat(sprintf("   [%d] %s\n", i, names(commands)[i]))
}

# Export reproducible script
# ExportCommands(ca, "my_analysis_pipeline.R")

# =============================================================================
# PART 2: ADVANCED STATISTICAL METHODS
# =============================================================================

cat("\n")
cat("===========================================================================\n")
cat("  PART 2: ADVANCED STATISTICAL METHODS\n")
cat("===========================================================================\n")

# -----------------------------------------------------------------------------
# 2.1 Hidden Markov Models (RunHMM)
# -----------------------------------------------------------------------------

cat("\n>>> 5. HIDDEN MARKOV MODELS\n")
cat("---------------------------------------------------------------------------\n")

ca <- RunHMM(ca, n_states = 3, n_iter = 50)
cat("   Fitted 3-state HMM to population activity\n")

hmm_result <- ca@misc$hmm
cat(sprintf("   Log-likelihood: %.2f\n", hmm_result$log_likelihood))
cat(sprintf("   Most common state: %d (%.1f%% of time)\n",
            which.max(table(hmm_result$states)),
            max(table(hmm_result$states)) / length(hmm_result$states) * 100))

# -----------------------------------------------------------------------------
# 2.2 Switching Linear Dynamical Systems (RunSLDS)
# -----------------------------------------------------------------------------

cat("\n>>> 6. SWITCHING LINEAR DYNAMICAL SYSTEMS\n")
cat("---------------------------------------------------------------------------\n")

ca <- RunSLDS(ca, n_states = 2, latent_dim = 3, n_iter = 30)
cat("   Fitted 2-state SLDS with 3D latent dynamics\n")

slds_result <- ca@misc$slds
cat(sprintf("   State 1 occupancy: %.1f%%\n",
            mean(slds_result$states == 1) * 100))
cat(sprintf("   Latent trajectory dims: %d x %d\n",
            nrow(slds_result$latent_states), ncol(slds_result$latent_states)))

# -----------------------------------------------------------------------------
# 2.3 Tensor Decomposition (RunTensorDecomp)
# -----------------------------------------------------------------------------

cat("\n>>> 7. TENSOR DECOMPOSITION\n")
cat("---------------------------------------------------------------------------\n")

# Create trial structure (split into 10 pseudo-trials)
ca <- RunTensorDecomp(ca, n_trials = 10, rank = 5, method = "cp")
cat("   Fitted rank-5 CP decomposition (neurons x time x trials)\n")

tensor_result <- ca@misc$tensor_decomp
cat(sprintf("   Neuron factors: %d x %d\n",
            nrow(tensor_result$neuron_factors), ncol(tensor_result$neuron_factors)))
cat(sprintf("   Temporal factors: %d x %d\n",
            nrow(tensor_result$temporal_factors), ncol(tensor_result$temporal_factors)))

# -----------------------------------------------------------------------------
# 2.4 Changepoint Detection (RunChangepoints)
# -----------------------------------------------------------------------------

cat("\n>>> 8. CHANGEPOINT DETECTION\n")
cat("---------------------------------------------------------------------------\n")

ca <- RunChangepoints(ca, method = "pelt", penalty = "BIC", max_changepoints = 10)
cat("   Detected changepoints in population dynamics\n")

cp_result <- ca@misc$changepoints
n_cp <- length(cp_result$changepoints)
cat(sprintf("   Changepoints found: %d\n", n_cp))
if (n_cp > 0) {
  cat(sprintf("   First changepoint at frame: %d (%.1f sec)\n",
              cp_result$changepoints[1], cp_result$changepoints[1] / frame_rate))
}

# -----------------------------------------------------------------------------
# 2.5 Topological Data Analysis (RunTDA)
# -----------------------------------------------------------------------------

cat("\n>>> 9. TOPOLOGICAL DATA ANALYSIS\n")
cat("---------------------------------------------------------------------------\n")

ca <- RunTDA(ca, max_dim = 2, n_points = 500)
cat("   Computed persistent homology of neural manifold\n")

tda_result <- ca@misc$tda
cat(sprintf("   H0 features (components): %d\n", tda_result$betti_numbers[1]))
cat(sprintf("   H1 features (loops): %d\n", tda_result$betti_numbers[2]))
if (length(tda_result$betti_numbers) > 2) {
  cat(sprintf("   H2 features (voids): %d\n", tda_result$betti_numbers[3]))
}

# -----------------------------------------------------------------------------
# 2.6 Recurrence Quantification Analysis (RunRQA)
# -----------------------------------------------------------------------------

cat("\n>>> 10. RECURRENCE QUANTIFICATION ANALYSIS\n")
cat("---------------------------------------------------------------------------\n")

ca <- RunRQA(ca, embedding_dim = 3, delay = 5, threshold = 0.1)
cat("   Computed recurrence plot and RQA metrics\n")

rqa_result <- ca@misc$rqa
cat(sprintf("   Recurrence rate: %.3f\n", rqa_result$recurrence_rate))
cat(sprintf("   Determinism: %.3f\n", rqa_result$determinism))
cat(sprintf("   Laminarity: %.3f\n", rqa_result$laminarity))

# -----------------------------------------------------------------------------
# 2.7 Causal Discovery (RunCausalDiscovery)
# -----------------------------------------------------------------------------

cat("\n>>> 11. CAUSAL DISCOVERY\n")
cat("---------------------------------------------------------------------------\n")

# Use subset for speed
ca <- RunCausalDiscovery(ca, method = "pc", alpha = 0.05, subset_cells = 1:15)
cat("   Inferred causal graph using PC algorithm\n")

causal_result <- ca@misc$causal_graph
n_edges <- sum(causal_result$adjacency_matrix) / 2
cat(sprintf("   Causal edges discovered: %d\n", n_edges))

# -----------------------------------------------------------------------------
# 2.8 Neural ODEs (RunNeuralODE)
# -----------------------------------------------------------------------------

cat("\n>>> 12. NEURAL ODEs\n")
cat("---------------------------------------------------------------------------\n")

ca <- RunNeuralODE(ca, latent_dim = 3, hidden_dim = 32, n_iter = 100)
cat("   Fitted Neural ODE for continuous-time latent dynamics\n")

node_result <- ca@misc$neural_ode
cat(sprintf("   Latent trajectories: %d timepoints x %d dims\n",
            nrow(node_result$latent_trajectory), ncol(node_result$latent_trajectory)))

# -----------------------------------------------------------------------------
# 2.9 Gaussian Process Smoothing (RunGPSmooth)
# -----------------------------------------------------------------------------

cat("\n>>> 13. GAUSSIAN PROCESS SMOOTHING\n")
cat("---------------------------------------------------------------------------\n")

ca <- RunGPSmooth(ca, length_scale = 10, subset_cells = 1:10)
cat("   Fitted GP regression for smooth interpolation\n")

gp_result <- ca@misc$gp_smooth
cat(sprintf("   Smoothed traces with uncertainty quantification\n"))
cat(sprintf("   Mean predictive std: %.4f\n", mean(gp_result$uncertainty)))

# -----------------------------------------------------------------------------
# 2.10 Optimal Transport Comparison (RunOTCompare)
# -----------------------------------------------------------------------------

cat("\n>>> 14. OPTIMAL TRANSPORT\n")
cat("---------------------------------------------------------------------------\n")

# Compare first and second half of recording
ca <- RunOTCompare(ca, split_point = 0.5, method = "sinkhorn")
cat("   Computed Wasserstein distance between recording halves\n")

ot_result <- ca@misc$ot_comparison
cat(sprintf("   Wasserstein distance: %.4f\n", ot_result$wasserstein_distance))

# =============================================================================
# PART 3: PUBLICATION-READY VISUALIZATIONS
# =============================================================================

cat("\n")
cat("===========================================================================\n")
cat("  PART 3: PUBLICATION-READY VISUALIZATIONS\n")
cat("===========================================================================\n")

if (interactive() && requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  cat("\n>>> 15. GENERATING PUBLICATION FIGURES\n")
  cat("---------------------------------------------------------------------------\n")

  # -------------------------------------------------------------------------
  # Neural Raster Plot
  # -------------------------------------------------------------------------

  cat("   [1/8] Neural raster plot...\n")
  p1 <- plot_neural_raster(
    GetSpikes(ca)$spike_predictions[1:20, 1:3000],
    frame_rate = frame_rate,
    title = "Population Activity"
  )

  # -------------------------------------------------------------------------
  # Connectivity Matrix
  # -------------------------------------------------------------------------

  cat("   [2/8] Connectivity matrix...\n")
  p2 <- plot_connectivity_matrix(
    GetGraph(ca, "connectivity_correlation")[1:25, 1:25],
    title = "Functional Connectivity"
  )

  # -------------------------------------------------------------------------
  # State Space Trajectory
  # -------------------------------------------------------------------------

  cat("   [3/8] GPFA trajectory...\n")
  pca_coords <- GetReduction(ca, "pca")$x[1:1000, 1:3]
  p3 <- plot_gpfa_trajectory(
    pca_coords,
    dims = c(1, 2, 3),
    color_by = "time",
    title = "Neural Trajectory"
  )

  # -------------------------------------------------------------------------
  # Assembly Timeline
  # -------------------------------------------------------------------------

  cat("   [4/8] Assembly timeline...\n")
  p4 <- plot_assembly_timeline(
    GetAssemblies(ca),
    frame_rate = frame_rate,
    title = "Assembly Activation"
  )

  # -------------------------------------------------------------------------
  # Circular Network
  # -------------------------------------------------------------------------

  cat("   [5/8] Circular network...\n")
  p5 <- plot_circular_network(
    GetGraph(ca, "connectivity_correlation")[1:20, 1:20],
    threshold = 0.3,
    title = "Network Structure"
  )

  # -------------------------------------------------------------------------
  # Changepoint Plot
  # -------------------------------------------------------------------------

  cat("   [6/8] Changepoint detection...\n")
  p6 <- plot_changepoints(
    ca@misc$changepoints,
    traces = colMeans(GetTraces(ca, "dff")),
    frame_rate = frame_rate,
    title = "Regime Changes"
  )

  # -------------------------------------------------------------------------
  # Recurrence Plot
  # -------------------------------------------------------------------------

  cat("   [7/8] Recurrence plot...\n")
  p7 <- plot_recurrence(
    ca@misc$rqa,
    title = "Recurrence Structure"
  )

  # -------------------------------------------------------------------------
  # Persistence Diagram
  # -------------------------------------------------------------------------

  cat("   [8/8] Persistence diagram...\n")
  p8 <- plot_persistence_diagram(
    ca@misc$tda,
    title = "Topological Features"
  )

  # -------------------------------------------------------------------------
  # Combine into multi-panel figure
  # -------------------------------------------------------------------------

  cat("\n   Combining into publication figure...\n")

  # Create 2x4 panel figure
  combined <- combine_panels(
    list(p1, p2, p3, p4, p5, p6, p7, p8),
    ncol = 4,
    labels = "AUTO"
  )

  # Save publication figure
  # save_pub_figure(combined, "calcium_analysis_figure.pdf",
  #                 width = 16, height = 8)

  print(combined)

  cat("\n   Publication figures generated!\n")
  cat("   Use save_pub_figure() to export at journal resolution.\n")

} else {
  cat("\n   [Visualization requires interactive session with ggplot2]\n")
}

# =============================================================================
# PART 4: MIXED-EFFECTS & CONFORMAL PREDICTION
# =============================================================================

cat("\n")
cat("===========================================================================\n")
cat("  PART 4: STATISTICAL INFERENCE\n")
cat("===========================================================================\n")

# -----------------------------------------------------------------------------
# 4.1 Mixed-Effects Models
# -----------------------------------------------------------------------------

cat("\n>>> 16. MIXED-EFFECTS MODELS\n")
cat("---------------------------------------------------------------------------\n")

# Prepare data for mixed model (simulating multiple subjects)
mixed_data <- prepare_mixed_data(
  traces = GetTraces(ca, "dff")[1:20, ],
  subject_ids = rep(1:4, each = 5),  # 4 subjects, 5 cells each
  condition = rep(c("A", "B"), each = 10)
)

# Fit linear mixed model
lmm_fit <- fit_lmm(
  data = mixed_data,
  formula = activity ~ condition + (1 | subject)
)
cat("   Fitted linear mixed model: activity ~ condition + (1|subject)\n")

# Compute ICC
icc <- compute_icc(lmm_fit)
cat(sprintf("   Intraclass correlation (ICC): %.3f\n", icc$icc))

# Test condition effect
condition_test <- test_condition_mixed(lmm_fit, "condition")
cat(sprintf("   Condition effect p-value: %.4f\n", condition_test$p_value))

# -----------------------------------------------------------------------------
# 4.2 Conformal Prediction
# -----------------------------------------------------------------------------

cat("\n>>> 17. CONFORMAL PREDICTION\n")
cat("---------------------------------------------------------------------------\n")

# Create conformal predictor for spike detection
train_traces <- GetTraces(ca, "dff")[, 1:6000]
test_traces <- GetTraces(ca, "dff")[, 6001:9000]

conformal <- conformal_predictor(
  train_data = train_traces,
  train_labels = GetSpikes(ca)$spike_predictions[, 1:6000],
  method = "split"
)
cat("   Created split conformal predictor\n")

# Get prediction intervals
predictions <- conformal_predict(
  conformal,
  newdata = test_traces,
  alpha = 0.1  # 90% coverage
)
cat(sprintf("   Prediction intervals at 90%% coverage\n"))
cat(sprintf("   Mean interval width: %.4f\n", mean(predictions$interval_width)))

# =============================================================================
# PART 5: PERFORMANCE FEATURES
# =============================================================================

cat("\n")
cat("===========================================================================\n")
cat("  PART 5: PERFORMANCE FEATURES\n")
cat("===========================================================================\n")

# -----------------------------------------------------------------------------
# 5.1 data.table Operations
# -----------------------------------------------------------------------------

cat("\n>>> 18. DATA.TABLE OPERATIONS\n")
cat("---------------------------------------------------------------------------\n")

if (requireNamespace("data.table", quietly = TRUE)) {
  # Convert to data.table for fast operations
  dt <- traces_to_dt(GetTraces(ca, "dff"))
  cat(sprintf("   Converted to data.table: %d rows\n", nrow(dt)))

  # Fast dF/F computation
  dt <- dt_compute_dff(dt, baseline_frames = 1:100)

  # Fast event detection
  dt <- dt_detect_events(dt, threshold = 2.5)

  # Cell summary statistics
  summary_dt <- dt_cell_summary(dt)
  cat(sprintf("   Cell summaries computed: %d cells\n", nrow(summary_dt)))
}

# -----------------------------------------------------------------------------
# 5.2 Parallel Processing
# -----------------------------------------------------------------------------

cat("\n>>> 19. PARALLEL PROCESSING\n")
cat("---------------------------------------------------------------------------\n")

cat("   Available parallel backends:\n")
cat("   - setup_parallel(workers = 4) for multi-core\n")
cat("   - parallel_spike_detection() for parallel spike inference\n")
cat("   - parallel_cross_correlation() for pairwise computations\n")

# -----------------------------------------------------------------------------
# 5.3 Arrow/DuckDB Backend
# -----------------------------------------------------------------------------

cat("\n>>> 20. OUT-OF-CORE PROCESSING\n")
cat("---------------------------------------------------------------------------\n")

cat("   For datasets exceeding RAM:\n")
cat("   - save_traces_parquet(traces, 'data.parquet')\n")
cat("   - query_traces('SELECT * FROM traces WHERE cell_id < 100')\n")
cat("   - stream_process_traces(parquet_file, chunk_size = 1000)\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("===========================================================================\n")
cat("  SHOWCASE COMPLETE\n")
cat("===========================================================================\n\n")

cat("CaExperiment object contains:\n")
cat(sprintf("  - %d cells, %d frames\n", ncells(ca), nframes(ca)))
cat(sprintf("  - %d assays (raw, corrected, dff)\n", length(ca@assays)))
cat(sprintf("  - %d reductions (pca)\n", length(ca@reductions)))
cat(sprintf("  - %d graphs (connectivity)\n", length(ca@graphs)))
cat(sprintf("  - %d commands logged\n", length(GetCommands(ca))))

cat("\nAdvanced analyses performed:\n")
cat("  - HMM (3 states)\n")
cat("  - SLDS (2 states, 3D latent)\n")
cat("  - Tensor decomposition (rank 5)\n")
cat("  - Changepoint detection (PELT)\n")
cat("  - TDA (persistent homology)\n")
cat("  - RQA (recurrence analysis)\n")
cat("  - Causal discovery (PC algorithm)\n")
cat("  - Neural ODE (continuous dynamics)\n")
cat("  - GP smoothing (uncertainty)\n")
cat("  - Optimal transport (distribution comparison)\n")

cat("\nFor more information:\n")
cat("  - help(CaExperiment)\n")
cat("  - help(RunHMM), help(RunSLDS), etc.\n")
cat("  - https://kevinj24fr.github.io/CaImagingAnalysisFr/\n\n")
