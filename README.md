# CaImagingAnalysisFr

> **A professional, robust, and state-of-the-art R package for advanced calcium imaging analysis—no Python required.**

---

## 🚀 Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Comprehensive Usage](#comprehensive-usage)
- [Vignettes & Tutorials](#vignettes--tutorials)
- [Segmentation & Quality Guidance](#segmentation--quality-guidance)
- [API Reference](#api-reference)
- [Contributing](#contributing)
- [Citing](#citing)
- [License](#license)
- [Acknowledgments](#acknowledgments)
- [Contact](#contact)

---

## ✨ Features
- 🧠 **Automated Cell Segmentation** (Suite2p, Cellpose, CaImAn, k-means, threshold) — all in base R
- 🔬 **Deep Learning Spike Inference** (base R, no Python required)
- 🧬 **Batch Effect Correction** (ComBat, MNN, deep harmonization)
- 🧹 **Advanced Denoising** (NMF, ICA, wavelet)
- 🕸️ **Dynamic Network & Causality Analysis**
- 🤖 **Unsupervised Learning & Anomaly Detection**
- 📊 **Bayesian Modeling & Model Comparison**
- 📋 **Automated Reporting & Interactive QC**
- 🗂️ **Data Curation & Metadata Handling**
- 🛠️ **End-to-End Professional Workflow**

---

## 🛠️ Installation

**From CRAN:**
```r
install.packages("CaImagingAnalysisFr")
```

---

## ⚡ Quick Start

```r
library(CaImagingAnalysisFr)
set.seed(123)

# Generate synthetic calcium imaging data
synthetic_data <- generate_synthetic_data(
  n_cells = 10,
  n_time = 500,
  spike_prob = 0.1
)

# Extract traces
traces <- synthetic_data[, grep("^Cell_", names(synthetic_data))]

# Basic calcium correction
corrected_traces <- calcium_correction(traces, method = "modern")

# Spike inference
spike_results <- infer_spikes(corrected_traces, method = "oasis")

# Plot results
plot_cell_trace(
  corrected_df = data.frame(Time = 1:ncol(traces), Cell_1 = traces[1, ]),
  cell = "Cell_1"
)

# Batch correction example
batch_labels <- rep(1:2, each = 5)
corrected_batch <- batch_correction(
  traces, 
  batch = batch_labels, 
  method = "combat"
)

# Network analysis
network_result <- functional_connectivity(
  corrected_traces, 
  method = "correlation", 
  threshold = 0.3
)

# Bayesian modeling
bayesian_result <- bayesian_spike_inference(
  traces[1, ], 
  model_type = "poisson", 
  n_samples = 100
)

# Generate quality control report
qc_report <- generate_qc_report(data = traces)

print("Quick start completed successfully!")
```

---

## 📚 Comprehensive Usage

```r
# Advanced segmentation
segmentation_result <- segment_cells(
  image_data = array(rnorm(100*100*50), dim = c(100, 100, 50)),
  method = "threshold"
)

# Deep learning spike inference
deep_spikes <- deep_spike_inference(
  traces[1, ], 
  model_type = "lstm"
)

# Denoising with NMF
nmf_result <- nmf_decompose(
  abs(traces), 
  n_components = 3
)

# Dynamic network analysis
dynamic_network <- time_varying_connectivity(
  traces, 
  window_size = 50, 
  step_size = 10
)

# Bayesian model comparison
model_comparison <- bayesian_model_comparison(
  list(model1 = bayesian_result, model2 = deep_spikes)
)
```

---

## 📖 Vignettes & Tutorials
- Getting Started: Basic and intermediate workflow
- Power Features: Advanced analytics and unique capabilities
- Basic Tutorial: Step-by-step guide for new users

---

## 🧩 Segmentation & Quality Guidance
- **Number of ROIs (`n_roi`)**: Should match expected cell count. Too high = over-segmentation; too low = missed cells.
- **Area**: Should match expected cell size. Large variance = inconsistent segmentation.
- **Eccentricity**: 0 (circle) to 1 (line). Most cells: 0.3–0.8. High = elongated/merged.
- **Solidity**: Area/convex hull area. Near 1 = compact/convex (good). Low = irregular/fragmented.
- **Jaccard**: Overlap with ground truth. Near 1 = ideal; <0.5 = poor.

Use these metrics to tune segmentation and assess reliability.

---

## 📚 API Reference
See the package documentation for complete function reference.

---

## 🤝 Contributing
We welcome contributions! Please see the contributing guidelines for code style and how to run tests. All contributors and feedback are appreciated.

---

## 📖 Citing
If you use CaImagingAnalysisFr in your research, please cite:

> Calcium Team (2024). CaImagingAnalysisFr: Professional Calcium Imaging Analysis in R.

---

## 📝 License
This project is licensed under the MIT License. See LICENSE file for details.

---

## 🙏 Acknowledgments
- Inspired by Suite2p, CaImAn, Cellpose, and the open-source neuroscience community.
- Thanks to all contributors and users!

---

## 📬 Contact
For questions, issues, or support, please contact the maintainer.
