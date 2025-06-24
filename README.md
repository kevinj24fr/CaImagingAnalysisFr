# CaImagingAnalysisFr <img src="https://img.shields.io/badge/R-Professional-blue" align="right" height="30"/>

[![CRAN Status](https://www.r-pkg.org/badges/version/CaImagingAnalysisFr)](https://cran.r-project.org/package=CaImagingAnalysisFr)
[![Build Status](https://github.com/yourusername/CaImagingAnalysisFr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yourusername/CaImagingAnalysisFr/actions)
[![Coverage Status](https://codecov.io/gh/yourusername/CaImagingAnalysisFr/branch/main/graph/badge.svg)](https://codecov.io/gh/yourusername/CaImagingAnalysisFr)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/CaImagingAnalysisFr)](https://cran.r-project.org/package=CaImagingAnalysisFr)

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

**From GitHub (latest):**
```r
# install.packages("devtools")
devtools::install_github("yourusername/CaImagingAnalysisFr")
```

---

## ⚡ Quick Start

```r
library(CaImagingAnalysisFr)
set.seed(123)

# Simulate data and run a basic workflow
raw_img <- array(rnorm(100*100*100), dim = c(100, 100, 100))
seg <- segment_cells(raw_img, method = "suite2p")
trace <- rnorm(500)
spikes <- infer_spikes(trace, method = "oasis")

# Visualize segmentation overlay and save as PNG
props <- region_properties(seg$rois, apply(raw_img, c(1,2), mean))
png("docs/segmentation_overlay_example.png", width = 600, height = 600)
plot_segmentation_overlay(apply(raw_img, c(1,2), mean), seg$rois, property = "area", props = props)
dev.off()

# Spike inference plot
png("docs/spike_inference_example.png", width = 600, height = 400)
plot(trace, type = "l", col = "gray", main = "Calcium Trace with Inferred Spikes", ylab = "Fluorescence")
points(which(spikes$spike > 0), trace[spikes$spike > 0], col = "red", pch = 19)
dev.off()

# Batch correction PCA plot
mat <- matrix(rnorm(1000), nrow = 50)
batch <- rep(1:2, each = 25)
corrected <- batch_correction(mat, batch = batch, method = "deep")
pca <- prcomp(t(cbind(mat, corrected)), scale. = TRUE)
batch_labels <- c(rep("Raw", ncol(mat)), rep("Corrected", ncol(corrected)))
png("docs/batch_correction_pca.png", width = 600, height = 400)
plot(pca$x[,1:2], col = ifelse(batch_labels == "Raw", "orange", "blue"), pch = 19,
     main = "Batch Correction: PCA Before/After", xlab = "PC1", ylab = "PC2")
legend("topright", legend = c("Raw", "Corrected"), col = c("orange", "blue"), pch = 19)
dev.off()

# Denoising example (NMF)
nmf_res <- nmf_decompose(abs(mat), n_components = 2)
png("docs/denoising_example.png", width = 600, height = 400)
plot(mat[,1], type = "l", col = "gray", main = "Raw vs. NMF Denoised Trace", ylab = "Signal")
lines(nmf_res$W[,1], col = "blue", lwd = 2)
legend("topright", legend = c("Raw", "NMF Denoised"), col = c("gray", "blue"), lwd = c(1,2))
dev.off()

# Network analysis plot
fc <- functional_connectivity(mat, method = "correlation", threshold = 0.3)
png("docs/network_analysis_example.png", width = 600, height = 600)
net <- network_visualization(fc$connectivity_matrix)
dev.off()

# Bayesian posterior plot
bayes <- bayesian_spike_inference(trace)
png("docs/bayesian_posterior_example.png", width = 600, height = 400)
hist(bayes$posterior_spikes, breaks = 30, col = "skyblue", main = "Posterior Spike Probability", xlab = "Probability")
dev.off()

# Benchmark segmentation quality
bench <- benchmark_segmentation_quality(props)
print(bench)
```

![Segmentation Overlay Example](docs/segmentation_overlay_example.png)
*Segmentation overlay: ROIs color-coded by area.*

![Spike Inference Example](docs/spike_inference_example.png)
*Calcium trace with inferred spikes (red dots).* 

![Batch Correction PCA](docs/batch_correction_pca.png)
*PCA plot: orange = raw, blue = batch-corrected.*

![Denoising Example](docs/denoising_example.png)
*Raw trace (gray) vs. NMF denoised (blue).* 

![Network Analysis Example](docs/network_analysis_example.png)
*Functional connectivity network plot.*

![Bayesian Posterior Example](docs/bayesian_posterior_example.png)
*Posterior spike probability distribution.*

*You can regenerate all output images above by running the code block in your R session.*

---

## 📚 Comprehensive Usage

```r
# Batch correction
mat <- matrix(rnorm(1000), nrow = 50)
batch <- rep(1:2, each = 25)
corrected <- batch_correction(mat, batch = batch, method = "deep")

# Denoising
nmf_res <- nmf_decompose(abs(mat), n_components = 2)

# Network analysis
fc <- functional_connectivity(mat, method = "correlation", threshold = 0.3)

# Bayesian modeling
bayes <- bayesian_spike_inference(trace)
```

---

## 📖 Vignettes & Tutorials
- [Getting Started](vignettes/getting-started.Rmd): Basic and intermediate workflow
- [Power Features](vignettes/power-features.Rmd): Advanced analytics and unique capabilities
- [Basic Tutorial](vignettes/tutorial-basic-workflow.Rmd): Step-by-step guide for new users

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
- [Function Reference (R Documentation)](https://yourusername.github.io/CaImagingAnalysisFr/reference/)

---

## 🤝 Contributing
We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines, code style, and how to run tests. All contributors and feedback are appreciated.

---

## 📖 Citing
If you use CaImagingAnalysisFr in your research, please cite:

> Your Name et al. (2024). CaImagingAnalysisFr: Professional Calcium Imaging Analysis in R. _GitHub_. https://github.com/yourusername/CaImagingAnalysisFr

---

## 📝 License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## 🙏 Acknowledgments
- Inspired by Suite2p, CaImAn, Cellpose, and the open-source neuroscience community.
- Thanks to all contributors and users!

---

## 📬 Contact
For questions, issues, or support, please open a [GitHub Issue](https://github.com/yourusername/CaImagingAnalysisFr/issues) or contact [your.email@domain.com](mailto:your.email@domain.com).
