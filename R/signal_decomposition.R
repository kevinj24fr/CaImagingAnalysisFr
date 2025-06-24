#' Non-negative Matrix Factorization (NMF) for Calcium Imaging
#'
#' Decompose a matrix of traces into non-negative components using NMF.
#'
#' @param data Matrix or data frame (cells x time)
#' @param n_components Number of components to extract
#' @param ... Additional arguments to NMF
#' @return List with W (basis) and H (coefficients)
#' @export
nmf_decompose <- function(data, n_components = 2, ...) {
  if (!requireNamespace("NMF", quietly = TRUE)) {
    stop("Package 'NMF' is required for NMF decomposition. Install with install.packages('NMF')")
  }
  nmf_result <- NMF::nmf(as.matrix(data), rank = n_components, ...)
  list(W = NMF::basis(nmf_result), H = NMF::coef(nmf_result))
}

#' Independent Component Analysis (ICA) for Calcium Imaging
#'
#' Decompose a matrix of traces into statistically independent components using ICA.
#'
#' @param data Matrix or data frame (cells x time)
#' @param n_components Number of components to extract
#' @return List with S (sources) and M (mixing matrix)
#' @export
ica_decompose <- function(data, n_components = 2) {
  if (!requireNamespace("fastICA", quietly = TRUE)) {
    stop("Package 'fastICA' is required for ICA decomposition. Install with install.packages('fastICA')")
  }
  ica_result <- fastICA::fastICA(as.matrix(data), n.comp = n_components)
  list(S = ica_result$S, M = ica_result$A)
}

#' Wavelet Denoising for Calcium Imaging Traces
#'
#' Denoise a numeric vector using wavelet shrinkage.
#'
#' @param trace Numeric vector
#' @param family Wavelet family (default: "haar")
#' @param verbose Show progress messages
#' @return Denoised vector
#' @export
wavelet_denoise <- function(trace, family = "haar", verbose = FALSE) {
  if (!requireNamespace("waveslim", quietly = TRUE)) {
    stop("Package 'waveslim' is required for wavelet denoising. Install with install.packages('waveslim')")
  }
  if (verbose) message("Performing wavelet denoising...")
  wd <- waveslim::dwt(trace, family)
  # Universal threshold
  sigma <- mad(wd$d1)
  threshold <- sigma * sqrt(2 * log(length(trace)))
  wd_thr <- lapply(wd, function(x) ifelse(abs(x) > threshold, x, 0))
  denoised <- waveslim::idwt(wd_thr, family)
  denoised[is.na(denoised)] <- trace[is.na(denoised)]
  denoised
}

#' Robust Principal Component Analysis (RPCA) via Python (reticulate)
#'
#' Decompose a matrix into low-rank and sparse components using RPCA (via Python's `r_pca` or `scikit-learn`)
#'
#' @param data Matrix or data frame (cells x time)
#' @param install_missing Whether to install missing Python dependencies
#' @param verbose Show progress messages
#' @return List with L (low-rank) and S (sparse)
#' @export
rpca_decompose <- function(data, install_missing = TRUE, verbose = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for RPCA via Python.")
  }
  # Check Python dependencies
  deps <- manage_python_dependencies(packages = c("r_pca"), install_missing = install_missing, verbose = verbose)
  if (!deps$r_pca$available) {
    stop("Python package 'r_pca' is not available and could not be installed.")
  }
  if (verbose) message("Running RPCA via Python...")
  r_pca <- reticulate::import("r_pca")
  rpca <- r_pca$R_pca(as.matrix(data))
  rpca$fit()
  list(L = rpca$L, S = rpca$S)
} 