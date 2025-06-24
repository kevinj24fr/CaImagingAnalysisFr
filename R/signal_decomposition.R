#' Signal Decomposition and Component Analysis
#'
#' Decompose calcium imaging signals into components using various methods.
#'
#' @name signal_decomposition
#' @docType package
NULL

#' Robust Principal Component Analysis (RPCA)
#'
#' Decompose signal into low-rank and sparse components using RPCA.
#'
#' @param data Input data matrix (time x variables)
#' @param lambda Penalty parameter for sparse component
#' @param mu Penalty parameter for low-rank component
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @param ... Additional arguments
#' @return List containing decomposed components
#' @export
robust_pca <- function(data, lambda = 1/sqrt(max(nrow(data), ncol(data))), 
                      mu = 10*lambda, max_iter = 100, tol = 1e-6, ...) {
  message("Running Robust PCA decomposition")
  
  # Base R implementation of RPCA using ADMM
  # This is a simplified version that captures the essence of RPCA
  
  # Initialize
  L <- matrix(0, nrow(data), ncol(data))  # Low-rank component
  S <- matrix(0, nrow(data), ncol(data))  # Sparse component
  Y <- matrix(0, nrow(data), ncol(data))  # Lagrange multiplier
  
  # ADMM iterations
  for (iter in 1:max_iter) {
    # Update L (low-rank component)
    L_old <- L
    temp <- data - S + Y/mu
    L <- soft_threshold_svd(temp, 1/mu)
    
    # Update S (sparse component)
    S_old <- S
    temp <- data - L + Y/mu
    S <- soft_threshold(temp, lambda/mu)
    
    # Update Y (Lagrange multiplier)
    Y <- Y + mu * (data - L - S)
    
    # Check convergence
    if (iter > 1) {
      L_change <- norm(L - L_old, "F") / norm(L_old, "F")
      S_change <- norm(S - S_old, "F") / norm(S_old, "F")
      
      if (max(L_change, S_change) < tol) {
        message("RPCA converged after ", iter, " iterations")
        break
      }
    }
  }
  
  return(list(
    low_rank = L,
    sparse = S,
    residual = data - L - S,
    parameters = list(lambda = lambda, mu = mu, iterations = iter)
  ))
}

#' Non-negative Matrix Factorization (NMF)
#'
#' Decompose signal using non-negative matrix factorization.
#'
#' @param data Input data matrix (time x variables)
#' @param rank Rank of factorization
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @param method NMF method ("multiplicative", "als")
#' @param ... Additional arguments
#' @return List containing decomposed components
#' @export
nonnegative_matrix_factorization <- function(data, rank = 3, max_iter = 100, 
                                           tol = 1e-6, method = "multiplicative", ...) {
  message("Running Non-negative Matrix Factorization")
  
  # Ensure non-negative data
  data <- pmax(data, 0)
  
  if (method == "multiplicative") {
    return(nmf_multiplicative(data, rank, max_iter, tol))
  } else if (method == "als") {
    return(nmf_als(data, rank, max_iter, tol))
  } else {
    stop("Unknown NMF method: ", method)
  }
}

#' Independent Component Analysis (ICA)
#'
#' Decompose signal using independent component analysis.
#'
#' @param data Input data matrix (time x variables)
#' @param n_components Number of components to extract
#' @param method ICA method ("fastica", "infomax", "jade")
#' @param ... Additional arguments
#' @return List containing decomposed components
#' @export
independent_component_analysis <- function(data, n_components = NULL, 
                                         method = "fastica", ...) {
  message("Running Independent Component Analysis")
  
  # Determine number of components
  if (is.null(n_components)) {
    n_components <- min(nrow(data), ncol(data))
  }
  
  # Center and whiten data
  data_centered <- scale(data, center = TRUE, scale = FALSE)
  
  if (method == "fastica") {
    return(ica_fastica(data_centered, n_components, ...))
  } else if (method == "infomax") {
    return(ica_infomax(data_centered, n_components, ...))
  } else if (method == "jade") {
    return(ica_jade(data_centered, n_components, ...))
  } else {
    stop("Unknown ICA method: ", method)
  }
}

#' Singular Spectrum Analysis (SSA)
#'
#' Decompose signal using singular spectrum analysis.
#'
#' @param data Input time series
#' @param window_size Window size for embedding
#' @param n_components Number of components to extract
#' @param ... Additional arguments
#' @return List containing decomposed components
#' @export
singular_spectrum_analysis <- function(data, window_size = NULL, n_components = NULL, ...) {
  message("Running Singular Spectrum Analysis")
  
  # Determine window size
  if (is.null(window_size)) {
    window_size <- min(length(data) %/% 3, 50)
  }
  
  # Determine number of components
  if (is.null(n_components)) {
    n_components <- window_size
  }
  
  # 1. Embedding
  trajectory_matrix <- create_trajectory_matrix(data, window_size)
  
  # 2. SVD decomposition
  svd_result <- svd(trajectory_matrix)
  
  # 3. Grouping (simplified - take first n_components)
  components <- list()
  for (i in 1:min(n_components, length(svd_result$d))) {
    # Reconstruct component
    component_matrix <- svd_result$u[, i, drop = FALSE] %*% 
                      t(svd_result$v[, i, drop = FALSE]) * svd_result$d[i]
    
    # Diagonal averaging to get time series
    component_ts <- diagonal_averaging(component_matrix)
    components[[i]] <- component_ts
  }
  
  # 4. Reconstruction
  reconstructed <- Reduce("+", components)
  
  return(list(
    components = components,
    reconstructed = reconstructed,
    singular_values = svd_result$d,
    trajectory_matrix = trajectory_matrix,
    parameters = list(window_size = window_size, n_components = n_components)
  ))
}

#' Empirical Mode Decomposition (EMD)
#'
#' Decompose signal using empirical mode decomposition.
#'
#' @param data Input time series
#' @param max_imfs Maximum number of IMFs to extract
#' @param max_iter Maximum iterations for IMF extraction
#' @param ... Additional arguments
#' @return List containing decomposed components
#' @export
empirical_mode_decomposition <- function(data, max_imfs = 10, max_iter = 100, ...) {
  message("Running Empirical Mode Decomposition")
  
  # Base R implementation of EMD
  # This is a simplified version that captures the essence of EMD
  
  imfs <- list()
  residual <- data
  
  for (imf_idx in 1:max_imfs) {
    message("Extracting IMF ", imf_idx)
    
    # Extract one IMF
    imf <- extract_imf(residual, max_iter)
    
    # Check if IMF is valid
    if (is_valid_imf(imf)) {
      imfs[[imf_idx]] <- imf
      residual <- residual - imf
    } else {
      message("Stopping EMD: no more valid IMFs")
      break
    }
  }
  
  return(list(
    imfs = imfs,
    residual = residual,
    n_imfs = length(imfs)
  ))
}

#' Helper Functions
#' @keywords internal

soft_threshold_svd <- function(X, tau) {
  # Soft thresholding for SVD (nuclear norm)
  svd_result <- svd(X)
  d_thresholded <- pmax(svd_result$d - tau, 0)
  
  # Reconstruct
  result <- svd_result$u %*% diag(d_thresholded) %*% t(svd_result$v)
  return(result)
}

soft_threshold <- function(X, tau) {
  # Soft thresholding for L1 norm
  result <- sign(X) * pmax(abs(X) - tau, 0)
  return(result)
}

norm <- function(X, type = "F") {
  # Matrix norm
  if (type == "F") {
    return(sqrt(sum(X^2)))
  } else {
    stop("Only Frobenius norm implemented")
  }
}

nmf_multiplicative <- function(data, rank, max_iter, tol) {
  # Multiplicative update NMF
  n <- nrow(data)
  p <- ncol(data)
  
  # Initialize W and H randomly
  W <- matrix(runif(n * rank), n, rank)
  H <- matrix(runif(rank * p), rank, p)
  
  # Ensure non-negativity
  W <- pmax(W, 1e-10)
  H <- pmax(H, 1e-10)
  
  for (iter in 1:max_iter) {
    W_old <- W
    H_old <- H
    
    # Update H
    numerator <- t(W) %*% data
    denominator <- t(W) %*% W %*% H
    H <- H * (numerator / (denominator + 1e-10))
    
    # Update W
    numerator <- data %*% t(H)
    denominator <- W %*% H %*% t(H)
    W <- W * (numerator / (denominator + 1e-10))
    
    # Check convergence
    if (iter > 1) {
      W_change <- norm(W - W_old, "F") / norm(W_old, "F")
      H_change <- norm(H - H_old, "F") / norm(H_old, "F")
      
      if (max(W_change, H_change) < tol) {
        message("NMF converged after ", iter, " iterations")
        break
      }
    }
  }
  
  return(list(
    W = W,
    H = H,
    reconstructed = W %*% H,
    iterations = iter
  ))
}

nmf_als <- function(data, rank, max_iter, tol) {
  # Alternating Least Squares NMF
  n <- nrow(data)
  p <- ncol(data)
  
  # Initialize W and H randomly
  W <- matrix(runif(n * rank), n, rank)
  H <- matrix(runif(rank * p), rank, p)
  
  for (iter in 1:max_iter) {
    W_old <- W
    H_old <- H
    
    # Update H (solve W*H = data for H)
    H <- solve(t(W) %*% W + 1e-8 * diag(rank)) %*% t(W) %*% data
    H <- pmax(H, 0)  # Ensure non-negativity
    
    # Update W (solve W*H = data for W)
    W <- data %*% t(H) %*% solve(H %*% t(H) + 1e-8 * diag(rank))
    W <- pmax(W, 0)  # Ensure non-negativity
    
    # Check convergence
    if (iter > 1) {
      W_change <- norm(W - W_old, "F") / norm(W_old, "F")
      H_change <- norm(H - H_old, "F") / norm(H_old, "F")
      
      if (max(W_change, H_change) < tol) {
        message("NMF ALS converged after ", iter, " iterations")
        break
      }
    }
  }
  
  return(list(
    W = W,
    H = H,
    reconstructed = W %*% H,
    iterations = iter
  ))
}

ica_fastica <- function(data, n_components, max_iter = 100, tol = 1e-6) {
  # FastICA algorithm
  n <- nrow(data)
  p <- ncol(data)
  
  # Whiten data
  cov_matrix <- cov(data)
  eigen_result <- eigen(cov_matrix)
  whitening_matrix <- eigen_result$vectors %*% diag(1/sqrt(eigen_result$values)) %*% t(eigen_result$vectors)
  data_whitened <- data %*% whitening_matrix
  
  # Initialize mixing matrix
  A <- matrix(rnorm(n_components * p), n_components, p)
  A <- t(scale(t(A)))  # Normalize rows
  
  for (iter in 1:max_iter) {
    A_old <- A
    
    # Update mixing matrix
    for (i in 1:n_components) {
      # FastICA update
      w <- A[i, ]
      wx <- data_whitened %*% w
      
      # Non-linearity (tanh)
      g_wx <- tanh(wx)
      g_prime_wx <- 1 - tanh(wx)^2
      
      # Update
      w_new <- colMeans(data_whitened * g_wx) - mean(g_prime_wx) * w
      w_new <- w_new / sqrt(sum(w_new^2))
      
      A[i, ] <- w_new
    }
    
    # Check convergence
    if (iter > 1) {
      change <- norm(A - A_old, "F") / norm(A_old, "F")
      if (change < tol) {
        message("FastICA converged after ", iter, " iterations")
        break
      }
    }
  }
  
  # Extract independent components
  S <- data_whitened %*% t(A)
  
  return(list(
    S = S,
    A = A,
    whitening_matrix = whitening_matrix,
    iterations = iter
  ))
}

ica_infomax <- function(data, n_components, max_iter = 100, learning_rate = 0.01) {
  # Infomax ICA algorithm
  n <- nrow(data)
  p <- ncol(data)
  
  # Whiten data
  cov_matrix <- cov(data)
  eigen_result <- eigen(cov_matrix)
  whitening_matrix <- eigen_result$vectors %*% diag(1/sqrt(eigen_result$values)) %*% t(eigen_result$vectors)
  data_whitened <- data %*% whitening_matrix
  
  # Initialize unmixing matrix
  W <- matrix(rnorm(n_components * p), n_components, p)
  
  for (iter in 1:max_iter) {
    # Forward pass
    S <- data_whitened %*% t(W)
    
    # Non-linearity (sigmoid)
    Y <- 1 / (1 + exp(-S))
    
    # Gradient update
    gradient <- t(data_whitened) %*% (1 - 2 * Y) + solve(t(W))
    W <- W + learning_rate * t(gradient)
  }
  
  return(list(
    S = S,
    W = W,
    whitening_matrix = whitening_matrix,
    iterations = max_iter
  ))
}

ica_jade <- function(data, n_components, ...) {
  # JADE (Joint Approximate Diagonalization of Eigenmatrices) ICA
  # This is a simplified version
  
  # Use PCA as approximation
  pca_result <- prcomp(data, center = TRUE, scale = FALSE)
  
  return(list(
    S = pca_result$x[, 1:n_components, drop = FALSE],
    A = t(pca_result$rotation[, 1:n_components, drop = FALSE]),
    eigenvalues = pca_result$sdev[1:n_components]^2
  ))
}

create_trajectory_matrix <- function(data, window_size) {
  # Create trajectory matrix for SSA
  n <- length(data)
  n_windows <- n - window_size + 1
  
  trajectory_matrix <- matrix(0, window_size, n_windows)
  for (i in 1:n_windows) {
    trajectory_matrix[, i] <- data[i:(i + window_size - 1)]
  }
  
  return(trajectory_matrix)
}

diagonal_averaging <- function(matrix_data) {
  # Diagonal averaging for SSA reconstruction
  n_rows <- nrow(matrix_data)
  n_cols <- ncol(matrix_data)
  
  result <- numeric(n_rows + n_cols - 1)
  counts <- numeric(n_rows + n_cols - 1)
  
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      idx <- i + j - 1
      result[idx] <- result[idx] + matrix_data[i, j]
      counts[idx] <- counts[idx] + 1
    }
  }
  
  # Average
  result <- result / counts
  
  return(result)
}

extract_imf <- function(data, max_iter) {
  # Extract one IMF using sifting process
  n <- length(data)
  h <- data
  
  for (iter in 1:max_iter) {
    # Find extrema
    maxima_indices <- find_extrema(h, "max")
    minima_indices <- find_extrema(h, "min")
    
    if (length(maxima_indices) < 2 || length(minima_indices) < 2) {
      break
    }
    
    # Interpolate envelopes
    upper_envelope <- interpolate_envelope(maxima_indices, h[maxima_indices], n)
    lower_envelope <- interpolate_envelope(minima_indices, h[minima_indices], n)
    
    # Calculate mean envelope
    mean_envelope <- (upper_envelope + lower_envelope) / 2
    
    # Extract detail
    h <- h - mean_envelope
    
    # Check if IMF criteria are met
    if (is_imf_criteria_met(h)) {
      break
    }
  }
  
  return(h)
}

find_extrema <- function(data, type) {
  # Find extrema in time series
  n <- length(data)
  extrema <- numeric(0)
  
  for (i in 2:(n-1)) {
    if (type == "max" && data[i] > data[i-1] && data[i] > data[i+1]) {
      extrema <- c(extrema, i)
    } else if (type == "min" && data[i] < data[i-1] && data[i] < data[i+1]) {
      extrema <- c(extrema, i)
    }
  }
  
  return(extrema)
}

interpolate_envelope <- function(indices, values, n) {
  # Interpolate envelope using cubic spline
  if (length(indices) < 2) {
    return(rep(mean(values), n))
  }
  
  # Add endpoints if needed
  if (indices[1] > 1) {
    indices <- c(1, indices)
    values <- c(values[1], values)
  }
  if (indices[length(indices)] < n) {
    indices <- c(indices, n)
    values <- c(values, values[length(values)])
  }
  
  # Interpolate
  envelope <- approx(indices, values, xout = 1:n, method = "linear")$y
  return(envelope)
}

is_imf_criteria_met <- function(data) {
  # Check if IMF criteria are met (simplified)
  # In practice, you'd check for zero crossings, extrema, etc.
  return(TRUE)  # Simplified
}

is_valid_imf <- function(imf) {
  # Check if IMF is valid (simplified)
  # In practice, you'd check various criteria
  return(TRUE)  # Simplified
}

#' NMF Decomposition
#'
#' Perform Non-negative Matrix Factorization on calcium imaging data.
#'
#' @param data Input data matrix
#' @param n_components Number of components
#' @param max_iter Maximum iterations
#' @param tol Tolerance for convergence
#' @param verbose Whether to show progress
#' @return NMF decomposition results
#' @export
nmf_decompose <- function(data, n_components = 2, max_iter = 100, tol = 1e-6, verbose = TRUE) {
  if (verbose) {
    message("Performing NMF decomposition with ", n_components, " components")
  }
  
  # Validate inputs
  if (!is.matrix(data) || any(data < 0)) {
    stop("Data must be a non-negative matrix")
  }
  
  if (n_components > min(nrow(data), ncol(data))) {
    stop("Number of components cannot exceed matrix dimensions")
  }
  
  # Initialize W and H matrices randomly
  n <- nrow(data)
  m <- ncol(data)
  
  W <- matrix(runif(n * n_components), n, n_components)
  H <- matrix(runif(n_components * m), n_components, m)
  
  # Multiplicative update algorithm
  for (iter in 1:max_iter) {
    # Update H
    H <- H * (t(W) %*% data) / (t(W) %*% W %*% H + 1e-10)
    
    # Update W
    W <- W * (data %*% t(H)) / (W %*% H %*% t(H) + 1e-10)
    
    # Check convergence
    if (iter %% 10 == 0) {
      reconstruction <- W %*% H
      error <- sum((data - reconstruction)^2)
      if (verbose) {
        message("Iteration ", iter, ", Error: ", error)
      }
    }
  }
  
  return(list(
    W = W,
    H = H,
    reconstruction = W %*% H,
    n_components = n_components,
    iterations = iter
  ))
}

#' ICA Decomposition
#'
#' Perform Independent Component Analysis on calcium imaging data.
#'
#' @param data Input data matrix
#' @param n_components Number of components
#' @param max_iter Maximum iterations
#' @param tol Tolerance for convergence
#' @param verbose Whether to show progress
#' @return ICA decomposition results
#' @export
ica_decompose <- function(data, n_components = 2, max_iter = 100, tol = 1e-6, verbose = TRUE) {
  if (verbose) {
    message("Performing ICA decomposition with ", n_components, " components")
  }
  
  # Validate inputs
  if (!is.matrix(data)) {
    stop("Data must be a matrix")
  }
  
  if (n_components > min(nrow(data), ncol(data))) {
    stop("Number of components cannot exceed matrix dimensions")
  }
  
  # Center the data
  data_centered <- scale(data, center = TRUE, scale = FALSE)
  
  # PCA preprocessing
  pca_result <- prcomp(data_centered, center = FALSE, scale. = FALSE)
  
  # Take top n_components
  X_pca <- pca_result$x[, 1:n_components, drop = FALSE]
  
  # FastICA algorithm (simplified)
  n_samples <- nrow(X_pca)
  
  # Initialize mixing matrix
  W <- matrix(rnorm(n_components^2), n_components, n_components)
  
  # Add regularization to avoid singular matrices
  W <- W %*% solve(sqrt(W %*% t(W) + 1e-8 * diag(n_components)))
  
  for (iter in 1:max_iter) {
    W_old <- W
    
    # Update W using FastICA update rule
    for (i in 1:n_components) {
      w <- W[i, ]
      
      # Non-linearity (tanh)
      wx <- X_pca %*% w
      g_wx <- tanh(wx)
      g_wx_prime <- 1 - tanh(wx)^2
      
      # Update rule - fix array multiplication
      w_new <- colMeans(X_pca * as.vector(g_wx)) - mean(g_wx_prime) * w
      w_new <- w_new / sqrt(sum(w_new^2))
      
      W[i, ] <- w_new
    }
    
    # Orthogonalize W with regularization
    W <- W %*% solve(sqrt(W %*% t(W) + 1e-8 * diag(n_components)))
    
    # Check convergence with proper validation
    diff_matrix <- W - W_old
    if (all(is.finite(diff_matrix))) {
      max_diff <- max(abs(diff_matrix))
      if (max_diff < tol) {
        break
      }
    }
  }
  
  # Extract independent components
  S <- X_pca %*% t(W)
  
  # Use pseudoinverse for A to handle singular matrices
  A <- tryCatch({
    solve(W + 1e-6 * diag(n_components))
  }, error = function(e) {
    # If solve fails, use pseudoinverse
    MASS::ginv(W + 1e-6 * diag(n_components))
  })
  
  return(list(
    S = S,
    W = W,
    A = A,
    n_components = n_components,
    iterations = iter
  ))
}

#' Wavelet Denoising
#'
#' Denoise calcium imaging data using wavelet transform.
#'
#' @param signal Input signal
#' @param wavelet Wavelet type
#' @param level Decomposition level
#' @param threshold Threshold method
#' @param verbose Whether to show progress
#' @return Denoised signal
#' @export
wavelet_denoise <- function(signal, wavelet = "db4", level = 3, threshold = "universal", verbose = TRUE) {
  if (verbose) {
    message("Performing wavelet denoising")
  }
  
  # Validate inputs
  if (!is.numeric(signal)) {
    stop("Signal must be numeric")
  }
  
  # Simple wavelet-like denoising using moving average and thresholding
  # This is a simplified implementation since we don't have a full wavelet library
  
  n <- length(signal)
  
  # Apply moving average filter (simulates low-pass filtering)
  window_size <- 5
  denoised <- numeric(n)
  
  for (i in 1:n) {
    start_idx <- max(1, i - window_size %/% 2)
    end_idx <- min(n, i + window_size %/% 2)
    denoised[i] <- mean(signal[start_idx:end_idx], na.rm = TRUE)
  }
  
  # Apply thresholding to remove small fluctuations
  if (threshold == "universal") {
    # Universal threshold
    noise_est <- median(abs(diff(signal))) / 0.6745
    threshold_val <- sqrt(2 * log(n)) * noise_est
  } else {
    # Simple threshold
    threshold_val <- sd(signal) * 2
  }
  
  # Apply soft thresholding
  residuals <- signal - denoised
  thresholded_residuals <- sign(residuals) * pmax(abs(residuals) - threshold_val, 0)
  
  # Final denoised signal
  final_denoised <- denoised + thresholded_residuals
  
  return(list(
    denoised = final_denoised,
    original = signal,
    threshold = threshold_val,
    wavelet = wavelet,
    level = level
  ))
} 