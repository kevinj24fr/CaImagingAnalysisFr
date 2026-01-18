#' Tensor Decomposition for Multi-way Neural Data
#'
#' CP (CANDECOMP/PARAFAC) and Tucker decomposition for analyzing
#' neurons x time x trials x conditions data structures.
#'
#' @name tensor_decomposition
NULL

# ============================================================================
# CP (CANDECOMP/PARAFAC) Decomposition
# ============================================================================

#' CP Tensor Decomposition
#'
#' Decomposes a tensor as a sum of rank-one tensors using Alternating Least Squares.
#' Natural for neurons x time x trials data to extract interpretable factors.
#'
#' @param X 3D or 4D array (neurons x time x trials [x conditions])
#' @param rank Number of components (rank of decomposition)
#' @param n_iter Maximum ALS iterations
#' @param tol Convergence tolerance
#' @param init Initialization method ("random", "svd", "hosvd")
#' @param normalize Normalize factor matrices at each iteration
#' @param verbose Print progress
#'
#' @return List with:
#'   - factors: List of factor matrices for each mode
#'   - weights: Component weights (lambdas)
#'   - fit: R-squared fit quality
#'   - reconstruction: Reconstructed tensor
#'   - iter: Number of iterations
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Decompose trial-structured data
#' # X is neurons x time x trials
#' cp <- cp_decompose(X, rank = 5)
#'
#' # Neuron factors show which cells participate
#' # Time factors show temporal patterns
#' # Trial factors show trial-to-trial variation
#' }
cp_decompose <- function(X,
                         rank = 3,
                         n_iter = 100,
                         tol = 1e-6,
                         init = c("svd", "random", "hosvd"),
                         normalize = TRUE,
                         verbose = FALSE) {

  init <- match.arg(init)

  # Ensure array
  if (!is.array(X)) {
    stop("X must be a 3D or 4D array")
  }

  dims <- dim(X)
  n_modes <- length(dims)

  if (n_modes < 3 || n_modes > 4) {
    stop("X must be 3D or 4D array")
  }

  # Initialize factor matrices
  factors <- initialize_cp(X, rank, init)

  # ALS iterations
  fit_history <- numeric(n_iter)
  X_norm_sq <- sum(X^2, na.rm = TRUE)

  for (iter in seq_len(n_iter)) {
    # Update each factor matrix
    for (mode in seq_len(n_modes)) {
      # Compute Khatri-Rao product of all other factors
      kr_prod <- khatri_rao_all(factors[-mode])

      # Unfold tensor along this mode
      X_unf <- unfold_tensor(X, mode)

      # Least squares update
      V <- compute_hadamard_gramians(factors[-mode])
      factors[[mode]] <- X_unf %*% kr_prod %*% solve(V + 1e-8 * diag(rank))
    }

    # Normalize columns (store norms in weights)
    if (normalize) {
      weights <- numeric(rank)
      for (r in seq_len(rank)) {
        col_norms <- sapply(factors, function(f) sqrt(sum(f[, r]^2)))
        weights[r] <- prod(col_norms)
        for (mode in seq_len(n_modes)) {
          if (col_norms[mode] > 0) {
            factors[[mode]][, r] <- factors[[mode]][, r] / col_norms[mode]
          }
        }
      }
    } else {
      weights <- rep(1, rank)
    }

    # Compute fit
    X_hat <- reconstruct_cp(factors, weights)
    resid_sq <- sum((X - X_hat)^2, na.rm = TRUE)
    fit <- 1 - resid_sq / X_norm_sq
    fit_history[iter] <- fit

    if (verbose && iter %% 10 == 0) {
      message(sprintf("Iteration %d: R² = %.4f", iter, fit))
    }

    # Check convergence
    if (iter > 1 && abs(fit_history[iter] - fit_history[iter - 1]) < tol) {
      if (verbose) message("Converged at iteration ", iter)
      break
    }
  }

  # Sort by weight magnitude
  ord <- order(abs(weights), decreasing = TRUE)
  weights <- weights[ord]
  factors <- lapply(factors, function(f) f[, ord, drop = FALSE])

  # Name dimensions
  mode_names <- c("neuron", "time", "trial", "condition")[1:n_modes]
  names(factors) <- mode_names

  structure(
    list(
      factors = factors,
      weights = weights,
      fit = fit,
      fit_history = fit_history[1:iter],
      dims = dims,
      rank = rank,
      n_modes = n_modes,
      iter = iter
    ),
    class = "cp_decomposition"
  )
}

#' Initialize CP Decomposition
#' @keywords internal
initialize_cp <- function(X, rank, method) {
  dims <- dim(X)
  n_modes <- length(dims)

  if (method == "random") {
    factors <- lapply(dims, function(d) matrix(rnorm(d * rank), d, rank))
  } else if (method == "svd") {
    # SVD-based initialization on unfoldings
    factors <- vector("list", n_modes)
    for (mode in seq_len(n_modes)) {
      X_unf <- unfold_tensor(X, mode)
      svd_result <- svd(X_unf, nu = rank, nv = 0)
      factors[[mode]] <- svd_result$u[, 1:min(rank, ncol(svd_result$u)), drop = FALSE]
      if (ncol(factors[[mode]]) < rank) {
        factors[[mode]] <- cbind(factors[[mode]],
                                  matrix(rnorm(dims[mode] * (rank - ncol(factors[[mode]]))),
                                         dims[mode], rank - ncol(factors[[mode]])))
      }
    }
  } else {
    # HOSVD initialization
    factors <- hosvd_init(X, rank)
  }

  factors
}

#' HOSVD Initialization
#' @keywords internal
hosvd_init <- function(X, rank) {
  dims <- dim(X)
  n_modes <- length(dims)
  factors <- vector("list", n_modes)

  for (mode in seq_len(n_modes)) {
    X_unf <- unfold_tensor(X, mode)
    svd_result <- svd(X_unf, nu = min(rank, dims[mode]), nv = 0)
    factors[[mode]] <- svd_result$u
    if (ncol(factors[[mode]]) < rank) {
      factors[[mode]] <- cbind(factors[[mode]],
                                matrix(rnorm(dims[mode] * (rank - ncol(factors[[mode]]))),
                                       dims[mode], rank - ncol(factors[[mode]])))
    }
  }

  factors
}

#' Unfold Tensor Along a Mode
#' @keywords internal
unfold_tensor <- function(X, mode) {
  dims <- dim(X)
  n_modes <- length(dims)

  # Permute so mode is first
  perm <- c(mode, setdiff(1:n_modes, mode))
  X_perm <- aperm(X, perm)

  # Reshape to matrix
  matrix(X_perm, nrow = dims[mode], ncol = prod(dims[-mode]))
}

#' Khatri-Rao Product of All Factor Matrices
#' @keywords internal
khatri_rao_all <- function(factors) {
  n_factors <- length(factors)
  if (n_factors == 1) return(factors[[1]])

  result <- factors[[n_factors]]
  for (i in (n_factors - 1):1) {
    result <- khatri_rao(factors[[i]], result)
  }
  result
}

#' Khatri-Rao Product (Column-wise Kronecker)
#' @keywords internal
khatri_rao <- function(A, B) {
  if (ncol(A) != ncol(B)) stop("Matrices must have same number of columns")

  R <- ncol(A)
  result <- matrix(0, nrow(A) * nrow(B), R)

  for (r in 1:R) {
    result[, r] <- as.vector(outer(A[, r], B[, r]))
  }

  result
}

#' Compute Hadamard Product of Gramian Matrices
#' @keywords internal
compute_hadamard_gramians <- function(factors) {
  rank <- ncol(factors[[1]])
  result <- matrix(1, rank, rank)

  for (f in factors) {
    result <- result * (t(f) %*% f)
  }

  result
}

#' Reconstruct Tensor from CP Decomposition
#' @keywords internal
reconstruct_cp <- function(factors, weights) {
  dims <- sapply(factors, nrow)
  n_modes <- length(factors)
  rank <- length(weights)

  # Build tensor as sum of rank-1 tensors
  result <- array(0, dims)

  for (r in seq_len(rank)) {
    outer_prod <- factors[[1]][, r]
    for (mode in 2:n_modes) {
      outer_prod <- outer(outer_prod, factors[[mode]][, r])
    }
    result <- result + weights[r] * outer_prod
  }

  result
}

#' @export
print.cp_decomposition <- function(x, ...) {
  cat("CP Tensor Decomposition\n")
  cat("=======================\n")
  cat(sprintf("Dimensions: %s\n", paste(x$dims, collapse = " x ")))
  cat(sprintf("Rank: %d\n", x$rank))
  cat(sprintf("R² fit: %.4f\n", x$fit))
  cat(sprintf("Iterations: %d\n", x$iter))
  cat("\nComponent weights:\n")
  cat(sprintf("  %s\n", paste(round(x$weights, 3), collapse = ", ")))
  invisible(x)
}

# ============================================================================
# Tucker Decomposition
# ============================================================================

#' Tucker Tensor Decomposition
#'
#' Decomposes a tensor into a core tensor multiplied by factor matrices.
#' More flexible than CP, allows different ranks per mode.
#'
#' @param X 3D or 4D array
#' @param ranks Vector of ranks for each mode (or single value for all)
#' @param n_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param init Initialization method ("hosvd", "random")
#' @param verbose Print progress
#'
#' @return List with:
#'   - core: Core tensor
#'   - factors: List of factor matrices
#'   - fit: R-squared fit quality
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Different compression per mode
#' tucker <- tucker_decompose(X, ranks = c(10, 5, 3))
#' }
tucker_decompose <- function(X,
                              ranks = 3,
                              n_iter = 50,
                              tol = 1e-6,
                              init = c("hosvd", "random"),
                              verbose = FALSE) {

  init <- match.arg(init)

  dims <- dim(X)
  n_modes <- length(dims)

  # Handle single rank input
  if (length(ranks) == 1) {
    ranks <- pmin(rep(ranks, n_modes), dims)
  }
  ranks <- pmin(ranks, dims)

  # Initialize factor matrices
  if (init == "hosvd") {
    factors <- vector("list", n_modes)
    for (mode in seq_len(n_modes)) {
      X_unf <- unfold_tensor(X, mode)
      svd_result <- svd(X_unf, nu = ranks[mode], nv = 0)
      factors[[mode]] <- svd_result$u[, 1:ranks[mode], drop = FALSE]
    }
  } else {
    factors <- lapply(seq_len(n_modes), function(m) {
      qr.Q(qr(matrix(rnorm(dims[m] * ranks[m]), dims[m], ranks[m])))
    })
  }

  # HOOI (Higher-Order Orthogonal Iteration)
  X_norm_sq <- sum(X^2, na.rm = TRUE)
  fit_history <- numeric(n_iter)

  for (iter in seq_len(n_iter)) {
    for (mode in seq_len(n_modes)) {
      # Compute n-mode product with all factors except mode
      Y <- X
      for (m in setdiff(seq_len(n_modes), mode)) {
        Y <- ttm(Y, t(factors[[m]]), m)
      }

      # Update factor via SVD of unfolding
      Y_unf <- unfold_tensor(Y, mode)
      svd_result <- svd(Y_unf, nu = ranks[mode], nv = 0)
      factors[[mode]] <- svd_result$u[, 1:ranks[mode], drop = FALSE]
    }

    # Compute core tensor
    core <- X
    for (m in seq_len(n_modes)) {
      core <- ttm(core, t(factors[[m]]), m)
    }

    # Compute fit
    X_hat <- core
    for (m in seq_len(n_modes)) {
      X_hat <- ttm(X_hat, factors[[m]], m)
    }

    resid_sq <- sum((X - X_hat)^2, na.rm = TRUE)
    fit <- 1 - resid_sq / X_norm_sq
    fit_history[iter] <- fit

    if (verbose && iter %% 10 == 0) {
      message(sprintf("Iteration %d: R² = %.4f", iter, fit))
    }

    if (iter > 1 && abs(fit_history[iter] - fit_history[iter - 1]) < tol) {
      if (verbose) message("Converged at iteration ", iter)
      break
    }
  }

  mode_names <- c("neuron", "time", "trial", "condition")[1:n_modes]
  names(factors) <- mode_names

  structure(
    list(
      core = core,
      factors = factors,
      ranks = ranks,
      fit = fit,
      fit_history = fit_history[1:iter],
      dims = dims,
      n_modes = n_modes,
      iter = iter
    ),
    class = "tucker_decomposition"
  )
}

#' Tensor Times Matrix (n-mode product)
#' @keywords internal
ttm <- function(X, M, mode) {
  dims <- dim(X)
  n_modes <- length(dims)

  # Unfold, multiply, refold
  X_unf <- unfold_tensor(X, mode)
  Y_unf <- M %*% X_unf

  # New dimensions
  new_dims <- dims
  new_dims[mode] <- nrow(M)

  # Refold
  perm <- c(mode, setdiff(1:n_modes, mode))
  inv_perm <- order(perm)

  Y <- array(Y_unf, dim = new_dims[perm])
  aperm(Y, inv_perm)
}

#' @export
print.tucker_decomposition <- function(x, ...) {
  cat("Tucker Tensor Decomposition\n")
  cat("===========================\n")
  cat(sprintf("Input dimensions: %s\n", paste(x$dims, collapse = " x ")))
  cat(sprintf("Core dimensions: %s\n", paste(x$ranks, collapse = " x ")))
  cat(sprintf("Compression: %.1f%%\n", 100 * (1 - prod(x$ranks) / prod(x$dims))))
  cat(sprintf("R² fit: %.4f\n", x$fit))
  cat(sprintf("Iterations: %d\n", x$iter))
  invisible(x)
}

# ============================================================================
# Tensor Utilities
# ============================================================================

#' Create Trial Tensor from Traces
#'
#' Reshape traces matrix and trial structure into 3D tensor.
#'
#' @param traces Matrix (cells x total_time)
#' @param trial_starts Vector of trial start indices
#' @param trial_length Length of each trial (frames)
#'
#' @return 3D array (cells x time x trials)
#' @export
create_trial_tensor <- function(traces, trial_starts, trial_length) {
  n_cells <- nrow(traces)
  n_trials <- length(trial_starts)

  X <- array(NA, dim = c(n_cells, trial_length, n_trials))

  for (t in seq_len(n_trials)) {
    start <- trial_starts[t]
    end <- min(start + trial_length - 1, ncol(traces))
    len <- end - start + 1

    if (len > 0) {
      X[, 1:len, t] <- traces[, start:end]
    }
  }

  X
}

#' Extract Temporal Factors from CP/Tucker
#'
#' Get the temporal components from tensor decomposition.
#'
#' @param decomp CP or Tucker decomposition result
#'
#' @return Matrix of temporal factors (time x components)
#' @export
get_temporal_factors <- function(decomp) {
  if (inherits(decomp, "cp_decomposition")) {
    decomp$factors$time
  } else if (inherits(decomp, "tucker_decomposition")) {
    decomp$factors$time
  } else {
    stop("decomp must be cp_decomposition or tucker_decomposition")
  }
}

#' Extract Neuron Factors from CP/Tucker
#'
#' Get the neuron loading factors from tensor decomposition.
#'
#' @param decomp CP or Tucker decomposition result
#'
#' @return Matrix of neuron factors (neurons x components)
#' @export
get_neuron_factors <- function(decomp) {
  if (inherits(decomp, "cp_decomposition")) {
    decomp$factors$neuron
  } else if (inherits(decomp, "tucker_decomposition")) {
    decomp$factors$neuron
  } else {
    stop("decomp must be cp_decomposition or tucker_decomposition")
  }
}

#' Extract Trial Factors from CP/Tucker
#'
#' Get the trial variation factors from tensor decomposition.
#'
#' @param decomp CP or Tucker decomposition result
#'
#' @return Matrix of trial factors (trials x components)
#' @export
get_trial_factors <- function(decomp) {
  if (inherits(decomp, "cp_decomposition")) {
    decomp$factors$trial
  } else if (inherits(decomp, "tucker_decomposition")) {
    decomp$factors$trial
  } else {
    stop("decomp must be cp_decomposition or tucker_decomposition")
  }
}

#' Plot Tensor Decomposition Results
#'
#' @param x CP or Tucker decomposition
#' @param which Which factors to plot ("neuron", "time", "trial", "all")
#' @param n_components Number of components to show
#' @param ... Additional arguments
#'
#' @export
plot_tensor <- function(x, which = "all", n_components = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  if (is.null(n_components)) {
    n_components <- min(5, x$rank)
  }

  plots <- list()

  if (which %in% c("all", "neuron") && "neuron" %in% names(x$factors)) {
    f <- x$factors$neuron[, 1:n_components, drop = FALSE]
    df <- data.frame(
      neuron = rep(1:nrow(f), n_components),
      component = rep(1:n_components, each = nrow(f)),
      value = as.vector(f)
    )
    plots$neuron <- ggplot2::ggplot(df, ggplot2::aes(x = neuron, y = value)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(~paste("Component", component), scales = "free_y") +
      ggplot2::labs(title = "Neuron Factors", x = "Neuron", y = "Loading") +
      ggplot2::theme_minimal()
  }

  if (which %in% c("all", "time") && "time" %in% names(x$factors)) {
    f <- x$factors$time[, 1:n_components, drop = FALSE]
    df <- data.frame(
      time = rep(1:nrow(f), n_components),
      component = rep(1:n_components, each = nrow(f)),
      value = as.vector(f)
    )
    plots$time <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = value, color = factor(component))) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(title = "Temporal Factors", x = "Time", y = "Value", color = "Component") +
      ggplot2::theme_minimal()
  }

  if (which %in% c("all", "trial") && "trial" %in% names(x$factors)) {
    f <- x$factors$trial[, 1:n_components, drop = FALSE]
    df <- data.frame(
      trial = rep(1:nrow(f), n_components),
      component = rep(1:n_components, each = nrow(f)),
      value = as.vector(f)
    )
    plots$trial <- ggplot2::ggplot(df, ggplot2::aes(x = trial, y = value)) +
      ggplot2::geom_point() +
      ggplot2::geom_line(alpha = 0.5) +
      ggplot2::facet_wrap(~paste("Component", component), scales = "free_y") +
      ggplot2::labs(title = "Trial Factors", x = "Trial", y = "Loading") +
      ggplot2::theme_minimal()
  }

  if (length(plots) == 1) {
    return(plots[[1]])
  }

  plots
}

#' Cross-validate Tensor Rank
#'
#' Use missing value imputation to select optimal rank.
#'
#' @param X Tensor array
#' @param ranks Vector of ranks to try
#' @param n_folds Number of CV folds
#' @param method "cp" or "tucker"
#'
#' @return List with optimal rank and CV errors
#' @export
cv_tensor_rank <- function(X, ranks = 1:10, n_folds = 5, method = c("cp", "tucker")) {
  method <- match.arg(method)

  n_elements <- length(X)
  fold_idx <- sample(rep(1:n_folds, length.out = n_elements))

  cv_errors <- matrix(NA, length(ranks), n_folds)

  for (r_idx in seq_along(ranks)) {
    rank <- ranks[r_idx]

    for (fold in 1:n_folds) {
      # Create missing mask
      X_train <- X
      test_mask <- fold_idx == fold
      X_train[test_mask] <- NA

      # Fit with missing values
      if (method == "cp") {
        fit <- cp_decompose(X_train, rank = rank, verbose = FALSE)
        X_hat <- reconstruct_cp(fit$factors, fit$weights)
      } else {
        fit <- tucker_decompose(X_train, ranks = rank, verbose = FALSE)
        X_hat <- fit$core
        for (m in seq_along(fit$factors)) {
          X_hat <- ttm(X_hat, fit$factors[[m]], m)
        }
      }

      # Error on held-out entries
      cv_errors[r_idx, fold] <- mean((X[test_mask] - X_hat[test_mask])^2, na.rm = TRUE)
    }
  }

  mean_errors <- rowMeans(cv_errors, na.rm = TRUE)
  se_errors <- apply(cv_errors, 1, sd, na.rm = TRUE) / sqrt(n_folds)

  optimal_rank <- ranks[which.min(mean_errors)]

  list(
    optimal_rank = optimal_rank,
    ranks = ranks,
    mean_error = mean_errors,
    se_error = se_errors,
    cv_errors = cv_errors
  )
}
