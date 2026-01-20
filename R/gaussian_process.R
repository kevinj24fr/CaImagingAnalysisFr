#' Gaussian Process Regression for Neural Data
#'
#' GP regression for smooth trajectory interpolation,
#' missing data imputation, and uncertainty quantification.
#'
#' @name gaussian_process
NULL

# ============================================================================
# Gaussian Process Regression
# ============================================================================

#' Fit Gaussian Process Regression
#'
#' Fit a GP to neural traces for smooth interpolation with uncertainty.
#'
#' @param x Input locations (e.g., time points)
#' @param y Observations (vector or matrix with observations in columns)
#' @param kernel Kernel type ("rbf", "matern32", "matern52", "periodic")
#' @param length_scale Initial length scale parameter
#' @param signal_var Signal variance
#' @param noise_var Observation noise variance
#' @param optimize Optimize hyperparameters via marginal likelihood
#' @param n_restarts Number of random restarts for optimization
#'
#' @return GP fit object with:
#'   - mean: Posterior mean function
#'   - var: Posterior variance function
#'   - hyperparams: Optimized hyperparameters
#'   - log_marginal_lik: Log marginal likelihood
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Smooth a noisy trace with uncertainty
#' time <- 1:100
#' trace <- sin(time/10) + rnorm(100, 0, 0.2)
#' gp <- fit_gp(time, trace, kernel = "rbf")
#'
#' # Predict at finer resolution
#' pred <- predict_gp(gp, seq(1, 100, by = 0.5))
#' }
fit_gp <- function(x, y,
                   kernel = c("rbf", "matern32", "matern52", "periodic"),
                   length_scale = NULL,
                   signal_var = NULL,
                   noise_var = NULL,
                   optimize = TRUE,
                   n_restarts = 3) {

  kernel <- match.arg(kernel)

  x <- as.numeric(x)
  y <- as.numeric(y)
  n <- length(x)

  # Initialize hyperparameters
  if (is.null(length_scale)) {
    length_scale <- diff(range(x)) / 10
  }
  if (is.null(signal_var)) {
    signal_var <- var(y)
  }
  if (is.null(noise_var)) {
    noise_var <- var(y) * 0.1
  }

  # Optimize hyperparameters
  if (optimize) {
    opt_result <- optimize_gp_hyperparams(x, y, kernel, length_scale,
                                           signal_var, noise_var, n_restarts)
    length_scale <- opt_result$length_scale
    signal_var <- opt_result$signal_var
    noise_var <- opt_result$noise_var
    log_marginal_lik <- opt_result$log_marginal_lik
  } else {
    log_marginal_lik <- gp_log_marginal_likelihood(
      x, y, kernel, length_scale, signal_var, noise_var)
  }

  # Compute kernel matrix
  K <- compute_kernel(x, x, kernel, length_scale, signal_var)
  K_noisy <- K + noise_var * diag(n)

  # Cholesky decomposition for stable inversion
  L <- chol(K_noisy + 1e-8 * diag(n))
  alpha <- backsolve(L, forwardsolve(t(L), y))

  structure(
    list(
      x = x,
      y = y,
      kernel = kernel,
      length_scale = length_scale,
      signal_var = signal_var,
      noise_var = noise_var,
      L = L,
      alpha = alpha,
      log_marginal_lik = log_marginal_lik
    ),
    class = "gp_fit"
  )
}

#' Compute Kernel Matrix
#' @keywords internal
compute_kernel <- function(x1, x2, kernel, length_scale, signal_var) {
  n1 <- length(x1)
  n2 <- length(x2)

  # Squared distance matrix
  D2 <- outer(x1, x2, function(a, b) (a - b)^2)

  if (kernel == "rbf") {
    # Radial basis function (squared exponential)
    K <- signal_var * exp(-D2 / (2 * length_scale^2))
  } else if (kernel == "matern32") {
    # Matern 3/2
    D <- sqrt(D2)
    r <- sqrt(3) * D / length_scale
    K <- signal_var * (1 + r) * exp(-r)
  } else if (kernel == "matern52") {
    # Matern 5/2
    D <- sqrt(D2)
    r <- sqrt(5) * D / length_scale
    K <- signal_var * (1 + r + r^2/3) * exp(-r)
  } else if (kernel == "periodic") {
    # Periodic kernel
    D <- sqrt(D2)
    K <- signal_var * exp(-2 * sin(pi * D / length_scale)^2 / length_scale^2)
  }

  K
}

#' GP Log Marginal Likelihood
#' @keywords internal
gp_log_marginal_likelihood <- function(x, y, kernel, length_scale, signal_var, noise_var) {
  n <- length(x)

  K <- compute_kernel(x, x, kernel, length_scale, signal_var)
  K_noisy <- K + noise_var * diag(n)

  # Use Cholesky for numerical stability
  tryCatch({
    L <- chol(K_noisy + 1e-8 * diag(n))
    alpha <- backsolve(L, forwardsolve(t(L), y))

    # Log marginal likelihood
    -0.5 * sum(y * alpha) - sum(log(diag(L))) - n/2 * log(2 * pi)
  }, error = function(e) {
    -Inf
  })
}

#' Optimize GP Hyperparameters
#' @keywords internal
optimize_gp_hyperparams <- function(x, y, kernel, length_scale, signal_var,
                                     noise_var, n_restarts) {

  # Objective function (negative log marginal likelihood)
  objective <- function(log_params) {
    ls <- exp(log_params[1])
    sv <- exp(log_params[2])
    nv <- exp(log_params[3])
    -gp_log_marginal_likelihood(x, y, kernel, ls, sv, nv)
  }

  best_value <- Inf
  best_params <- c(length_scale, signal_var, noise_var)

  for (restart in seq_len(n_restarts)) {
    if (restart == 1) {
      init <- log(c(length_scale, signal_var, noise_var))
    } else {
      init <- log(c(
        length_scale * exp(rnorm(1, 0, 0.5)),
        signal_var * exp(rnorm(1, 0, 0.5)),
        noise_var * exp(rnorm(1, 0, 0.5))
      ))
    }

    opt <- tryCatch({
      optim(init, objective, method = "L-BFGS-B",
            lower = log(c(1e-5, 1e-5, 1e-5)),
            upper = log(c(1e5, 1e5, 1e5)))
    }, error = function(e) {
      list(value = Inf, par = init)
    })

    if (opt$value < best_value) {
      best_value <- opt$value
      best_params <- exp(opt$par)
    }
  }

  list(
    length_scale = best_params[1],
    signal_var = best_params[2],
    noise_var = best_params[3],
    log_marginal_lik = -best_value
  )
}

#' Predict with Gaussian Process
#'
#' @param object GP fit object
#' @param newx New input locations for prediction
#' @param return_var Return variance in addition to mean
#' @param return_cov Return full covariance matrix
#' @param ... Additional arguments (ignored)
#'
#' @return List with mean, variance, and optionally covariance
#' @export
predict_gp <- function(object, newx, return_var = TRUE, return_cov = FALSE, ...) {
  if (!inherits(object, "gp_fit")) {
    stop("object must be a gp_fit")
  }

  newx <- as.numeric(newx)
  n_new <- length(newx)

  # Kernel between new and training points
  K_star <- compute_kernel(newx, object$x, object$kernel,
                           object$length_scale, object$signal_var)

  # Posterior mean
  mean <- as.vector(K_star %*% object$alpha)

  result <- list(mean = mean)

  if (return_var || return_cov) {
    # K_star_star
    K_ss <- compute_kernel(newx, newx, object$kernel,
                           object$length_scale, object$signal_var)

    # v = L \ K_star^T
    v <- forwardsolve(t(object$L), t(K_star))

    # Posterior covariance
    cov <- K_ss - t(v) %*% v

    if (return_var) {
      result$var <- pmax(0, diag(cov))
      result$sd <- sqrt(result$var)
    }

    if (return_cov) {
      result$cov <- cov
    }
  }

  result
}

#' @export
print.gp_fit <- function(x, ...) {
  cat("Gaussian Process Regression\n")
  cat("===========================\n")
  cat(sprintf("Kernel: %s\n", x$kernel))
  cat(sprintf("Training points: %d\n", length(x$x)))
  cat(sprintf("Length scale: %.4f\n", x$length_scale))
  cat(sprintf("Signal variance: %.4f\n", x$signal_var))
  cat(sprintf("Noise variance: %.4f\n", x$noise_var))
  cat(sprintf("Log marginal likelihood: %.2f\n", x$log_marginal_lik))
  invisible(x)
}

#' Plot GP Fit
#'
#' @param x GP fit object
#' @param n_points Number of prediction points
#' @param conf_level Confidence level for bands
#' @param ... Additional arguments
#'
#' @export
plot_gp <- function(x, n_points = 200, conf_level = 0.95, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  # Prediction grid
  x_new <- seq(min(x$x), max(x$x), length.out = n_points)
  pred <- predict_gp(x, x_new, return_var = TRUE)

  z <- qnorm(1 - (1 - conf_level) / 2)

  df_pred <- data.frame(
    x = x_new,
    mean = pred$mean,
    lower = pred$mean - z * pred$sd,
    upper = pred$mean + z * pred$sd
  )

  df_obs <- data.frame(x = x$x, y = x$y)

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = df_pred,
                         ggplot2::aes(x = x, ymin = lower, ymax = upper),
                         fill = "steelblue", alpha = 0.3) +
    ggplot2::geom_line(data = df_pred,
                       ggplot2::aes(x = x, y = mean),
                       color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(data = df_obs,
                        ggplot2::aes(x = x, y = y),
                        size = 2, alpha = 0.6) +
    ggplot2::labs(title = "Gaussian Process Regression",
                  x = "Input", y = "Output") +
    ggplot2::theme_minimal()
}

# ============================================================================
# Multi-output GP
# ============================================================================

#' Fit Multi-output GP for Multiple Traces
#'
#' Fit separate GPs to each cell trace with shared hyperparameters.
#'
#' @param time Time vector
#' @param traces Matrix of traces (cells x time)
#' @param kernel Kernel type
#' @param share_hyperparams Share hyperparameters across cells
#' @param length_scale Fixed length scale (NULL for automatic optimization)
#' @param ... Additional arguments to fit_gp
#'
#' @return List of GP fits
#' @export
fit_gp_traces <- function(time, traces, kernel = "rbf", share_hyperparams = TRUE,
                          length_scale = NULL, ...) {
  n_cells <- nrow(traces)

  # If length_scale is provided, use it directly
  if (!is.null(length_scale)) {
    gp_list <- lapply(seq_len(n_cells), function(i) {
      fit_gp(time, traces[i, ], kernel = kernel,
             length_scale = length_scale,
             optimize = TRUE, ...)
    })
  } else if (share_hyperparams) {
    # Fit one GP to get shared hyperparams
    sample_trace <- traces[1, ]
    gp_sample <- fit_gp(time, sample_trace, kernel = kernel, optimize = TRUE, ...)

    # Apply to all traces
    gp_list <- lapply(seq_len(n_cells), function(i) {
      fit_gp(time, traces[i, ], kernel = kernel,
             length_scale = gp_sample$length_scale,
             signal_var = gp_sample$signal_var,
             noise_var = gp_sample$noise_var,
             optimize = FALSE)
    })
  } else {
    gp_list <- lapply(seq_len(n_cells), function(i) {
      fit_gp(time, traces[i, ], kernel = kernel, optimize = TRUE, ...)
    })
  }

  structure(
    list(
      gp_list = gp_list,
      time = time,
      n_cells = n_cells,
      kernel = kernel,
      length_scale = length_scale
    ),
    class = "gp_traces"
  )
}

#' Predict with Multi-output GP
#'
#' @param object gp_traces object
#' @param newtime New time points
#' @param return_var Return variance
#' @param ... Additional arguments
#'
#' @return List with mean matrix and optionally variance matrix
#' @export
predict_gp_traces <- function(object, newtime, return_var = TRUE, ...) {
  n_cells <- object$n_cells
  n_time <- length(newtime)

  mean_mat <- matrix(0, n_cells, n_time)
  var_mat <- if (return_var) matrix(0, n_cells, n_time) else NULL

  for (i in seq_len(n_cells)) {
    pred <- predict_gp(object$gp_list[[i]], newtime, return_var = return_var)
    mean_mat[i, ] <- pred$mean
    if (return_var) {
      var_mat[i, ] <- pred$var
    }
  }

  list(mean = mean_mat, var = var_mat)
}

#' Impute Missing Values with GP
#'
#' @param time Time vector
#' @param trace Trace with NAs
#' @param kernel Kernel type
#'
#' @return List with imputed trace and uncertainty
#' @export
impute_missing_gp <- function(time, trace, kernel = "rbf") {
  observed <- !is.na(trace)

  if (all(observed)) {
    return(list(imputed = trace, sd = rep(0, length(trace))))
  }

  if (sum(observed) < 3) {
    warning("Too few observations for GP imputation")
    imputed <- trace
    imputed[!observed] <- mean(trace, na.rm = TRUE)
    return(list(imputed = imputed, sd = rep(NA, length(trace))))
  }

  # Fit GP to observed data
  gp <- fit_gp(time[observed], trace[observed], kernel = kernel)

  # Predict at all time points
  pred <- predict_gp(gp, time, return_var = TRUE)

  # Use predictions for missing values
  imputed <- trace
  imputed[!observed] <- pred$mean[!observed]

  list(imputed = imputed, sd = pred$sd, observed = observed)
}

# ============================================================================
# Sparse GP for Large Data
# ============================================================================

#' Sparse Gaussian Process Regression
#'
#' Sparse GP using inducing points for large datasets.
#'
#' @param x Input locations
#' @param y Observations
#' @param n_inducing Number of inducing points
#' @param kernel Kernel type
#' @param optimize Optimize hyperparameters and inducing locations
#'
#' @return Sparse GP fit
#' @export
fit_sparse_gp <- function(x, y, n_inducing = 20, kernel = "rbf", optimize = TRUE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  n <- length(x)

  # Initialize inducing points (spread across input range)
  z <- seq(min(x), max(x), length.out = n_inducing)

  # Initialize hyperparameters
  length_scale <- diff(range(x)) / 5
  signal_var <- var(y)
  noise_var <- var(y) * 0.1

  if (optimize) {
    # Optimize via variational lower bound (FITC approximation)
    # Simplified optimization
    gp_full <- fit_gp(x, y, kernel = kernel, optimize = TRUE)
    length_scale <- gp_full$length_scale
    signal_var <- gp_full$signal_var
    noise_var <- gp_full$noise_var
  }

  # Compute required matrices
  Kzz <- compute_kernel(z, z, kernel, length_scale, signal_var)
  Kzx <- compute_kernel(z, x, kernel, length_scale, signal_var)
  Kxz <- t(Kzx)

  # FITC approximation
  Lzz <- chol(Kzz + 1e-8 * diag(n_inducing))
  V <- forwardsolve(t(Lzz), Kzx)

  # Diagonal of Kxx (for noise correction)
  Kxx_diag <- rep(signal_var, n)
  Qxx_diag <- colSums(V^2)
  Lambda <- pmax(Kxx_diag - Qxx_diag + noise_var, noise_var)

  # Woodbury inversion
  B <- diag(n_inducing) + V %*% diag(1/Lambda) %*% t(V)
  Lb <- chol(B + 1e-8 * diag(n_inducing))

  beta <- forwardsolve(t(Lzz), Kzx %*% (y / Lambda))
  alpha <- backsolve(Lb, forwardsolve(t(Lb), beta))
  alpha <- forwardsolve(t(Lzz), alpha)

  structure(
    list(
      x = x,
      y = y,
      z = z,
      kernel = kernel,
      length_scale = length_scale,
      signal_var = signal_var,
      noise_var = noise_var,
      Lzz = Lzz,
      alpha = alpha,
      n_inducing = n_inducing
    ),
    class = "sparse_gp_fit"
  )
}

#' Predict with Sparse GP
#'
#' @param object Sparse GP fit
#' @param newx New input locations
#' @param return_var Return variance
#' @param ... Additional arguments
#'
#' @return Predictions
#' @export
predict.sparse_gp_fit <- function(object, newx, return_var = TRUE, ...) {
  newx <- as.numeric(newx)

  Kxz <- compute_kernel(newx, object$z, object$kernel,
                        object$length_scale, object$signal_var)

  mean <- as.vector(Kxz %*% object$alpha)

  result <- list(mean = mean)

  if (return_var) {
    v <- forwardsolve(t(object$Lzz), t(Kxz))
    var <- object$signal_var - colSums(v^2)
    result$var <- pmax(0, var)
    result$sd <- sqrt(result$var)
  }

  result
}
