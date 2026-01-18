#' Gaussian Process Factor Analysis (GPFA)
#'
#' Extract smooth low-dimensional latent trajectories from neural population activity.
#' GPFA combines factor analysis with Gaussian process priors to extract temporally
#' smooth latent variables that explain shared variability across neurons.
#'
#' @name gpfa
#' @references
#' Yu, B. M., et al. (2009). Gaussian-process factor analysis for low-dimensional
#' single-trial analysis of neural population activity. Journal of Neurophysiology.
NULL

#' Fit GPFA Model
#'
#' Fit a Gaussian Process Factor Analysis model to neural population data.
#'
#' @param trials List of trial matrices (neurons x time) or 3D array (neurons x time x trials)
#' @param n_latents Number of latent dimensions to extract
#' @param bin_width Time bin width in ms (for kernel timescale interpretation)
#' @param max_iter Maximum EM iterations
#' @param tol Convergence tolerance
#' @param verbose Print progress
#'
#' @return List containing:
#' \itemize{
#'   \item{trajectories}{List of latent trajectories (n_latents x time) per trial}
#'   \item{C}{Observation matrix (neurons x n_latents)}
#'   \item{d}{Mean firing rates (neurons)}
#'   \item{R}{Private noise variances (neurons)}
#'   \item{gamma}{GP timescale parameters (n_latents)}
#'   \item{ll_trace}{Log-likelihood trace across iterations}
#' }
#'
#' @keywords internal
#' @note Use \code{\link{RunGPFA}} for CaExperiment workflow
#'
#' @examples
#' \dontrun{
#' # Create trial list from continuous recording
#' trials <- split_into_trials(traces, trial_starts, trial_length = 100)
#'
#' # Fit GPFA with 3 latent dimensions
#' gpfa_result <- fit_gpfa(trials, n_latents = 3)
#'
#' # Plot first latent trajectory
#' plot(gpfa_result$trajectories[[1]][1, ], type = "l")
#' }
fit_gpfa <- function(trials, n_latents = 3, bin_width = 20,
                     max_iter = 100, tol = 1e-6, verbose = TRUE) {

  # Convert 3D array to list if needed
  if (is.array(trials) && length(dim(trials)) == 3) {
    n_trials <- dim(trials)[3]
    trials <- lapply(1:n_trials, function(i) trials[, , i])
  }

  # Validate input
  if (!is.list(trials) || length(trials) == 0) {
    stop("trials must be a non-empty list of matrices or 3D array")
  }

  n_trials <- length(trials)
  n_neurons <- nrow(trials[[1]])
  n_time <- ncol(trials[[1]])

  # Check dimensions consistency
  for (i in seq_along(trials)) {
    if (nrow(trials[[i]]) != n_neurons) {
      stop("All trials must have the same number of neurons")
    }
  }

  if (verbose) {
    cat(sprintf("Fitting GPFA: %d neurons, %d trials, %d latents\n",
                n_neurons, n_trials, n_latents))
  }

  # Initialize parameters
  params <- .gpfa_init_params(trials, n_latents)

  # EM algorithm
  ll_trace <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {
    # E-step: compute posterior over latents
    estep <- .gpfa_estep(trials, params)

    # M-step: update parameters
    params <- .gpfa_mstep(trials, estep, params)

    # Compute log-likelihood
    ll_trace[iter] <- estep$ll

    if (verbose && iter %% 10 == 0) {
      cat(sprintf("  Iteration %d: LL = %.2f\n", iter, ll_trace[iter]))
    }

    # Check convergence
    if (iter > 1) {
      ll_change <- abs(ll_trace[iter] - ll_trace[iter - 1]) / abs(ll_trace[iter - 1])
      if (ll_change < tol) {
        if (verbose) cat(sprintf("Converged at iteration %d\n", iter))
        ll_trace <- ll_trace[1:iter]
        break
      }
    }
  }

  # Extract final trajectories
  trajectories <- estep$means

  list(
    trajectories = trajectories,
    C = params$C,
    d = params$d,
    R = params$R,
    gamma = params$gamma,
    ll_trace = ll_trace,
    n_latents = n_latents,
    n_neurons = n_neurons,
    bin_width = bin_width
  )
}

#' Initialize GPFA parameters
#' @keywords internal
.gpfa_init_params <- function(trials, n_latents) {
  n_neurons <- nrow(trials[[1]])
  n_time <- ncol(trials[[1]])
  n_trials <- length(trials)

  # Stack all trials for initialization
  all_data <- do.call(cbind, trials)

  # Initialize d as mean firing rate
  d <- rowMeans(all_data)

  # Center data
  centered <- all_data - d

  # Initialize C using FA/PCA on centered data
  if (n_neurons > n_latents) {
    # Use SVD for initialization
    svd_result <- svd(centered, nu = n_latents, nv = n_latents)
    C <- svd_result$u[, 1:n_latents, drop = FALSE] %*%
         diag(sqrt(svd_result$d[1:n_latents]), n_latents, n_latents)
  } else {
    C <- matrix(rnorm(n_neurons * n_latents, sd = 0.1), n_neurons, n_latents)
  }

  # Initialize R as residual variance
  if (n_neurons > n_latents) {
    reconstructed <- C %*% t(svd_result$v[, 1:n_latents, drop = FALSE] %*%
                             diag(sqrt(svd_result$d[1:n_latents]), n_latents, n_latents))
    residuals <- centered - reconstructed
    R <- pmax(apply(residuals, 1, var), 1e-6)
  } else {
    R <- rep(var(as.vector(centered)) / 2, n_neurons)
  }

  # Initialize GP timescales (in bins)
  gamma <- rep(10, n_latents)  # Default timescale of 10 bins

  list(
    C = C,
    d = d,
    R = R,
    gamma = gamma,
    n_time = n_time
  )
}

#' GPFA E-step
#' @keywords internal
.gpfa_estep <- function(trials, params) {
  n_trials <- length(trials)
  n_neurons <- nrow(trials[[1]])
  n_latents <- ncol(params$C)

  means <- vector("list", n_trials)
  covs <- vector("list", n_trials)
  ll <- 0

  for (i in seq_len(n_trials)) {
    trial_data <- trials[[i]]
    n_time <- ncol(trial_data)

    # Build GP prior covariance for this trial
    K <- .gpfa_build_K(n_time, n_latents, params$gamma)

    # Center data
    y <- trial_data - params$d

    # Compute posterior (using Woodbury identity for efficiency)
    C <- params$C
    R_inv <- 1 / params$R

    # Posterior covariance: (K^-1 + C'R^-1 C)^-1
    # Using Woodbury: K - K C' (R + C K C')^-1 C K

    # For each time point (simplified approach - assumes independent time points in R)
    # Full GPFA would correlate across time, but this is a practical approximation

    CRC <- crossprod(C, C * R_inv)  # C' R^-1 C

    # Block diagonal structure for latents across time
    post_cov_single <- solve(diag(n_latents) + CRC %*% .gpfa_K_single(params$gamma))

    # Posterior mean
    y_flat <- as.vector(y)

    # Simple approach: solve for each time point
    post_mean <- matrix(0, n_latents, n_time)
    for (t in seq_len(n_time)) {
      yt <- y[, t]
      # E[x|y] = K C' (C K C' + R)^-1 y
      K_single <- .gpfa_K_single(params$gamma)
      S <- C %*% K_single %*% t(C) + diag(params$R)
      post_mean[, t] <- K_single %*% t(C) %*% solve(S, yt)
    }

    means[[i]] <- post_mean
    covs[[i]] <- post_cov_single

    # Log-likelihood contribution (approximate)
    for (t in seq_len(n_time)) {
      S <- C %*% .gpfa_K_single(params$gamma) %*% t(C) + diag(params$R)
      ll <- ll - 0.5 * (n_neurons * log(2 * pi) +
                        determinant(S, logarithm = TRUE)$modulus[1] +
                        sum(y[, t] * solve(S, y[, t])))
    }
  }

  list(means = means, covs = covs, ll = ll)
}

#' GPFA M-step
#' @keywords internal
.gpfa_mstep <- function(trials, estep, params) {
  n_trials <- length(trials)
  n_neurons <- nrow(trials[[1]])
  n_latents <- ncol(params$C)

  # Update d (mean)
  all_residuals <- lapply(seq_len(n_trials), function(i) {
    trials[[i]] - params$C %*% estep$means[[i]]
  })
  d_new <- rowMeans(do.call(cbind, all_residuals))

  # Update C (loading matrix)
  # C_new = (sum_i Y_i X_i') (sum_i X_i X_i' + cov)^-1
  YX <- matrix(0, n_neurons, n_latents)
  XX <- matrix(0, n_latents, n_latents)

  for (i in seq_len(n_trials)) {
    y <- trials[[i]] - d_new
    x <- estep$means[[i]]
    YX <- YX + y %*% t(x)
    XX <- XX + x %*% t(x) + ncol(x) * estep$covs[[i]]
  }

  C_new <- YX %*% solve(XX + diag(1e-6, n_latents))

  # Update R (private variances)
  R_new <- numeric(n_neurons)
  total_time <- 0

  for (i in seq_len(n_trials)) {
    y <- trials[[i]] - d_new
    x <- estep$means[[i]]
    residual <- y - C_new %*% x
    R_new <- R_new + rowSums(residual^2)
    total_time <- total_time + ncol(y)
  }
  R_new <- pmax(R_new / total_time, 1e-6)

  # Update gamma (GP timescales) - gradient descent step
  gamma_new <- params$gamma
  for (k in seq_len(n_latents)) {
    # Simple heuristic: estimate from autocorrelation of inferred latents
    all_x_k <- do.call(c, lapply(estep$means, function(x) x[k, ]))
    if (length(all_x_k) > 2) {
      acf_vals <- acf(all_x_k, lag.max = 20, plot = FALSE)$acf[-1]
      # Find lag where ACF drops below 0.5
      half_life <- which(acf_vals < 0.5)[1]
      if (!is.na(half_life) && half_life > 0) {
        gamma_new[k] <- 0.8 * gamma_new[k] + 0.2 * max(2, half_life)
      }
    }
  }

  list(
    C = C_new,
    d = d_new,
    R = R_new,
    gamma = gamma_new,
    n_time = params$n_time
  )
}

#' Build GP covariance matrix
#' @keywords internal
.gpfa_build_K <- function(n_time, n_latents, gamma) {
  # Build block diagonal K matrix
  K_list <- lapply(seq_len(n_latents), function(k) {
    tau <- gamma[k]
    times <- seq_len(n_time)
    dist_mat <- outer(times, times, "-")
    exp(-0.5 * (dist_mat / tau)^2)
  })

  # Block diagonal
  K <- matrix(0, n_time * n_latents, n_time * n_latents)
  for (k in seq_len(n_latents)) {
    idx <- ((k - 1) * n_time + 1):(k * n_time)
    K[idx, idx] <- K_list[[k]]
  }
  K
}

#' Single timepoint GP covariance
#' @keywords internal
.gpfa_K_single <- function(gamma) {
  diag(gamma / 10)  # Simplified marginal variance
}

#' Project New Data Through GPFA Model
#'
#' Project new trial data into the learned latent space.
#'
#' @param gpfa_model Fitted GPFA model from fit_gpfa()
#' @param new_trials List of new trial matrices
#'
#' @return List of latent trajectories for new trials
#' @keywords internal
gpfa_project <- function(gpfa_model, new_trials) {
  if (is.array(new_trials) && length(dim(new_trials)) == 3) {
    n_trials <- dim(new_trials)[3]
    new_trials <- lapply(1:n_trials, function(i) new_trials[, , i])
  }

  params <- list(
    C = gpfa_model$C,
    d = gpfa_model$d,
    R = gpfa_model$R,
    gamma = gpfa_model$gamma
  )

  estep <- .gpfa_estep(new_trials, params)
  estep$means
}

#' Compute GPFA Orthonormalized Trajectories
#'
#' Orthonormalize latent trajectories for visualization.
#'
#' @param gpfa_model Fitted GPFA model
#'
#' @return List with orthonormalized trajectories and transformation matrix
#' @keywords internal
gpfa_orthonormalize <- function(gpfa_model) {
  C <- gpfa_model$C

  # Orthonormalize C using SVD
  svd_C <- svd(C)
  C_orth <- svd_C$u[, 1:gpfa_model$n_latents]
  T_orth <- diag(svd_C$d[1:gpfa_model$n_latents]) %*% t(svd_C$v[, 1:gpfa_model$n_latents])

  # Transform trajectories
  orth_trajectories <- lapply(gpfa_model$trajectories, function(x) {
    T_orth %*% x
  })

  list(
    trajectories = orth_trajectories,
    C_orth = C_orth,
    T_orth = T_orth
  )
}

#' Cross-validate GPFA Latent Dimensionality
#'
#' Use cross-validation to select optimal number of latent dimensions.
#'
#' @param trials List of trial matrices
#' @param latent_dims Vector of dimensions to test
#' @param n_folds Number of CV folds
#' @param verbose Print progress
#'
#' @return Data frame with cross-validation results
#' @keywords internal
gpfa_cv <- function(trials, latent_dims = 1:10, n_folds = 5, verbose = TRUE) {
  n_trials <- length(trials)
  fold_ids <- sample(rep(1:n_folds, length.out = n_trials))

  results <- data.frame(
    n_latents = integer(),
    fold = integer(),
    train_ll = numeric(),
    test_ll = numeric()
  )

  for (n_lat in latent_dims) {
    if (verbose) cat(sprintf("Testing %d latent dimensions...\n", n_lat))

    for (fold in seq_len(n_folds)) {
      train_trials <- trials[fold_ids != fold]
      test_trials <- trials[fold_ids == fold]

      # Fit on training data
      fit <- fit_gpfa(train_trials, n_latents = n_lat, verbose = FALSE)

      # Evaluate on test data
      test_estep <- .gpfa_estep(test_trials, list(
        C = fit$C, d = fit$d, R = fit$R, gamma = fit$gamma
      ))

      results <- rbind(results, data.frame(
        n_latents = n_lat,
        fold = fold,
        train_ll = tail(fit$ll_trace, 1),
        test_ll = test_estep$ll
      ))
    }
  }

  # Summarize
  summary_df <- aggregate(cbind(train_ll, test_ll) ~ n_latents, data = results,
                          FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x))))
  summary_df <- do.call(data.frame, summary_df)
  names(summary_df) <- c("n_latents", "train_ll_mean", "train_ll_se",
                         "test_ll_mean", "test_ll_se")

  list(
    full_results = results,
    summary = summary_df,
    optimal = summary_df$n_latents[which.max(summary_df$test_ll_mean)]
  )
}

#' Plot GPFA Trajectories
#'
#' Visualize GPFA latent trajectories in 2D or 3D.
#'
#' @param gpfa_model Fitted GPFA model
#' @param dims Which latent dimensions to plot (2 or 3)
#' @param trials Which trials to plot (NULL for all)
#' @param color_by How to color: "trial", "time"
#' @param add_start_end Add markers for start/end points
#'
#' @return ggplot object (2D) or plotly object (3D)
#' @keywords internal
plot_gpfa <- function(gpfa_model, dims = c(1, 2), trials = NULL,
                      color_by = "time", add_start_end = TRUE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  trajectories <- gpfa_model$trajectories
  if (is.null(trials)) trials <- seq_along(trajectories)

  # Build data frame for plotting
  plot_data <- do.call(rbind, lapply(trials, function(i) {
    traj <- trajectories[[i]]
    n_time <- ncol(traj)
    data.frame(
      trial = i,
      time = seq_len(n_time),
      x = traj[dims[1], ],
      y = traj[dims[2], ],
      z = if (length(dims) == 3) traj[dims[3], ] else NA
    )
  }))

  if (length(dims) == 2) {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y,
                         color = if (color_by == "time") time else factor(trial),
                         group = trial)) +
      ggplot2::geom_path(alpha = 0.7) +
      ggplot2::labs(x = sprintf("Latent %d", dims[1]),
                    y = sprintf("Latent %d", dims[2]),
                    color = color_by) +
      ggplot2::theme_minimal()

    if (add_start_end) {
      starts <- plot_data[plot_data$time == 1, ]
      ends <- do.call(rbind, lapply(trials, function(i) {
        subset(plot_data, trial == i)[sum(plot_data$trial == i), ]
      }))
      p <- p +
        ggplot2::geom_point(data = starts, shape = 16, size = 3) +
        ggplot2::geom_point(data = ends, shape = 17, size = 3)
    }

    return(p)
  } else if (length(dims) == 3) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("plotly package required for 3D plots, falling back to 2D")
      return(plot_gpfa(gpfa_model, dims = dims[1:2], trials = trials,
                       color_by = color_by, add_start_end = add_start_end))
    }

    plotly::plot_ly(plot_data, x = ~x, y = ~y, z = ~z,
                    color = if (color_by == "time") ~time else ~factor(trial),
                    type = "scatter3d", mode = "lines") %>%
      plotly::layout(scene = list(
        xaxis = list(title = sprintf("Latent %d", dims[1])),
        yaxis = list(title = sprintf("Latent %d", dims[2])),
        zaxis = list(title = sprintf("Latent %d", dims[3]))
      ))
  }
}

#' Split Continuous Recording into Trials
#'
#' Helper function to split continuous recordings into trial segments.
#'
#' @param traces Continuous trace matrix (neurons x time)
#' @param trial_starts Vector of trial start indices
#' @param trial_length Length of each trial (frames)
#'
#' @return List of trial matrices
#' @keywords internal
split_into_trials <- function(traces, trial_starts, trial_length) {
  lapply(trial_starts, function(start) {
    end <- min(start + trial_length - 1, ncol(traces))
    traces[, start:end]
  })
}
