#' State Space Models: HMM and Switching Linear Dynamical Systems
#'
#' Hidden Markov Models and Switching Linear Dynamical Systems for
#' identifying discrete neural/behavioral states with distinct dynamics.
#'
#' @name state_space_models
NULL

# ============================================================================
# Hidden Markov Models (HMM)
# ============================================================================

#' Fit Hidden Markov Model to Neural Data
#'
#' Fits an HMM to identify discrete latent states from continuous observations.
#' Uses the Baum-Welch (EM) algorithm for parameter estimation.
#'
#' @param traces Matrix of traces (cells x time) or vector for single cell
#' @param n_states Number of hidden states
#' @param n_iter Maximum EM iterations
#' @param tol Convergence tolerance
#' @param init_method Initialization method ("kmeans", "random", "uniform")
#' @param emission Type of emission distribution ("gaussian", "poisson")
#' @param verbose Print progress
#'
#' @return List with:
#'   - states: Most likely state sequence (Viterbi)
#'   - state_probs: Posterior state probabilities (n_states x time)
#'   - transition: Transition probability matrix
#'   - initial: Initial state distribution
#'   - emission_params: Emission distribution parameters
#'   - log_likelihood: Final log-likelihood
#'   - convergence: Convergence info
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit 3-state HMM to population activity
#' pop_activity <- colMeans(traces)
#' hmm <- fit_hmm(pop_activity, n_states = 3)
#' plot(pop_activity, col = hmm$states)
#' }
fit_hmm <- function(traces,
                    n_states = 3,
                    n_iter = 100,
                    tol = 1e-6,
                    init_method = c("kmeans", "random", "uniform"),
                    emission = c("gaussian", "poisson"),
                    verbose = FALSE) {

  init_method <- match.arg(init_method)
  emission <- match.arg(emission)

 # Handle matrix input - use first PC or mean
  if (is.matrix(traces)) {
    if (nrow(traces) > 1) {
      # Use first PC for dimensionality reduction
      pca <- prcomp(t(traces), center = TRUE, scale. = FALSE)
      obs <- pca$x[, 1]
    } else {
      obs <- as.numeric(traces)
    }
  } else {
    obs <- as.numeric(traces)
  }

  n_obs <- length(obs)

  # Initialize parameters
  init_params <- initialize_hmm(obs, n_states, init_method, emission)
  A <- init_params$A  # Transition matrix
  pi <- init_params$pi  # Initial distribution
  emit <- init_params$emit  # Emission parameters

  # EM algorithm (Baum-Welch)
  log_lik_history <- numeric(n_iter)

  for (iter in seq_len(n_iter)) {
    # E-step: Forward-Backward
    fb <- forward_backward(obs, A, pi, emit, emission)

    log_lik_history[iter] <- fb$log_likelihood

    if (verbose && iter %% 10 == 0) {
      message(sprintf("Iteration %d: log-likelihood = %.4f", iter, fb$log_likelihood))
    }

    # Check convergence
    if (iter > 1) {
      if (abs(log_lik_history[iter] - log_lik_history[iter - 1]) < tol) {
        if (verbose) message("Converged at iteration ", iter)
        break
      }
    }

    # M-step: Update parameters
    # Update transition matrix
    A_new <- matrix(0, n_states, n_states)
    for (i in seq_len(n_states)) {
      for (j in seq_len(n_states)) {
        A_new[i, j] <- sum(fb$xi[i, j, ], na.rm = TRUE)
      }
      A_new[i, ] <- A_new[i, ] / sum(A_new[i, ])
    }
    A <- A_new

    # Update initial distribution
    pi <- fb$gamma[, 1]
    pi <- pi / sum(pi)

    # Update emission parameters
    emit <- update_emission(obs, fb$gamma, emission)
  }

  # Viterbi decoding for most likely state sequence
  states <- viterbi_decode(obs, A, pi, emit, emission)

  # Compute state occupancy and dwell times
  state_stats <- compute_state_statistics(states, n_states)

  structure(
    list(
      states = states,
      state_probs = fb$gamma,
      transition = A,
      initial = pi,
      emission_params = emit,
      emission_type = emission,
      log_likelihood = log_lik_history[iter],
      log_lik_history = log_lik_history[1:iter],
      n_states = n_states,
      state_occupancy = state_stats$occupancy,
      dwell_times = state_stats$dwell_times,
      n_transitions = state_stats$n_transitions
    ),
    class = "hmm_fit"
  )
}

#' Initialize HMM Parameters
#' @keywords internal
initialize_hmm <- function(obs, n_states, method, emission) {
  n_obs <- length(obs)

  # Initialize transition matrix (slightly diagonal-heavy for stability)
  A <- matrix(1 / n_states, n_states, n_states)
  diag(A) <- 0.7 + 0.3 / n_states
  A <- A / rowSums(A)

  # Initial distribution
  pi <- rep(1 / n_states, n_states)

  # Initialize emission parameters
  if (method == "kmeans") {
    km <- kmeans(obs, centers = n_states, nstart = 10)
    state_assign <- km$cluster

    if (emission == "gaussian") {
      means <- tapply(obs, state_assign, mean)
      sds <- tapply(obs, state_assign, sd)
      sds[sds < 0.01 * sd(obs)] <- 0.01 * sd(obs)  # Floor on variance
      emit <- list(means = as.numeric(means), sds = as.numeric(sds))
    } else {
      # Poisson
      lambdas <- pmax(0.1, tapply(obs, state_assign, mean))
      emit <- list(lambdas = as.numeric(lambdas))
    }
  } else if (method == "random") {
    if (emission == "gaussian") {
      means <- sort(runif(n_states, min(obs), max(obs)))
      sds <- rep(sd(obs) / sqrt(n_states), n_states)
      emit <- list(means = means, sds = sds)
    } else {
      lambdas <- sort(runif(n_states, 0.1, max(obs)))
      emit <- list(lambdas = lambdas)
    }
  } else {
    # Uniform initialization
    quantiles <- quantile(obs, probs = seq(0, 1, length.out = n_states + 1))
    if (emission == "gaussian") {
      means <- (quantiles[-1] + quantiles[-(n_states + 1)]) / 2
      sds <- rep(sd(obs) / sqrt(n_states), n_states)
      emit <- list(means = means, sds = sds)
    } else {
      lambdas <- pmax(0.1, (quantiles[-1] + quantiles[-(n_states + 1)]) / 2)
      emit <- list(lambdas = lambdas)
    }
  }

  list(A = A, pi = pi, emit = emit)
}

#' Forward-Backward Algorithm
#' @keywords internal
forward_backward <- function(obs, A, pi, emit, emission) {
  n_obs <- length(obs)
  n_states <- length(pi)

  # Compute emission probabilities
  B <- compute_emission_probs(obs, emit, emission)

  # Forward pass (scaled for numerical stability)
  alpha <- matrix(0, n_states, n_obs)
  scale <- numeric(n_obs)

  alpha[, 1] <- pi * B[, 1]
  scale[1] <- sum(alpha[, 1])
  alpha[, 1] <- alpha[, 1] / scale[1]

  for (t in 2:n_obs) {
    for (j in seq_len(n_states)) {
      alpha[j, t] <- sum(alpha[, t - 1] * A[, j]) * B[j, t]
    }
    scale[t] <- sum(alpha[, t])
    if (scale[t] > 0) {
      alpha[, t] <- alpha[, t] / scale[t]
    }
  }

  # Backward pass
  beta <- matrix(0, n_states, n_obs)
  beta[, n_obs] <- 1

  for (t in (n_obs - 1):1) {
    for (i in seq_len(n_states)) {
      beta[i, t] <- sum(A[i, ] * B[, t + 1] * beta[, t + 1])
    }
    if (scale[t + 1] > 0) {
      beta[, t] <- beta[, t] / scale[t + 1]
    }
  }

  # Posterior state probabilities (gamma)
  gamma <- alpha * beta
  gamma <- gamma / (colSums(gamma) + 1e-300)

  # Transition posteriors (xi)
  xi <- array(0, dim = c(n_states, n_states, n_obs - 1))
  for (t in 1:(n_obs - 1)) {
    denom <- sum(alpha[, t] %o% (A * outer(rep(1, n_states), B[, t + 1] * beta[, t + 1])))
    if (denom > 0) {
      for (i in seq_len(n_states)) {
        for (j in seq_len(n_states)) {
          xi[i, j, t] <- alpha[i, t] * A[i, j] * B[j, t + 1] * beta[j, t + 1] / denom
        }
      }
    }
  }

  # Log-likelihood
  log_likelihood <- sum(log(scale + 1e-300))

  list(
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    xi = xi,
    scale = scale,
    log_likelihood = log_likelihood
  )
}

#' Compute Emission Probabilities
#' @keywords internal
compute_emission_probs <- function(obs, emit, emission) {
  n_states <- if (emission == "gaussian") length(emit$means) else length(emit$lambdas)
  n_obs <- length(obs)

  B <- matrix(0, n_states, n_obs)

  if (emission == "gaussian") {
    for (k in seq_len(n_states)) {
      B[k, ] <- dnorm(obs, mean = emit$means[k], sd = emit$sds[k])
    }
  } else {
    # Poisson
    for (k in seq_len(n_states)) {
      B[k, ] <- dpois(round(pmax(0, obs)), lambda = emit$lambdas[k])
    }
  }

  # Floor for numerical stability
  B[B < 1e-300] <- 1e-300

  B
}

#' Update Emission Parameters
#' @keywords internal
update_emission <- function(obs, gamma, emission) {
  n_states <- nrow(gamma)

  if (emission == "gaussian") {
    means <- numeric(n_states)
    sds <- numeric(n_states)

    for (k in seq_len(n_states)) {
      weights <- gamma[k, ]
      weights <- weights / sum(weights)
      means[k] <- sum(weights * obs)
      sds[k] <- sqrt(sum(weights * (obs - means[k])^2))
    }

    sds[sds < 0.01 * sd(obs)] <- 0.01 * sd(obs)

    list(means = means, sds = sds)
  } else {
    lambdas <- numeric(n_states)
    for (k in seq_len(n_states)) {
      weights <- gamma[k, ]
      lambdas[k] <- sum(weights * obs) / sum(weights)
    }
    lambdas <- pmax(0.1, lambdas)
    list(lambdas = lambdas)
  }
}

#' Viterbi Decoding
#' @keywords internal
viterbi_decode <- function(obs, A, pi, emit, emission) {
  n_obs <- length(obs)
  n_states <- length(pi)

  B <- compute_emission_probs(obs, emit, emission)

  # Log scale for numerical stability
  log_A <- log(A + 1e-300)
  log_pi <- log(pi + 1e-300)
  log_B <- log(B)

  # Viterbi algorithm
  delta <- matrix(-Inf, n_states, n_obs)
  psi <- matrix(0L, n_states, n_obs)

  delta[, 1] <- log_pi + log_B[, 1]

  for (t in 2:n_obs) {
    for (j in seq_len(n_states)) {
      temp <- delta[, t - 1] + log_A[, j]
      psi[j, t] <- which.max(temp)
      delta[j, t] <- max(temp) + log_B[j, t]
    }
  }

  # Backtrack
  states <- integer(n_obs)
  states[n_obs] <- which.max(delta[, n_obs])

  for (t in (n_obs - 1):1) {
    states[t] <- psi[states[t + 1], t + 1]
  }

  states
}

#' Compute State Statistics
#' @keywords internal
compute_state_statistics <- function(states, n_states) {
  n_obs <- length(states)

  # Occupancy
  occupancy <- table(factor(states, levels = 1:n_states)) / n_obs

  # Dwell times
  rle_states <- rle(states)
  dwell_times <- lapply(1:n_states, function(k) {
    rle_states$lengths[rle_states$values == k]
  })

  # Number of transitions
  n_transitions <- sum(diff(states) != 0)

  list(
    occupancy = as.numeric(occupancy),
    dwell_times = dwell_times,
    n_transitions = n_transitions
  )
}

#' @export
print.hmm_fit <- function(x, ...) {
  cat("Hidden Markov Model Fit\n")
  cat("=======================\n")
  cat(sprintf("States: %d\n", x$n_states))
  cat(sprintf("Emission: %s\n", x$emission_type))
  cat(sprintf("Log-likelihood: %.2f\n", x$log_likelihood))
  cat(sprintf("Transitions: %d\n", x$n_transitions))
  cat("\nState occupancy:\n")
  for (k in seq_len(x$n_states)) {
    cat(sprintf("  State %d: %.1f%%\n", k, 100 * x$state_occupancy[k]))
  }
  invisible(x)
}

# ============================================================================
# Switching Linear Dynamical Systems (SLDS)
# ============================================================================

#' Fit Switching Linear Dynamical System
#'
#' Fits an SLDS where discrete states govern different linear dynamics.
#' Combines HMM for discrete states with LDS for continuous latent dynamics.
#'
#' @param traces Matrix of traces (cells x time)
#' @param n_states Number of discrete states
#' @param latent_dim Dimension of continuous latent state
#' @param n_iter Maximum EM iterations
#' @param tol Convergence tolerance
#' @param verbose Print progress
#'
#' @return List with:
#'   - discrete_states: Inferred discrete state sequence
#'   - latent_states: Continuous latent trajectories (latent_dim x time)
#'   - dynamics: List of dynamics matrices (A) for each state
#'   - observation: Observation matrix (C)
#'   - transition: Discrete state transition matrix
#'   - state_noise: Process noise covariance per state
#'   - obs_noise: Observation noise covariance
#'
#' @export
#'
#' @examples
#' \dontrun{
#' slds <- fit_slds(traces, n_states = 2, latent_dim = 3)
#' plot(slds$latent_states[1, ], col = slds$discrete_states)
#' }
fit_slds <- function(traces,
                     n_states = 2,
                     latent_dim = 3,
                     n_iter = 50,
                     tol = 1e-4,
                     verbose = FALSE) {

  if (!is.matrix(traces)) {
    traces <- matrix(traces, nrow = 1)
  }

  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Initialize with PCA for latent states
  pca <- prcomp(t(traces), center = TRUE, scale. = FALSE)
  z_init <- t(pca$x[, 1:min(latent_dim, ncol(pca$x))])
  if (nrow(z_init) < latent_dim) {
    z_init <- rbind(z_init, matrix(0, latent_dim - nrow(z_init), n_time))
  }

  # Initialize discrete states with k-means on latent
  km <- kmeans(t(z_init), centers = n_states, nstart = 10)
  s_init <- km$cluster

  # Initialize parameters
  # Observation matrix
  C <- pca$rotation[, 1:latent_dim]

  # Mean observation
  d <- pca$center

  # Observation noise
  reconstructed <- C %*% z_init
  residuals <- traces - sweep(reconstructed, 1, d, "+")
  R <- diag(pmax(0.01, apply(residuals, 1, var)))

  # Dynamics matrices (one per discrete state) - start near identity
  A_list <- lapply(1:n_states, function(k) {
    0.95 * diag(latent_dim) + 0.05 * matrix(rnorm(latent_dim^2, 0, 0.1), latent_dim, latent_dim)
  })

  # State noise covariance per state
  Q_list <- lapply(1:n_states, function(k) {
    0.1 * diag(latent_dim)
  })

  # Transition matrix for discrete states
  Pi <- matrix(0.1 / (n_states - 1), n_states, n_states)
  diag(Pi) <- 0.9

  # Initial discrete state distribution
  pi0 <- rep(1 / n_states, n_states)

  # EM iterations
  z <- z_init
  s <- s_init
  log_lik_history <- numeric(n_iter)

  for (iter in seq_len(n_iter)) {
    # E-step: Infer latent states given discrete states
    z_new <- slds_infer_latent(traces, s, A_list, C, d, Q_list, R)

    # E-step: Infer discrete states given latent states
    s_probs <- slds_infer_discrete(z_new, A_list, Q_list, Pi, pi0)
    s_new <- apply(s_probs, 2, which.max)

    # Compute log-likelihood
    log_lik <- slds_log_likelihood(traces, z_new, s_new, A_list, C, d, Q_list, R)
    log_lik_history[iter] <- log_lik

    if (verbose && iter %% 5 == 0) {
      message(sprintf("Iteration %d: log-likelihood = %.2f", iter, log_lik))
    }

    # Check convergence
    if (iter > 1 && abs(log_lik_history[iter] - log_lik_history[iter - 1]) < tol) {
      if (verbose) message("Converged at iteration ", iter)
      break
    }

    # M-step: Update dynamics matrices
    for (k in 1:n_states) {
      idx <- which(s_new[-1] == k)
      if (length(idx) >= latent_dim) {
        Z_prev <- z_new[, idx]
        Z_curr <- z_new[, idx + 1]
        # Ridge regression for stability
        lambda <- 0.01
        A_list[[k]] <- Z_curr %*% t(Z_prev) %*%
          solve(Z_prev %*% t(Z_prev) + lambda * diag(latent_dim))

        # Update state noise
        resid <- Z_curr - A_list[[k]] %*% Z_prev
        Q_list[[k]] <- resid %*% t(resid) / length(idx)
        Q_list[[k]] <- (Q_list[[k]] + t(Q_list[[k]])) / 2 + 0.001 * diag(latent_dim)
      }
    }

    # M-step: Update observation matrix
    C <- traces %*% t(z_new) %*% solve(z_new %*% t(z_new) + 0.01 * diag(latent_dim))

    # M-step: Update observation noise
    resid <- traces - C %*% z_new
    R <- diag(pmax(0.01, rowMeans(resid^2)))

    # M-step: Update transition matrix
    for (i in 1:n_states) {
      for (j in 1:n_states) {
        Pi[i, j] <- sum(s_new[-1] == j & s_new[-n_time] == i)
      }
      if (sum(Pi[i, ]) > 0) Pi[i, ] <- Pi[i, ] / sum(Pi[i, ])
    }

    z <- z_new
    s <- s_new
  }

  # Compute state statistics
  state_stats <- compute_state_statistics(s, n_states)

  structure(
    list(
      discrete_states = s,
      state_probs = s_probs,
      latent_states = z,
      dynamics = A_list,
      observation = C,
      obs_mean = d,
      transition = Pi,
      state_noise = Q_list,
      obs_noise = R,
      n_states = n_states,
      latent_dim = latent_dim,
      log_likelihood = log_lik_history[iter],
      log_lik_history = log_lik_history[1:iter],
      state_occupancy = state_stats$occupancy,
      dwell_times = state_stats$dwell_times
    ),
    class = "slds_fit"
  )
}

#' Infer Latent States for SLDS (Kalman-like)
#' @keywords internal
slds_infer_latent <- function(y, s, A_list, C, d, Q_list, R) {
  n_cells <- nrow(y)
  n_time <- ncol(y)
  latent_dim <- ncol(C)

  z <- matrix(0, latent_dim, n_time)
  P <- array(0, dim = c(latent_dim, latent_dim, n_time))

  # Initialize
  z[, 1] <- solve(t(C) %*% solve(R) %*% C + diag(latent_dim)) %*% t(C) %*% solve(R) %*% (y[, 1] - d)
  P[, , 1] <- diag(latent_dim)

  # Forward pass (Kalman filter)
  for (t in 2:n_time) {
    k <- s[t]
    A <- A_list[[k]]
    Q <- Q_list[[k]]

    # Predict
    z_pred <- A %*% z[, t - 1]
    P_pred <- A %*% P[, , t - 1] %*% t(A) + Q

    # Update
    y_pred <- C %*% z_pred + d
    S <- C %*% P_pred %*% t(C) + R
    K <- P_pred %*% t(C) %*% solve(S)

    z[, t] <- z_pred + K %*% (y[, t] - y_pred)
    P[, , t] <- (diag(latent_dim) - K %*% C) %*% P_pred
  }

  z
}

#' Infer Discrete States for SLDS
#' @keywords internal
slds_infer_discrete <- function(z, A_list, Q_list, Pi, pi0) {
  n_time <- ncol(z)
  n_states <- length(A_list)
  latent_dim <- nrow(z)

  # Compute emission likelihoods (how well each state's dynamics explains transitions)
  log_emit <- matrix(0, n_states, n_time)

  for (t in 2:n_time) {
    z_curr <- z[, t]
    z_prev <- z[, t - 1]

    for (k in 1:n_states) {
      z_pred <- A_list[[k]] %*% z_prev
      resid <- z_curr - z_pred
      # Multivariate normal log-likelihood
      log_emit[k, t] <- -0.5 * (t(resid) %*% solve(Q_list[[k]]) %*% resid +
                                   log(det(Q_list[[k]])) + latent_dim * log(2 * pi))
    }
  }
  log_emit[, 1] <- 0  # No dynamics constraint for first time point

  # Forward-backward on discrete states
  log_alpha <- matrix(-Inf, n_states, n_time)
  log_alpha[, 1] <- log(pi0 + 1e-300) + log_emit[, 1]

  for (t in 2:n_time) {
    for (j in 1:n_states) {
      log_alpha[j, t] <- log_sum_exp(log_alpha[, t - 1] + log(Pi[, j] + 1e-300)) + log_emit[j, t]
    }
  }

  # Backward
  log_beta <- matrix(0, n_states, n_time)

  for (t in (n_time - 1):1) {
    for (i in 1:n_states) {
      log_beta[i, t] <- log_sum_exp(log(Pi[i, ] + 1e-300) + log_emit[, t + 1] + log_beta[, t + 1])
    }
  }

  # Posterior
  log_gamma <- log_alpha + log_beta
  log_gamma <- log_gamma - apply(log_gamma, 2, log_sum_exp)

  exp(log_gamma)
}

#' Log-Sum-Exp for Numerical Stability
#' @keywords internal
log_sum_exp <- function(x) {
  max_x <- max(x)
  if (is.infinite(max_x)) return(max_x)
  max_x + log(sum(exp(x - max_x)))
}

#' SLDS Log-Likelihood
#' @keywords internal
slds_log_likelihood <- function(y, z, s, A_list, C, d, Q_list, R) {
  n_time <- ncol(y)

  # Observation likelihood
  resid_obs <- y - C %*% z - d
  log_lik_obs <- -0.5 * sum(diag(t(resid_obs) %*% solve(R) %*% resid_obs))

  # Dynamics likelihood
  log_lik_dyn <- 0
  for (t in 2:n_time) {
    k <- s[t]
    z_pred <- A_list[[k]] %*% z[, t - 1]
    resid <- z[, t] - z_pred
    log_lik_dyn <- log_lik_dyn - 0.5 * t(resid) %*% solve(Q_list[[k]]) %*% resid
  }

  as.numeric(log_lik_obs + log_lik_dyn)
}

#' @export
print.slds_fit <- function(x, ...) {
  cat("Switching Linear Dynamical System Fit\n")
  cat("=====================================\n")
  cat(sprintf("Discrete states: %d\n", x$n_states))
  cat(sprintf("Latent dimension: %d\n", x$latent_dim))
  cat(sprintf("Log-likelihood: %.2f\n", x$log_likelihood))
  cat("\nState occupancy:\n")
  for (k in seq_len(x$n_states)) {
    cat(sprintf("  State %d: %.1f%%\n", k, 100 * x$state_occupancy[k]))
  }
  cat("\nDynamics eigenvalues per state:\n")
  for (k in seq_len(x$n_states)) {
    eig <- eigen(x$dynamics[[k]])$values
    cat(sprintf("  State %d: max |eig| = %.3f\n", k, max(abs(eig))))
  }
  invisible(x)
}

#' Plot HMM or SLDS Results
#'
#' @param x HMM or SLDS fit object
#' @param trace Optional trace to overlay
#' @param ... Additional arguments
#'
#' @export
plot_state_model <- function(x, trace = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required for plotting")
  }

  if (inherits(x, "hmm_fit")) {
    n_time <- length(x$states)
    df <- data.frame(
      time = 1:n_time,
      state = factor(x$states)
    )

    if (!is.null(trace)) {
      df$value <- if (is.matrix(trace)) colMeans(trace) else trace

      p <- ggplot2::ggplot(df, ggplot2::aes(x = time)) +
        ggplot2::geom_line(ggplot2::aes(y = value), alpha = 0.7) +
        ggplot2::geom_point(ggplot2::aes(y = value, color = state), size = 0.5) +
        ggplot2::labs(title = "HMM State Sequence", x = "Time", y = "Activity") +
        ggplot2::theme_minimal()
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = 1, fill = state)) +
        ggplot2::geom_tile() +
        ggplot2::labs(title = "HMM State Sequence", x = "Time", y = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.y = ggplot2::element_blank())
    }
  } else if (inherits(x, "slds_fit")) {
    n_time <- ncol(x$latent_states)

    df <- data.frame(
      time = rep(1:n_time, x$latent_dim),
      dim = rep(1:x$latent_dim, each = n_time),
      value = as.vector(t(x$latent_states)),
      state = factor(rep(x$discrete_states, x$latent_dim))
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = value, color = state)) +
      ggplot2::geom_line(alpha = 0.8) +
      ggplot2::facet_wrap(~paste("Latent", dim), scales = "free_y", ncol = 1) +
      ggplot2::labs(title = "SLDS Latent Trajectories", x = "Time", y = "Latent State") +
      ggplot2::theme_minimal()
  }

  p
}

#' Predict States for New Data
#'
#' @param object Fitted HMM or SLDS object
#' @param newdata New observations
#' @param ... Additional arguments
#'
#' @return Predicted state sequence
#' @export
predict_states <- function(object, newdata, ...) {
  UseMethod("predict_states")
}

#' @export
predict_states.hmm_fit <- function(object, newdata, ...) {
  if (is.matrix(newdata) && nrow(newdata) > 1) {
    pca <- prcomp(t(newdata), center = TRUE, scale. = FALSE)
    obs <- pca$x[, 1]
  } else {
    obs <- as.numeric(newdata)
  }

  viterbi_decode(obs, object$transition, object$initial,
                 object$emission_params, object$emission_type)
}
