#' Neural Ordinary Differential Equations
#'
#' Continuous-time latent dynamics models using neural networks
#' to parameterize the derivative function.
#'
#' @name neural_ode
NULL

# ============================================================================
# Neural ODE Implementation
# ============================================================================

#' Fit Neural ODE to Neural Trajectories
#'
#' Learn continuous-time dynamics from discrete observations.
#' Uses a neural network to parameterize dx/dt = f_theta(x, t).
#'
#' @param trajectories Matrix of trajectories (time x dimensions)
#' @param time Time vector
#' @param latent_dim Dimension of latent state
#' @param hidden_dims Hidden layer dimensions for dynamics network
#' @param n_epochs Training epochs
#' @param learning_rate Learning rate
#' @param solver ODE solver ("euler", "rk4", "dopri5")
#' @param dt Integration time step
#' @param verbose Print training progress
#'
#' @return Neural ODE model object
#' @export
#'
#' @examples
#' \dontrun{
#' # Learn continuous dynamics from discrete observations
#' node <- fit_neural_ode(trajectories, time, latent_dim = 3)
#'
#' # Generate continuous trajectories
#' t_fine <- seq(min(time), max(time), length.out = 500)
#' pred <- predict_neural_ode(node, trajectories[1, ], t_fine)
#' }
fit_neural_ode <- function(trajectories,
                            time,
                            latent_dim = NULL,
                            hidden_dims = c(32, 32),
                            n_epochs = 100,
                            learning_rate = 0.01,
                            solver = c("euler", "rk4", "dopri5"),
                            dt = NULL,
                            verbose = FALSE) {

  solver <- match.arg(solver)

  if (!is.matrix(trajectories)) {
    trajectories <- matrix(trajectories, ncol = 1)
  }

  n_time <- nrow(trajectories)
  obs_dim <- ncol(trajectories)

  if (is.null(latent_dim)) {
    latent_dim <- obs_dim
  }

  if (is.null(dt)) {
    dt <- min(diff(time)) / 2
  }

  # Initialize neural network parameters
  # Simple MLP: input -> hidden -> ... -> output
  dims <- c(latent_dim, hidden_dims, latent_dim)
  n_layers <- length(dims) - 1

  params <- list()
  for (l in seq_len(n_layers)) {
    # Xavier initialization
    fan_in <- dims[l]
    fan_out <- dims[l + 1]
    scale <- sqrt(2 / (fan_in + fan_out))

    params[[paste0("W", l)]] <- matrix(rnorm(fan_in * fan_out, 0, scale), fan_in, fan_out)
    params[[paste0("b", l)]] <- rep(0, fan_out)
  }

  # Encoder (if latent_dim != obs_dim)
  if (latent_dim != obs_dim) {
    params$encoder_W <- matrix(rnorm(obs_dim * latent_dim, 0, 0.1), obs_dim, latent_dim)
    params$encoder_b <- rep(0, latent_dim)
    params$decoder_W <- matrix(rnorm(latent_dim * obs_dim, 0, 0.1), latent_dim, obs_dim)
    params$decoder_b <- rep(0, obs_dim)
  }

  # Training loop
  loss_history <- numeric(n_epochs)

  for (epoch in seq_len(n_epochs)) {
    # Forward pass
    if (latent_dim != obs_dim) {
      z0 <- sweep(trajectories[1, , drop = FALSE] %*% params$encoder_W, 2, params$encoder_b, "+")
    } else {
      z0 <- trajectories[1, , drop = FALSE]
    }

    # Solve ODE
    z_traj <- solve_ode(z0, time, params, dims, solver, dt)

    # Decode if necessary
    if (latent_dim != obs_dim) {
      pred <- sweep(z_traj %*% params$decoder_W, 2, params$decoder_b, "+")
    } else {
      pred <- z_traj
    }

    # Loss (MSE)
    loss <- mean((pred - trajectories)^2)
    loss_history[epoch] <- loss

    if (verbose && epoch %% 10 == 0) {
      message(sprintf("Epoch %d: loss = %.6f", epoch, loss))
    }

    # Backward pass (gradient computation via finite differences for simplicity)
    # In practice, would use automatic differentiation
    grads <- compute_node_gradients(trajectories, z_traj, time, params, dims,
                                     latent_dim, obs_dim, solver, dt)

    # Update parameters
    for (name in names(params)) {
      if (name %in% names(grads)) {
        params[[name]] <- params[[name]] - learning_rate * grads[[name]]
      }
    }

    # Learning rate decay
    if (epoch %% 50 == 0) {
      learning_rate <- learning_rate * 0.9
    }
  }

  structure(
    list(
      params = params,
      dims = dims,
      latent_dim = latent_dim,
      obs_dim = obs_dim,
      solver = solver,
      dt = dt,
      time = time,
      loss_history = loss_history
    ),
    class = "neural_ode"
  )
}

#' Solve ODE with Neural Network Dynamics
#' @keywords internal
solve_ode <- function(z0, time, params, dims, solver, dt) {
  n_time <- length(time)
  latent_dim <- ncol(z0)

  z_traj <- matrix(0, n_time, latent_dim)
  z_traj[1, ] <- z0

  z <- as.vector(z0)

  for (i in 2:n_time) {
    t_span <- time[i] - time[i - 1]
    n_steps <- ceiling(t_span / dt)
    h <- t_span / n_steps

    for (step in seq_len(n_steps)) {
      if (solver == "euler") {
        dz <- neural_dynamics(z, params, dims)
        z <- z + h * dz
      } else if (solver == "rk4") {
        k1 <- neural_dynamics(z, params, dims)
        k2 <- neural_dynamics(z + h/2 * k1, params, dims)
        k3 <- neural_dynamics(z + h/2 * k2, params, dims)
        k4 <- neural_dynamics(z + h * k3, params, dims)
        z <- z + h/6 * (k1 + 2*k2 + 2*k3 + k4)
      } else {
        # Adaptive dopri5 (simplified fixed step)
        k1 <- neural_dynamics(z, params, dims)
        k2 <- neural_dynamics(z + h/5 * k1, params, dims)
        k3 <- neural_dynamics(z + h * (3/40 * k1 + 9/40 * k2), params, dims)
        k4 <- neural_dynamics(z + h * (44/45 * k1 - 56/15 * k2 + 32/9 * k3), params, dims)
        k5 <- neural_dynamics(z + h * (19372/6561 * k1 - 25360/2187 * k2 + 64448/6561 * k3 - 212/729 * k4), params, dims)
        k6 <- neural_dynamics(z + h * (9017/3168 * k1 - 355/33 * k2 + 46732/5247 * k3 + 49/176 * k4 - 5103/18656 * k5), params, dims)
        z <- z + h * (35/384 * k1 + 500/1113 * k3 + 125/192 * k4 - 2187/6784 * k5 + 11/84 * k6)
      }
    }

    z_traj[i, ] <- z
  }

  z_traj
}

#' Neural Network Dynamics Function
#' @keywords internal
neural_dynamics <- function(z, params, dims) {
  n_layers <- length(dims) - 1
  x <- z

  for (l in seq_len(n_layers)) {
    W <- params[[paste0("W", l)]]
    b <- params[[paste0("b", l)]]

    x <- x %*% W + b

    # Activation (tanh for all but last layer)
    if (l < n_layers) {
      x <- tanh(x)
    }
  }

  as.vector(x)
}

#' Compute Gradients for Neural ODE
#' @keywords internal
compute_node_gradients <- function(trajectories, z_traj, time, params, dims,
                                    latent_dim, obs_dim, solver, dt) {
  eps <- 1e-5
  grads <- list()

  for (name in names(params)) {
    param <- params[[name]]
    grad <- array(0, dim = dim(param) %||% length(param))

    # Finite difference
    for (idx in seq_along(param)) {
      params_plus <- params
      params_plus[[name]][idx] <- param[idx] + eps

      z0 <- if (latent_dim != obs_dim) {
        sweep(trajectories[1, , drop = FALSE] %*% params_plus$encoder_W, 2, params_plus$encoder_b, "+")
      } else {
        trajectories[1, , drop = FALSE]
      }

      z_plus <- solve_ode(z0, time, params_plus, dims, solver, dt)

      if (latent_dim != obs_dim) {
        pred_plus <- sweep(z_plus %*% params_plus$decoder_W, 2, params_plus$decoder_b, "+")
      } else {
        pred_plus <- z_plus
      }

      loss_plus <- mean((pred_plus - trajectories)^2)

      params_minus <- params
      params_minus[[name]][idx] <- param[idx] - eps

      z0 <- if (latent_dim != obs_dim) {
        sweep(trajectories[1, , drop = FALSE] %*% params_minus$encoder_W, 2, params_minus$encoder_b, "+")
      } else {
        trajectories[1, , drop = FALSE]
      }

      z_minus <- solve_ode(z0, time, params_minus, dims, solver, dt)

      if (latent_dim != obs_dim) {
        pred_minus <- sweep(z_minus %*% params_minus$decoder_W, 2, params_minus$decoder_b, "+")
      } else {
        pred_minus <- z_minus
      }

      loss_minus <- mean((pred_minus - trajectories)^2)

      grad[idx] <- (loss_plus - loss_minus) / (2 * eps)
    }

    grads[[name]] <- grad
  }

  grads
}

#' Predict with Neural ODE
#'
#' @param object Neural ODE model
#' @param z0 Initial state
#' @param time Time points for prediction
#' @param ... Additional arguments
#'
#' @return Predicted trajectories
#' @export
predict_neural_ode <- function(object, z0, time, ...) {
  if (!inherits(object, "neural_ode")) {
    stop("object must be a neural_ode")
  }

  # Encode if necessary
  if (object$latent_dim != object$obs_dim) {
    z0_latent <- as.vector(z0 %*% object$params$encoder_W + object$params$encoder_b)
    z0 <- matrix(z0_latent, nrow = 1)
  } else {
    z0 <- matrix(z0, nrow = 1)
  }

  # Solve ODE
  z_traj <- solve_ode(z0, time, object$params, object$dims, object$solver, object$dt)

  # Decode if necessary
  if (object$latent_dim != object$obs_dim) {
    pred <- sweep(z_traj %*% object$params$decoder_W, 2, object$params$decoder_b, "+")
  } else {
    pred <- z_traj
  }

  pred
}

#' @export
print.neural_ode <- function(x, ...) {
  cat("Neural ODE Model\n")
  cat("================\n")
  cat(sprintf("Observation dimension: %d\n", x$obs_dim))
  cat(sprintf("Latent dimension: %d\n", x$latent_dim))
  cat(sprintf("Hidden layers: %s\n", paste(x$dims[2:(length(x$dims)-1)], collapse = " -> ")))
  cat(sprintf("Solver: %s (dt = %.4f)\n", x$solver, x$dt))
  cat(sprintf("Final loss: %.6f\n", tail(x$loss_history, 1)))
  invisible(x)
}

#' Plot Neural ODE Training
#'
#' @param x Neural ODE model
#' @param ... Additional arguments
#'
#' @export
plot.neural_ode <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  df <- data.frame(
    epoch = seq_along(x$loss_history),
    loss = x$loss_history
  )

  ggplot2::ggplot(df, ggplot2::aes(x = epoch, y = loss)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Neural ODE Training", x = "Epoch", y = "Loss (log scale)") +
    ggplot2::theme_minimal()
}

# ============================================================================
# Latent Neural ODE for Population Dynamics
# ============================================================================

#' Fit Latent Neural ODE for Population Activity
#'
#' Learn continuous latent dynamics from high-dimensional neural recordings.
#'
#' @param traces Neural traces (cells x time)
#' @param time Time vector
#' @param latent_dim Latent dimension
#' @param hidden_dims Hidden layer sizes
#' @param n_epochs Training epochs
#' @param verbose Print progress
#'
#' @return Latent Neural ODE model
#' @export
fit_latent_node <- function(traces, time, latent_dim = 3, hidden_dims = c(32, 32),
                             n_epochs = 100, verbose = FALSE) {

  # PCA preprocessing
  pca <- prcomp(t(traces), center = TRUE, scale. = FALSE)
  trajectories <- pca$x[, 1:min(latent_dim * 2, ncol(pca$x))]

  # Fit Neural ODE
  node <- fit_neural_ode(
    trajectories = trajectories,
    time = time,
    latent_dim = latent_dim,
    hidden_dims = hidden_dims,
    n_epochs = n_epochs,
    verbose = verbose
  )

  node$pca <- pca
  node$traces_dim <- nrow(traces)

  class(node) <- c("latent_node", class(node))
  node
}

#' Generate Trajectories from Latent Neural ODE
#'
#' @param object Latent Neural ODE
#' @param t_fine Fine time grid for continuous trajectories
#' @param initial_time Start time index (default: 1)
#'
#' @return List with latent and reconstructed trajectories
#' @export
generate_trajectories <- function(object, t_fine, initial_time = 1) {
  if (!inherits(object, "latent_node")) {
    stop("object must be a latent_node")
  }

  # Initial condition from PCA
  z0 <- object$pca$x[initial_time, 1:min(object$latent_dim * 2, ncol(object$pca$x))]

  # Predict latent trajectory
  latent_traj <- predict_neural_ode(object, z0, t_fine)

  # Reconstruct (inverse PCA)
  # Note: This is approximate since we only use top PCs
  reconstructed <- latent_traj %*% t(object$pca$rotation[, 1:ncol(latent_traj)])
  reconstructed <- sweep(reconstructed, 2, object$pca$center, "+")

  list(
    time = t_fine,
    latent = latent_traj,
    reconstructed = t(reconstructed)  # cells x time
  )
}

#' Plot Latent Neural ODE Phase Space
#'
#' @param object Latent Neural ODE
#' @param dims Dimensions to plot (default: 1:3 for 3D)
#'
#' @export
plot_node_phase_space <- function(object, dims = 1:3) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  # Generate fine trajectory
  t_fine <- seq(min(object$time), max(object$time), length.out = 500)
  traj <- generate_trajectories(object, t_fine)

  latent <- traj$latent
  dims <- dims[dims <= ncol(latent)]

  if (length(dims) >= 3) {
    # 3D plot (simplified 2D projection)
    df <- data.frame(
      x = latent[, dims[1]],
      y = latent[, dims[2]],
      z = latent[, dims[3]],
      time = t_fine
    )

    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = time)) +
      ggplot2::geom_path(linewidth = 0.8) +
      ggplot2::scale_color_viridis_c() +
      ggplot2::labs(title = "Neural ODE Latent Dynamics",
                    x = paste("Dim", dims[1]),
                    y = paste("Dim", dims[2])) +
      ggplot2::theme_minimal()
  } else {
    df <- data.frame(
      x = latent[, dims[1]],
      y = if (length(dims) > 1) latent[, dims[2]] else 0,
      time = t_fine
    )

    ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = time)) +
      ggplot2::geom_path(linewidth = 0.8) +
      ggplot2::scale_color_viridis_c() +
      ggplot2::theme_minimal()
  }
}
