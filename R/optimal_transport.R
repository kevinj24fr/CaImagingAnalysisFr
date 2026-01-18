#' Optimal Transport Methods for Neural Data
#'
#' Wasserstein distance and Gromov-Wasserstein optimal transport
#' for comparing neural activity distributions across conditions.
#'
#' @name optimal_transport
NULL

# ============================================================================
# Wasserstein Distance
# ============================================================================

#' Compute Wasserstein Distance Between Distributions
#'
#' Earth Mover's Distance / Wasserstein-1 distance between activity distributions.
#' Superior to correlation for comparing neural representations.
#'
#' @param x First distribution (vector or matrix of samples)
#' @param y Second distribution (same size as x)
#' @param p Order of Wasserstein distance (1 or 2)
#' @param method Computation method ("exact", "sinkhorn", "sliced")
#' @param reg Regularization for Sinkhorn (entropy parameter)
#' @param n_projections Number of projections for sliced Wasserstein
#'
#' @return Wasserstein distance value
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare activity distributions between conditions
#' baseline_activity <- colMeans(traces[, 1:100])
#' stim_activity <- colMeans(traces[, 101:200])
#' dist <- wasserstein_distance(baseline_activity, stim_activity)
#' }
wasserstein_distance <- function(x, y,
                                  p = 1,
                                  method = c("exact", "sinkhorn", "sliced"),
                                  reg = 0.01,
                                  n_projections = 100) {
  method <- match.arg(method)

  # Ensure numeric vectors
  x <- as.numeric(x)
  y <- as.numeric(y)

  if (method == "exact") {
    wasserstein_exact(x, y, p)
  } else if (method == "sinkhorn") {
    wasserstein_sinkhorn(x, y, p, reg)
  } else {
    wasserstein_sliced(x, y, p, n_projections)
  }
}

#' Exact Wasserstein Distance (1D)
#' @keywords internal
wasserstein_exact <- function(x, y, p = 1) {
  n <- length(x)
  m <- length(y)

  # For 1D, exact solution is sorting-based
  x_sorted <- sort(x)
  y_sorted <- sort(y)

  if (n == m) {
    # Same size: direct comparison
    mean(abs(x_sorted - y_sorted)^p)^(1/p)
  } else {
    # Different sizes: interpolate to common grid
    n_grid <- max(n, m)
    x_interp <- approx(seq_along(x_sorted), x_sorted, n = n_grid)$y
    y_interp <- approx(seq_along(y_sorted), y_sorted, n = n_grid)$y
    mean(abs(x_interp - y_interp)^p)^(1/p)
  }
}

#' Sinkhorn-Knopp Algorithm for Entropic OT
#' @keywords internal
wasserstein_sinkhorn <- function(x, y, p = 2, reg = 0.01, n_iter = 100, tol = 1e-6) {
  n <- length(x)
  m <- length(y)

  # Cost matrix
  C <- outer(x, y, function(a, b) abs(a - b)^p)
  K <- exp(-C / reg)

  # Marginals (uniform)
  a <- rep(1/n, n)
  b <- rep(1/m, m)

  # Sinkhorn iterations
  u <- rep(1, n)
  v <- rep(1, m)

  for (iter in seq_len(n_iter)) {
    u_old <- u
    u <- a / (K %*% v)
    v <- b / (t(K) %*% u)

    if (max(abs(u - u_old)) < tol) break
  }

  # Transport plan
  P <- diag(u) %*% K %*% diag(v)

  # Wasserstein distance
  sum(P * C)^(1/p)
}

#' Sliced Wasserstein Distance (for high-D data)
#' @keywords internal
wasserstein_sliced <- function(x, y, p = 2, n_projections = 100) {
  if (is.vector(x)) x <- matrix(x, ncol = 1)
  if (is.vector(y)) y <- matrix(y, ncol = 1)

  d <- ncol(x)
  distances <- numeric(n_projections)

  for (i in seq_len(n_projections)) {
    # Random projection direction
    theta <- rnorm(d)
    theta <- theta / sqrt(sum(theta^2))

    # Project data
    x_proj <- x %*% theta
    y_proj <- y %*% theta

    # 1D Wasserstein
    distances[i] <- wasserstein_exact(x_proj, y_proj, p)^p
  }

  mean(distances)^(1/p)
}

#' Wasserstein Distance Matrix Between Multiple Conditions
#'
#' @param data_list List of matrices/vectors, one per condition
#' @param p Wasserstein order
#' @param method Computation method
#'
#' @return Symmetric distance matrix
#' @export
wasserstein_matrix <- function(data_list, p = 1, method = "exact") {
  n <- length(data_list)
  D <- matrix(0, n, n)

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      D[i, j] <- wasserstein_distance(data_list[[i]], data_list[[j]], p, method)
      D[j, i] <- D[i, j]
    }
  }

  if (!is.null(names(data_list))) {
    rownames(D) <- colnames(D) <- names(data_list)
  }

  D
}

# ============================================================================
# Gromov-Wasserstein Optimal Transport
# ============================================================================

#' Gromov-Wasserstein Distance
#'
#' Unsupervised alignment between neural representations via structure-preserving
#' optimal transport. Works across different neurons/animals by comparing
#' internal relationship structures.
#'
#' @param X First dataset (samples x features) or distance matrix
#' @param Y Second dataset (samples x features) or distance matrix
#' @param p Distance power for Gromov-Wasserstein
#' @param reg Entropic regularization
#' @param n_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param is_distance Are X and Y already distance matrices?
#'
#' @return List with:
#'   - distance: GW distance
#'   - coupling: Optimal transport plan
#'   - log: Convergence log
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Align neural representations across animals
#' gw <- gromov_wasserstein(animal1_activity, animal2_activity)
#' alignment <- gw$coupling
#' }
gromov_wasserstein <- function(X, Y,
                                p = 2,
                                reg = 0.01,
                                n_iter = 100,
                                tol = 1e-6,
                                is_distance = FALSE) {

  # Compute distance matrices if needed
  if (!is_distance) {
    C1 <- as.matrix(dist(X))
    C2 <- as.matrix(dist(Y))
  } else {
    C1 <- as.matrix(X)
    C2 <- as.matrix(Y)
  }

  n <- nrow(C1)
  m <- nrow(C2)

  # Normalize distance matrices
  C1 <- C1 / max(C1)
  C2 <- C2 / max(C2)

  # Uniform marginals
  a <- rep(1/n, n)
  b <- rep(1/m, m)

  # Initialize coupling as outer product of marginals
  T <- outer(a, b)

  log <- numeric(n_iter)

  for (iter in seq_len(n_iter)) {
    # Compute GW cost tensor (quadratic)
    # L(C1, C2, T) = sum_{i,j,k,l} |C1_{ik} - C2_{jl}|^p * T_{ij} * T_{kl}

    # Efficient computation via matrix operations
    constC1 <- C1^p %*% (T %*% matrix(1, m, m))
    constC2 <- matrix(1, n, n) %*% (T %*% t(C2^p))

    # Cross term
    hC1 <- C1^(p/2)
    hC2 <- C2^(p/2)
    cross <- hC1 %*% T %*% t(hC2)

    # Gradient (cost matrix for Sinkhorn)
    G <- constC1 + constC2 - 2 * cross

    # Sinkhorn step
    K <- exp(-G / reg)
    u <- rep(1, n)
    v <- rep(1, m)

    for (sink_iter in 1:50) {
      u <- a / (K %*% v + 1e-300)
      v <- b / (t(K) %*% u + 1e-300)
    }

    T_new <- diag(as.vector(u)) %*% K %*% diag(as.vector(v))

    # Compute GW distance
    gw_dist <- sum(G * T_new)
    log[iter] <- gw_dist

    # Check convergence
    if (iter > 1 && abs(log[iter] - log[iter-1]) < tol) {
      break
    }

    T <- T_new
  }

  structure(
    list(
      distance = gw_dist,
      coupling = T,
      log = log[1:iter],
      iterations = iter
    ),
    class = "gromov_wasserstein"
  )
}

#' @export
print.gromov_wasserstein <- function(x, ...) {
  cat("Gromov-Wasserstein Optimal Transport\n")
  cat("====================================\n")
  cat(sprintf("GW distance: %.4f\n", x$distance))
  cat(sprintf("Coupling dimensions: %d x %d\n", nrow(x$coupling), ncol(x$coupling)))
  cat(sprintf("Iterations: %d\n", x$iterations))
  invisible(x)
}

#' Wasserstein Barycenter
#'
#' Compute the Wasserstein barycenter (Frechet mean) of multiple distributions.
#'
#' @param distributions List of distribution samples
#' @param weights Optional weights for each distribution
#' @param n_support Size of barycenter support
#' @param n_iter Maximum iterations
#' @param reg Sinkhorn regularization
#'
#' @return Vector representing the barycenter distribution
#' @export
wasserstein_barycenter <- function(distributions,
                                    weights = NULL,
                                    n_support = NULL,
                                    n_iter = 100,
                                    reg = 0.01) {

  K <- length(distributions)

  if (is.null(weights)) {
    weights <- rep(1/K, K)
  }
  weights <- weights / sum(weights)

  if (is.null(n_support)) {
    n_support <- median(sapply(distributions, length))
  }

  # Initialize barycenter as weighted average of sorted distributions
  bary <- numeric(n_support)
  for (k in seq_len(K)) {
    interp <- approx(seq_along(distributions[[k]]),
                     sort(distributions[[k]]),
                     n = n_support)$y
    bary <- bary + weights[k] * interp
  }

  # Fixed-point iterations
  for (iter in seq_len(n_iter)) {
    bary_old <- bary

    # Update barycenter
    updates <- matrix(0, n_support, K)
    for (k in seq_len(K)) {
      # Compute transport from distribution k to barycenter
      d_k <- sort(distributions[[k]])
      d_k_interp <- approx(seq_along(d_k), d_k, n = n_support)$y
      updates[, k] <- d_k_interp
    }

    bary <- rowSums(sweep(updates, 2, weights, "*"))

    if (max(abs(bary - bary_old)) < 1e-6) break
  }

  bary
}

# ============================================================================
# Applications to Neural Data
# ============================================================================

#' Compare Neural Representations via Optimal Transport
#'
#' Compare population activity patterns between conditions using OT.
#'
#' @param traces1 Traces from condition 1 (cells x time)
#' @param traces2 Traces from condition 2 (cells x time)
#' @param method "wasserstein" or "gromov_wasserstein"
#' @param by Comparison level ("timepoint", "cell", "population")
#'
#' @return Distance value or distance time series
#' @export
compare_representations_ot <- function(traces1, traces2,
                                        method = c("wasserstein", "gromov_wasserstein"),
                                        by = c("timepoint", "cell", "population")) {
  method <- match.arg(method)
  by <- match.arg(by)

  if (by == "timepoint") {
    # Distance at each time point (across cells)
    n_time <- min(ncol(traces1), ncol(traces2))
    distances <- numeric(n_time)

    for (t in seq_len(n_time)) {
      if (method == "wasserstein") {
        distances[t] <- wasserstein_distance(traces1[, t], traces2[, t])
      } else {
        gw <- gromov_wasserstein(matrix(traces1[, t], ncol = 1),
                                  matrix(traces2[, t], ncol = 1))
        distances[t] <- gw$distance
      }
    }
    distances

  } else if (by == "cell") {
    # Distance for each cell (across time)
    n_cells <- min(nrow(traces1), nrow(traces2))
    distances <- numeric(n_cells)

    for (c in seq_len(n_cells)) {
      if (method == "wasserstein") {
        distances[c] <- wasserstein_distance(traces1[c, ], traces2[c, ])
      } else {
        gw <- gromov_wasserstein(matrix(traces1[c, ], ncol = 1),
                                  matrix(traces2[c, ], ncol = 1))
        distances[c] <- gw$distance
      }
    }
    distances

  } else {
    # Single population-level distance
    if (method == "wasserstein") {
      # Flatten and compare
      wasserstein_distance(as.vector(traces1), as.vector(traces2))
    } else {
      gw <- gromov_wasserstein(t(traces1), t(traces2))
      gw$distance
    }
  }
}

#' Optimal Transport Alignment Between Sessions
#'
#' Align neural activity across recording sessions using GW-OT.
#'
#' @param session1 Activity matrix from session 1 (cells x time)
#' @param session2 Activity matrix from session 2 (cells x time)
#' @param reg Regularization parameter
#'
#' @return List with:
#'   - alignment: Soft cell-to-cell correspondence
#'   - distance: GW distance
#'   - hard_assignment: Hard 1-to-1 mapping
#'
#' @export
align_sessions_ot <- function(session1, session2, reg = 0.01) {

  # Compute cell-cell similarity within each session
  sim1 <- cor(t(session1))
  sim2 <- cor(t(session2))

  # Handle NAs
  sim1[is.na(sim1)] <- 0
  sim2[is.na(sim2)] <- 0

  # Convert to distance
  D1 <- 1 - sim1
  D2 <- 1 - sim2

  # GW alignment
  gw <- gromov_wasserstein(D1, D2, is_distance = TRUE, reg = reg)

  # Hard assignment (Hungarian-like from coupling)
  coupling <- gw$coupling
  hard_assignment <- apply(coupling, 1, which.max)

  list(
    alignment = coupling,
    distance = gw$distance,
    hard_assignment = hard_assignment,
    confidence = apply(coupling, 1, max) / rowSums(coupling)
  )
}

#' Earth Mover's Distance for Place Cell Comparison
#'
#' Compare place field distributions using EMD.
#'
#' @param field1 First place field (1D or 2D)
#' @param field2 Second place field (same size)
#' @param positions Optional position coordinates
#'
#' @return EMD value
#' @export
place_field_emd <- function(field1, field2, positions = NULL) {

  # Normalize to probability distributions
  field1 <- field1 / sum(field1 + 1e-10)
  field2 <- field2 / sum(field2 + 1e-10)

  # Flatten if 2D
  if (is.matrix(field1)) {
    if (is.null(positions)) {
      nr <- nrow(field1)
      nc <- ncol(field1)
      positions <- expand.grid(x = 1:nc, y = 1:nr)
    }
    field1 <- as.vector(field1)
    field2 <- as.vector(field2)
  } else {
    if (is.null(positions)) {
      positions <- data.frame(x = seq_along(field1))
    }
  }

  # Cost matrix (Euclidean distance between positions)
  pos_mat <- as.matrix(positions)
  C <- as.matrix(dist(pos_mat))

  # EMD via Sinkhorn
  wasserstein_sinkhorn_weighted(field1, field2, C)
}

#' Weighted Sinkhorn for Non-uniform Marginals
#' @keywords internal
wasserstein_sinkhorn_weighted <- function(a, b, C, reg = 0.01, n_iter = 100) {
  n <- length(a)
  m <- length(b)

  K <- exp(-C / reg)

  u <- rep(1, n)
  v <- rep(1, m)

  for (iter in seq_len(n_iter)) {
    u <- a / (K %*% v + 1e-300)
    v <- b / (t(K) %*% u + 1e-300)
  }

  P <- diag(as.vector(u)) %*% K %*% diag(as.vector(v))
  sum(P * C)
}

#' Plot Optimal Transport Coupling
#'
#' @param coupling Transport coupling matrix from GW/Wasserstein
#' @param labels1 Labels for first distribution
#' @param labels2 Labels for second distribution
#'
#' @export
plot_ot_coupling <- function(coupling, labels1 = NULL, labels2 = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  n1 <- nrow(coupling)
  n2 <- ncol(coupling)

  if (is.null(labels1)) labels1 <- paste0("S1_", 1:n1)
  if (is.null(labels2)) labels2 <- paste0("S2_", 1:n2)

  df <- expand.grid(from = labels1, to = labels2)
  df$weight <- as.vector(coupling)

  ggplot2::ggplot(df, ggplot2::aes(x = from, y = to, fill = weight)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "plasma") +
    ggplot2::labs(title = "Optimal Transport Coupling",
                  x = "Source", y = "Target", fill = "Coupling") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
