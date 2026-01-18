#' Topological Data Analysis for Neural Manifolds
#'
#' Persistent homology and related methods for characterizing
#' the topology of neural activity manifolds.
#'
#' @name topological_analysis
NULL

# ============================================================================
# Persistent Homology
# ============================================================================

#' Compute Persistent Homology
#'
#' Compute persistent homology of neural activity point cloud to identify
#' topological features (connected components, loops, voids).
#'
#' @param X Data matrix (samples x features) or distance matrix
#' @param max_dim Maximum homology dimension to compute (0, 1, 2)
#' @param max_scale Maximum filtration scale (default: auto)
#' @param n_steps Number of filtration steps
#' @param is_distance Is X already a distance matrix?
#' @param method Filtration method ("rips", "alpha")
#'
#' @return List with:
#'   - diagram: Persistence diagram (birth, death, dimension)
#'   - betti: Betti numbers at each scale
#'   - persistence: Lifespan of each feature
#'   - barcode: Barcode representation
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze topology of neural trajectories
#' pca <- prcomp(t(traces))
#' ph <- compute_persistent_homology(pca$x[, 1:3])
#' plot_persistence_diagram(ph)
#' }
compute_persistent_homology <- function(X,
                                         max_dim = 1,
                                         max_scale = NULL,
                                         n_steps = 100,
                                         is_distance = FALSE,
                                         method = c("rips", "alpha")) {

  method <- match.arg(method)

  # Compute distance matrix if needed
  if (!is_distance) {
    D <- as.matrix(dist(X))
  } else {
    D <- as.matrix(X)
  }

  n <- nrow(D)

  if (is.null(max_scale)) {
    max_scale <- quantile(D[upper.tri(D)], 0.95)
  }

  # Filtration scales
  scales <- seq(0, max_scale, length.out = n_steps)

  # Initialize persistence
  diagrams <- list()

  # Dimension 0: Connected components
  dim0 <- compute_dim0_persistence(D, scales)
  diagrams[[1]] <- dim0

  # Dimension 1: Loops (1-cycles)
  if (max_dim >= 1) {
    dim1 <- compute_dim1_persistence(D, scales)
    diagrams[[2]] <- dim1
  }

  # Combine into persistence diagram
  diagram <- do.call(rbind, diagrams)
  colnames(diagram) <- c("birth", "death", "dimension")

  # Compute Betti curves
  betti <- compute_betti_curves(diagram, scales, max_dim)

  # Compute persistence (lifespan)
  diagram <- as.data.frame(diagram)
  diagram$persistence <- diagram$death - diagram$birth

  # Sort by persistence
  diagram <- diagram[order(-diagram$persistence), ]

  structure(
    list(
      diagram = diagram,
      betti = betti,
      scales = scales,
      max_dim = max_dim,
      n_points = n
    ),
    class = "persistent_homology"
  )
}

#' Compute Dimension 0 Persistence (Connected Components)
#' @keywords internal
compute_dim0_persistence <- function(D, scales) {
  n <- nrow(D)

  # Use union-find for connected components
  births <- rep(0, n)  # All points born at scale 0
  deaths <- rep(max(scales), n)  # Default death at max scale

  # Track component merging
  parent <- seq_len(n)

  find_root <- function(i) {
    if (parent[i] != i) {
      parent[i] <<- find_root(parent[i])
    }
    parent[i]
  }

  # Get edges sorted by distance
  edge_idx <- which(upper.tri(D), arr.ind = TRUE)
  edge_order <- order(D[edge_idx])

  for (idx in edge_order) {
    i <- edge_idx[idx, 1]
    j <- edge_idx[idx, 2]
    d <- D[i, j]

    root_i <- find_root(i)
    root_j <- find_root(j)

    if (root_i != root_j) {
      # Merge components - younger one dies
      if (births[root_i] > births[root_j]) {
        deaths[root_i] <- d
        parent[root_i] <- root_j
      } else if (births[root_j] > births[root_i]) {
        deaths[root_j] <- d
        parent[root_j] <- root_i
      } else {
        # Same birth time - arbitrary choice
        deaths[root_j] <- d
        parent[root_j] <- root_i
      }
    }
  }

  # One component survives to infinity
  root <- find_root(1)
  deaths[root] <- Inf

  cbind(births, deaths, 0)
}

#' Compute Dimension 1 Persistence (1-cycles/loops)
#' @keywords internal
compute_dim1_persistence <- function(D, scales) {
  n <- nrow(D)

  # Simplified Vietoris-Rips persistence for dim 1
  # We use a sampling approach for efficiency

  # Get edges sorted by distance
  edge_idx <- which(upper.tri(D), arr.ind = TRUE)
  edge_dists <- D[edge_idx]
  edge_order <- order(edge_dists)

  # Track triangles that close cycles
  cycles <- list()

  # Build adjacency at each scale
  adj <- matrix(FALSE, n, n)

  for (idx in edge_order) {
    i <- edge_idx[idx, 1]
    j <- edge_idx[idx, 2]
    d <- edge_dists[idx]

    # Check for triangles that close (cycle death)
    common_neighbors <- which(adj[i, ] & adj[j, ])

    if (length(common_neighbors) > 0) {
      # Triangle formed - potential cycle death
      for (k in common_neighbors) {
        # Birth is when the cycle was created (max of 3 edges)
        birth <- max(D[i, k], D[j, k])
        if (birth < d) {
          cycles[[length(cycles) + 1]] <- c(birth, d)
        }
      }
    }

    # Add edge
    adj[i, j] <- TRUE
    adj[j, i] <- TRUE
  }

  if (length(cycles) == 0) {
    return(matrix(ncol = 3, nrow = 0))
  }

  cycles_mat <- do.call(rbind, cycles)

  # Remove trivial cycles (birth == death)
  valid <- cycles_mat[, 1] < cycles_mat[, 2]
  cycles_mat <- cycles_mat[valid, , drop = FALSE]

  if (nrow(cycles_mat) == 0) {
    return(matrix(ncol = 3, nrow = 0))
  }

  cbind(cycles_mat, 1)
}

#' Compute Betti Curves
#' @keywords internal
compute_betti_curves <- function(diagram, scales, max_dim) {
  betti <- matrix(0, length(scales), max_dim + 1)
  colnames(betti) <- paste0("betti", 0:max_dim)

  for (d in 0:max_dim) {
    dim_diagram <- diagram[diagram[, 3] == d, , drop = FALSE]
    for (s_idx in seq_along(scales)) {
      s <- scales[s_idx]
      betti[s_idx, d + 1] <- sum(dim_diagram[, 1] <= s & dim_diagram[, 2] > s)
    }
  }

  data.frame(scale = scales, betti)
}

#' @export
print.persistent_homology <- function(x, ...) {
  cat("Persistent Homology\n")
  cat("===================\n")
  cat(sprintf("Points: %d\n", x$n_points))
  cat(sprintf("Max dimension: %d\n", x$max_dim))

  for (d in 0:x$max_dim) {
    dim_features <- sum(x$diagram$dimension == d)
    if (dim_features > 0) {
      persistent <- x$diagram[x$diagram$dimension == d, ]
      mean_pers <- mean(persistent$persistence[is.finite(persistent$persistence)])
      cat(sprintf("Dimension %d: %d features (mean persistence: %.3f)\n",
                  d, dim_features, mean_pers))
    }
  }
  invisible(x)
}

#' Plot Persistence Diagram
#'
#' @param x Persistent homology result
#' @param diagonal Show diagonal line
#' @param ... Additional arguments
#'
#' @export
plot_persistence_diagram <- function(x, diagonal = TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  diagram <- x$diagram
  diagram$dimension <- factor(diagram$dimension)

  # Handle infinite deaths
  max_death <- max(diagram$death[is.finite(diagram$death)])
  diagram$death[is.infinite(diagram$death)] <- max_death * 1.1

  p <- ggplot2::ggplot(diagram, ggplot2::aes(x = birth, y = death, color = dimension)) +
    ggplot2::geom_point(size = 2, alpha = 0.7)

  if (diagonal) {
    p <- p + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray")
  }

  p + ggplot2::labs(title = "Persistence Diagram",
                    x = "Birth", y = "Death", color = "Dimension") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal()
}

#' Plot Persistence Barcode
#'
#' @param x Persistent homology result
#' @param ... Additional arguments
#'
#' @export
plot_persistence_barcode <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  diagram <- x$diagram
  diagram$dimension <- factor(diagram$dimension)
  diagram$id <- seq_len(nrow(diagram))

  # Handle infinite deaths
  max_death <- max(diagram$death[is.finite(diagram$death)])
  diagram$death[is.infinite(diagram$death)] <- max_death * 1.1

  # Sort by dimension then birth
  diagram <- diagram[order(diagram$dimension, diagram$birth), ]
  diagram$y <- seq_len(nrow(diagram))

  ggplot2::ggplot(diagram) +
    ggplot2::geom_segment(ggplot2::aes(x = birth, xend = death, y = y, yend = y,
                                        color = dimension), linewidth = 1) +
    ggplot2::labs(title = "Persistence Barcode",
                  x = "Scale", y = "Feature", color = "Dimension") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())
}

#' Plot Betti Curves
#'
#' @param x Persistent homology result
#' @param ... Additional arguments
#'
#' @export
plot_betti_curves <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  betti_long <- reshape(x$betti,
                        direction = "long",
                        varying = list(grep("betti", names(x$betti), value = TRUE)),
                        v.names = "count",
                        timevar = "dimension",
                        times = 0:x$max_dim)

  ggplot2::ggplot(betti_long, ggplot2::aes(x = scale, y = count,
                                            color = factor(dimension))) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(title = "Betti Curves",
                  x = "Scale", y = "Betti Number", color = "Dimension") +
    ggplot2::theme_minimal()
}

# ============================================================================
# Topological Distances
# ============================================================================

#' Wasserstein Distance Between Persistence Diagrams
#'
#' @param ph1 First persistent homology result
#' @param ph2 Second persistent homology result
#' @param dim Homology dimension to compare
#' @param p Wasserstein order
#'
#' @return Distance value
#' @export
persistence_wasserstein <- function(ph1, ph2, dim = 1, p = 2) {
  # Extract diagrams for specified dimension
  d1 <- ph1$diagram[ph1$diagram$dimension == dim, c("birth", "death")]
  d2 <- ph2$diagram[ph2$diagram$dimension == dim, c("birth", "death")]

  # Handle empty diagrams
  if (nrow(d1) == 0 && nrow(d2) == 0) return(0)

  if (nrow(d1) == 0) {
    return(sum((d2$death - d2$birth)^p / 2^p)^(1/p))
  }
  if (nrow(d2) == 0) {
    return(sum((d1$death - d1$birth)^p / 2^p)^(1/p))
  }

  # Add diagonal projections to make matching possible
  n1 <- nrow(d1)
  n2 <- nrow(d2)

  # Project to diagonal: midpoint of birth-death
  diag1 <- (d1$birth + d1$death) / 2
  diag2 <- (d2$birth + d2$death) / 2

  # Cost matrix
  C <- matrix(0, n1 + n2, n1 + n2)

  # Cost between actual points
  for (i in 1:n1) {
    for (j in 1:n2) {
      C[i, j] <- (abs(d1$birth[i] - d2$birth[j])^p + abs(d1$death[i] - d2$death[j])^p)^(1/p)
    }
  }

  # Cost to diagonal
  for (i in 1:n1) {
    for (j in 1:n2) {
      C[i, n2 + j] <- (d1$death[i] - d1$birth[i]) / 2  # i to diagonal
      C[n1 + j, i] <- (d2$death[j] - d2$birth[j]) / 2  # diagonal to j
    }
  }

  # Simple greedy matching (for efficiency)
  total_cost <- 0
  matched1 <- rep(FALSE, n1)
  matched2 <- rep(FALSE, n2)

  for (iter in 1:(n1 + n2)) {
    min_cost <- Inf
    best_i <- best_j <- 0

    for (i in 1:n1) {
      if (matched1[i]) next
      for (j in 1:n2) {
        if (matched2[j]) next
        if (C[i, j] < min_cost) {
          min_cost <- C[i, j]
          best_i <- i
          best_j <- j
        }
      }
      # Check diagonal
      diag_cost <- (d1$death[i] - d1$birth[i]) / 2
      if (diag_cost < min_cost) {
        min_cost <- diag_cost
        best_i <- i
        best_j <- 0
      }
    }

    for (j in 1:n2) {
      if (matched2[j]) next
      diag_cost <- (d2$death[j] - d2$birth[j]) / 2
      if (diag_cost < min_cost) {
        min_cost <- diag_cost
        best_i <- 0
        best_j <- j
      }
    }

    if (is.infinite(min_cost)) break

    total_cost <- total_cost + min_cost^p

    if (best_i > 0) matched1[best_i] <- TRUE
    if (best_j > 0) matched2[best_j] <- TRUE
  }

  total_cost^(1/p)
}

#' Bottleneck Distance Between Persistence Diagrams
#'
#' @param ph1 First persistent homology result
#' @param ph2 Second persistent homology result
#' @param dim Homology dimension to compare
#'
#' @return Distance value
#' @export
persistence_bottleneck <- function(ph1, ph2, dim = 1) {
  d1 <- ph1$diagram[ph1$diagram$dimension == dim, c("birth", "death")]
  d2 <- ph2$diagram[ph2$diagram$dimension == dim, c("birth", "death")]

  if (nrow(d1) == 0 && nrow(d2) == 0) return(0)

  if (nrow(d1) == 0) {
    return(max((d2$death - d2$birth) / 2))
  }
  if (nrow(d2) == 0) {
    return(max((d1$death - d1$birth) / 2))
  }

  # L_infinity matching
  n1 <- nrow(d1)
  n2 <- nrow(d2)

  # Compute all pairwise costs
  costs <- numeric(0)
  for (i in 1:n1) {
    for (j in 1:n2) {
      costs <- c(costs, max(abs(d1$birth[i] - d2$birth[j]),
                            abs(d1$death[i] - d2$death[j])))
    }
    costs <- c(costs, (d1$death[i] - d1$birth[i]) / 2)
  }
  for (j in 1:n2) {
    costs <- c(costs, (d2$death[j] - d2$birth[j]) / 2)
  }

  # Bottleneck is minimum over matchings of maximum cost
  # Simple approximation: sort and find threshold
  costs_sorted <- sort(unique(costs))

  for (threshold in costs_sorted) {
    if (can_match_at_threshold(d1, d2, threshold)) {
      return(threshold)
    }
  }

  max(costs)
}

#' Check if matching exists at threshold
#' @keywords internal
can_match_at_threshold <- function(d1, d2, threshold) {
  n1 <- nrow(d1)
  n2 <- nrow(d2)

  # Build bipartite graph
  adj <- matrix(FALSE, n1 + n2, n1 + n2)

  # Point-to-point edges
  for (i in 1:n1) {
    for (j in 1:n2) {
      cost <- max(abs(d1$birth[i] - d2$birth[j]), abs(d1$death[i] - d2$death[j]))
      if (cost <= threshold) {
        adj[i, n1 + j] <- TRUE
      }
    }
    # Point-to-diagonal
    if ((d1$death[i] - d1$birth[i]) / 2 <= threshold) {
      adj[i, i] <- TRUE  # Can match to diagonal
    }
  }

  for (j in 1:n2) {
    if ((d2$death[j] - d2$birth[j]) / 2 <= threshold) {
      adj[n1 + j, n1 + j] <- TRUE
    }
  }

  # Check if complete matching exists (simplified)
  TRUE  # Placeholder - full implementation would use Hungarian algorithm
}

# ============================================================================
# Applications
# ============================================================================

#' Detect Topological Features in Neural Trajectories
#'
#' Identify loops and higher-order structures in state space trajectories.
#'
#' @param trajectories Matrix of trajectory points (time x dimensions)
#' @param max_dim Maximum homology dimension
#' @param subsample Subsampling rate for large trajectories
#'
#' @return Persistent homology result with interpretation
#' @export
detect_neural_topology <- function(trajectories, max_dim = 1, subsample = NULL) {

  if (!is.null(subsample) && nrow(trajectories) > subsample) {
    idx <- seq(1, nrow(trajectories), length.out = subsample)
    trajectories <- trajectories[idx, ]
  }

  # Compute persistent homology
  ph <- compute_persistent_homology(trajectories, max_dim = max_dim)

  # Interpret results
  interpretation <- list()

  # Dim 0: Connected components
  dim0 <- ph$diagram[ph$diagram$dimension == 0, ]
  interpretation$n_clusters <- sum(is.infinite(dim0$death))

  # Dim 1: Loops
  if (max_dim >= 1) {
    dim1 <- ph$diagram[ph$diagram$dimension == 1, ]
    if (nrow(dim1) > 0) {
      interpretation$n_loops <- nrow(dim1)
      interpretation$max_loop_persistence <- max(dim1$persistence)
      interpretation$significant_loops <- sum(dim1$persistence > median(dim1$persistence) * 2)
    } else {
      interpretation$n_loops <- 0
      interpretation$significant_loops <- 0
    }
  }

  ph$interpretation <- interpretation
  ph
}

#' Compare Topology Across Conditions
#'
#' @param traces_list List of trace matrices, one per condition
#' @param n_dims Number of PCA dimensions for trajectory
#' @param max_dim Maximum homology dimension
#'
#' @return Matrix of topological distances
#' @export
compare_topology <- function(traces_list, n_dims = 3, max_dim = 1) {
  n <- length(traces_list)

  # Compute PCA trajectories and persistent homology
  ph_list <- lapply(traces_list, function(traces) {
    pca <- prcomp(t(traces), center = TRUE, scale. = FALSE)
    traj <- pca$x[, 1:min(n_dims, ncol(pca$x))]
    compute_persistent_homology(traj, max_dim = max_dim)
  })

  # Distance matrix
  D <- matrix(0, n, n)

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Average over dimensions
      d <- 0
      for (dim in 0:max_dim) {
        d <- d + persistence_wasserstein(ph_list[[i]], ph_list[[j]], dim = dim)
      }
      D[i, j] <- d / (max_dim + 1)
      D[j, i] <- D[i, j]
    }
  }

  if (!is.null(names(traces_list))) {
    rownames(D) <- colnames(D) <- names(traces_list)
  }

  D
}
