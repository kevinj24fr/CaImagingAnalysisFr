#' State Space Analysis
#'
#' Population-level trajectory analysis using dimensionality reduction.
#' Visualize and quantify how neural population activity evolves over time.
#'
#' @name state_space
#' @keywords internal
NULL

#' Compute population trajectories
#'
#' Project neural population activity into low-dimensional state space.
#'
#' @param traces Matrix of traces (cells x time)
#' @param method Dimensionality reduction: "pca", "umap", "tsne", "fa"
#' @param n_dims Number of dimensions for projection
#' @param ... Additional arguments passed to reduction method
#'
#' @return List with trajectories and reduction details
#' @export
#'
#' @examples
#' \dontrun{
#' traj <- compute_trajectories(traces, method = "pca", n_dims = 3)
#' plot_trajectory(traj)
#' }
compute_trajectories <- function(traces, method = "pca", n_dims = 3, ...) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Transpose: time points as observations, cells as features
  data_matrix <- t(traces)

  # Z-score normalize
  data_scaled <- scale(data_matrix)

  # Apply dimensionality reduction
  result <- switch(method,
    pca = {
      pca <- prcomp(data_scaled, center = FALSE, scale. = FALSE)
      list(
        coords = pca$x[, 1:n_dims, drop = FALSE],
        loadings = pca$rotation[, 1:n_dims, drop = FALSE],
        variance_explained = pca$sdev^2 / sum(pca$sdev^2),
        model = pca
      )
    },
    umap = {
      if (!requireNamespace("uwot", quietly = TRUE)) {
        stop("uwot package required for UMAP")
      }
      umap_result <- uwot::umap(data_scaled, n_components = n_dims, ...)
      list(
        coords = umap_result,
        model = attr(umap_result, "model")
      )
    },
    tsne = {
      if (!requireNamespace("Rtsne", quietly = TRUE)) {
        stop("Rtsne package required for t-SNE")
      }
      tsne_result <- Rtsne::Rtsne(data_scaled, dims = n_dims, ...)
      list(
        coords = tsne_result$Y
      )
    },
    fa = {
      fa_result <- factanal(data_scaled, factors = n_dims, scores = "regression", ...)
      list(
        coords = fa_result$scores,
        loadings = fa_result$loadings[, 1:n_dims],
        model = fa_result
      )
    },
    stop("Unknown method: ", method)
  )

  # Add time index
  result$time <- seq_len(n_time)

  structure(
    c(result, list(
      method = method,
      n_dims = n_dims,
      n_cells = n_cells,
      n_time = n_time
    )),
    class = "state_trajectory"
  )
}

#' Compute trajectory velocity
#'
#' Calculate speed of movement through state space.
#'
#' @param trajectory State trajectory object
#'
#' @return Vector of velocities at each time point
#' @export
trajectory_velocity <- function(trajectory) {
  coords <- trajectory$coords
  n_time <- nrow(coords)

  # Euclidean distance between consecutive points
  diffs <- diff(coords)
  velocities <- sqrt(rowSums(diffs^2))

  c(NA, velocities)
}

#' Find fixed points in trajectory
#'
#' Identify stable states where the system tends to dwell.
#'
#' @param trajectory State trajectory object
#' @param n_clusters Number of fixed points to find
#' @param min_dwell Minimum frames to be considered dwelling
#'
#' @return List with fixed point locations and dwell times
#' @export
find_fixed_points <- function(trajectory, n_clusters = 5, min_dwell = 10) {
  coords <- trajectory$coords
  velocity <- trajectory_velocity(trajectory)

  # Cluster low-velocity points
  slow_idx <- which(velocity < quantile(velocity, 0.3, na.rm = TRUE))

  if (length(slow_idx) < n_clusters * 2) {
    n_clusters <- max(2, length(slow_idx) %/% 3)
  }

  if (length(slow_idx) < 3) {
    return(list(
      centers = matrix(nrow = 0, ncol = ncol(coords)),
      dwell_times = numeric(),
      assignments = integer()
    ))
  }

  km <- kmeans(coords[slow_idx, , drop = FALSE], centers = n_clusters, nstart = 10)

  # Compute dwell times for each cluster
  assignments <- rep(NA, nrow(coords))
  assignments[slow_idx] <- km$cluster

  dwell_times <- sapply(seq_len(n_clusters), function(k) {
    in_cluster <- which(assignments == k)
    if (length(in_cluster) == 0) return(0)

    # Find contiguous runs
    runs <- rle(diff(in_cluster) == 1)
    max(runs$lengths[runs$values], 0, na.rm = TRUE)
  })

  list(
    centers = km$centers,
    dwell_times = dwell_times,
    assignments = assignments,
    n_fixed_points = sum(dwell_times >= min_dwell)
  )
}

#' Compute trajectory curvature
#'
#' Measure how curved/straight the trajectory is.
#'
#' @param trajectory State trajectory object
#' @param window Smoothing window
#'
#' @return Vector of curvature values
#' @export
trajectory_curvature <- function(trajectory, window = 5) {
  coords <- trajectory$coords
  n_time <- nrow(coords)

  if (ncol(coords) < 2) {
    return(rep(NA, n_time))
  }

  # First and second derivatives
  d1 <- apply(coords, 2, function(x) {
    c(NA, diff(x))
  })

  d2 <- apply(d1, 2, function(x) {
    c(NA, diff(x))
  })

  # Curvature in 2D/3D
  if (ncol(coords) == 2) {
    curvature <- abs(d1[, 1] * d2[, 2] - d1[, 2] * d2[, 1]) /
      (sqrt(d1[, 1]^2 + d1[, 2]^2)^3 + 1e-10)
  } else {
    # Approximate curvature for higher dimensions
    speed <- sqrt(rowSums(d1^2, na.rm = TRUE))
    accel <- sqrt(rowSums(d2^2, na.rm = TRUE))
    curvature <- accel / (speed^2 + 1e-10)
  }

  # Smooth
  curvature <- stats::filter(curvature, rep(1/window, window), sides = 2)

  as.numeric(curvature)
}

#' Compare trajectories between conditions
#'
#' Quantify differences in state space trajectories.
#'
#' @param traj1 First trajectory
#' @param traj2 Second trajectory
#' @param method Comparison method: "frechet", "hausdorff", "dtw", "procrustes"
#'
#' @return Distance or similarity measure
#' @export
compare_trajectories <- function(traj1, traj2, method = "procrustes") {
  coords1 <- traj1$coords
  coords2 <- traj2$coords

  # Ensure same dimensionality
  n_dims <- min(ncol(coords1), ncol(coords2))
  coords1 <- coords1[, 1:n_dims, drop = FALSE]
  coords2 <- coords2[, 1:n_dims, drop = FALSE]

  result <- switch(method,
    procrustes = {
      # Procrustes analysis
      proc <- .procrustes_align(coords1, coords2)
      list(
        distance = proc$distance,
        rotation = proc$rotation,
        scale = proc$scale,
        aligned_coords = proc$aligned
      )
    },
    hausdorff = {
      # Hausdorff distance
      d <- .hausdorff_distance(coords1, coords2)
      list(distance = d)
    },
    dtw = {
      # Dynamic time warping
      if (!requireNamespace("dtw", quietly = TRUE)) {
        stop("dtw package required")
      }
      alignment <- dtw::dtw(coords1, coords2)
      list(
        distance = alignment$distance,
        normalized_distance = alignment$normalizedDistance,
        alignment = alignment
      )
    },
    frechet = {
      # Frechet distance (discrete approximation)
      d <- .frechet_distance(coords1, coords2)
      list(distance = d)
    }
  )

  c(result, list(method = method))
}

.procrustes_align <- function(X, Y) {
  n <- nrow(X)

  # Center
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  Y_centered <- scale(Y, center = TRUE, scale = FALSE)

  # SVD of cross-covariance
  svd_result <- svd(t(X_centered) %*% Y_centered)

  # Rotation
  R <- svd_result$v %*% t(svd_result$u)

  # Scale
  s <- sum(diag(t(X_centered) %*% Y_centered %*% R)) / sum(X_centered^2)

  # Transform Y
  Y_aligned <- s * Y_centered %*% R

  # Distance
  distance <- sqrt(sum((X_centered - Y_aligned)^2) / n)

  list(
    aligned = Y_aligned,
    rotation = R,
    scale = s,
    distance = distance
  )
}

.hausdorff_distance <- function(X, Y) {
  # Compute all pairwise distances
  dist_XY <- as.matrix(dist(rbind(X, Y)))
  n1 <- nrow(X)
  n2 <- nrow(Y)

  dist_XY <- dist_XY[1:n1, (n1 + 1):(n1 + n2)]

  # Hausdorff: max(max_x min_y d(x,y), max_y min_x d(x,y))
  max(max(apply(dist_XY, 1, min)), max(apply(dist_XY, 2, min)))
}

.frechet_distance <- function(X, Y) {
  n <- nrow(X)
  m <- nrow(Y)

  # Distance matrix
  ca <- matrix(Inf, n, m)

  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      d <- sqrt(sum((X[i, ] - Y[j, ])^2))
      if (i == 1 && j == 1) {
        ca[i, j] <- d
      } else if (i == 1) {
        ca[i, j] <- max(ca[i, j - 1], d)
      } else if (j == 1) {
        ca[i, j] <- max(ca[i - 1, j], d)
      } else {
        ca[i, j] <- max(min(ca[i - 1, j], ca[i, j - 1], ca[i - 1, j - 1]), d)
      }
    }
  }

  ca[n, m]
}

#' Plot state space trajectory
#'
#' @param trajectory State trajectory object
#' @param color_by What to color by: "time", "velocity", "condition", or vector
#' @param dims Which dimensions to plot (default: first 2 or 3)
#' @param show_points Show individual points
#'
#' @return ggplot or plotly object
#' @export
plot_trajectory <- function(trajectory, color_by = "time", dims = NULL,
                            show_points = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  coords <- trajectory$coords

  if (is.null(dims)) {
    dims <- seq_len(min(3, ncol(coords)))
  }

  df <- as.data.frame(coords[, dims, drop = FALSE])
  names(df) <- paste0("dim", dims)
  df$time <- seq_len(nrow(coords))

  # Set up color variable
  if (is.character(color_by) && length(color_by) == 1) {
    if (color_by == "time") {
      df$color <- df$time
    } else if (color_by == "velocity") {
      df$color <- trajectory_velocity(trajectory)
    }
  } else {
    df$color <- color_by
  }

  if (length(dims) == 2) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, color = color)) +
      ggplot2::geom_path(linewidth = 0.5) +
      ggplot2::scale_color_viridis_c() +
      ggplot2::labs(
        title = sprintf("State Space Trajectory (%s)", trajectory$method),
        x = sprintf("Dimension %d", dims[1]),
        y = sprintf("Dimension %d", dims[2])
      ) +
      ggplot2::theme_minimal()

    if (show_points) {
      p <- p + ggplot2::geom_point(size = 0.5)
    }

    # Mark start and end
    p <- p +
      ggplot2::geom_point(data = df[1, ], color = "green", size = 3, shape = 17) +
      ggplot2::geom_point(data = df[nrow(df), ], color = "red", size = 3, shape = 15)

  } else if (length(dims) >= 3 && requireNamespace("plotly", quietly = TRUE)) {
    p <- plotly::plot_ly(df, x = ~dim1, y = ~dim2, z = ~dim3,
                         color = ~color, colors = viridisLite::viridis(100),
                         type = "scatter3d", mode = "lines")
  } else {
    # Fall back to 2D
    p <- plot_trajectory(trajectory, color_by, dims = dims[1:2], show_points)
  }

  p
}

#' Compute neural dimensionality
#'
#' Estimate the intrinsic dimensionality of neural activity.
#'
#' @param traces Matrix of traces (cells x time)
#' @param method Method: "participation_ratio", "broken_stick", "eigenvalue"
#'
#' @return Estimated dimensionality
#' @export
estimate_dimensionality <- function(traces, method = "participation_ratio") {
  # PCA
  pca <- prcomp(t(traces), center = TRUE, scale. = TRUE)
  eigenvalues <- pca$sdev^2
  total_var <- sum(eigenvalues)

  if (method == "participation_ratio") {
    # Participation ratio: (sum(lambda))^2 / sum(lambda^2)
    pr <- sum(eigenvalues)^2 / sum(eigenvalues^2)
    dimensionality <- pr

  } else if (method == "broken_stick") {
    # Compare to broken stick distribution
    n <- length(eigenvalues)
    expected <- rev(cumsum(1 / (n:1))) / n * total_var

    dimensionality <- sum(eigenvalues > expected)

  } else if (method == "eigenvalue") {
    # Count eigenvalues > 1 (Kaiser criterion)
    dimensionality <- sum(eigenvalues > 1)

  } else if (method == "variance_threshold") {
    # Number of PCs for 95% variance
    cum_var <- cumsum(eigenvalues) / total_var
    dimensionality <- which(cum_var >= 0.95)[1]
  }

  list(
    dimensionality = dimensionality,
    method = method,
    eigenvalues = eigenvalues,
    variance_explained = eigenvalues / total_var
  )
}
