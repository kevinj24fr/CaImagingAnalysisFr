#' Neural Assembly Detection
#'
#' Methods for detecting recurring co-activation patterns (neural assemblies)
#' in calcium imaging data. Identifies groups of neurons that fire together.
#'
#' @name neural_assemblies
#' @keywords internal
NULL

#' Detect neural assemblies using PCA/ICA
#'
#' Identifies cell assemblies as groups of neurons with correlated activity
#' using dimensionality reduction methods.
#'
#' @param traces Matrix of traces (cells x time) or spike matrix
#' @param method Detection method: "pca", "ica", "nmf", "svd"
#' @param n_assemblies Number of assemblies to detect (NULL for automatic)
#' @param threshold Correlation threshold for assembly membership
#' @param min_cells Minimum cells per assembly
#'
#' @return List with assembly memberships, weights, and activation time courses
#' @export
#'
#' @examples
#' \dontrun{
#' assemblies <- detect_assemblies(traces, method = "ica", n_assemblies = 10)
#' plot_assemblies(assemblies)
#' }
detect_assemblies <- function(traces, method = "pca", n_assemblies = NULL,
                               threshold = 0.3, min_cells = 3) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Z-score normalize
  traces_z <- t(scale(t(traces)))

  # Estimate number of assemblies if not provided
  if (is.null(n_assemblies)) {
    n_assemblies <- .estimate_n_assemblies(traces_z)
  }

  # Detect using specified method
  result <- switch(method,
    pca = .detect_assemblies_pca(traces_z, n_assemblies),
    ica = .detect_assemblies_ica(traces_z, n_assemblies),
    nmf = .detect_assemblies_nmf(traces, n_assemblies),
    svd = .detect_assemblies_svd(traces_z, n_assemblies),
    stop("Unknown method: ", method)
  )

  # Extract assembly memberships
  weights <- result$weights
  activations <- result$activations

  # Threshold weights to get discrete memberships
  memberships <- lapply(seq_len(n_assemblies), function(i) {
    members <- which(abs(weights[, i]) > threshold)
    if (length(members) >= min_cells) {
      list(
        cells = members,
        weights = weights[members, i],
        strength = mean(abs(weights[members, i]))
      )
    } else {
      NULL
    }
  })

  # Remove empty assemblies
  memberships <- Filter(Negate(is.null), memberships)

  structure(
    list(
      memberships = memberships,
      weights = weights,
      activations = activations,
      method = method,
      n_assemblies = length(memberships),
      n_cells = n_cells,
      threshold = threshold
    ),
    class = "neural_assemblies"
  )
}

.estimate_n_assemblies <- function(traces_z) {
  # Use parallel analysis or eigenvalue > 1 rule
  pca <- prcomp(t(traces_z), center = FALSE, scale. = FALSE)
  eigenvalues <- pca$sdev^2

  # Marcenko-Pastur threshold for random matrix
  n <- nrow(traces_z)
  p <- ncol(traces_z)
  lambda_max <- (1 + sqrt(n/p))^2

  sum(eigenvalues > lambda_max)
}

.detect_assemblies_pca <- function(traces_z, n_assemblies) {
  pca <- prcomp(t(traces_z), center = FALSE, scale. = FALSE)

  list(
    weights = pca$rotation[, 1:n_assemblies, drop = FALSE],
    activations = t(pca$x[, 1:n_assemblies, drop = FALSE])
  )
}

.detect_assemblies_ica <- function(traces_z, n_assemblies) {
  if (!requireNamespace("fastICA", quietly = TRUE)) {
    warning("fastICA not available, falling back to PCA")
    return(.detect_assemblies_pca(traces_z, n_assemblies))
  }

  ica <- fastICA::fastICA(t(traces_z), n.comp = n_assemblies, method = "C")

  list(
    weights = ica$S,
    activations = t(ica$A)
  )
}

.detect_assemblies_nmf <- function(traces, n_assemblies) {
  # Simple NMF using alternating least squares
  traces_pos <- pmax(traces, 0)

  n_cells <- nrow(traces_pos)
  n_time <- ncol(traces_pos)

  # Initialize
  W <- matrix(runif(n_cells * n_assemblies), n_cells, n_assemblies)
  H <- matrix(runif(n_assemblies * n_time), n_assemblies, n_time)

  # Iterate

  for (iter in 1:100) {
    # Update H
    H <- H * (t(W) %*% traces_pos) / (t(W) %*% W %*% H + 1e-10)
    # Update W
    W <- W * (traces_pos %*% t(H)) / (W %*% H %*% t(H) + 1e-10)
  }

  list(
    weights = W,
    activations = H
  )
}

.detect_assemblies_svd <- function(traces_z, n_assemblies) {
  svd_result <- svd(traces_z, nu = n_assemblies, nv = n_assemblies)

  list(
    weights = svd_result$u,
    activations = diag(svd_result$d[1:n_assemblies]) %*% t(svd_result$v)
  )
}

#' Detect assemblies using correlation clustering
#'
#' Groups cells by their pairwise correlations using hierarchical clustering.
#'
#' @param traces Matrix of traces (cells x time)
#' @param method Clustering method: "hierarchical", "kmeans", "spectral"
#' @param n_assemblies Number of assemblies (NULL for automatic)
#' @param cor_method Correlation method
#'
#' @return List with assembly assignments
#' @export
detect_assemblies_correlation <- function(traces, method = "hierarchical",
                                           n_assemblies = NULL,
                                           cor_method = "pearson") {
  # Compute correlation matrix
  cor_mat <- cor(t(traces), method = cor_method)

  if (method == "hierarchical") {
    # Distance from correlation
    dist_mat <- as.dist(1 - cor_mat)
    hc <- hclust(dist_mat, method = "ward.D2")

    if (is.null(n_assemblies)) {
      # Use silhouette to find optimal k
      n_assemblies <- .find_optimal_k(dist_mat, hc, max_k = 20)
    }

    assignments <- cutree(hc, k = n_assemblies)

  } else if (method == "kmeans") {
    if (is.null(n_assemblies)) n_assemblies <- 5

    # K-means on correlation features
    km <- kmeans(cor_mat, centers = n_assemblies, nstart = 25)
    assignments <- km$cluster

  } else if (method == "spectral") {
    assignments <- .spectral_clustering(cor_mat, n_assemblies)
  }

  # Create membership lists
  memberships <- lapply(unique(assignments), function(k) {
    cells <- which(assignments == k)
    list(
      cells = cells,
      mean_correlation = mean(cor_mat[cells, cells][upper.tri(cor_mat[cells, cells])])
    )
  })

  structure(
    list(
      memberships = memberships,
      assignments = assignments,
      correlation_matrix = cor_mat,
      method = method,
      n_assemblies = length(unique(assignments))
    ),
    class = "neural_assemblies"
  )
}

.find_optimal_k <- function(dist_mat, hc, max_k = 20) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    return(5)  # Default
  }

  sil_widths <- sapply(2:max_k, function(k) {
    clusters <- cutree(hc, k = k)
    mean(cluster::silhouette(clusters, dist_mat)[, 3])
  })

  which.max(sil_widths) + 1
}

.spectral_clustering <- function(similarity, k) {
  if (is.null(k)) k <- 5

  # Degree matrix
  D <- diag(rowSums(similarity))

  # Laplacian

  L <- D - similarity

  # Normalized Laplacian
  D_inv_sqrt <- diag(1 / sqrt(diag(D) + 1e-10))
  L_norm <- D_inv_sqrt %*% L %*% D_inv_sqrt

  # Eigen decomposition
  eig <- eigen(L_norm, symmetric = TRUE)

  # Use k smallest eigenvectors (excluding first)
  features <- eig$vectors[, (ncol(eig$vectors) - k + 1):ncol(eig$vectors)]

  # K-means on eigenvector features
  km <- kmeans(features, centers = k, nstart = 25)
  km$cluster
}

#' Compute assembly activation strength over time
#'
#' @param traces Matrix of traces (cells x time)
#' @param assemblies Neural assemblies object from detect_assemblies
#'
#' @return Matrix of assembly activations (assemblies x time)
#' @export
compute_assembly_activation <- function(traces, assemblies) {
  n_assemblies <- length(assemblies$memberships)
  n_time <- ncol(traces)

  activations <- matrix(0, n_assemblies, n_time)

  for (i in seq_len(n_assemblies)) {
    members <- assemblies$memberships[[i]]$cells
    weights <- assemblies$memberships[[i]]$weights

    if (is.null(weights)) {
      # Equal weighting
      activations[i, ] <- colMeans(traces[members, , drop = FALSE])
    } else {
      # Weighted average
      activations[i, ] <- colSums(traces[members, , drop = FALSE] * abs(weights)) /
        sum(abs(weights))
    }
  }

  activations
}

#' Detect assembly activation events
#'
#' Find time points where assemblies are significantly activated.
#'
#' @param activations Assembly activation matrix
#' @param threshold_sd Threshold in standard deviations
#' @param min_duration Minimum event duration in frames
#'
#' @return List of activation events per assembly
#' @export
detect_assembly_events <- function(activations, threshold_sd = 2,
                                    min_duration = 3) {
  n_assemblies <- nrow(activations)

  events <- lapply(seq_len(n_assemblies), function(i) {
    act <- activations[i, ]
    threshold <- mean(act) + threshold_sd * sd(act)

    # Find threshold crossings
    above <- act > threshold

    # Find event boundaries
    starts <- which(diff(c(FALSE, above)) == 1)
    ends <- which(diff(c(above, FALSE)) == -1)

    # Filter by duration
    durations <- ends - starts + 1
    valid <- durations >= min_duration

    if (sum(valid) == 0) {
      return(data.frame(
        start = integer(),
        end = integer(),
        duration = integer(),
        peak_activation = numeric()
      ))
    }

    data.frame(
      start = starts[valid],
      end = ends[valid],
      duration = durations[valid],
      peak_activation = sapply(which(valid), function(j) {
        max(act[starts[j]:ends[j]])
      })
    )
  })

  names(events) <- paste0("assembly_", seq_len(n_assemblies))
  events
}

#' Test assembly coactivation significance
#'
#' Test whether assemblies coactivate more than expected by chance.
#'
#' @param activations Assembly activation matrix
#' @param n_shuffles Number of shuffles for null distribution
#' @param method Test method: "correlation", "coincidence"
#'
#' @return Matrix of p-values for pairwise coactivation
#' @export
test_assembly_coactivation <- function(activations, n_shuffles = 1000,
                                        method = "correlation") {
  n_assemblies <- nrow(activations)

  # Observed coactivation
  if (method == "correlation") {
    observed <- cor(t(activations))
  } else {
    # Coincidence: fraction of time both active
    binary <- activations > (rowMeans(activations) + rowSds(activations))
    observed <- (binary %*% t(binary)) / ncol(activations)
  }

  # Null distribution by circular shuffling
  null_dist <- array(0, dim = c(n_assemblies, n_assemblies, n_shuffles))

  for (s in seq_len(n_shuffles)) {
    shuffled <- activations
    for (i in seq_len(n_assemblies)) {
      shift <- sample(ncol(activations), 1)
      shuffled[i, ] <- c(activations[i, (shift + 1):ncol(activations)],
                         activations[i, 1:shift])
    }

    if (method == "correlation") {
      null_dist[, , s] <- cor(t(shuffled))
    } else {
      binary <- shuffled > (rowMeans(shuffled) + rowSds(shuffled))
      null_dist[, , s] <- (binary %*% t(binary)) / ncol(shuffled)
    }
  }

  # Compute p-values
  p_values <- matrix(0, n_assemblies, n_assemblies)
  for (i in seq_len(n_assemblies)) {
    for (j in seq_len(n_assemblies)) {
      p_values[i, j] <- mean(null_dist[i, j, ] >= observed[i, j])
    }
  }

  list(
    observed = observed,
    p_values = p_values,
    significant = p_values < 0.05,
    n_shuffles = n_shuffles
  )
}

#' Plot neural assemblies
#'
#' @param assemblies Neural assemblies object
#' @param traces Optional traces for activation plot
#' @param type Plot type: "weights", "correlation", "activation"
#'
#' @return ggplot object
#' @export
plot_assemblies <- function(assemblies, traces = NULL, type = "weights") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  if (type == "weights" && !is.null(assemblies$weights)) {
    # Heatmap of assembly weights
    weights <- assemblies$weights
    df <- expand.grid(
      cell = seq_len(nrow(weights)),
      assembly = seq_len(ncol(weights))
    )
    df$weight <- as.vector(weights)

    ggplot2::ggplot(df, ggplot2::aes(x = assembly, y = cell, fill = weight)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
      ggplot2::labs(title = "Assembly Weights", x = "Assembly", y = "Cell") +
      ggplot2::theme_minimal()

  } else if (type == "correlation" && !is.null(assemblies$correlation_matrix)) {
    # Correlation matrix ordered by assembly
    cor_mat <- assemblies$correlation_matrix
    order_idx <- order(assemblies$assignments)
    cor_ordered <- cor_mat[order_idx, order_idx]

    df <- expand.grid(cell1 = seq_len(nrow(cor_ordered)),
                      cell2 = seq_len(ncol(cor_ordered)))
    df$correlation <- as.vector(cor_ordered)

    ggplot2::ggplot(df, ggplot2::aes(x = cell1, y = cell2, fill = correlation)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                    limits = c(-1, 1)) +
      ggplot2::labs(title = "Correlation Matrix (ordered by assembly)") +
      ggplot2::theme_minimal() +
      ggplot2::coord_fixed()

  } else if (type == "activation" && !is.null(traces)) {
    # Activation time courses
    activations <- compute_assembly_activation(traces, assemblies)

    df <- data.frame(
      time = rep(seq_len(ncol(activations)), nrow(activations)),
      activation = as.vector(t(activations)),
      assembly = factor(rep(seq_len(nrow(activations)), each = ncol(activations)))
    )

    ggplot2::ggplot(df, ggplot2::aes(x = time, y = activation, color = assembly)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ assembly, scales = "free_y", ncol = 1) +
      ggplot2::labs(title = "Assembly Activations", x = "Time", y = "Activation") +
      ggplot2::theme_minimal()
  }
}

# Helper: row standard deviations
rowSds <- function(x) {
  apply(x, 1, sd)
}
