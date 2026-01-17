#' Information Theory Measures
#'
#' Information-theoretic analysis of calcium imaging data.
#' Mutual information, transfer entropy, and related measures.
#'
#' @name information_theory
#' @keywords internal
NULL

#' Compute mutual information between cells
#'
#' Estimate mutual information from calcium traces.
#'
#' @param traces Matrix of traces (cells x time)
#' @param method Estimation method: "binned", "knn", "kernel"
#' @param n_bins Number of bins for histogram method
#' @param k Number of neighbors for knn method
#'
#' @return Matrix of pairwise mutual information
#' @export
#'
#' @examples
#' \dontrun{
#' mi_matrix <- mutual_information(traces, method = "binned")
#' }
mutual_information <- function(traces, method = "binned", n_bins = 10, k = 3) {
  n_cells <- nrow(traces)

  mi_matrix <- matrix(0, n_cells, n_cells)

  for (i in seq_len(n_cells)) {
    for (j in i:n_cells) {
      mi <- switch(method,
        binned = .mi_binned(traces[i, ], traces[j, ], n_bins),
        knn = .mi_knn(traces[i, ], traces[j, ], k),
        kernel = .mi_kernel(traces[i, ], traces[j, ])
      )
      mi_matrix[i, j] <- mi
      mi_matrix[j, i] <- mi
    }
  }

  mi_matrix
}

.mi_binned <- function(x, y, n_bins) {
  # Discretize
  x_bins <- cut(x, breaks = n_bins, labels = FALSE)
  y_bins <- cut(y, breaks = n_bins, labels = FALSE)

  # Joint and marginal frequencies
  joint <- table(x_bins, y_bins) / length(x)
  px <- table(x_bins) / length(x)
  py <- table(y_bins) / length(y)

  # MI = sum p(x,y) * log(p(x,y) / (p(x)*p(y)))
  mi <- 0
  for (i in seq_along(px)) {
    for (j in seq_along(py)) {
      if (joint[i, j] > 0) {
        mi <- mi + joint[i, j] * log2(joint[i, j] / (px[i] * py[j]))
      }
    }
  }

  max(0, mi)  # MI should be non-negative
}

.mi_knn <- function(x, y, k) {
  n <- length(x)

  # Standardize
  x <- (x - mean(x)) / sd(x)
  y <- (y - mean(y)) / sd(y)

  # Combined space
  xy <- cbind(x, y)

  # K-nearest neighbor distances
  dist_xy <- as.matrix(dist(xy))
  dist_x <- as.matrix(dist(x))
  dist_y <- as.matrix(dist(y))

  diag(dist_xy) <- Inf
  diag(dist_x) <- Inf
  diag(dist_y) <- Inf

  mi <- 0
  for (i in seq_len(n)) {
    # Distance to k-th neighbor in joint space
    eps <- sort(dist_xy[i, ])[k]

    # Count neighbors within eps in marginals
    nx <- sum(dist_x[i, ] < eps)
    ny <- sum(dist_y[i, ] < eps)

    mi <- mi + digamma(k) - digamma(nx + 1) - digamma(ny + 1)
  }

  mi <- mi / n + digamma(n)
  max(0, mi / log(2))  # Convert to bits
}

.mi_kernel <- function(x, y) {
  # Kernel density estimation
  n <- length(x)

  # Standardize
  x <- (x - mean(x)) / sd(x)
  y <- (y - mean(y)) / sd(y)

  # Bandwidth (Silverman's rule)
  bw_x <- 1.06 * sd(x) * n^(-1/5)
  bw_y <- 1.06 * sd(y) * n^(-1/5)

  # Estimate densities at each point
  log_pxy <- log_px <- log_py <- numeric(n)

  for (i in seq_len(n)) {
    # Marginals
    px_i <- mean(dnorm(x[i], x, bw_x))
    py_i <- mean(dnorm(y[i], y, bw_y))

    # Joint (product kernel)
    pxy_i <- mean(dnorm(x[i], x, bw_x) * dnorm(y[i], y, bw_y))

    log_px[i] <- log(px_i + 1e-10)
    log_py[i] <- log(py_i + 1e-10)
    log_pxy[i] <- log(pxy_i + 1e-10)
  }

  mi <- mean(log_pxy - log_px - log_py) / log(2)
  max(0, mi)
}

#' Compute transfer entropy
#'
#' Estimate directed information transfer between cells.
#'
#' @param source Source cell trace
#' @param target Target cell trace
#' @param lag Time lag for transfer entropy
#' @param n_bins Number of bins for discretization
#' @param history History length
#'
#' @return Transfer entropy value
#' @export
transfer_entropy <- function(source, target, lag = 1, n_bins = 5, history = 1) {
  n <- length(source)

  # Discretize
  source_d <- cut(source, breaks = n_bins, labels = FALSE)
  target_d <- cut(target, breaks = n_bins, labels = FALSE)

  # Create lagged versions
  # TE = I(target_future ; source_past | target_past)
  valid_idx <- (max(lag, history) + 1):n

  target_future <- target_d[valid_idx]
  source_past <- source_d[valid_idx - lag]
  target_past <- target_d[valid_idx - history]

  # Compute conditional MI
  # TE = H(target_future | target_past) - H(target_future | target_past, source_past)

  # H(X|Y) = H(X,Y) - H(Y)
  h_tf_tp <- .joint_entropy(target_future, target_past, n_bins)
  h_tp <- .entropy(target_past, n_bins)
  h_cond1 <- h_tf_tp - h_tp

  # H(X|Y,Z)
  h_tf_tp_sp <- .joint_entropy_3(target_future, target_past, source_past, n_bins)
  h_tp_sp <- .joint_entropy(target_past, source_past, n_bins)
  h_cond2 <- h_tf_tp_sp - h_tp_sp

  te <- h_cond1 - h_cond2
  max(0, te)
}

.entropy <- function(x, n_bins) {
  p <- table(x) / length(x)
  -sum(p * log2(p + 1e-10))
}

.joint_entropy <- function(x, y, n_bins) {
  p <- table(x, y) / length(x)
  -sum(p * log2(p + 1e-10))
}

.joint_entropy_3 <- function(x, y, z, n_bins) {
  df <- data.frame(x = x, y = y, z = z)
  counts <- table(df)
  p <- counts / sum(counts)
  -sum(p * log2(p + 1e-10))
}

#' Compute transfer entropy matrix
#'
#' Pairwise transfer entropy between all cells.
#'
#' @param traces Matrix of traces (cells x time)
#' @param lag Time lag
#' @param n_bins Number of bins
#'
#' @return Matrix of directed transfer entropy
#' @export
transfer_entropy_matrix <- function(traces, lag = 1, n_bins = 5) {
  n_cells <- nrow(traces)

  te_matrix <- matrix(0, n_cells, n_cells)

  for (i in seq_len(n_cells)) {
    for (j in seq_len(n_cells)) {
      if (i != j) {
        te_matrix[i, j] <- transfer_entropy(traces[i, ], traces[j, ], lag, n_bins)
      }
    }
  }

  te_matrix
}

#' Test transfer entropy significance
#'
#' @param source Source trace
#' @param target Target trace
#' @param lag Time lag
#' @param n_shuffles Number of shuffles
#' @param n_bins Number of bins
#'
#' @return List with TE value and p-value
#' @export
test_transfer_entropy <- function(source, target, lag = 1, n_shuffles = 100,
                                   n_bins = 5) {
  observed_te <- transfer_entropy(source, target, lag, n_bins)

  # Null distribution by shuffling source
  null_te <- replicate(n_shuffles, {
    shuffled <- sample(source)
    transfer_entropy(shuffled, target, lag, n_bins)
  })

  p_value <- mean(null_te >= observed_te)

  list(
    te = observed_te,
    p_value = p_value,
    null_mean = mean(null_te),
    null_sd = sd(null_te),
    z_score = (observed_te - mean(null_te)) / sd(null_te)
  )
}

#' Compute active information storage
#'
#' How much information the cell stores about its own past.
#'
#' @param trace Single cell trace
#' @param history History length
#' @param n_bins Number of bins
#'
#' @return Active information storage value
#' @export
active_information_storage <- function(trace, history = 1, n_bins = 10) {
  n <- length(trace)

  # Discretize
  trace_d <- cut(trace, breaks = n_bins, labels = FALSE)

  valid_idx <- (history + 1):n
  current <- trace_d[valid_idx]
  past <- trace_d[valid_idx - history]

  # AIS = I(current ; past)
  .mi_binned(current, past, n_bins)
}

#' Compute entropy rate
#'
#' Rate at which the system produces new information.
#'
#' @param trace Single cell trace
#' @param max_history Maximum history to consider
#' @param n_bins Number of bins
#'
#' @return Entropy rate estimate
#' @export
entropy_rate <- function(trace, max_history = 5, n_bins = 10) {
  # Discretize
  trace_d <- cut(trace, breaks = n_bins, labels = FALSE)

  # H(X_n | X_{n-1}, ..., X_{n-k}) for increasing k
  cond_entropies <- sapply(1:max_history, function(k) {
    n <- length(trace_d)
    valid_idx <- (k + 1):n

    # Create history matrix
    current <- trace_d[valid_idx]

    # Joint with history
    if (k == 1) {
      past <- trace_d[valid_idx - 1]
      h_joint <- .joint_entropy(current, past, n_bins)
      h_past <- .entropy(past, n_bins)
    } else {
      # Multi-step history
      hist_cols <- sapply(1:k, function(i) trace_d[valid_idx - i])
      hist_string <- apply(hist_cols, 1, paste, collapse = "_")

      joint_string <- paste(current, hist_string, sep = "_")
      h_joint <- .entropy(joint_string, length(unique(joint_string)))
      h_past <- .entropy(hist_string, length(unique(hist_string)))
    }

    h_joint - h_past
  })

  list(
    entropy_rate = cond_entropies[max_history],
    conditional_entropies = cond_entropies
  )
}

#' Information-theoretic network analysis
#'
#' Build directed network based on transfer entropy.
#'
#' @param traces Matrix of traces (cells x time)
#' @param lag Time lag
#' @param threshold_percentile Percentile threshold for edges
#' @param n_bins Number of bins
#'
#' @return igraph network object
#' @export
information_network <- function(traces, lag = 1, threshold_percentile = 95,
                                 n_bins = 5) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package required")
  }

  te_matrix <- transfer_entropy_matrix(traces, lag, n_bins)

  # Threshold
  threshold <- quantile(te_matrix[te_matrix > 0], threshold_percentile / 100)
  adj_matrix <- te_matrix
  adj_matrix[adj_matrix < threshold] <- 0

  # Create directed graph
  g <- igraph::graph_from_adjacency_matrix(
    adj_matrix,
    mode = "directed",
    weighted = TRUE
  )

  # Add node attributes
  igraph::V(g)$mean_activity <- rowMeans(traces)

  # Compute network metrics
  igraph::V(g)$in_strength <- igraph::strength(g, mode = "in")
  igraph::V(g)$out_strength <- igraph::strength(g, mode = "out")

  g
}

#' Compute redundancy and synergy
#'
#' Partial information decomposition for triplets of cells.
#'
#' @param traces Matrix of traces (cells x time)
#' @param target Target cell index
#' @param sources Vector of two source cell indices
#' @param n_bins Number of bins
#'
#' @return List with redundancy, synergy, and unique information
#' @export
partial_information <- function(traces, target, sources, n_bins = 5) {
  t_trace <- traces[target, ]
  s1_trace <- traces[sources[1], ]
  s2_trace <- traces[sources[2], ]

  # I(T; S1)
  mi_t_s1 <- .mi_binned(t_trace, s1_trace, n_bins)

  # I(T; S2)
  mi_t_s2 <- .mi_binned(t_trace, s2_trace, n_bins)

  # I(T; S1, S2) - joint MI
  # Approximate by combining sources
  s_combined <- paste(
    cut(s1_trace, n_bins, labels = FALSE),
    cut(s2_trace, n_bins, labels = FALSE),
    sep = "_"
  )
  t_d <- cut(t_trace, n_bins, labels = FALSE)
  mi_t_s12 <- .mi_binned(as.numeric(factor(s_combined)), t_d, n_bins^2)

  # Redundancy (minimum MI)
  redundancy <- min(mi_t_s1, mi_t_s2)

  # Unique information
  unique_s1 <- mi_t_s1 - redundancy
  unique_s2 <- mi_t_s2 - redundancy

  # Synergy
  synergy <- mi_t_s12 - mi_t_s1 - mi_t_s2 + redundancy

  list(
    redundancy = redundancy,
    unique_s1 = unique_s1,
    unique_s2 = unique_s2,
    synergy = synergy,
    total = mi_t_s12
  )
}

#' Plot information matrix
#'
#' @param info_matrix MI or TE matrix
#' @param title Plot title
#' @param directed Is the matrix directed (asymmetric)?
#'
#' @return ggplot object
#' @export
plot_information_matrix <- function(info_matrix, title = "Information Matrix",
                                     directed = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  n <- nrow(info_matrix)
  df <- expand.grid(source = 1:n, target = 1:n)
  df$value <- as.vector(info_matrix)

  ggplot2::ggplot(df, ggplot2::aes(x = target, y = source, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(title = title, x = "Target Cell", y = "Source Cell",
                  fill = if (directed) "TE (bits)" else "MI (bits)") +
    ggplot2::theme_minimal() +
    ggplot2::coord_fixed()
}
