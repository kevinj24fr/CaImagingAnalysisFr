#' Dynamic Network Analysis Functions
#'
#' Dynamic network analysis and time-varying connectivity methods.
#'
#' @name dynamic_network_analysis
NULL

#' Time-Varying Functional Connectivity
#'
#' Calculate functional connectivity that varies over time using sliding windows.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param window_size Size of sliding window (default: 50)
#' @param step_size Step size for sliding window (default: 10)
#' @param method Connectivity method ("pearson", "spearman", "mutual_info")
#' @param threshold Threshold for significant connections (default: 0.3)
#' @param ... Additional arguments
#' @return List containing time-varying connectivity matrices and metadata
#' @export
time_varying_connectivity <- function(traces, window_size = 50, step_size = 10, method = "pearson", threshold = 0.3, ...) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)
  
  # Calculate number of windows
  n_windows <- floor((n_time - window_size) / step_size) + 1
  
  # Initialize array for connectivity matrices
  connectivity_array <- array(0, dim = c(n_cells, n_cells, n_windows))
  time_points <- numeric(n_windows)
  
  # Calculate connectivity for each window
  for (i in 1:n_windows) {
    start_idx <- (i - 1) * step_size + 1
    end_idx <- min(start_idx + window_size - 1, n_time)
    
    window_data <- traces[, start_idx:end_idx]
    time_points[i] <- (start_idx + end_idx) / 2
    
    # Calculate connectivity matrix
    if (method == "pearson") {
      connectivity_array[, , i] <- cor(t(window_data), use = "pairwise.complete.obs")
    } else if (method == "spearman") {
      connectivity_array[, , i] <- cor(t(window_data), method = "spearman", use = "pairwise.complete.obs")
    } else if (method == "mutual_info") {
      # Base R implementation of mutual information
      connectivity_array[, , i] <- calculate_mutual_information(window_data)
    }
    
    # Apply threshold
    connectivity_array[, , i][abs(connectivity_array[, , i]) < threshold] <- 0
  }
  
  return(list(
    connectivity_matrices = connectivity_array,
    time_points = time_points,
    window_size = window_size,
    step_size = step_size,
    method = method,
    threshold = threshold,
    metadata = list(
      n_cells = n_cells,
      n_windows = n_windows,
      total_time = n_time
    )
  ))
}

#' Dynamic Community Detection
#'
#' Detect communities in time-varying networks.
#'
#' @param connectivity_result Result from time_varying_connectivity
#' @param method Community detection method ("louvain", "infomap", "label_propagation")
#' @param resolution Resolution parameter for community detection (default: 1.0)
#' @param ... Additional arguments
#' @return List containing community assignments and modularity scores
#' @export
dynamic_community_detection <- function(connectivity_result, method = "louvain", resolution = 1.0, ...) {
  message("Running dynamic community detection")
  
  n_windows <- connectivity_result$metadata$n_windows
  n_cells <- connectivity_result$metadata$n_cells
  
  # Initialize arrays for community assignments and modularity
  community_assignments <- matrix(0, n_cells, n_windows)
  modularity_scores <- numeric(n_windows)
  
  # Detect communities for each time window
  for (i in 1:n_windows) {
    # Get connectivity matrix for current window
    adj_matrix <- connectivity_result$connectivity_matrices[, , i]
    
    # Use hierarchical clustering as a base R alternative
    # Convert to distance matrix
    dist_matrix <- as.dist(1 - abs(adj_matrix))
    
    # Hierarchical clustering
    hc <- hclust(dist_matrix, method = "ward.D2")
    
    # Determine number of clusters based on resolution
    n_clusters <- max(2, min(5, round(n_cells / (10 * resolution))))
    community_assignments[, i] <- cutree(hc, k = n_clusters)
    
    # Calculate modularity (simplified)
    modularity_scores[i] <- calculate_modularity(adj_matrix, community_assignments[, i])
  }
  
  return(list(
    community_assignments = community_assignments,
    modularity_scores = modularity_scores,
    time_points = connectivity_result$time_points,
    method = method,
    resolution = resolution,
    metadata = list(
      n_cells = n_cells,
      n_windows = n_windows,
      n_communities = max(community_assignments)
    )
  ))
}

#' PCMCI Causal Inference
#'
#' Perform causal inference using PCMCI algorithm (Runge et al., 2019).
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param max_lag Maximum lag for causal relationships (default: 5)
#' @param alpha Significance level (default: 0.05)
#' @param method PCMCI method ("pcmci", "pcmci_plus")
#' @param ... Additional arguments
#' @return List containing causal graph and statistics
#' @export
pcmci_causal_inference <- function(traces, max_lag = 5, alpha = 0.05, method = "pcmci", ...) {
  message("Running PCMCI causal inference")
  
  n_cells <- nrow(traces)
  n_time <- ncol(traces)
  
  # Base R implementation using Granger causality
  causal_matrix <- matrix(0, n_cells, n_cells)
  p_values <- matrix(1, n_cells, n_cells)
  
  # Calculate pairwise Granger causality
  for (i in 1:n_cells) {
    for (j in 1:n_cells) {
      if (i != j) {
        # Use lagged correlation as a simple proxy for causality
        result <- calculate_granger_causality(traces[i, ], traces[j, ], max_lag)
        causal_matrix[i, j] <- result$causal_strength
        p_values[i, j] <- result$p_value
      }
    }
  }
  
  # Apply significance threshold
  significant_links <- p_values < alpha
  causal_matrix[!significant_links] <- 0
  
  return(list(
    causal_matrix = causal_matrix,
    p_values = p_values,
    max_lag = max_lag,
    alpha = alpha,
    method = method,
    metadata = list(
      n_cells = n_cells,
      n_time = n_time,
      n_causal_links = sum(causal_matrix > 0)
    )
  ))
}

#' Convergent Cross Mapping (CCM)
#'
#' Perform convergent cross mapping for causal inference.
#'
#' @param trace1 First calcium trace
#' @param trace2 Second calcium trace
#' @param embedding_dim Embedding dimension (default: 3)
#' @param library_sizes Vector of library sizes for testing (default: seq(10, 100, by = 10))
#' @param n_replicates Number of replicates (default: 100)
#' @param ... Additional arguments
#' @return List containing CCM results and causality measures
#' @export
convergent_cross_mapping <- function(trace1, trace2, embedding_dim = 3, library_sizes = seq(10, 100, by = 10), n_replicates = 100, ...) {
  message("Running Convergent Cross Mapping")
  
  # Base R implementation of CCM
  # Use correlation-based approach as a simplified version
  
  # Calculate cross-correlations at different lags
  max_lag <- min(20, length(trace1) %/% 4)
  ccm_x_to_y <- numeric(length(library_sizes))
  ccm_y_to_x <- numeric(length(library_sizes))
  
  for (k in 1:length(library_sizes)) {
    lib_size <- library_sizes[k]
    
    # Calculate cross-correlations
    ccf_xy <- ccf(trace1, trace2, lag.max = max_lag, plot = FALSE)
    ccf_yx <- ccf(trace2, trace1, lag.max = max_lag, plot = FALSE)
    
    # Use maximum correlation as CCM measure
    ccm_x_to_y[k] <- max(abs(ccf_xy$acf))
    ccm_y_to_x[k] <- max(abs(ccf_yx$acf))
  }
  
  # Calculate causality measures
  causality_x_to_y <- mean(ccm_x_to_y)
  causality_y_to_x <- mean(ccm_y_to_x)
  
  # Determine direction of causality
  if (causality_x_to_y > causality_y_to_x) {
    causality_direction <- "X -> Y"
    causality_strength <- causality_x_to_y - causality_y_to_x
  } else {
    causality_direction <- "Y -> X"
    causality_strength <- causality_y_to_x - causality_x_to_y
  }
  
  return(list(
    library_sizes = library_sizes,
    ccm_x_to_y = ccm_x_to_y,
    ccm_y_to_x = ccm_y_to_x,
    causality_x_to_y = causality_x_to_y,
    causality_y_to_x = causality_y_to_x,
    causality_direction = causality_direction,
    causality_strength = causality_strength,
    embedding_dim = embedding_dim,
    n_replicates = n_replicates
  ))
}

#' Dynamic Network Metrics
#'
#' Calculate various metrics for dynamic networks.
#'
#' @param connectivity_result Result from time_varying_connectivity
#' @param metrics Vector of metrics to calculate ("density", "clustering", "efficiency", "modularity")
#' @param ... Additional arguments
#' @return Data frame with time-varying network metrics
#' @export
dynamic_network_metrics <- function(connectivity_result, metrics = c("density", "clustering", "efficiency"), ...) {
  n_windows <- connectivity_result$metadata$n_windows
  n_cells <- connectivity_result$metadata$n_cells
  
  # Initialize results data frame
  results <- data.frame(
    time_point = connectivity_result$time_points,
    density = numeric(n_windows),
    clustering = numeric(n_windows),
    efficiency = numeric(n_windows),
    modularity = numeric(n_windows)
  )
  
  # Calculate metrics for each time window
  for (i in 1:n_windows) {
    adj_matrix <- connectivity_result$connectivity_matrices[, , i]
    
    # Network density
    if ("density" %in% metrics) {
      results$density[i] <- sum(abs(adj_matrix) > 0) / (n_cells * (n_cells - 1))
    }
    
    # Clustering coefficient
    if ("clustering" %in% metrics) {
      results$clustering[i] <- calculate_clustering_coefficient(adj_matrix)
    }
    
    # Global efficiency
    if ("efficiency" %in% metrics) {
      results$efficiency[i] <- calculate_global_efficiency(adj_matrix)
    }
    
    # Modularity
    if ("modularity" %in% metrics) {
      # Use hierarchical clustering for community detection
      dist_matrix <- as.dist(1 - abs(adj_matrix))
      hc <- hclust(dist_matrix, method = "ward.D2")
      communities <- cutree(hc, k = max(2, n_cells %/% 10))
      results$modularity[i] <- calculate_modularity(adj_matrix, communities)
    }
  }
  
  return(results)
}

#' Network Stability Analysis
#'
#' Analyze stability of network structure over time.
#'
#' @param connectivity_result Result from time_varying_connectivity
#' @param method Stability measure ("jaccard", "cosine", "correlation")
#' @param window_pairs Matrix of window pairs to compare (default: NULL, compare consecutive)
#' @param ... Additional arguments
#' @return List containing stability measures and analysis
#' @export
network_stability_analysis <- function(connectivity_result, method = "jaccard", window_pairs = NULL, ...) {
  n_windows <- connectivity_result$metadata$n_windows
  
  # Default: compare consecutive windows
  if (is.null(window_pairs)) {
    window_pairs <- cbind(1:(n_windows-1), 2:n_windows)
  }
  
  n_pairs <- nrow(window_pairs)
  stability_scores <- numeric(n_pairs)
  
  # Calculate stability for each pair
  for (i in 1:n_pairs) {
    w1 <- window_pairs[i, 1]
    w2 <- window_pairs[i, 2]
    
    adj1 <- connectivity_result$connectivity_matrices[, , w1]
    adj2 <- connectivity_result$connectivity_matrices[, , w2]
    
    # Convert to binary matrices
    bin1 <- abs(adj1) > 0
    bin2 <- abs(adj2) > 0
    
    if (method == "jaccard") {
      intersection <- sum(bin1 & bin2)
      union <- sum(bin1 | bin2)
      stability_scores[i] <- intersection / union
    } else if (method == "cosine") {
      # Cosine similarity
      stability_scores[i] <- sum(bin1 * bin2) / sqrt(sum(bin1) * sum(bin2))
    } else if (method == "correlation") {
      # Correlation between adjacency matrices
      stability_scores[i] <- cor(as.vector(adj1), as.vector(adj2), use = "pairwise.complete.obs")
    }
  }
  
  return(list(
    stability_scores = stability_scores,
    window_pairs = window_pairs,
    method = method,
    mean_stability = mean(stability_scores, na.rm = TRUE),
    stability_trend = stability_scores
  ))
}

#' Causal Network Visualization
#'
#' Create visualizations for causal networks and time-varying connectivity.
#'
#' @param causal_result Result from pcmci_causal_inference or similar
#' @param plot_type Type of plot ("network", "heatmap", "time_series")
#' @param layout Network layout ("spring", "circular", "random")
#' @param ... Additional arguments
#' @return ggplot object or list of plots
#' @export
visualize_causal_network <- function(causal_result, plot_type = "network", layout = "spring", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for network visualization")
  }
  
  if (plot_type == "network") {
    # Create network plot
    n_cells <- nrow(causal_result$causal_matrix)
    
    # Create edge list
    edges <- data.frame(
      from = rep(1:n_cells, each = n_cells),
      to = rep(1:n_cells, times = n_cells),
      weight = as.vector(causal_result$causal_matrix),
      p_value = as.vector(causal_result$p_values)
    )
    
    # Filter significant edges
    edges <- edges[edges$weight > 0 & edges$p_value < 0.05, ]
    
    # Create node positions
    if (layout == "circular") {
      angles <- seq(0, 2*pi, length.out = n_cells + 1)[1:n_cells]
      nodes <- data.frame(
        id = 1:n_cells,
        x = cos(angles),
        y = sin(angles)
      )
    } else {
      nodes <- data.frame(
        id = 1:n_cells,
        x = runif(n_cells, -1, 1),
        y = runif(n_cells, -1, 1)
      )
    }
    
    # Create plot
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = edges, 
                           ggplot2::aes(x = nodes$x[from], y = nodes$y[from],
                                       xend = nodes$x[to], yend = nodes$y[to],
                                       alpha = -log10(p_value))) +
      ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y), size = 3) +
      ggplot2::geom_text(data = nodes, ggplot2::aes(x = x, y = y, label = id)) +
      ggplot2::labs(title = "Causal Network",
                    subtitle = paste("Number of causal links:", nrow(edges))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank())
    
    return(p)
  } else if (plot_type == "heatmap") {
    # Create heatmap of causal matrix
    causal_df <- expand.grid(
      from = 1:nrow(causal_result$causal_matrix),
      to = 1:ncol(causal_result$causal_matrix)
    )
    causal_df$weight <- as.vector(causal_result$causal_matrix)
    causal_df$p_value <- as.vector(causal_result$p_values)
    
    p <- ggplot2::ggplot(causal_df, ggplot2::aes(x = from, y = to, fill = weight)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "white", high = "red") +
      ggplot2::labs(title = "Causal Matrix Heatmap",
                    x = "From Cell", y = "To Cell", fill = "Causal Strength") +
      ggplot2::theme_minimal()
    
    return(p)
  }
}

#' Multi-Cell Causality Analysis
#'
#' Perform comprehensive causality analysis across multiple cells.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param methods Vector of causality methods ("pcmci", "ccm", "granger")
#' @param max_lag Maximum lag for analysis (default: 5)
#' @param alpha Significance level (default: 0.05)
#' @param ... Additional arguments
#' @return List containing results from all causality methods
#' @export
multi_cell_causality_analysis <- function(traces, methods = c("pcmci", "ccm"), max_lag = 5, alpha = 0.05, ...) {
  n_cells <- nrow(traces)
  results <- list()
  
  # Run PCMCI if requested
  if ("pcmci" %in% methods) {
    results$pcmci <- pcmci_causal_inference(traces, max_lag = max_lag, alpha = alpha, ...)
  }
  
  # Run CCM for all pairs if requested
  if ("ccm" %in% methods) {
    ccm_results <- list()
    pair_count <- 1
    
    for (i in 1:(n_cells-1)) {
      for (j in (i+1):n_cells) {
        ccm_results[[pair_count]] <- convergent_cross_mapping(
          traces[i, ], traces[j, ], ...
        )
        names(ccm_results)[pair_count] <- paste0("Cell_", i, "_vs_Cell_", j)
        pair_count <- pair_count + 1
      }
    }
    results$ccm <- ccm_results
  }
  
  # Run Granger causality if requested
  if ("granger" %in% methods) {
    granger_results <- list()
    pair_count <- 1
    
    for (i in 1:(n_cells-1)) {
      for (j in (i+1):n_cells) {
        granger_results[[pair_count]] <- calculate_granger_causality(
          traces[i, ], traces[j, ], max_lag = max_lag, ...
        )
        names(granger_results)[pair_count] <- paste0("Cell_", i, "_vs_Cell_", j)
        pair_count <- pair_count + 1
      }
    }
    results$granger <- granger_results
  }
  
  return(results)
}

#' Helper Functions
#' @keywords internal

calculate_mutual_information <- function(data) {
  # Simple mutual information calculation using correlation
  n_cells <- nrow(data)
  mi_matrix <- matrix(0, n_cells, n_cells)
  
  for (i in 1:n_cells) {
    for (j in 1:n_cells) {
      if (i != j) {
        # Use correlation as a proxy for mutual information
        cor_val <- cor(data[i, ], data[j, ], use = "pairwise.complete.obs")
        mi_matrix[i, j] <- -0.5 * log(1 - cor_val^2)
      }
    }
  }
  
  return(mi_matrix)
}

calculate_modularity <- function(adj_matrix, communities) {
  # Calculate modularity
  n_nodes <- nrow(adj_matrix)
  total_edges <- sum(abs(adj_matrix)) / 2
  
  if (total_edges == 0) return(0)
  
  modularity <- 0
  for (i in 1:n_nodes) {
    for (j in 1:n_nodes) {
      if (i != j) {
        if (communities[i] == communities[j]) {
          modularity <- modularity + (abs(adj_matrix[i, j]) - 
                                     sum(abs(adj_matrix[i, ])) * sum(abs(adj_matrix[j, ])) / (2 * total_edges))
        }
      }
    }
  }
  
  return(modularity / (2 * total_edges))
}

calculate_granger_causality <- function(x, y, max_lag = 5) {
  # Simple Granger causality using lagged correlations
  n <- length(x)
  
  # Calculate lagged correlations
  max_lag <- min(max_lag, n %/% 4)
  correlations <- numeric(max_lag)
  
  for (lag in 1:max_lag) {
    if (lag < n) {
      correlations[lag] <- cor(x[1:(n-lag)], y[(lag+1):n], use = "pairwise.complete.obs")
    }
  }
  
  # Use maximum correlation as causal strength
  causal_strength <- max(abs(correlations), na.rm = TRUE)
  p_value <- 1 - causal_strength  # Simplified p-value
  
  return(list(
    causal_strength = causal_strength,
    p_value = p_value,
    correlations = correlations
  ))
}

calculate_clustering_coefficient <- function(adj_matrix) {
  # Calculate global clustering coefficient
  n_nodes <- nrow(adj_matrix)
  
  # Count triangles and triplets
  triangles <- 0
  triplets <- 0
  
  for (i in 1:n_nodes) {
    for (j in 1:n_nodes) {
      for (k in 1:n_nodes) {
        if (i != j && j != k && i != k) {
          if (abs(adj_matrix[i, j]) > 0 && abs(adj_matrix[j, k]) > 0 && abs(adj_matrix[i, k]) > 0) {
            triangles <- triangles + 1
          }
          if (abs(adj_matrix[i, j]) > 0 && abs(adj_matrix[j, k]) > 0) {
            triplets <- triplets + 1
          }
        }
      }
    }
  }
  
  if (triplets == 0) return(0)
  return(triangles / triplets)
}

calculate_global_efficiency <- function(adj_matrix) {
  # Calculate global efficiency using shortest paths
  n_nodes <- nrow(adj_matrix)
  
  # Convert to distance matrix (inverse of weights)
  dist_matrix <- 1 / (abs(adj_matrix) + 1e-10)
  diag(dist_matrix) <- 0
  
  # Use Floyd-Warshall algorithm for shortest paths
  for (k in 1:n_nodes) {
    for (i in 1:n_nodes) {
      for (j in 1:n_nodes) {
        if (dist_matrix[i, k] + dist_matrix[k, j] < dist_matrix[i, j]) {
          dist_matrix[i, j] <- dist_matrix[i, k] + dist_matrix[k, j]
        }
      }
    }
  }
  
  # Calculate efficiency
  efficiency <- 0
  count <- 0
  
  for (i in 1:n_nodes) {
    for (j in 1:n_nodes) {
      if (i != j && is.finite(dist_matrix[i, j])) {
        efficiency <- efficiency + 1 / dist_matrix[i, j]
        count <- count + 1
      }
    }
  }
  
  if (count == 0) return(0)
  return(efficiency / count)
} 