#' Network Analysis Functions
#'
#' Network analysis and connectivity methods for calcium imaging data.
#'
#' @name network_analysis
NULL

#' Functional Connectivity Analysis
#'
#' Compute functional connectivity between calcium imaging traces.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param method Connectivity method ("correlation", "mutual_info", "granger", "transfer_entropy")
#' @param threshold_method Thresholding method ("absolute", "percentile", "density")
#' @param threshold_value Threshold value
#' @param threshold Alternative parameter name for threshold_value (for backward compatibility)
#' @param ... Additional arguments
#' @return List containing connectivity matrix and network properties
#' @export
functional_connectivity <- function(traces, method = "correlation", 
                                  threshold_method = "percentile", 
                                  threshold_value = 0.95, 
                                  threshold = NULL, ...) {
  message("Computing functional connectivity")
  
  # Handle backward compatibility for threshold parameter
  if (!is.null(threshold)) {
    threshold_value <- threshold
  }
  
  # Base R implementation of functional connectivity analysis
  
  if (method == "correlation") {
    connectivity_matrix <- correlation_connectivity(traces, ...)
  } else if (method == "mutual_info") {
    connectivity_matrix <- mutual_info_connectivity(traces, ...)
  } else if (method == "granger") {
    connectivity_matrix <- granger_causality(traces, ...)
  } else if (method == "transfer_entropy") {
    connectivity_matrix <- transfer_entropy_connectivity(traces, ...)
  } else {
    stop("Unknown connectivity method: ", method)
  }
  
  # Apply thresholding
  thresholded_matrix <- threshold_connectivity(connectivity_matrix, 
                                             threshold_method, 
                                             threshold_value)
  
  # Compute network properties
  network_properties <- compute_network_properties(thresholded_matrix)
  
  # Extract connections (edge list)
  connections <- extract_connections(thresholded_matrix)
  
  return(list(
    connectivity_matrix = connectivity_matrix,
    thresholded_matrix = thresholded_matrix,
    network_properties = network_properties,
    connections = connections,
    method = method,
    parameters = list(threshold_method = threshold_method, 
                     threshold_value = threshold_value)
  ))
}

#' Correlation-based Connectivity
#'
#' Compute correlation-based functional connectivity.
#'
#' @param traces Matrix of calcium traces
#' @param correlation_type Type of correlation ("pearson", "spearman", "kendall")
#' @param ... Additional arguments
#' @return Correlation matrix
#' @keywords internal
correlation_connectivity <- function(traces, correlation_type = "pearson", ...) {
  # Compute correlation matrix
  if (correlation_type == "pearson") {
    cor_matrix <- cor(t(traces), method = "pearson")
  } else if (correlation_type == "spearman") {
    cor_matrix <- cor(t(traces), method = "spearman")
  } else if (correlation_type == "kendall") {
    cor_matrix <- cor(t(traces), method = "kendall")
  } else {
    stop("Unknown correlation type: ", correlation_type)
  }
  
  # Remove diagonal (self-connections)
  diag(cor_matrix) <- 0
  
  return(cor_matrix)
}

#' Mutual Information Connectivity
#'
#' Compute mutual information-based functional connectivity.
#'
#' @param traces Matrix of calcium traces
#' @param n_bins Number of bins for discretization
#' @param ... Additional arguments
#' @return Mutual information matrix
#' @keywords internal
mutual_info_connectivity <- function(traces, n_bins = 10, ...) {
  n_cells <- nrow(traces)
  mi_matrix <- matrix(0, n_cells, n_cells)
  
  for (i in 1:n_cells) {
    for (j in 1:n_cells) {
      if (i != j) {
        mi_matrix[i, j] <- compute_mutual_information(traces[i, ], traces[j, ], n_bins)
      }
    }
  }
  
  return(mi_matrix)
}

#' Granger Causality
#'
#' Compute Granger causality between time series.
#'
#' @param traces Matrix of calcium traces OR first time series vector
#' @param max_lag Maximum lag for causality analysis OR second time series vector
#' @param ... Additional arguments
#' @return Granger causality matrix or pairwise result
#' @export
granger_causality <- function(traces, max_lag = 5, ...) {
  # Check if first argument is a vector (pairwise analysis)
  if (is.vector(traces) && is.vector(max_lag)) {
    # Pairwise Granger causality between two vectors
    x <- traces
    y <- max_lag
    
    # Extract max_lag from additional arguments
    args <- list(...)
    max_lag_val <- if ("max_lag" %in% names(args)) args$max_lag else 5
    
    return(compute_granger_causality(x, y, max_lag_val))
  }
  
  # Matrix-based analysis (original functionality)
  n_cells <- nrow(traces)
  gc_matrix <- matrix(0, n_cells, n_cells)
  
  for (i in 1:n_cells) {
    for (j in 1:n_cells) {
      if (i != j) {
        gc_matrix[i, j] <- compute_granger_causality(traces[i, ], traces[j, ], max_lag)
      }
    }
  }
  
  return(gc_matrix)
}

#' Transfer Entropy Connectivity
#'
#' Compute transfer entropy-based functional connectivity.
#'
#' @param traces Matrix of calcium traces
#' @param n_bins Number of bins for discretization
#' @param max_lag Maximum lag for transfer entropy
#' @param ... Additional arguments
#' @return Transfer entropy matrix
#' @keywords internal
transfer_entropy_connectivity <- function(traces, n_bins = 5, max_lag = 3, ...) {
  n_cells <- nrow(traces)
  te_matrix <- matrix(0, n_cells, n_cells)
  
  for (i in 1:n_cells) {
    for (j in 1:n_cells) {
      if (i != j) {
        te_matrix[i, j] <- compute_transfer_entropy(traces[i, ], traces[j, ], n_bins, max_lag)
      }
    }
  }
  
  return(te_matrix)
}

#' Threshold Connectivity Matrix
#'
#' Apply thresholding to connectivity matrix.
#'
#' @param connectivity_matrix Input connectivity matrix
#' @param method Thresholding method
#' @param value Threshold value
#' @return Thresholded matrix
#' @keywords internal
threshold_connectivity <- function(connectivity_matrix, method, value) {
  if (method == "absolute") {
    thresholded <- connectivity_matrix
    thresholded[abs(thresholded) < value] <- 0
  } else if (method == "percentile") {
    threshold <- quantile(abs(connectivity_matrix), value, na.rm = TRUE)
    thresholded <- connectivity_matrix
    thresholded[abs(thresholded) < threshold] <- 0
  } else if (method == "density") {
    # Keep top connections to achieve desired density
    n_connections <- round(value * length(connectivity_matrix))
    sorted_values <- sort(abs(connectivity_matrix), decreasing = TRUE)
    threshold <- sorted_values[n_connections]
    thresholded <- connectivity_matrix
    thresholded[abs(thresholded) < threshold] <- 0
  } else {
    stop("Unknown thresholding method: ", method)
  }
  
  return(thresholded)
}

#' Compute Network Properties
#'
#' Compute various network properties.
#'
#' @param adjacency_matrix Adjacency matrix
#' @return List of network properties
#' @keywords internal
compute_network_properties <- function(adjacency_matrix) {
  # Convert to igraph object
  if (!requireNamespace("igraph", quietly = TRUE)) {
    warning("Package 'igraph' not available. Returning basic properties.")
    return(compute_basic_properties(adjacency_matrix))
  }
  
  # Create undirected graph (take absolute values)
  abs_matrix <- abs(adjacency_matrix)
  g <- igraph::graph_from_adjacency_matrix(abs_matrix, mode = "undirected", weighted = TRUE)
  
  # Compute properties
  properties <- list(
    n_nodes = igraph::vcount(g),
    n_edges = igraph::ecount(g),
    density = igraph::edge_density(g),
    average_degree = mean(igraph::degree(g)),
    clustering_coefficient = igraph::transitivity(g, type = "global"),
    average_path_length = igraph::mean_distance(g),
    diameter = igraph::diameter(g),
    modularity = compute_modularity(g),
    degree_distribution = igraph::degree_distribution(g),
    betweenness_centrality = igraph::betweenness(g),
    closeness_centrality = igraph::closeness(g),
    eigenvector_centrality = igraph::eigen_centrality(g)$vector
  )
  
  return(properties)
}

#' Compute Basic Network Properties
#'
#' Compute basic network properties without igraph dependency.
#'
#' @param adjacency_matrix Adjacency matrix
#' @return List of basic properties
#' @keywords internal
compute_basic_properties <- function(adjacency_matrix) {
  n_nodes <- nrow(adjacency_matrix)
  
  # Basic properties
  n_edges <- sum(adjacency_matrix != 0) / 2  # Undirected graph
  density <- n_edges / (n_nodes * (n_nodes - 1) / 2)
  
  # Degree distribution
  degrees <- rowSums(adjacency_matrix != 0)
  average_degree <- mean(degrees)
  
  # Clustering coefficient (simplified)
  clustering_coefficient <- compute_simple_clustering(adjacency_matrix)
  
  return(list(
    n_nodes = n_nodes,
    n_edges = n_edges,
    density = density,
    average_degree = average_degree,
    clustering_coefficient = clustering_coefficient,
    degree_distribution = degrees
  ))
}

#' Compute Simple Clustering Coefficient
#'
#' Compute clustering coefficient without igraph.
#'
#' @param adjacency_matrix Adjacency matrix
#' @return Clustering coefficient
#' @keywords internal
compute_simple_clustering <- function(adjacency_matrix) {
  n_nodes <- nrow(adjacency_matrix)
  total_triangles <- 0
  total_triplets <- 0
  
  for (i in 1:n_nodes) {
    for (j in 1:n_nodes) {
      for (k in 1:n_nodes) {
        if (i != j && j != k && i != k) {
          # Check if i-j-k forms a triplet
          if (adjacency_matrix[i, j] != 0 && adjacency_matrix[j, k] != 0) {
            total_triplets <- total_triplets + 1
            
            # Check if i-k also exists (forming a triangle)
            if (adjacency_matrix[i, k] != 0) {
              total_triangles <- total_triangles + 1
            }
          }
        }
      }
    }
  }
  
  if (total_triplets == 0) {
    return(0)
  } else {
    return(total_triangles / total_triplets)
  }
}

#' Compute Modularity
#'
#' Compute modularity of network.
#'
#' @param graph igraph object
#' @return Modularity value
#' @keywords internal
compute_modularity <- function(graph) {
  # Use Louvain community detection
  communities <- igraph::cluster_louvain(graph)
  modularity <- igraph::modularity(graph, igraph::membership(communities))
  return(modularity)
}

#' Compute Mutual Information
#'
#' Compute mutual information between two variables.
#'
#' @param x First variable
#' @param y Second variable
#' @param n_bins Number of bins for discretization
#' @return Mutual information value
#' @keywords internal
compute_mutual_information <- function(x, y, n_bins) {
  # Discretize variables
  x_binned <- cut(x, breaks = n_bins, labels = FALSE)
  y_binned <- cut(y, breaks = n_bins, labels = FALSE)
  
  # Compute joint and marginal distributions
  joint_dist <- table(x_binned, y_binned) / length(x)
  x_dist <- table(x_binned) / length(x)
  y_dist <- table(y_binned) / length(y)
  
  # Compute mutual information
  mi <- 0
  for (i in 1:nrow(joint_dist)) {
    for (j in 1:ncol(joint_dist)) {
      if (joint_dist[i, j] > 0) {
        mi <- mi + joint_dist[i, j] * log2(joint_dist[i, j] / (x_dist[i] * y_dist[j]))
      }
    }
  }
  
  return(mi)
}

#' Compute Granger Causality
#'
#' Compute Granger causality between two time series.
#'
#' @param x First time series
#' @param y Second time series
#' @param max_lag Maximum lag
#' @return Granger causality value
#' @keywords internal
compute_granger_causality <- function(x, y, max_lag) {
  n <- length(x)
  
  # Ensure max_lag is valid
  if (max_lag >= n) {
    max_lag <- max(1, n - 1)
  }
  
  # Fit unrestricted model (y ~ y_lags + x_lags)
  y_lags <- create_lags(y, max_lag)
  x_lags <- create_lags(x, max_lag)
  
  # Ensure both matrices have the same dimensions
  min_rows <- min(nrow(y_lags), nrow(x_lags))
  if (min_rows < 1) {
    return(0)
  }
  
  y_lags <- y_lags[1:min_rows, , drop = FALSE]
  x_lags <- x_lags[1:min_rows, , drop = FALSE]
  
  # Remove NA values
  valid_indices <- complete.cases(y_lags, x_lags)
  y_valid <- y[(max_lag + 1):n][1:min_rows][valid_indices]
  y_lags_valid <- y_lags[valid_indices, , drop = FALSE]
  x_lags_valid <- x_lags[valid_indices, , drop = FALSE]
  
  if (length(y_valid) < max_lag + 1) {
    return(0)
  }
  
  # Fit models
  unrestricted_model <- lm(y_valid ~ y_lags_valid + x_lags_valid)
  restricted_model <- lm(y_valid ~ y_lags_valid)
  
  # Compute F-statistic
  rss_unrestricted <- sum(residuals(unrestricted_model)^2)
  rss_restricted <- sum(residuals(restricted_model)^2)
  
  if (rss_restricted == 0) {
    return(0)
  }
  
  f_stat <- ((rss_restricted - rss_unrestricted) / max_lag) / (rss_unrestricted / (length(y_valid) - 2 * max_lag - 1))
  
  return(f_stat)
}

#' Compute Transfer Entropy
#'
#' Compute transfer entropy between two time series.
#'
#' @param x Source time series
#' @param y Target time series
#' @param n_bins Number of bins
#' @param max_lag Maximum lag
#' @return Transfer entropy value
#' @keywords internal
compute_transfer_entropy <- function(x, y, n_bins, max_lag) {
  n <- length(x)
  
  # Discretize time series
  x_binned <- cut(x, breaks = n_bins, labels = FALSE)
  y_binned <- cut(y, breaks = n_bins, labels = FALSE)
  
  # Create lagged versions
  y_lags <- create_lags(y_binned, max_lag)
  x_lags <- create_lags(x_binned, max_lag)
  
  # Remove NA values
  valid_indices <- complete.cases(y_binned[(max_lag + 1):n], y_lags, x_lags)
  
  if (sum(valid_indices) < 10) {
    return(0)
  }
  
  # Compute conditional entropies
  y_future <- y_binned[(max_lag + 1):n][valid_indices]
  y_past <- y_lags[valid_indices, ]
  x_past <- x_lags[valid_indices, ]
  
  # H(Y_future | Y_past)
  h_y_future_given_y_past <- compute_conditional_entropy(y_future, y_past)
  
  # H(Y_future | Y_past, X_past)
  h_y_future_given_both <- compute_conditional_entropy(y_future, cbind(y_past, x_past))
  
  # Transfer entropy
  te <- h_y_future_given_y_past - h_y_future_given_both
  
  return(max(0, te))  # Transfer entropy should be non-negative
}

#' Create Lagged Variables
#'
#' Create lagged versions of a time series.
#'
#' @param x Time series
#' @param max_lag Maximum lag
#' @return Matrix of lagged variables
#' @keywords internal
create_lags <- function(x, max_lag) {
  n <- length(x)
  
  # Ensure max_lag is valid
  if (max_lag >= n) {
    max_lag <- max(1, n - 1)
  }
  
  lags <- matrix(NA, n - max_lag, max_lag)
  
  for (i in 1:max_lag) {
    lags[, i] <- x[i:(n - max_lag + i - 1)]
  }
  
  return(lags)
}

#' Compute Conditional Entropy
#'
#' Compute conditional entropy H(X|Y).
#'
#' @param x Variable X
#' @param y Variable Y (can be matrix)
#' @return Conditional entropy
#' @keywords internal
compute_conditional_entropy <- function(x, y) {
  # Create contingency table
  if (is.matrix(y)) {
    # For multivariate Y, combine into single factor
    y_combined <- apply(y, 1, paste, collapse = "_")
  } else {
    y_combined <- y
  }
  
  # Compute joint distribution
  joint_dist <- table(x, y_combined) / length(x)
  
  # Compute conditional entropy
  ce <- 0
  for (j in 1:ncol(joint_dist)) {
    p_y <- sum(joint_dist[, j])
    if (p_y > 0) {
      for (i in 1:nrow(joint_dist)) {
        if (joint_dist[i, j] > 0) {
          ce <- ce - joint_dist[i, j] * log2(joint_dist[i, j] / p_y)
        }
      }
    }
  }
  
  return(ce)
}

#' Community Detection
#'
#' Detect communities in functional connectivity network.
#'
#' @param adjacency_matrix Adjacency matrix
#' @param method Community detection method ("louvain", "girvan_newman", "label_propagation")
#' @param ... Additional arguments
#' @return Community detection results
#' @export
community_detection <- function(adjacency_matrix, method = "louvain", ...) {
  message("Detecting communities")
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for community detection")
  }
  
  # Create graph
  abs_matrix <- abs(adjacency_matrix)
  g <- igraph::graph_from_adjacency_matrix(abs_matrix, mode = "undirected", weighted = TRUE)
  
  # Detect communities
  if (method == "louvain") {
    communities <- igraph::cluster_louvain(g)
  } else if (method == "girvan_newman") {
    communities <- igraph::cluster_edge_betweenness(g)
  } else if (method == "label_propagation") {
    communities <- igraph::cluster_label_prop(g)
  } else {
    stop("Unknown community detection method: ", method)
  }
  
  return(list(
    communities = communities,
    membership = igraph::membership(communities),
    modularity = igraph::modularity(g, igraph::membership(communities)),
    method = method
  ))
}

#' Network Visualization
#'
#' Create network visualization.
#'
#' @param adjacency_matrix Adjacency matrix
#' @param layout Layout method ("fr", "kk", "circle", "random")
#' @param vertex_size Vertex size
#' @param edge_width Edge width scaling
#' @param ... Additional arguments
#' @return Network plot
#' @export
network_visualization <- function(adjacency_matrix, layout = "fr", 
                                vertex_size = 5, edge_width = 1, ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network visualization")
  }
  
  # Create graph
  abs_matrix <- abs(adjacency_matrix)
  g <- igraph::graph_from_adjacency_matrix(abs_matrix, mode = "undirected", weighted = TRUE)
  
  # Set layout
  if (layout == "fr") {
    coords <- igraph::layout_with_fr(g)
  } else if (layout == "kk") {
    coords <- igraph::layout_with_kk(g)
  } else if (layout == "circle") {
    coords <- igraph::layout_in_circle(g)
  } else if (layout == "random") {
    coords <- igraph::layout_randomly(g)
  } else {
    stop("Unknown layout method: ", layout)
  }
  
  # Create plot
  plot(g, layout = coords, vertex.size = vertex_size, 
       edge.width = edge_width * abs(igraph::E(g)$weight), ...)
  
  return(list(graph = g, layout = coords))
}

#' Graph Metrics
#'
#' Compute various graph metrics for network analysis.
#'
#' @param adjacency_matrix Adjacency matrix
#' @param ... Additional arguments
#' @return List of graph metrics
#' @export
graph_metrics <- function(adjacency_matrix, ...) {
  if (is.list(adjacency_matrix) && "connectivity_matrix" %in% names(adjacency_matrix)) {
    adjacency_matrix <- adjacency_matrix$connectivity_matrix
  }
  
  n_nodes <- nrow(adjacency_matrix)
  
  # Basic metrics
  degrees <- rowSums(adjacency_matrix != 0)
  density <- sum(adjacency_matrix != 0) / (n_nodes * (n_nodes - 1))
  
  # Centrality measures
  betweenness <- compute_betweenness_centrality(adjacency_matrix)
  closeness <- compute_closeness_centrality(adjacency_matrix)
  
  # Clustering coefficient
  clustering <- compute_simple_clustering(adjacency_matrix)
  
  return(list(
    degree = degrees,
    density = density,
    betweenness = betweenness,
    closeness = closeness,
    clustering = clustering,
    n_nodes = n_nodes
  ))
}

#' Compute Betweenness Centrality
#'
#' @param adjacency_matrix Adjacency matrix
#' @return Betweenness centrality values
#' @keywords internal
compute_betweenness_centrality <- function(adjacency_matrix) {
  n_nodes <- nrow(adjacency_matrix)
  betweenness <- numeric(n_nodes)
  
  # Simplified betweenness centrality
  for (i in 1:n_nodes) {
    centrality <- 0
    for (j in 1:n_nodes) {
      for (k in 1:n_nodes) {
        if (i != j && j != k && i != k) {
          # Check if i is on shortest path between j and k
          if (adjacency_matrix[j, i] != 0 && adjacency_matrix[i, k] != 0) {
            centrality <- centrality + 1
          }
        }
      }
    }
    betweenness[i] <- centrality
  }
  
  return(betweenness)
}

#' Compute Closeness Centrality
#'
#' @param adjacency_matrix Adjacency matrix
#' @return Closeness centrality values
#' @keywords internal
compute_closeness_centrality <- function(adjacency_matrix) {
  n_nodes <- nrow(adjacency_matrix)
  closeness <- numeric(n_nodes)
  
  for (i in 1:n_nodes) {
    # Compute distances from node i to all other nodes
    distances <- rep(Inf, n_nodes)
    distances[i] <- 0
    
    # Simple breadth-first search
    queue <- i
    visited <- logical(n_nodes)
    visited[i] <- TRUE
    
    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]
      
      for (j in 1:n_nodes) {
        if (!visited[j] && adjacency_matrix[current, j] != 0) {
          distances[j] <- distances[current] + 1
          visited[j] <- TRUE
          queue <- c(queue, j)
        }
      }
    }
    
    # Closeness centrality is inverse of average distance
    avg_distance <- mean(distances[distances != Inf])
    closeness[i] <- if (avg_distance == Inf) 0 else 1 / avg_distance
  }
  
  return(closeness)
}

#' Extract Connections from Thresholded Matrix
#'
#' Convert thresholded adjacency matrix to edge list.
#'
#' @param thresholded_matrix Thresholded adjacency matrix
#' @return Data frame with columns: from, to, weight
#' @keywords internal
extract_connections <- function(thresholded_matrix) {
  n_nodes <- nrow(thresholded_matrix)
  
  # Find non-zero entries (upper triangle only to avoid duplicates)
  connections <- data.frame(
    from = integer(0),
    to = integer(0),
    weight = numeric(0)
  )
  
  for (i in 1:(n_nodes-1)) {
    for (j in (i+1):n_nodes) {
      if (thresholded_matrix[i, j] != 0) {
        connections <- rbind(connections, data.frame(
          from = i,
          to = j,
          weight = thresholded_matrix[i, j]
        ))
      }
    }
  }
  
  return(connections)
} 