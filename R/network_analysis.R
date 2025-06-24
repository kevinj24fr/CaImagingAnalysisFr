#' Functional Connectivity Analysis for Calcium Imaging
#'
#' Calculate functional connectivity between calcium traces using correlation.
#'
#' @param data Matrix or data frame (cells x time)
#' @param method Correlation method: "pearson", "spearman", or "kendall"
#' @param threshold Threshold for significant connections (default: 0.5)
#' @return Matrix of correlation coefficients
#' @export
functional_connectivity <- function(data, method = "pearson", threshold = 0.5) {
  cor_matrix <- cor(t(as.matrix(data)), method = method)
  # Apply threshold
  cor_matrix[abs(cor_matrix) < threshold] <- 0
  cor_matrix
}

#' Granger Causality Analysis
#'
#' Test for Granger causality between pairs of time series.
#'
#' @param x First time series
#' @param y Second time series
#' @param max_lag Maximum lag to test
#' @param verbose Show progress messages
#' @return List with test statistics and p-values
#' @export
granger_causality <- function(x, y, max_lag = 5, verbose = FALSE) {
  if (!requireNamespace("lmtest", quietly = TRUE)) {
    stop("Package 'lmtest' is required for Granger causality. Install with install.packages('lmtest')")
  }
  if (verbose) message("Testing Granger causality...")
  
  # Create lagged variables
  n <- length(x)
  results <- list()
  
  for (lag in 1:max_lag) {
    # Model 1: y ~ lagged y
    y_lagged <- embed(y, lag + 1)
    model1 <- lm(y_lagged[, 1] ~ y_lagged[, -1])
    
    # Model 2: y ~ lagged y + lagged x
    x_lagged <- embed(x, lag + 1)
    if (nrow(x_lagged) < nrow(y_lagged)) {
      x_lagged <- x_lagged[1:nrow(y_lagged), ]
    }
    model2 <- lm(y_lagged[, 1] ~ y_lagged[, -1] + x_lagged[, -1])
    
    # Granger test
    test_result <- lmtest::waldtest(model1, model2)
    results[[paste0("lag_", lag)]] <- list(
      f_stat = test_result$F[2],
      p_value = test_result$`Pr(>F)`[2]
    )
  }
  results
}

#' Transfer Entropy Analysis
#'
#' Calculate transfer entropy between time series using Python's `pyitlib`.
#'
#' @param x First time series
#' @param y Second time series
#' @param k History length for x
#' @param l History length for y
#' @param install_missing Whether to install missing Python dependencies
#' @param verbose Show progress messages
#' @return Transfer entropy value
#' @export
transfer_entropy <- function(x, y, k = 1, l = 1, install_missing = TRUE, verbose = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for transfer entropy.")
  }
  # Check Python dependencies
  deps <- manage_python_dependencies(packages = c("pyitlib"), install_missing = install_missing, verbose = verbose)
  if (!deps$pyitlib$available) {
    stop("Python package 'pyitlib' is not available and could not be installed.")
  }
  if (verbose) message("Calculating transfer entropy...")
  pyitlib <- reticulate::import("pyitlib")
  te <- pyitlib$information_transfer_entropy(x, y, k, l)
  te
}

#' Graph Metrics for Functional Networks
#'
#' Calculate graph-theoretic metrics for a functional connectivity matrix.
#'
#' @param adj_matrix Adjacency matrix (from functional_connectivity)
#' @param metrics Vector of metrics to calculate
#' @return List of graph metrics
#' @export
graph_metrics <- function(adj_matrix, metrics = c("degree", "clustering", "betweenness")) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for graph metrics. Install with install.packages('igraph')")
  }
  # Create undirected graph
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
  results <- list()
  
  if ("degree" %in% metrics) {
    results$degree <- igraph::degree(g)
  }
  if ("clustering" %in% metrics) {
    results$clustering <- igraph::transitivity(g, type = "local")
  }
  if ("betweenness" %in% metrics) {
    results$betweenness <- igraph::betweenness(g)
  }
  if ("eigenvector" %in% metrics) {
    results$eigenvector <- igraph::eigen_centrality(g)$vector
  }
  results
}

#' Network Visualization
#'
#' Create a network plot from functional connectivity matrix.
#'
#' @param adj_matrix Adjacency matrix
#' @param layout Layout algorithm (default: "fr")
#' @param vertex_size Vertex size
#' @param edge_width Edge width scaling
#' @return ggplot2 network plot
#' @export
plot_network <- function(adj_matrix, layout = "fr", vertex_size = 3, edge_width = 1) {
  if (!requireNamespace("igraph", quietly = TRUE) || !requireNamespace("ggraph", quietly = TRUE)) {
    stop("Packages 'igraph' and 'ggraph' are required for network visualization.")
  }
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
  ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(aes(alpha = weight), width = edge_width) +
    ggraph::geom_node_point(size = vertex_size) +
    ggraph::theme_graph()
} 