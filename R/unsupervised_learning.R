#' K-means Clustering for Calcium Imaging Data
#'
#' Cluster cells or time points using k-means.
#'
#' @param data Matrix or data frame (cells x features)
#' @param centers Number of clusters
#' @param ... Additional arguments to kmeans
#' @return List with cluster assignments and centers
#' @export
kmeans_clustering <- function(data, centers = 3, ...) {
  km <- kmeans(as.matrix(data), centers = centers, ...)
  list(cluster = km$cluster, centers = km$centers)
}

#' Hierarchical Clustering for Calcium Imaging Data
#'
#' Cluster cells or time points using hierarchical clustering.
#'
#' @param data Matrix or data frame (cells x features)
#' @param method Linkage method (default: "ward.D2")
#' @return hclust object
#' @export
hierarchical_clustering <- function(data, method = "ward.D2") {
  d <- dist(as.matrix(data))
  hc <- hclust(d, method = method)
  hc
}

#' UMAP Dimensionality Reduction
#'
#' Reduce dimensionality using UMAP.
#'
#' @param data Matrix or data frame (cells x features)
#' @param n_neighbors Number of neighbors
#' @param min_dist Minimum distance
#' @param n_components Number of output dimensions
#' @param ... Additional arguments to umap
#' @return Matrix of reduced dimensions
#' @export
umap_reduce <- function(data, n_neighbors = 15, min_dist = 0.1, n_components = 2, ...) {
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for UMAP. Install with install.packages('uwot')")
  }
  uwot::umap(as.matrix(data), n_neighbors = n_neighbors, min_dist = min_dist, n_components = n_components, ...)
}

#' t-SNE Dimensionality Reduction
#'
#' Reduce dimensionality using t-SNE.
#'
#' @param data Matrix or data frame (cells x features)
#' @param dims Number of output dimensions
#' @param perplexity Perplexity parameter
#' @param ... Additional arguments to Rtsne
#' @return Matrix of reduced dimensions
#' @export
tsne_reduce <- function(data, dims = 2, perplexity = 30, ...) {
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Package 'Rtsne' is required for t-SNE. Install with install.packages('Rtsne')")
  }
  tsne <- Rtsne::Rtsne(as.matrix(data), dims = dims, perplexity = perplexity, ...)
  tsne$Y
}

#' Autoencoder Dimensionality Reduction (Python Keras)
#'
#' Reduce dimensionality using a simple autoencoder (requires Python keras).
#'
#' @param data Matrix or data frame (cells x features)
#' @param encoding_dim Size of encoding layer
#' @param epochs Number of training epochs
#' @param install_missing Whether to install missing Python dependencies
#' @param verbose Show progress messages
#' @return Matrix of encoded representations
#' @export
autoencoder_reduce <- function(data, encoding_dim = 2, epochs = 20, install_missing = TRUE, verbose = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for autoencoder reduction.")
  }
  # Check Python dependencies
  deps <- manage_python_dependencies(packages = c("keras"), install_missing = install_missing, verbose = verbose)
  if (!deps$keras$available) {
    stop("Python package 'keras' is not available and could not be installed.")
  }
  keras <- reticulate::import("keras")
  np <- reticulate::import("numpy")
  x <- as.matrix(data)
  input_dim <- ncol(x)
  model <- keras$Sequential()
  model$add(keras$layers$InputLayer(input_shape = as.integer(input_dim)))
  model$add(keras$layers$Dense(as.integer(encoding_dim), activation = 'relu'))
  model$add(keras$layers$Dense(as.integer(input_dim), activation = 'linear'))
  model$compile(optimizer = 'adam', loss = 'mse')
  model$fit(x, x, epochs = as.integer(epochs), verbose = as.integer(verbose))
  encoder <- keras$Model(inputs = model$input, outputs = model$layers[[2]]$output)
  encoded <- encoder$predict(x)
  encoded
}

#' Anomaly Detection via Isolation Forest (Python)
#'
#' Detect anomalies using isolation forest (requires Python scikit-learn).
#'
#' @param data Matrix or data frame (cells x features)
#' @param install_missing Whether to install missing Python dependencies
#' @param verbose Show progress messages
#' @return Vector of anomaly scores
#' @export
anomaly_detection <- function(data, install_missing = TRUE, verbose = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for anomaly detection.")
  }
  # Check Python dependencies
  deps <- manage_python_dependencies(packages = c("scikit_learn"), install_missing = install_missing, verbose = verbose)
  if (!deps$scikit_learn$available) {
    stop("Python package 'scikit-learn' is not available and could not be installed.")
  }
  # Note: scikit-learn is installed as scikit-learn but imported as sklearn in Python.
  sklearn <- reticulate::import("sklearn.ensemble")
  x <- as.matrix(data)
  clf <- sklearn$IsolationForest()
  clf$fit(x)
  scores <- clf$decision_function(x)
  scores
} 