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
  
  # Adjust parameters for small datasets
  n_samples <- nrow(data)
  if (n_samples < n_neighbors) {
    n_neighbors <- max(2, n_samples - 1)
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
  
  # Adjust perplexity for small datasets
  n_samples <- nrow(data)
  if (n_samples < perplexity * 3) {
    perplexity <- max(1, floor(n_samples / 3))
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

#' Anomaly Detection via Statistical Methods
#'
#' Detect anomalies using statistical outlier detection methods.
#'
#' @param data Matrix or data frame (cells x features)
#' @param method Detection method ("mahalanobis", "iqr", "zscore")
#' @param threshold Threshold for anomaly detection (default: 0.05)
#' @param ... Additional arguments
#' @return Vector of anomaly scores
#' @export
anomaly_detection <- function(data, method = "mahalanobis", threshold = 0.05, ...) {
  x <- as.matrix(data)
  n_samples <- nrow(x)
  
  if (method == "mahalanobis") {
    # Mahalanobis distance-based anomaly detection
    if (n_samples > ncol(x)) {
      # Calculate Mahalanobis distance
      mu <- colMeans(x, na.rm = TRUE)
      sigma <- cov(x, use = "pairwise.complete.obs")
      
      # Handle singular covariance matrix
      if (any(is.na(sigma)) || det(sigma) == 0) {
        # Fallback to Euclidean distance
        scores <- apply(x, 1, function(row) {
          sqrt(sum((row - mu)^2, na.rm = TRUE))
        })
      } else {
        scores <- apply(x, 1, function(row) {
          mahalanobis(row, mu, sigma)
        })
      }
    } else {
      # Fallback to Euclidean distance for small sample sizes
      mu <- colMeans(x, na.rm = TRUE)
      scores <- apply(x, 1, function(row) {
        sqrt(sum((row - mu)^2, na.rm = TRUE))
      })
    }
    
  } else if (method == "iqr") {
    # IQR-based outlier detection
    scores <- numeric(n_samples)
    for (i in 1:n_samples) {
      sample_scores <- numeric(ncol(x))
      for (j in 1:ncol(x)) {
        col_data <- x[, j]
        q1 <- quantile(col_data, 0.25, na.rm = TRUE)
        q3 <- quantile(col_data, 0.75, na.rm = TRUE)
        iqr <- q3 - q1
        if (iqr > 0) {
          sample_scores[j] <- abs(x[i, j] - median(col_data, na.rm = TRUE)) / iqr
        }
      }
      scores[i] <- mean(sample_scores, na.rm = TRUE)
    }
    
  } else if (method == "zscore") {
    # Z-score based anomaly detection
    scores <- numeric(n_samples)
    for (i in 1:n_samples) {
      sample_scores <- numeric(ncol(x))
      for (j in 1:ncol(x)) {
        col_data <- x[, j]
        mean_val <- mean(col_data, na.rm = TRUE)
        sd_val <- sd(col_data, na.rm = TRUE)
        if (sd_val > 0) {
          sample_scores[j] <- abs(x[i, j] - mean_val) / sd_val
        }
      }
      scores[i] <- mean(sample_scores, na.rm = TRUE)
    }
    
  } else {
    stop("Unknown method. Use 'mahalanobis', 'iqr', or 'zscore'")
  }
  
  # Normalize scores to [0, 1] range
  if (max(scores, na.rm = TRUE) > min(scores, na.rm = TRUE)) {
    scores <- (scores - min(scores, na.rm = TRUE)) / 
              (max(scores, na.rm = TRUE) - min(scores, na.rm = TRUE))
  }
  
  return(scores)
} 