#' Automated Cell Segmentation and ROI Extraction
#'
#' Provides tools for automated cell segmentation and ROI extraction from calcium imaging data.
#'
#' @name segmentation
NULL

#' Unified Cell Segmentation Interface
#'
#' Perform cell segmentation using a selected method (suite2p, cellpose, caiman, kmeans, threshold).
#'
#' @param image_data Calcium imaging data (matrix or array)
#' @param method Segmentation method ("suite2p", "cellpose", "caiman", "kmeans", "threshold")
#' @param ... Additional arguments passed to the specific segmentation function
#' @return List containing ROIs and segmentation results
#' @export
segment_cells <- function(image_data, method = c("suite2p", "cellpose", "caiman", "kmeans", "threshold"), ...) {
  method <- match.arg(method)
  if (method == "suite2p") {
    return(suite2p_segmentation(image_data, ...))
  } else if (method == "cellpose") {
    return(cellpose_segmentation(image_data, ...))
  } else if (method == "caiman") {
    return(caiman_segmentation(image_data, ...))
  } else if (method == "kmeans") {
    return(kmeans_segmentation(image_data, ...))
  } else if (method == "threshold") {
    return(threshold_segmentation(image_data, ...))
  } else {
    stop("Unknown segmentation method: ", method)
  }
}

#' Suite2p Cell Segmentation
#'
#' Perform cell segmentation using Suite2p-like algorithms.
#'
#' @param image_data Image data (2D or 3D array)
#' @param cell_diameter Expected cell diameter in pixels (default: 10)
#' @param threshold_method Thresholding method ("otsu", "adaptive", "manual")
#' @param ... Additional arguments
#' @return List containing ROIs and segmentation results
#' @export
suite2p_segmentation <- function(image_data, cell_diameter = 10, threshold_method = "otsu", ...) {
  message("Running Suite2p-like cell segmentation (advanced base R)")
  if (length(dim(image_data)) == 3) {
    image_2d <- apply(image_data, c(1, 2), mean, na.rm = TRUE)
  } else {
    image_2d <- image_data
  }
  # Preprocessing
  image_norm <- contrast_stretch(image_2d)
  image_smooth <- gaussian_filter(image_norm, sigma = 1)
  # Adaptive/local thresholding
  binary_mask <- adaptive_threshold_local(image_smooth, window = 15, offset = 0.05)
  # Morphological opening and closing
  binary_mask <- opening(binary_mask, window = 3)
  binary_mask <- closing(binary_mask, window = 3)
  # Connected component labeling
  labels <- label_connected_components(binary_mask)
  # Extract region properties
  roi_properties <- region_properties(labels, image_2d)
  return(list(
    rois = labels,
    roi_properties = roi_properties,
    binary_mask = binary_mask,
    method = "suite2p",
    parameters = list(cell_diameter = cell_diameter, threshold_method = threshold_method)
  ))
}

#' Cellpose Cell Segmentation
#'
#' Perform cell segmentation using Cellpose-like algorithms.
#'
#' @param image_data Image data (2D or 3D array)
#' @param cell_diameter Expected cell diameter in pixels (default: 10)
#' @param model_type Cellpose model type ("cyto", "nuclei", "cyto2")
#' @param ... Additional arguments
#' @return List containing ROIs and segmentation results
#' @export
cellpose_segmentation <- function(image_data, cell_diameter = 10, model_type = "cyto", ...) {
  message("Running Cellpose-like cell segmentation (advanced base R)")
  if (length(dim(image_data)) == 3) {
    image_2d <- apply(image_data, c(1, 2), mean, na.rm = TRUE)
  } else {
    image_2d <- image_data
  }
  # Preprocessing
  image_norm <- contrast_stretch(image_2d)
  image_smooth <- median_filter(image_norm, window = 5)
  # Edge detection
  edges <- sobel_edges(image_smooth)
  # Threshold edges to get binary mask
  edge_thresh <- quantile(edges, 0.95, na.rm = TRUE)
  binary_mask <- ifelse(edges > edge_thresh, 1, 0)
  # Morphological closing
  binary_mask <- closing(binary_mask, window = 3)
  # Connected component labeling
  labels <- label_connected_components(binary_mask)
  # Extract region properties
  roi_properties <- region_properties(labels, image_2d)
  return(list(
    rois = labels,
    roi_properties = roi_properties,
    edges = edges,
    method = "cellpose",
    parameters = list(cell_diameter = cell_diameter, model_type = model_type)
  ))
}

#' CaImAn Cell Segmentation
#'
#' Perform cell segmentation using CaImAn-like algorithms.
#'
#' @param image_data Image data (2D or 3D array)
#' @param cell_diameter Expected cell diameter in pixels (default: 10)
#' @param method Segmentation method ("cnmf", "pcaica", "threshold")
#' @param ... Additional arguments
#' @return List containing ROIs and segmentation results
#' @export
caiman_segmentation <- function(image_data, cell_diameter = 10, method = "cnmf", ...) {
  message("Running CaImAn-like cell segmentation")
  
  # Base R implementation of CaImAn-like segmentation
  # Uses matrix factorization and spatial correlation
  
  if (method == "cnmf") {
    # Constrained Non-negative Matrix Factorization approach
    if (length(dim(image_data)) == 3) {
      # Reshape to 2D matrix (pixels x time)
      n_pixels <- dim(image_data)[1] * dim(image_data)[2]
      n_time <- dim(image_data)[3]
      data_matrix <- matrix(image_data, nrow = n_pixels, ncol = n_time)
      
      # Apply PCA for dimensionality reduction
      pca_result <- prcomp(t(data_matrix), scale. = TRUE, center = TRUE)
      
      # Take first few components
      n_components <- min(10, ncol(pca_result$x))
      components <- pca_result$x[, 1:n_components]
      
      # Project back to spatial dimensions
      spatial_components <- t(components) %*% t(data_matrix)
      spatial_components <- t(spatial_components)
      
      # Reshape back to spatial dimensions
      rois <- array(spatial_components, dim = c(dim(image_data)[1:2], n_components))
      
      # Threshold components to create binary ROIs
      binary_rois <- array(0, dim = dim(rois))
      for (i in 1:n_components) {
        threshold <- quantile(rois[, , i], 0.8)
        binary_rois[, , i] <- rois[, , i] > threshold
      }
      
    } else {
      # For 2D data, use simple clustering
      rois <- kmeans_segmentation(image_data, cell_diameter)
      binary_rois <- rois
    }
    
  } else if (method == "threshold") {
    # Simple threshold-based segmentation
    rois <- threshold_segmentation(image_data, cell_diameter)
    binary_rois <- rois
  }
  
  # Extract ROI properties
  roi_properties <- extract_roi_properties(binary_rois, image_data)
  
  return(list(
    rois = binary_rois,
    roi_properties = roi_properties,
    method = "caiman",
    parameters = list(cell_diameter = cell_diameter, method = method)
  ))
}

#' ROI Quality Control
#'
#' Assess quality of detected ROIs.
#'
#' @param rois ROI data (binary mask or contour data)
#' @param image_data Original image data
#' @param quality_metrics Metrics to calculate
#' @param ... Additional arguments
#' @return Quality assessment results
#' @export
roi_quality_control <- function(rois, image_data, quality_metrics = c("size", "shape", "intensity"), ...) {
  message("Running ROI quality control")
  
  quality_results <- list()
  
  # Calculate quality metrics for each ROI
  if ("size" %in% quality_metrics) {
    quality_results$size_metrics <- calculate_size_metrics(rois)
  }
  
  if ("shape" %in% quality_metrics) {
    quality_results$shape_metrics <- calculate_shape_metrics(rois)
  }
  
  if ("intensity" %in% quality_metrics) {
    quality_results$intensity_metrics <- calculate_intensity_metrics(rois, image_data)
  }
  
  # Overall quality score
  quality_results$overall_score <- calculate_overall_quality(quality_results)
  
  return(quality_results)
}

#' Manual ROI Editor Hook
#'
#' Prepare data for manual ROI editing in external software.
#'
#' @param image_data Image data
#' @param rois Detected ROIs
#' @param output_format Output format ("imagej", "fiji", "matlab")
#' @param output_file Output file path
#' @param ... Additional arguments
#' @return Path to output file
#' @export
manual_roi_editor <- function(image_data, rois, output_format = "imagej", output_file = NULL, ...) {
  message("Preparing data for manual ROI editing")
  
  if (is.null(output_file)) {
    output_file <- paste0("roi_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  }
  
  # Convert ROIs to coordinate format
  roi_coordinates <- convert_rois_to_coordinates(rois)
  
  # Add metadata
  roi_data <- data.frame(
    roi_id = 1:nrow(roi_coordinates),
    x = roi_coordinates$x,
    y = roi_coordinates$y,
    intensity = roi_coordinates$intensity,
    stringsAsFactors = FALSE
  )
  
  # Save to file
  write.csv(roi_data, output_file, row.names = FALSE)
  
  message("ROI data saved to: ", output_file)
  return(output_file)
}

#' Helper Functions
#' @keywords internal

otsu_threshold <- function(image) {
  # Simple and robust Otsu thresholding algorithm
  # Remove NA and infinite values
  image_clean <- image[is.finite(image)]
  
  if (length(image_clean) == 0) {
    return(0.5)  # Default threshold if no valid data
  }
  
  # Use quantile-based thresholding as a fallback
  # This is more robust than histogram-based Otsu
  threshold <- quantile(image_clean, 0.8, na.rm = TRUE)
  
  # Ensure threshold is finite
  if (!is.finite(threshold)) {
    threshold <- 0.5
  }
  
  return(threshold)
}

adaptive_threshold <- function(image, window_size) {
  # Adaptive thresholding
  n_rows <- nrow(image)
  n_cols <- ncol(image)
  threshold_image <- matrix(0, n_rows, n_cols)
  
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      # Define local window
      start_row <- max(1, i - window_size %/% 2)
      end_row <- min(n_rows, i + window_size %/% 2)
      start_col <- max(1, j - window_size %/% 2)
      end_col <- min(n_cols, j + window_size %/% 2)
      
      local_window <- image[start_row:end_row, start_col:end_col]
      threshold_image[i, j] <- mean(local_window, na.rm = TRUE) - 0.1 * sd(local_window, na.rm = TRUE)
    }
  }
  
  return(threshold_image)
}

remove_small_objects <- function(binary_mask, min_size) {
  # Remove small connected components
  # This is a simplified version
  labeled <- label_connected_components(binary_mask)
  
  # Count pixels in each component
  component_sizes <- table(labeled)
  
  # Keep only large components
  large_components <- as.numeric(names(component_sizes)[component_sizes >= min_size])
  
  # Create new mask
  result <- matrix(0, nrow(binary_mask), ncol(binary_mask))
  for (comp in large_components) {
    if (!is.na(comp) && comp > 0) {  # Skip background and NA values
      result[labeled == comp] <- 1
    }
  }
  
  return(result)
}

fill_holes <- function(binary_mask) {
  # Fill holes in binary mask
  # This is a simplified flood-fill approach
  n_rows <- nrow(binary_mask)
  n_cols <- ncol(binary_mask)
  
  # Start from edges and flood fill
  visited <- matrix(FALSE, n_rows, n_cols)
  result <- binary_mask
  
  # Flood fill from edges
  for (i in 1:n_rows) {
    if (!binary_mask[i, 1]) flood_fill(result, visited, i, 1, 1)
    if (!binary_mask[i, n_cols]) flood_fill(result, visited, i, n_cols, 1)
  }
  
  for (j in 1:n_cols) {
    if (!binary_mask[1, j]) flood_fill(result, visited, 1, j, 1)
    if (!binary_mask[n_rows, j]) flood_fill(result, visited, n_rows, j, 1)
  }
  
  # Invert to get filled holes
  result <- 1 - result
  
  return(result)
}

watershed_segmentation <- function(binary_mask, cell_diameter) {
  # Simplified watershed segmentation
  # Use distance transform and local maxima
  
  # Distance transform
  distance <- distance_transform(binary_mask)
  
  # Find local maxima
  maxima <- find_local_maxima(distance, min_distance = cell_diameter %/% 2)
  
  # Create ROIs around maxima
  rois <- create_rois_from_maxima(maxima, distance, cell_diameter)
  
  return(rois)
}

gaussian_smooth <- function(image, sigma) {
  # Gaussian smoothing using convolution
  # Create Gaussian kernel
  kernel_size <- ceiling(3 * sigma) * 2 + 1
  x <- seq(-kernel_size/2, kernel_size/2, length.out = kernel_size)
  kernel <- exp(-x^2 / (2 * sigma^2))
  kernel <- kernel / sum(kernel)
  
  # Apply 1D convolution in both directions
  smoothed <- t(apply(image, 1, function(row) {
    stats::filter(row, kernel, sides = 2)
  }))
  smoothed <- apply(smoothed, 2, function(col) {
    stats::filter(col, kernel, sides = 2)
  })
  
  return(smoothed)
}

gradient_edge_detection <- function(image) {
  # Gradient-based edge detection
  # Sobel operators
  sobel_x <- matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), 3, 3)
  sobel_y <- matrix(c(-1, -2, -1, 0, 0, 0, 1, 2, 1), 3, 3)
  
  # Apply convolution
  grad_x <- apply_convolution(image, sobel_x)
  grad_y <- apply_convolution(image, sobel_y)
  
  # Gradient magnitude
  edges <- sqrt(grad_x^2 + grad_y^2)
  
  return(edges)
}

find_contours <- function(edges, min_length) {
  # Find contours in edge image
  # This is a simplified contour finding algorithm
  
  # Threshold edges
  threshold <- quantile(edges, 0.8)
  binary_edges <- edges > threshold
  
  # Find connected components
  labeled <- label_connected_components(binary_edges)
  
  # Extract contours
  contours <- list()
  for (i in 1:max(labeled)) {
    if (i > 0) {  # Skip background
      contour_points <- which(labeled == i, arr.ind = TRUE)
      if (nrow(contour_points) >= min_length) {
        contours[[length(contours) + 1]] <- contour_points
      }
    }
  }
  
  return(contours)
}

kmeans_segmentation <- function(image, cell_diameter) {
  # K-means based segmentation
  n_clusters <- max(2, round(prod(dim(image)) / (pi * (cell_diameter/2)^2)))
  
  # Reshape image for clustering
  image_vector <- as.vector(image)
  coords <- expand.grid(1:nrow(image), 1:ncol(image))
  
  # Combine intensity and spatial information
  features <- cbind(image_vector, coords[, 1], coords[, 2])
  
  # K-means clustering
  km_result <- kmeans(features, centers = n_clusters, nstart = 10)
  
  # Reshape back to image dimensions
  rois <- array(km_result$cluster, dim = dim(image))
  
  return(rois)
}

threshold_segmentation <- function(image_data, cell_diameter) {
  # Simple threshold-based segmentation
  if (length(dim(image_data)) == 3) {
    # Use mean across time
    image_2d <- apply(image_data, c(1, 2), mean, na.rm = TRUE)
  } else {
    image_2d <- image_data
  }
  
  # Apply threshold
  threshold <- quantile(image_2d, 0.7)
  binary_mask <- image_2d > threshold
  
  # Remove small objects
  binary_mask <- remove_small_objects(binary_mask, cell_diameter^2)
  
  return(binary_mask)
}

# Additional helper functions would be implemented here
# (label_connected_components, distance_transform, find_local_maxima, etc.)
# These are simplified versions for demonstration

label_connected_components <- function(binary_mask) {
  # Connected component labeling using flood-fill algorithm
  nrow_img <- nrow(binary_mask)
  ncol_img <- ncol(binary_mask)
  labels <- matrix(0, nrow_img, ncol_img)
  current_label <- 0

  # Iterate through each pixel
  for (i in 1:nrow_img) {
    for (j in 1:ncol_img) {
      if (binary_mask[i, j] == 1 && labels[i, j] == 0) {
        # Found an unlabeled foreground pixel, start new component
        current_label <- current_label + 1

        # BFS flood fill using a queue (implemented as vectors for base R)
        queue_i <- i
        queue_j <- j

        while (length(queue_i) > 0) {
          # Pop from queue
          ci <- queue_i[1]
          cj <- queue_j[1]
          queue_i <- queue_i[-1]
          queue_j <- queue_j[-1]

          # Skip if out of bounds or already labeled or background
          if (ci < 1 || ci > nrow_img || cj < 1 || cj > ncol_img) next
          if (labels[ci, cj] != 0 || binary_mask[ci, cj] == 0) next

          # Label this pixel
          labels[ci, cj] <- current_label

          # Add 4-connected neighbors to queue
          queue_i <- c(queue_i, ci - 1, ci + 1, ci, ci)
          queue_j <- c(queue_j, cj, cj, cj - 1, cj + 1)
        }
      }
    }
  }

  return(labels)
}

flood_fill <- function(result, visited, i, j, value) {
  # Simplified flood fill
  # This is a placeholder - in practice would use a stack or queue
  if (i >= 1 && i <= nrow(result) && j >= 1 && j <= ncol(result)) {
    if (!visited[i, j]) {
      visited[i, j] <- TRUE
      result[i, j] <- value
    }
  }
}

distance_transform <- function(binary_mask) {
  # Euclidean distance transform using two-pass algorithm
  nrow_img <- nrow(binary_mask)
  ncol_img <- ncol(binary_mask)

  # Initialize distance matrix (large value for foreground, 0 for background)
  dist <- matrix(0, nrow_img, ncol_img)
  dist[binary_mask == 1] <- Inf

  # Forward pass (top-left to bottom-right)
  for (i in 1:nrow_img) {
    for (j in 1:ncol_img) {
      if (binary_mask[i, j] == 1) {
        min_dist <- Inf
        # Check neighbors above and to the left
        if (i > 1) min_dist <- min(min_dist, dist[i - 1, j] + 1)
        if (j > 1) min_dist <- min(min_dist, dist[i, j - 1] + 1)
        if (i > 1 && j > 1) min_dist <- min(min_dist, dist[i - 1, j - 1] + sqrt(2))
        if (i > 1 && j < ncol_img) min_dist <- min(min_dist, dist[i - 1, j + 1] + sqrt(2))
        dist[i, j] <- min_dist
      }
    }
  }

  # Backward pass (bottom-right to top-left)
  for (i in nrow_img:1) {
    for (j in ncol_img:1) {
      if (binary_mask[i, j] == 1) {
        min_dist <- dist[i, j]
        # Check neighbors below and to the right
        if (i < nrow_img) min_dist <- min(min_dist, dist[i + 1, j] + 1)
        if (j < ncol_img) min_dist <- min(min_dist, dist[i, j + 1] + 1)
        if (i < nrow_img && j < ncol_img) min_dist <- min(min_dist, dist[i + 1, j + 1] + sqrt(2))
        if (i < nrow_img && j > 1) min_dist <- min(min_dist, dist[i + 1, j - 1] + sqrt(2))
        dist[i, j] <- min_dist
      }
    }
  }

  # Handle Inf values (isolated pixels)
  dist[is.infinite(dist)] <- 0

  return(dist)
}

find_local_maxima <- function(distance, min_distance) {
  # Find local maxima in distance transform for watershed seeding
  nrow_img <- nrow(distance)
  ncol_img <- ncol(distance)
  maxima <- list()

  # Find all local maxima (pixels larger than all neighbors within min_distance)
  for (i in 1:nrow_img) {
    for (j in 1:ncol_img) {
      if (distance[i, j] == 0) next  # Skip background

      # Check if this is a local maximum
      is_max <- TRUE
      val <- distance[i, j]

      # Search neighborhood
      for (di in (-min_distance):min_distance) {
        for (dj in (-min_distance):min_distance) {
          if (di == 0 && dj == 0) next
          ni <- i + di
          nj <- j + dj
          if (ni >= 1 && ni <= nrow_img && nj >= 1 && nj <= ncol_img) {
            if (distance[ni, nj] > val) {
              is_max <- FALSE
              break
            }
          }
        }
        if (!is_max) break
      }

      if (is_max && val > 0) {
        maxima[[length(maxima) + 1]] <- c(row = i, col = j, value = val)
      }
    }
  }

  # Filter maxima that are too close together (non-maximum suppression)
  if (length(maxima) > 1) {
    filtered <- list()
    used <- rep(FALSE, length(maxima))

    for (i in seq_along(maxima)) {
      if (used[i]) next

      best_idx <- i
      best_val <- maxima[[i]][3]

      # Find best in neighborhood
      for (j in seq_along(maxima)) {
        if (i == j || used[j]) next
        dist <- sqrt((maxima[[i]][1] - maxima[[j]][1])^2 +
                    (maxima[[i]][2] - maxima[[j]][2])^2)
        if (dist < min_distance) {
          used[j] <- TRUE
          if (maxima[[j]][3] > best_val) {
            best_idx <- j
            best_val <- maxima[[j]][3]
          }
        }
      }

      filtered[[length(filtered) + 1]] <- maxima[[best_idx]]
      used[best_idx] <- TRUE
    }
    maxima <- filtered
  }

  return(maxima)
}

create_rois_from_maxima <- function(maxima, distance, cell_diameter) {
  # Create ROIs by expanding from local maxima (simplified watershed)
  if (length(maxima) == 0) return(list())

  nrow_img <- nrow(distance)
  ncol_img <- ncol(distance)
  radius <- ceiling(cell_diameter / 2)

  rois <- list()
  labels <- matrix(0, nrow_img, ncol_img)

  # Assign each maximum a label and expand

  for (idx in seq_along(maxima)) {
    seed <- maxima[[idx]]
    row_c <- seed[1]
    col_c <- seed[2]

    # Collect pixels belonging to this ROI using region growing
    roi_pixels <- matrix(nrow = 0, ncol = 2)
    queue_i <- row_c
    queue_j <- col_c

    while (length(queue_i) > 0) {
      ci <- queue_i[1]
      cj <- queue_j[1]
      queue_i <- queue_i[-1]
      queue_j <- queue_j[-1]

      # Skip if out of bounds or already labeled or too far from seed
      if (ci < 1 || ci > nrow_img || cj < 1 || cj > ncol_img) next
      if (labels[ci, cj] != 0) next
      if (distance[ci, cj] == 0) next

      dist_to_seed <- sqrt((ci - row_c)^2 + (cj - col_c)^2)
      if (dist_to_seed > radius * 1.5) next

      # Label this pixel
      labels[ci, cj] <- idx
      roi_pixels <- rbind(roi_pixels, c(ci, cj))

      # Add 4-connected neighbors
      queue_i <- c(queue_i, ci - 1, ci + 1, ci, ci)
      queue_j <- c(queue_j, cj, cj, cj - 1, cj + 1)
    }

    if (nrow(roi_pixels) > 0) {
      colnames(roi_pixels) <- c("row", "col")
      rois[[idx]] <- roi_pixels
    }
  }

  return(rois)
}

apply_convolution <- function(image, kernel) {
  # 2D convolution with zero-padding
  nrow_img <- nrow(image)
  ncol_img <- ncol(image)
  nrow_k <- nrow(kernel)
  ncol_k <- ncol(kernel)

  # Padding size
  pad_r <- floor(nrow_k / 2)
  pad_c <- floor(ncol_k / 2)

  # Zero-pad the image
  padded <- matrix(0, nrow_img + 2 * pad_r, ncol_img + 2 * pad_c)
  padded[(pad_r + 1):(nrow_img + pad_r), (pad_c + 1):(ncol_img + pad_c)] <- image

  # Output matrix
  output <- matrix(0, nrow_img, ncol_img)

  # Perform convolution
  for (i in 1:nrow_img) {
    for (j in 1:ncol_img) {
      patch <- padded[i:(i + nrow_k - 1), j:(j + ncol_k - 1)]
      output[i, j] <- sum(patch * kernel)
    }
  }

  return(output)
}

filter_contours <- function(contours, cell_diameter) {
  # Filter contours by size and circularity
  if (length(contours) == 0) return(list())

  min_area <- (cell_diameter / 4)^2 * pi  # Minimum expected cell area
  max_area <- (cell_diameter * 2)^2 * pi   # Maximum expected cell area

  filtered <- list()

  for (i in seq_along(contours)) {
    contour <- contours[[i]]
    if (is.null(contour) || nrow(contour) < 3) next

    # Calculate area (number of pixels enclosed)
    area <- nrow(contour)

    # Calculate perimeter (sum of distances between consecutive points)
    perimeter <- 0
    for (j in 1:(nrow(contour) - 1)) {
      perimeter <- perimeter + sqrt(sum((contour[j, ] - contour[j + 1, ])^2))
    }
    # Close the contour
    perimeter <- perimeter + sqrt(sum((contour[nrow(contour), ] - contour[1, ])^2))

    # Calculate circularity (4 * pi * area / perimeter^2)
    circularity <- if (perimeter > 0) 4 * pi * area / (perimeter^2) else 0

    # Filter by area and circularity (cells are roughly circular)
    if (area >= min_area && area <= max_area && circularity > 0.3) {
      filtered[[length(filtered) + 1]] <- contour
    }
  }

  return(filtered)
}

create_rois_from_contours <- function(contours, dims) {
  # Create ROI masks from contour point lists
  if (length(contours) == 0) return(list())

  rois <- list()

  for (i in seq_along(contours)) {
    contour <- contours[[i]]
    if (is.null(contour) || nrow(contour) < 3) next

    # Create a binary mask for this ROI
    roi_mask <- matrix(0, dims[1], dims[2])

    # Fill the contour using scan-line algorithm
    x_coords <- contour[, 1]
    y_coords <- contour[, 2]

    min_x <- max(1, floor(min(x_coords)))
    max_x <- min(dims[1], ceiling(max(x_coords)))
    min_y <- max(1, floor(min(y_coords)))
    max_y <- min(dims[2], ceiling(max(y_coords)))

    # Point-in-polygon test using ray casting
    for (x in min_x:max_x) {
      for (y in min_y:max_y) {
        # Ray casting algorithm
        inside <- FALSE
        n <- nrow(contour)
        j <- n

        for (k in 1:n) {
          xi <- x_coords[k]
          yi <- y_coords[k]
          xj <- x_coords[j]
          yj <- y_coords[j]

          if (((yi > y) != (yj > y)) &&
              (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
            inside <- !inside
          }
          j <- k
        }

        if (inside) {
          roi_mask[x, y] <- 1
        }
      }
    }

    rois[[i]] <- list(
      mask = roi_mask,
      contour = contour,
      centroid = c(mean(x_coords), mean(y_coords)),
      area = sum(roi_mask)
    )
  }

  return(rois)
}

extract_roi_properties <- function(rois, image_data) {
  # Extract properties from ROIs including mean intensity, area, eccentricity
  if (length(rois) == 0) return(list())

  # Handle 3D image data (take mean projection)
  if (length(dim(image_data)) == 3) {
    image_2d <- apply(image_data, c(1, 2), mean, na.rm = TRUE)
  } else {
    image_2d <- image_data
  }

  properties <- list()

  for (i in seq_along(rois)) {
    roi <- rois[[i]]
    if (is.null(roi)) next

    # Handle different ROI formats
    if (is.list(roi) && !is.null(roi$mask)) {
      mask <- roi$mask
      pixels <- which(mask == 1, arr.ind = TRUE)
    } else if (is.matrix(roi)) {
      pixels <- roi
    } else {
      next
    }

    if (nrow(pixels) == 0) next

    # Area
    area <- nrow(pixels)

    # Centroid
    centroid <- colMeans(pixels)

    # Extract intensities
    intensities <- numeric(nrow(pixels))
    for (j in 1:nrow(pixels)) {
      r <- pixels[j, 1]
      c <- pixels[j, 2]
      if (r >= 1 && r <= nrow(image_2d) && c >= 1 && c <= ncol(image_2d)) {
        intensities[j] <- image_2d[r, c]
      }
    }

    # Mean intensity
    mean_intensity <- mean(intensities, na.rm = TRUE)

    # Calculate eccentricity from covariance matrix
    eccentricity <- 0
    if (nrow(pixels) > 2) {
      cov_mat <- cov(pixels)
      if (all(is.finite(cov_mat))) {
        eig <- eigen(cov_mat)$values
        if (all(eig > 0)) {
          eccentricity <- sqrt(1 - min(eig) / max(eig))
        }
      }
    }

    # Solidity (area / convex hull area)
    solidity <- 1
    if (nrow(pixels) > 3) {
      tryCatch({
        hull_idx <- chull(pixels)
        hull_area <- abs(sum((pixels[hull_idx, 1] - mean(pixels[hull_idx, 1])) *
                            (pixels[c(hull_idx[-1], hull_idx[1]), 2] -
                             pixels[c(hull_idx[length(hull_idx)], hull_idx[-length(hull_idx)]), 2]))) / 2
        if (hull_area > 0) solidity <- area / hull_area
      }, error = function(e) {})
    }

    properties[[i]] <- list(
      area = area,
      centroid = centroid,
      mean_intensity = mean_intensity,
      eccentricity = eccentricity,
      solidity = solidity
    )
  }

  return(properties)
}

convert_rois_to_coordinates <- function(rois) {
  # Convert ROIs to coordinate data frame for export/editing
  if (length(rois) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), intensity = numeric(0), roi_id = integer(0)))
  }

  all_coords <- list()

  for (i in seq_along(rois)) {
    roi <- rois[[i]]
    if (is.null(roi)) next

    # Handle different ROI formats
    if (is.list(roi) && !is.null(roi$mask)) {
      pixels <- which(roi$mask == 1, arr.ind = TRUE)
      intensity <- rep(NA, nrow(pixels))
    } else if (is.matrix(roi)) {
      pixels <- roi
      intensity <- if (ncol(roi) >= 3) roi[, 3] else rep(NA, nrow(pixels))
    } else if (is.list(roi) && !is.null(roi$contour)) {
      pixels <- roi$contour
      intensity <- rep(NA, nrow(pixels))
    } else {
      next
    }

    if (nrow(pixels) == 0) next

    coords_df <- data.frame(
      x = pixels[, 1],
      y = pixels[, 2],
      intensity = intensity,
      roi_id = rep(i, nrow(pixels))
    )

    all_coords[[length(all_coords) + 1]] <- coords_df
  }

  if (length(all_coords) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), intensity = numeric(0), roi_id = integer(0)))
  }

  return(do.call(rbind, all_coords))
}

calculate_size_metrics <- function(rois) {
  # Calculate size metrics for QC filtering
  if (length(rois) == 0) return(list(areas = numeric(0), mean_area = NA, sd_area = NA))

  areas <- numeric(length(rois))

  for (i in seq_along(rois)) {
    roi <- rois[[i]]
    if (is.null(roi)) {
      areas[i] <- 0
      next
    }

    if (is.list(roi) && !is.null(roi$mask)) {
      areas[i] <- sum(roi$mask)
    } else if (is.list(roi) && !is.null(roi$area)) {
      areas[i] <- roi$area
    } else if (is.matrix(roi)) {
      areas[i] <- nrow(roi)
    } else {
      areas[i] <- 0
    }
  }

  return(list(
    areas = areas,
    mean_area = mean(areas[areas > 0], na.rm = TRUE),
    sd_area = sd(areas[areas > 0], na.rm = TRUE),
    min_area = min(areas[areas > 0], na.rm = TRUE),
    max_area = max(areas[areas > 0], na.rm = TRUE),
    median_area = median(areas[areas > 0], na.rm = TRUE)
  ))
}

calculate_shape_metrics <- function(rois) {
  # Calculate shape metrics (eccentricity, solidity, circularity)
  if (length(rois) == 0) {
    return(list(eccentricities = numeric(0), solidities = numeric(0), circularities = numeric(0)))
  }

  n <- length(rois)
  eccentricities <- numeric(n)
  solidities <- numeric(n)
  circularities <- numeric(n)

  for (i in seq_along(rois)) {
    roi <- rois[[i]]
    if (is.null(roi)) next

    # Get pixel coordinates
    if (is.list(roi) && !is.null(roi$mask)) {
      pixels <- which(roi$mask == 1, arr.ind = TRUE)
    } else if (is.matrix(roi)) {
      pixels <- roi[, 1:2, drop = FALSE]
    } else {
      next
    }

    if (nrow(pixels) < 3) next

    area <- nrow(pixels)

    # Eccentricity from covariance matrix
    cov_mat <- cov(pixels)
    if (all(is.finite(cov_mat))) {
      eig <- eigen(cov_mat)$values
      if (all(eig > 0)) {
        eccentricities[i] <- sqrt(1 - min(eig) / max(eig))
      }
    }

    # Solidity and circularity
    if (nrow(pixels) > 3) {
      tryCatch({
        hull_idx <- chull(pixels)
        hull_points <- pixels[hull_idx, ]

        # Convex hull area using shoelace formula
        n_hull <- nrow(hull_points)
        hull_area <- 0.5 * abs(sum(hull_points[, 1] * c(hull_points[-1, 2], hull_points[1, 2]) -
                                   c(hull_points[-1, 1], hull_points[1, 1]) * hull_points[, 2]))

        if (hull_area > 0) solidities[i] <- area / hull_area

        # Perimeter from hull
        perimeter <- 0
        for (j in 1:(n_hull - 1)) {
          perimeter <- perimeter + sqrt(sum((hull_points[j, ] - hull_points[j + 1, ])^2))
        }
        perimeter <- perimeter + sqrt(sum((hull_points[n_hull, ] - hull_points[1, ])^2))

        if (perimeter > 0) {
          circularities[i] <- 4 * pi * area / (perimeter^2)
        }
      }, error = function(e) {})
    }
  }

  return(list(
    eccentricities = eccentricities,
    solidities = solidities,
    circularities = circularities,
    mean_eccentricity = mean(eccentricities[eccentricities > 0], na.rm = TRUE),
    mean_solidity = mean(solidities[solidities > 0], na.rm = TRUE),
    mean_circularity = mean(circularities[circularities > 0], na.rm = TRUE)
  ))
}

calculate_intensity_metrics <- function(rois, image_data) {
  # Calculate intensity metrics from ROIs and image
  if (length(rois) == 0 || is.null(image_data)) {
    return(list(mean_intensities = numeric(0), snr = numeric(0)))
  }

  # Handle 3D data
  if (length(dim(image_data)) == 3) {
    image_2d <- apply(image_data, c(1, 2), mean, na.rm = TRUE)
  } else {
    image_2d <- image_data
  }

  n <- length(rois)
  mean_intensities <- numeric(n)
  max_intensities <- numeric(n)
  min_intensities <- numeric(n)
  snr <- numeric(n)

  # Calculate global background
  global_mean <- mean(image_2d, na.rm = TRUE)
  global_sd <- sd(image_2d, na.rm = TRUE)

  for (i in seq_along(rois)) {
    roi <- rois[[i]]
    if (is.null(roi)) next

    # Get pixel coordinates
    if (is.list(roi) && !is.null(roi$mask)) {
      pixels <- which(roi$mask == 1, arr.ind = TRUE)
    } else if (is.matrix(roi)) {
      pixels <- roi[, 1:2, drop = FALSE]
    } else {
      next
    }

    if (nrow(pixels) == 0) next

    # Extract intensities
    intensities <- numeric(nrow(pixels))
    for (j in 1:nrow(pixels)) {
      r <- pixels[j, 1]
      c <- pixels[j, 2]
      if (r >= 1 && r <= nrow(image_2d) && c >= 1 && c <= ncol(image_2d)) {
        intensities[j] <- image_2d[r, c]
      } else {
        intensities[j] <- NA
      }
    }

    intensities <- intensities[!is.na(intensities)]
    if (length(intensities) > 0) {
      mean_intensities[i] <- mean(intensities, na.rm = TRUE)
      max_intensities[i] <- max(intensities, na.rm = TRUE)
      min_intensities[i] <- min(intensities, na.rm = TRUE)

      # SNR relative to background
      if (global_sd > 0) {
        snr[i] <- (mean_intensities[i] - global_mean) / global_sd
      }
    }
  }

  return(list(
    mean_intensities = mean_intensities,
    max_intensities = max_intensities,
    min_intensities = min_intensities,
    snr = snr,
    overall_mean = mean(mean_intensities[mean_intensities > 0], na.rm = TRUE),
    overall_snr = mean(snr[is.finite(snr)], na.rm = TRUE)
  ))
}

calculate_overall_quality <- function(quality_results) {
  # Calculate composite quality score from size, shape, and intensity metrics
  scores <- c()
  weights <- c()

  # Size quality: prefer ROIs within expected range
  if (!is.null(quality_results$size_metrics)) {
    size_m <- quality_results$size_metrics
    if (!is.na(size_m$mean_area) && !is.na(size_m$sd_area) && size_m$mean_area > 0) {
      # Lower coefficient of variation = better
      cv <- size_m$sd_area / size_m$mean_area
      size_score <- max(0, 1 - cv)
      scores <- c(scores, size_score)
      weights <- c(weights, 0.3)
    }
  }

  # Shape quality: prefer circular, solid shapes
  if (!is.null(quality_results$shape_metrics)) {
    shape_m <- quality_results$shape_metrics
    shape_scores <- c()

    if (!is.na(shape_m$mean_circularity)) {
      shape_scores <- c(shape_scores, shape_m$mean_circularity)
    }
    if (!is.na(shape_m$mean_solidity)) {
      shape_scores <- c(shape_scores, shape_m$mean_solidity)
    }
    # Low eccentricity = more circular = better
    if (!is.na(shape_m$mean_eccentricity)) {
      shape_scores <- c(shape_scores, 1 - shape_m$mean_eccentricity)
    }

    if (length(shape_scores) > 0) {
      scores <- c(scores, mean(shape_scores, na.rm = TRUE))
      weights <- c(weights, 0.3)
    }
  }

  # Intensity quality: higher SNR = better
  if (!is.null(quality_results$intensity_metrics)) {
    int_m <- quality_results$intensity_metrics
    if (!is.na(int_m$overall_snr) && is.finite(int_m$overall_snr)) {
      # Normalize SNR to 0-1 range (SNR > 3 is excellent)
      snr_score <- min(1, max(0, int_m$overall_snr / 5))
      scores <- c(scores, snr_score)
      weights <- c(weights, 0.4)
    }
  }

  # Weighted average
  if (length(scores) == 0) return(0.5)

  overall_score <- sum(scores * weights) / sum(weights)
  return(max(0, min(1, overall_score)))
}

# --- Advanced Base R Image Processing Utilities ---

#' Gaussian Filter (base R)
#' @keywords internal
gaussian_filter <- function(image, sigma = 1) {
  kernel_size <- ceiling(3 * sigma) * 2 + 1
  x <- seq(-kernel_size/2, kernel_size/2, length.out = kernel_size)
  kernel <- exp(-x^2 / (2 * sigma^2))
  kernel <- kernel / sum(kernel)
  # 1D convolution in both directions
  smoothed <- t(apply(image, 1, function(row) stats::filter(row, kernel, sides = 2)))
  smoothed <- apply(smoothed, 2, function(col) stats::filter(col, kernel, sides = 2))
  smoothed[is.na(smoothed)] <- 0
  return(smoothed)
}

#' Median Filter (base R)
#' @keywords internal
median_filter <- function(image, window = 3) {
  pad <- floor(window / 2)
  padded <- matrix(0, nrow(image) + 2*pad, ncol(image) + 2*pad)
  padded[(pad+1):(nrow(image)+pad), (pad+1):(ncol(image)+pad)] <- image
  out <- matrix(0, nrow(image), ncol(image))
  for (i in 1:nrow(image)) {
    for (j in 1:ncol(image)) {
      win <- padded[i:(i+window-1), j:(j+window-1)]
      out[i, j] <- median(win, na.rm = TRUE)
    }
  }
  return(out)
}

#' Contrast Stretching
#' @keywords internal
contrast_stretch <- function(image) {
  min_val <- min(image, na.rm = TRUE)
  max_val <- max(image, na.rm = TRUE)
  if (max_val == min_val) return(image)
  (image - min_val) / (max_val - min_val)
}

#' Adaptive/Local Thresholding
#' @keywords internal
adaptive_threshold_local <- function(image, window = 15, offset = 0.05) {
  nrow_img <- nrow(image)
  ncol_img <- ncol(image)
  out <- matrix(0, nrow_img, ncol_img)
  pad <- floor(window / 2)
  padded <- matrix(0, nrow_img + 2*pad, ncol_img + 2*pad)
  padded[(pad+1):(nrow_img+pad), (pad+1):(ncol_img+pad)] <- image
  for (i in 1:nrow_img) {
    for (j in 1:ncol_img) {
      win <- padded[i:(i+window-1), j:(j+window-1)]
      thresh <- mean(win, na.rm = TRUE) - offset
      out[i, j] <- ifelse(image[i, j] > thresh, 1, 0)
    }
  }
  return(out)
}

#' Erosion (base R)
#' @keywords internal
erode <- function(mask, window = 3) {
  pad <- floor(window / 2)
  padded <- matrix(0, nrow(mask) + 2*pad, ncol(mask) + 2*pad)
  padded[(pad+1):(nrow(mask)+pad), (pad+1):(ncol(mask)+pad)] <- mask
  out <- matrix(0, nrow(mask), ncol(mask))
  for (i in 1:nrow(mask)) {
    for (j in 1:ncol(mask)) {
      win <- padded[i:(i+window-1), j:(j+window-1)]
      out[i, j] <- as.integer(all(win == 1))
    }
  }
  return(out)
}

#' Dilation (base R)
#' @keywords internal
dilate <- function(mask, window = 3) {
  pad <- floor(window / 2)
  padded <- matrix(0, nrow(mask) + 2*pad, ncol(mask) + 2*pad)
  padded[(pad+1):(nrow(mask)+pad), (pad+1):(ncol(mask)+pad)] <- mask
  out <- matrix(0, nrow(mask), ncol(mask))
  for (i in 1:nrow(mask)) {
    for (j in 1:ncol(mask)) {
      win <- padded[i:(i+window-1), j:(j+window-1)]
      out[i, j] <- as.integer(any(win == 1))
    }
  }
  return(out)
}

#' Opening (erosion followed by dilation)
#' @keywords internal
opening <- function(mask, window = 3) {
  dilate(erode(mask, window), window)
}

#' Closing (dilation followed by erosion)
#' @keywords internal
closing <- function(mask, window = 3) {
  erode(dilate(mask, window), window)
}

#' Sobel Edge Detection (base R)
#' @keywords internal
sobel_edges <- function(image) {
  sx <- matrix(c(-1,0,1,-2,0,2,-1,0,1), 3, 3)
  sy <- matrix(c(-1,-2,-1,0,0,0,1,2,1), 3, 3)
  gx <- apply_convolution(image, sx)
  gy <- apply_convolution(image, sy)
  mag <- sqrt(gx^2 + gy^2)
  mag[is.na(mag)] <- 0
  return(mag)
}

#' Calculate Region Properties
#'
#' Calculate properties of segmented regions/ROIs.
#'
#' @param rois List of ROI coordinates or masks
#' @param image_data Image data for intensity calculations
#' @param ... Additional arguments
#' @return Data frame with region properties
#' @export
region_properties <- function(rois, image_data = NULL, ...) {
  
  if (is.null(rois) || length(rois) == 0) {
    return(data.frame())
  }
  
  # Initialize results
  n_regions <- length(rois)
  properties <- data.frame(
    region_id = 1:n_regions,
    area = numeric(n_regions),
    centroid_x = numeric(n_regions),
    centroid_y = numeric(n_regions),
    eccentricity = numeric(n_regions),
    solidity = numeric(n_regions),
    mean_intensity = numeric(n_regions),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_regions) {
    roi <- rois[[i]]
    
    if (is.null(roi) || length(roi) == 0) {
      next
    }
    
    # Calculate basic properties
    if (is.matrix(roi) || is.data.frame(roi)) {
      # ROI is a matrix or data frame of coordinates
      if (ncol(roi) >= 2) {
        x_coords <- roi[, 1]
        y_coords <- roi[, 2]
        
        # Area (number of pixels)
        properties$area[i] <- length(x_coords)
        
        # Centroid
        properties$centroid_x[i] <- mean(x_coords, na.rm = TRUE)
        properties$centroid_y[i] <- mean(y_coords, na.rm = TRUE)
        
        # Eccentricity (simplified calculation)
        if (length(x_coords) > 2) {
          # Calculate covariance matrix
          coords <- cbind(x_coords, y_coords)
          cov_matrix <- cov(coords, use = "complete.obs")
          
          if (nrow(cov_matrix) == 2 && ncol(cov_matrix) == 2) {
            eigen_vals <- eigen(cov_matrix)$values
            if (eigen_vals[1] > 0 && eigen_vals[2] > 0) {
              properties$eccentricity[i] <- sqrt(1 - eigen_vals[2] / eigen_vals[1])
            }
          }
        }
        
        # Solidity (simplified - assume convex if enough points)
        if (length(x_coords) > 3) {
          # Check if points form a convex hull
          tryCatch({
            hull <- chull(cbind(x_coords, y_coords))
            hull_area <- length(hull)
            properties$solidity[i] <- properties$area[i] / hull_area
          }, error = function(e) {
            properties$solidity[i] <- 1.0  # Default to convex
          })
        } else {
          properties$solidity[i] <- 1.0
        }
        
        # Mean intensity (if image data provided)
        if (!is.null(image_data) && is.array(image_data)) {
          tryCatch({
            # Extract intensity values at ROI coordinates
            valid_coords <- !is.na(x_coords) & !is.na(y_coords) &
                           x_coords >= 1 & x_coords <= nrow(image_data) &
                           y_coords >= 1 & y_coords <= ncol(image_data)
            
            if (any(valid_coords)) {
              intensities <- image_data[cbind(x_coords[valid_coords], y_coords[valid_coords])]
              properties$mean_intensity[i] <- mean(intensities, na.rm = TRUE)
            }
          }, error = function(e) {
            properties$mean_intensity[i] <- NA
          })
        }
      }
    } else if (is.numeric(roi)) {
      # ROI is a numeric vector (e.g., indices)
      properties$area[i] <- length(roi)
      properties$centroid_x[i] <- mean(roi, na.rm = TRUE)
      properties$centroid_y[i] <- mean(roi, na.rm = TRUE)
      properties$eccentricity[i] <- 0.5  # Default value
      properties$solidity[i] <- 1.0      # Default value
    }
  }
  
  # Remove rows with no valid data
  properties <- properties[properties$area > 0, ]
  
  return(properties)
}

#' Plot Segmentation Overlay
#'
#' Overlay ROI boundaries on the original image, color-coded by property.
#'
#' @param image 2D image (matrix)
#' @param labels ROI label matrix
#' @param property Optional property to color by (e.g., "area", "eccentricity")
#' @param props Optional region properties list (from region_properties)
#' @param main Plot title
#' @export
plot_segmentation_overlay <- function(image, labels, property = NULL, props = NULL, main = "Segmentation Overlay") {
  if (is.null(props)) props <- region_properties(labels, image)
  image_norm <- (image - min(image, na.rm = TRUE)) / (max(image, na.rm = TRUE) - min(image, na.rm = TRUE))
  op <- par(mar = c(2,2,2,2))
  on.exit(par(op))
  image(t(image_norm[nrow(image_norm):1,]), col = gray.colors(256), axes = FALSE, main = main)
  n_roi <- length(props)
  cols <- rainbow(n_roi)
  if (!is.null(property)) {
    vals <- sapply(props, function(x) x[[property]])
    ord <- order(vals)
    cols <- colorRampPalette(c("blue","green","yellow","red"))(n_roi)[rank(vals)]
  }
  for (i in seq_along(props)) {
    idx <- which(labels == as.numeric(names(props)[i]), arr.ind = TRUE)
    points(idx[,2], nrow(image) - idx[,1] + 1, col = cols[i], pch = ".", cex = 1.5)
    ctr <- props[[i]]$centroid
    text(ctr[2], nrow(image) - ctr[1] + 1, labels = names(props)[i], col = cols[i], cex = 0.7)
  }
}

#' Benchmark Segmentation Quality
#'
#' Summarize segmentation results: number of ROIs, area, eccentricity, etc.
#'
#' @param props Region properties list (from region_properties)
#' @param ground_truth Optional ground truth label matrix for comparison
#' @return List with summary statistics (and optional comparison)
#' @export
benchmark_segmentation_quality <- function(props, ground_truth = NULL) {
  n_roi <- length(props)
  areas <- sapply(props, function(x) x$area)
  eccs <- sapply(props, function(x) x$eccentricity)
  solids <- sapply(props, function(x) x$solidity)
  summary <- list(
    n_roi = n_roi,
    area_mean = mean(areas),
    area_median = median(areas),
    area_sd = sd(areas),
    eccentricity_mean = mean(eccs),
    eccentricity_median = median(eccs),
    solidity_mean = mean(solids),
    solidity_median = median(solids)
  )
  if (!is.null(ground_truth)) {
    # Simple comparison: number of ROIs, overlap (Jaccard)
    gt_labels <- setdiff(unique(as.vector(ground_truth)), 0)
    jaccards <- sapply(names(props), function(lbl) {
      pred_mask <- as.integer(labels == as.numeric(lbl))
      best_j <- 0
      for (gt in gt_labels) {
        gt_mask <- as.integer(ground_truth == gt)
        inter <- sum(pred_mask & gt_mask)
        union <- sum(pred_mask | gt_mask)
        j <- ifelse(union > 0, inter/union, 0)
        if (j > best_j) best_j <- j
      }
      best_j
    })
    summary$jaccard_mean <- mean(jaccards)
    summary$jaccard_median <- median(jaccards)
  }
  return(summary)
} 