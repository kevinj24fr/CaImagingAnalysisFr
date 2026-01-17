#' Rcpp Function Wrappers
#'
#' R wrappers for C++ implementations of performance-critical functions.
#' These provide fallback implementations when Rcpp is not available.
#'
#' @name rcpp_wrappers
#' @keywords internal
NULL

#' Check if compiled code is available
#'
#' @return Logical indicating if Rcpp functions are loaded
#' @keywords internal
.rcpp_available <- function() {
  # Check if the compiled functions are available
  tryCatch({
    exists("shift_image_cpp", mode = "function") ||
      is.loaded("shift_image_cpp")
  }, error = function(e) FALSE)
}

#' Shift image (with Rcpp acceleration)
#'
#' Shift an image by the specified displacement using bilinear interpolation.
#' Uses C++ implementation if available, falls back to R otherwise.
#'
#' @param image Numeric matrix
#' @param shift Numeric vector of length 2 (x, y shifts)
#'
#' @return Shifted image matrix
#' @export
shift_image_fast <- function(image, shift) {
  if (.rcpp_available() && exists("shift_image_cpp")) {
    return(shift_image_cpp(image, shift[1], shift[2]))
  }

  # R fallback
  shift_image(image, shift)
}

#' Phase correlation (with Rcpp acceleration)
#'
#' Compute shift between reference and target images.
#' Uses C++ implementation if available.
#'
#' @param ref Reference image matrix
#' @param target Target image matrix
#' @param max_shift Maximum allowed shift
#'
#' @return Numeric vector with shift (x, y)
#' @export
phase_correlation_fast <- function(ref, target, max_shift = 20) {
  if (.rcpp_available() && exists("phase_correlation_cpp")) {
    result <- phase_correlation_cpp(ref, target, as.integer(max_shift))
    return(result[1:2])
  }

  # R fallback
  phase_correlation(ref, target, max_shift = max_shift)
}

#' Connected components labeling (with Rcpp acceleration)
#'
#' Label connected components in binary image.
#' Uses C++ implementation if available.
#'
#' @param binary_image Binary matrix (0/1)
#'
#' @return Integer matrix with component labels
#' @export
connected_components_fast <- function(binary_image) {
  if (.rcpp_available() && exists("connected_components_cpp")) {
    return(connected_components_cpp(as.integer(binary_image)))
  }

  # R fallback using basic flood fill
  .connected_components_r(binary_image)
}

#' R fallback for connected components
#' @keywords internal
.connected_components_r <- function(binary_image) {
  nrow <- nrow(binary_image)
  ncol <- ncol(binary_image)
  labels <- matrix(0L, nrow, ncol)
  current_label <- 0L

  for (i in seq_len(nrow)) {
    for (j in seq_len(ncol)) {
      if (binary_image[i, j] > 0 && labels[i, j] == 0L) {
        current_label <- current_label + 1L
        # BFS
        queue <- list(c(i, j))
        labels[i, j] <- current_label

        while (length(queue) > 0) {
          p <- queue[[1]]
          queue <- queue[-1]
          pi <- p[1]
          pj <- p[2]

          # Check 8-neighbors
          for (di in -1:1) {
            for (dj in -1:1) {
              if (di == 0 && dj == 0) next
              ni <- pi + di
              nj <- pj + dj
              if (ni >= 1 && ni <= nrow && nj >= 1 && nj <= ncol) {
                if (binary_image[ni, nj] > 0 && labels[ni, nj] == 0L) {
                  labels[ni, nj] <- current_label
                  queue <- c(queue, list(c(ni, nj)))
                }
              }
            }
          }
        }
      }
    }
  }

  attr(labels, "n_components") <- current_label
  labels
}

#' Distance transform (with Rcpp acceleration)
#'
#' Compute Euclidean distance transform.
#' Uses C++ implementation if available.
#'
#' @param binary_image Binary matrix
#'
#' @return Numeric matrix with distances
#' @export
distance_transform_fast <- function(binary_image) {
  if (.rcpp_available() && exists("distance_transform_cpp")) {
    return(distance_transform_cpp(as.integer(binary_image)))
  }

  # R fallback (simplified)
  .distance_transform_r(binary_image)
}

#' R fallback for distance transform
#' @keywords internal
.distance_transform_r <- function(binary_image) {
  nrow <- nrow(binary_image)
  ncol <- ncol(binary_image)
  INF <- nrow + ncol

  dist <- matrix(0, nrow, ncol)
  dist[binary_image > 0] <- INF

  # Forward pass
  for (i in seq_len(nrow)) {
    for (j in seq_len(ncol)) {
      if (dist[i, j] > 0) {
        val <- INF
        if (i > 1) val <- min(val, dist[i-1, j] + 1)
        if (j > 1) val <- min(val, dist[i, j-1] + 1)
        if (i > 1 && j > 1) val <- min(val, dist[i-1, j-1] + 1.414)
        if (i > 1 && j < ncol) val <- min(val, dist[i-1, j+1] + 1.414)
        dist[i, j] <- min(dist[i, j], val)
      }
    }
  }

  # Backward pass
  for (i in seq(nrow, 1)) {
    for (j in seq(ncol, 1)) {
      if (dist[i, j] > 0) {
        val <- dist[i, j]
        if (i < nrow) val <- min(val, dist[i+1, j] + 1)
        if (j < ncol) val <- min(val, dist[i, j+1] + 1)
        if (i < nrow && j < ncol) val <- min(val, dist[i+1, j+1] + 1.414)
        if (i < nrow && j > 1) val <- min(val, dist[i+1, j-1] + 1.414)
        dist[i, j] <- val
      }
    }
  }

  dist
}

#' Gaussian blur (with Rcpp acceleration)
#'
#' Gaussian smoothing of image.
#' Uses C++ implementation if available.
#'
#' @param image Numeric matrix
#' @param sigma Standard deviation
#'
#' @return Smoothed image matrix
#' @export
gaussian_blur_fast <- function(image, sigma = 1) {
  if (.rcpp_available() && exists("gaussian_blur_cpp")) {
    return(gaussian_blur_cpp(image, sigma))
  }

  # R fallback using stats::filter
  ksize <- ceiling(6 * sigma)
  if (ksize %% 2 == 0) ksize <- ksize + 1
  half <- ksize %/% 2

  kernel <- dnorm(seq(-half, half), sd = sigma)
  kernel <- kernel / sum(kernel)

  # Apply separable filter
  temp <- t(apply(image, 1, function(row) {
    stats::filter(row, kernel, sides = 2)
  }))
  temp[is.na(temp)] <- image[is.na(temp)]

  result <- apply(temp, 2, function(col) {
    stats::filter(col, kernel, sides = 2)
  })
  result[is.na(result)] <- temp[is.na(result)]

  as.matrix(result)
}

#' Extract trace from ROI (with Rcpp acceleration)
#'
#' Extract fluorescence trace from ROI mask.
#' Uses C++ implementation if available.
#'
#' @param movie 3D array (height x width x frames)
#' @param mask Binary mask matrix
#' @param method Extraction method
#'
#' @return Numeric vector of trace values
#' @export
extract_trace_fast <- function(movie, mask, method = "mean") {
  if (.rcpp_available() && exists("extract_trace_cpp")) {
    return(extract_trace_cpp(movie, mask, method))
  }

  # R fallback
  mask_idx <- which(mask > 0)
  n_frames <- dim(movie)[3]

  vapply(seq_len(n_frames), function(f) {
    vals <- movie[, , f][mask_idx]
    if (method == "mean") {
      mean(vals, na.rm = TRUE)
    } else if (method == "sum") {
      sum(vals, na.rm = TRUE)
    } else if (method == "weighted") {
      weights <- mask[mask_idx]
      sum(vals * weights) / sum(weights)
    } else {
      median(vals, na.rm = TRUE)
    }
  }, numeric(1))
}

#' Local correlation image (with Rcpp acceleration)
#'
#' Compute pixel-wise local correlation for cell detection.
#' Uses C++ implementation if available.
#'
#' @param movie 3D array
#' @param radius Neighborhood radius
#'
#' @return Correlation image matrix
#' @export
local_correlation_fast <- function(movie, radius = 1) {
  if (.rcpp_available() && exists("local_correlation_cpp")) {
    return(local_correlation_cpp(movie, as.integer(radius)))
  }

  # R fallback (simplified, slower)
  dims <- dim(movie)
  nrow <- dims[1]
  ncol <- dims[2]
  n_frames <- dims[3]

  # Flatten time series
  pixels <- matrix(movie, nrow = nrow * ncol, ncol = n_frames)

  # Compute mean correlation with neighbors for each pixel
  corr_img <- matrix(0, nrow, ncol)

  for (i in (radius + 1):(nrow - radius)) {
    for (j in (radius + 1):(ncol - radius)) {
      idx <- (j - 1) * nrow + i
      center <- pixels[idx, ]

      corrs <- numeric()
      for (di in -radius:radius) {
        for (dj in -radius:radius) {
          if (di == 0 && dj == 0) next
          ni <- i + di
          nj <- j + dj
          nidx <- (nj - 1) * nrow + ni
          neighbor <- pixels[nidx, ]
          corrs <- c(corrs, cor(center, neighbor))
        }
      }

      corr_img[i, j] <- mean(corrs, na.rm = TRUE)
    }
  }

  corr_img
}

#' Motion correct movie (with Rcpp acceleration)
#'
#' Full motion correction pipeline.
#' Uses C++ implementation if available.
#'
#' @param movie 3D array
#' @param template Reference template or "mean"
#' @param max_shift Maximum shift
#'
#' @return List with corrected movie and shifts
#' @export
motion_correct_fast <- function(movie, template = "mean", max_shift = 20) {
  dims <- dim(movie)

  # Compute template if needed
  if (identical(template, "mean")) {
    if (.rcpp_available() && exists("compute_template_cpp")) {
      template <- compute_template_cpp(movie, min(100, dims[3]))
    } else {
      template <- apply(movie[, , 1:min(100, dims[3])], c(1, 2), mean)
    }
  }

  # Use C++ if available
  if (.rcpp_available() && exists("motion_correct_movie_cpp")) {
    return(motion_correct_movie_cpp(movie, template, as.integer(max_shift)))
  }

  # R fallback
  motion_correct(movie, template = template, max_shift = max_shift)
}

#' Median filter (with Rcpp acceleration)
#'
#' @param image Numeric matrix
#' @param radius Filter radius
#'
#' @return Filtered image matrix
#' @export
median_filter_fast <- function(image, radius = 1) {
  if (.rcpp_available() && exists("median_filter_cpp")) {
    return(median_filter_cpp(image, as.integer(radius)))
  }

  # R fallback
  nrow <- nrow(image)
  ncol <- ncol(image)
  result <- matrix(0, nrow, ncol)

  for (i in seq_len(nrow)) {
    for (j in seq_len(ncol)) {
      i_range <- max(1, i - radius):min(nrow, i + radius)
      j_range <- max(1, j - radius):min(ncol, j + radius)
      result[i, j] <- median(image[i_range, j_range])
    }
  }

  result
}

#' Threshold image (with Rcpp acceleration)
#'
#' @param image Numeric matrix
#' @param threshold Threshold value (ignored for "otsu")
#' @param method "binary" or "otsu"
#'
#' @return Binary integer matrix
#' @export
threshold_image_fast <- function(image, threshold = 0, method = "binary") {
  if (.rcpp_available() && exists("threshold_image_cpp")) {
    return(threshold_image_cpp(image, threshold, method))
  }

  # R fallback
  if (method == "otsu") {
    # Simple Otsu implementation
    hist_result <- hist(image, breaks = 256, plot = FALSE)
    counts <- hist_result$counts
    mids <- hist_result$mids

    total <- sum(counts)
    sum_total <- sum(counts * mids)

    best_thresh <- mids[1]
    best_var <- 0
    sum_bg <- 0
    weight_bg <- 0

    for (i in seq_along(counts)) {
      weight_bg <- weight_bg + counts[i]
      if (weight_bg == 0) next

      weight_fg <- total - weight_bg
      if (weight_fg == 0) break

      sum_bg <- sum_bg + counts[i] * mids[i]

      mean_bg <- sum_bg / weight_bg
      mean_fg <- (sum_total - sum_bg) / weight_fg

      var_between <- weight_bg * weight_fg * (mean_bg - mean_fg)^2
      if (var_between > best_var) {
        best_var <- var_between
        best_thresh <- mids[i]
      }
    }

    threshold <- best_thresh
  }

  result <- matrix(0L, nrow(image), ncol(image))
  result[image > threshold] <- 1L
  attr(result, "threshold") <- threshold
  result
}
