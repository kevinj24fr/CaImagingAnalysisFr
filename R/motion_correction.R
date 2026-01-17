#' Motion Correction Functions
#'
#' Functions for correcting motion artifacts in calcium imaging movies.
#'
#' @name motion_correction
NULL

#' Motion Correction
#'
#' Correct motion artifacts in calcium imaging data using template-based
#' or frame-by-frame registration.
#'
#' @param movie 3D array [height x width x frames] or path to TIFF file
#' @param method Registration method ("template", "piecewise", "optical_flow")
#' @param template Reference template (NULL = compute from data)
#' @param template_frames Frames to use for template computation
#' @param max_shift Maximum allowed shift in pixels
#' @param upsample_factor Subpixel precision (1 = pixel, 10 = 0.1 pixel)
#' @param batch_size Process frames in batches (memory management)
#' @param verbose Show progress messages
#' @return List with corrected movie and motion parameters
#'
#' @details
#' Available methods:
#' \describe{
#'   \item{template}{Rigid registration to a template image using cross-correlation}
#'   \item{piecewise}{Non-rigid correction by dividing frame into patches}
#'   \item{optical_flow}{Dense optical flow estimation (requires OpenCV)}
#' }
#'
#' @examples
#' \dontrun{
#' # Load movie
#' movie <- read_tiff_stack("recording.tif")
#'
#' # Basic motion correction
#' result <- motion_correct(movie)
#'
#' # With custom template
#' template <- apply(movie[,,1:100], c(1,2), mean)
#' result <- motion_correct(movie, template = template)
#'
#' # Access results
#' corrected_movie <- result$corrected
#' shifts <- result$shifts
#' }
#'
#' @export
motion_correct <- function(movie,
                           method = c("template", "piecewise", "optical_flow"),
                           template = NULL,
                           template_frames = NULL,
                           max_shift = 20,
                           upsample_factor = 10,
                           batch_size = 500,
                           verbose = TRUE) {

  method <- match.arg(method)

  # Load movie if path provided
  if (is.character(movie)) {
    if (verbose) message("Loading movie from: ", movie)
    movie <- read_tiff_stack(movie, verbose = verbose)
  }

  if (length(dim(movie)) != 3) {
    stop("Movie must be a 3D array [height x width x frames]")
  }

  height <- dim(movie)[1]
  width <- dim(movie)[2]
  n_frames <- dim(movie)[3]

  if (verbose) {
    message("Motion correction: ", height, "x", width, " x ", n_frames, " frames")
    message("Method: ", method)
  }

  # Compute template if not provided
  if (is.null(template)) {
    if (is.null(template_frames)) {
      # Use first 100 frames or 10% of movie
      template_frames <- 1:min(100, ceiling(n_frames * 0.1))
    }
    if (verbose) message("Computing template from frames ", min(template_frames), "-", max(template_frames))
    template <- compute_template(movie, template_frames)
  }

  # Apply motion correction based on method
  result <- switch(method,
    "template" = motion_correct_template(movie, template, max_shift, upsample_factor,
                                         batch_size, verbose),
    "piecewise" = motion_correct_piecewise(movie, template, max_shift, upsample_factor,
                                           batch_size, verbose),
    "optical_flow" = motion_correct_optical_flow(movie, template, verbose)
  )

  result$template <- template
  result$method <- method
  result$parameters <- list(
    max_shift = max_shift,
    upsample_factor = upsample_factor,
    template_frames = template_frames
  )

  return(result)
}

#' Compute Template
#'
#' Compute a reference template from a subset of frames.
#'
#' @param movie 3D array
#' @param frames Frame indices to use
#' @param method "mean" or "median"
#' @return 2D template matrix
#'
#' @export
compute_template <- function(movie, frames = NULL, method = "mean") {
  if (is.null(frames)) {
    frames <- seq_len(dim(movie)[3])
  }

  subset <- movie[, , frames, drop = FALSE]

  if (method == "mean") {
    template <- apply(subset, c(1, 2), mean, na.rm = TRUE)
  } else {
    template <- apply(subset, c(1, 2), median, na.rm = TRUE)
  }

  return(template)
}

#' Template-based Motion Correction
#' @keywords internal
motion_correct_template <- function(movie, template, max_shift, upsample_factor,
                                    batch_size, verbose) {
  height <- dim(movie)[1]
  width <- dim(movie)[2]
  n_frames <- dim(movie)[3]

  # Initialize output
  corrected <- array(0, dim = dim(movie))
  shifts <- matrix(0, nrow = n_frames, ncol = 2)
  colnames(shifts) <- c("dy", "dx")
  correlations <- numeric(n_frames)

  # Pre-compute template FFT
  template_padded <- pad_image(template, max_shift)
  template_fft <- fft(template_padded)

  # Process in batches
  n_batches <- ceiling(n_frames / batch_size)

  for (batch in seq_len(n_batches)) {
    start_frame <- (batch - 1) * batch_size + 1
    end_frame <- min(batch * batch_size, n_frames)

    if (verbose) message("  Processing frames ", start_frame, "-", end_frame)

    for (i in start_frame:end_frame) {
      frame <- movie[, , i]

      # Compute shift using phase correlation
      shift_result <- phase_correlation(frame, template, template_fft,
                                        max_shift, upsample_factor)

      shifts[i, ] <- c(shift_result$dy, shift_result$dx)
      correlations[i] <- shift_result$correlation

      # Apply shift
      corrected[, , i] <- shift_image(frame, shift_result$dy, shift_result$dx)
    }
  }

  if (verbose) {
    message("Motion correction complete")
    message("  Mean shift: dy=", round(mean(abs(shifts[, 1])), 2),
            " dx=", round(mean(abs(shifts[, 2])), 2), " pixels")
    message("  Max shift: dy=", round(max(abs(shifts[, 1])), 2),
            " dx=", round(max(abs(shifts[, 2])), 2), " pixels")
  }

  list(
    corrected = corrected,
    shifts = shifts,
    correlations = correlations
  )
}

#' Phase Correlation for Shift Detection
#'
#' Compute the shift between two images using phase correlation.
#'
#' @param image Target image
#' @param template Reference image
#' @param template_fft Pre-computed FFT of template (optional)
#' @param max_shift Maximum shift to search
#' @param upsample_factor Subpixel precision
#' @return List with dy, dx, and correlation
#'
#' @export
phase_correlation <- function(image, template, template_fft = NULL,
                              max_shift = 20, upsample_factor = 1) {
  # Pad images
  padded_image <- pad_image(image, max_shift)
  padded_template <- pad_image(template, max_shift)

  # Compute FFTs
  image_fft <- fft(padded_image)
  if (is.null(template_fft)) {
    template_fft <- fft(padded_template)
  }

  # Cross-power spectrum
  cross_power <- image_fft * Conj(template_fft)
  cross_power <- cross_power / (Mod(cross_power) + 1e-10)

  # Inverse FFT to get correlation
  correlation <- Re(fft(cross_power, inverse = TRUE)) / length(cross_power)

  # Find peak
  peak_idx <- which.max(correlation)
  dims <- dim(correlation)
  peak_row <- ((peak_idx - 1) %% dims[1]) + 1
  peak_col <- ((peak_idx - 1) %/% dims[1]) + 1

  # Convert to shifts (handle wrap-around)
  dy <- peak_row - 1
  dx <- peak_col - 1
  if (dy > dims[1] / 2) dy <- dy - dims[1]
  if (dx > dims[2] / 2) dx <- dx - dims[2]

  # Subpixel refinement if requested
  if (upsample_factor > 1) {
    subpixel <- subpixel_shift(correlation, peak_row, peak_col, upsample_factor)
    dy <- dy + subpixel$dy
    dx <- dx + subpixel$dx
  }

  # Limit to max_shift
  dy <- max(-max_shift, min(max_shift, dy))
  dx <- max(-max_shift, min(max_shift, dx))

  list(
    dy = dy,
    dx = dx,
    correlation = correlation[peak_row, peak_col]
  )
}

#' Subpixel Shift Refinement
#' @keywords internal
subpixel_shift <- function(correlation, peak_row, peak_col, upsample_factor) {
  # Fit parabola to peak neighborhood for subpixel precision
  dims <- dim(correlation)

  # Get 3x3 neighborhood
  rows <- (peak_row - 1):(peak_row + 1)
  cols <- (peak_col - 1):(peak_col + 1)

  # Handle boundaries
  rows <- ((rows - 1) %% dims[1]) + 1
  cols <- ((cols - 1) %% dims[2]) + 1

  neighborhood <- correlation[rows, cols]

  # Parabolic fit in y direction
  dy <- 0
  if (neighborhood[1, 2] != neighborhood[3, 2]) {
    dy <- 0.5 * (neighborhood[1, 2] - neighborhood[3, 2]) /
      (neighborhood[1, 2] - 2 * neighborhood[2, 2] + neighborhood[3, 2])
  }

  # Parabolic fit in x direction
  dx <- 0
  if (neighborhood[2, 1] != neighborhood[2, 3]) {
    dx <- 0.5 * (neighborhood[2, 1] - neighborhood[2, 3]) /
      (neighborhood[2, 1] - 2 * neighborhood[2, 2] + neighborhood[2, 3])
  }

  # Limit to reasonable range
  dy <- max(-0.5, min(0.5, dy))
  dx <- max(-0.5, min(0.5, dx))

  list(dy = dy, dx = dx)
}

#' Pad Image
#' @keywords internal
pad_image <- function(image, pad_size) {
  height <- nrow(image)
  width <- ncol(image)

  # Find next power of 2 for efficient FFT
  new_height <- 2^ceiling(log2(height + 2 * pad_size))
  new_width <- 2^ceiling(log2(width + 2 * pad_size))

  padded <- matrix(mean(image, na.rm = TRUE), nrow = new_height, ncol = new_width)

  # Center the original image
  start_row <- floor((new_height - height) / 2) + 1
  start_col <- floor((new_width - width) / 2) + 1

  padded[start_row:(start_row + height - 1), start_col:(start_col + width - 1)] <- image

  return(padded)
}

#' Shift Image
#'
#' Apply a subpixel shift to an image using bilinear interpolation.
#'
#' @param image 2D matrix
#' @param dy Vertical shift (positive = down)
#' @param dx Horizontal shift (positive = right)
#' @return Shifted image
#'
#' @export
shift_image <- function(image, dy, dx) {
  height <- nrow(image)
  width <- ncol(image)

  # Integer and fractional parts
  dy_int <- floor(dy)
  dx_int <- floor(dx)
  dy_frac <- dy - dy_int
  dx_frac <- dx - dx_int

  # Create output
  shifted <- matrix(mean(image, na.rm = TRUE), nrow = height, ncol = width)

  # Source coordinates
  src_rows <- 1:height - dy_int
  src_cols <- 1:width - dx_int

  # Valid indices
  valid_rows <- which(src_rows >= 1 & src_rows < height)
  valid_cols <- which(src_cols >= 1 & src_cols < width)

  # Bilinear interpolation weights
  w00 <- (1 - dy_frac) * (1 - dx_frac)
  w01 <- (1 - dy_frac) * dx_frac
  w10 <- dy_frac * (1 - dx_frac)
  w11 <- dy_frac * dx_frac

  for (r in valid_rows) {
    src_r <- src_rows[r]
    for (c in valid_cols) {
      src_c <- src_cols[c]

      # Bilinear interpolation
      val <- w00 * image[src_r, src_c]
      if (src_c + 1 <= width) val <- val + w01 * image[src_r, src_c + 1]
      if (src_r + 1 <= height) val <- val + w10 * image[src_r + 1, src_c]
      if (src_r + 1 <= height && src_c + 1 <= width) {
        val <- val + w11 * image[src_r + 1, src_c + 1]
      }

      shifted[r, c] <- val
    }
  }

  return(shifted)
}

#' Piecewise Motion Correction
#' @keywords internal
motion_correct_piecewise <- function(movie, template, max_shift, upsample_factor,
                                     batch_size, verbose) {
  height <- dim(movie)[1]
  width <- dim(movie)[2]
  n_frames <- dim(movie)[3]

  # Define patches (e.g., 4x4 grid)
  n_patches_y <- 4
  n_patches_x <- 4
  patch_height <- ceiling(height / n_patches_y)
  patch_width <- ceiling(width / n_patches_x)

  if (verbose) message("  Using ", n_patches_y, "x", n_patches_x, " patches")

  # Initialize output
  corrected <- array(0, dim = dim(movie))
  all_shifts <- list()

  for (py in seq_len(n_patches_y)) {
    for (px in seq_len(n_patches_x)) {
      # Define patch boundaries
      row_start <- (py - 1) * patch_height + 1
      row_end <- min(py * patch_height, height)
      col_start <- (px - 1) * patch_width + 1
      col_end <- min(px * patch_width, width)

      # Extract patch from template
      template_patch <- template[row_start:row_end, col_start:col_end]

      # Process each frame
      patch_shifts <- matrix(0, nrow = n_frames, ncol = 2)

      for (i in seq_len(n_frames)) {
        frame_patch <- movie[row_start:row_end, col_start:col_end, i]

        # Compute shift for this patch
        shift_result <- phase_correlation(frame_patch, template_patch,
                                          max_shift = max_shift / 2,
                                          upsample_factor = upsample_factor)

        patch_shifts[i, ] <- c(shift_result$dy, shift_result$dx)

        # Apply shift to patch
        corrected_patch <- shift_image(frame_patch, shift_result$dy, shift_result$dx)
        corrected[row_start:row_end, col_start:col_end, i] <- corrected_patch
      }

      all_shifts[[paste0("patch_", py, "_", px)]] <- patch_shifts
    }
  }

  if (verbose) message("Piecewise motion correction complete")

  list(
    corrected = corrected,
    shifts = all_shifts,
    n_patches = c(n_patches_y, n_patches_x)
  )
}

#' Optical Flow Motion Correction
#' @keywords internal
motion_correct_optical_flow <- function(movie, template, verbose) {
  # This would require OpenCV integration
  # For now, fall back to template-based
  warning("Optical flow method requires OpenCV. Falling back to template-based correction.")
  motion_correct_template(movie, template, max_shift = 20, upsample_factor = 10,
                          batch_size = 500, verbose = verbose)
}

#' Assess Motion Correction Quality
#'
#' Compute quality metrics for motion correction.
#'
#' @param original Original movie (3D array)
#' @param corrected Corrected movie (3D array)
#' @param template Reference template
#' @return List of quality metrics
#'
#' @export
assess_motion_correction <- function(original, corrected, template = NULL) {
  n_frames <- dim(original)[3]

  # Compute mean images
  mean_original <- apply(original, c(1, 2), mean, na.rm = TRUE)
  mean_corrected <- apply(corrected, c(1, 2), mean, na.rm = TRUE)

  # Compute sharpness (variance of Laplacian)
  laplacian_kernel <- matrix(c(0, 1, 0, 1, -4, 1, 0, 1, 0), nrow = 3)
  sharpness_original <- var(as.vector(convolve2d(mean_original, laplacian_kernel)), na.rm = TRUE)
  sharpness_corrected <- var(as.vector(convolve2d(mean_corrected, laplacian_kernel)), na.rm = TRUE)

  # Compute frame-to-frame correlation improvement
  cors_original <- numeric(n_frames - 1)
  cors_corrected <- numeric(n_frames - 1)

  for (i in 1:(n_frames - 1)) {
    cors_original[i] <- cor(as.vector(original[, , i]), as.vector(original[, , i + 1]))
    cors_corrected[i] <- cor(as.vector(corrected[, , i]), as.vector(corrected[, , i + 1]))
  }

  # Correlation with template
  template_cor_original <- NULL
  template_cor_corrected <- NULL
  if (!is.null(template)) {
    template_cor_original <- mean(sapply(seq_len(n_frames), function(i) {
      cor(as.vector(original[, , i]), as.vector(template))
    }))
    template_cor_corrected <- mean(sapply(seq_len(n_frames), function(i) {
      cor(as.vector(corrected[, , i]), as.vector(template))
    }))
  }

  list(
    sharpness_improvement = sharpness_corrected / sharpness_original,
    mean_correlation_original = mean(cors_original, na.rm = TRUE),
    mean_correlation_corrected = mean(cors_corrected, na.rm = TRUE),
    correlation_improvement = mean(cors_corrected, na.rm = TRUE) / mean(cors_original, na.rm = TRUE),
    template_correlation_original = template_cor_original,
    template_correlation_corrected = template_cor_corrected,
    n_frames = n_frames
  )
}

#' 2D Convolution
#' @keywords internal
convolve2d <- function(image, kernel) {
  # Simple 2D convolution
  k_height <- nrow(kernel)
  k_width <- ncol(kernel)
  pad_h <- floor(k_height / 2)
  pad_w <- floor(k_width / 2)

  height <- nrow(image)
  width <- ncol(image)

  # Pad image
  padded <- matrix(0, height + 2 * pad_h, width + 2 * pad_w)
  padded[(pad_h + 1):(pad_h + height), (pad_w + 1):(pad_w + width)] <- image

  # Convolve
  result <- matrix(0, height, width)
  for (i in 1:height) {
    for (j in 1:width) {
      patch <- padded[i:(i + k_height - 1), j:(j + k_width - 1)]
      result[i, j] <- sum(patch * kernel)
    }
  }

  return(result)
}

#' Plot Motion Correction Results
#'
#' Visualize motion correction shifts and quality.
#'
#' @param motion_result Result from motion_correct()
#' @param type Plot type ("shifts", "trajectory", "heatmap")
#' @return ggplot object
#'
#' @export
plot_motion_correction <- function(motion_result, type = "shifts") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }

  shifts <- motion_result$shifts

  if (is.matrix(shifts)) {
    df <- data.frame(
      frame = seq_len(nrow(shifts)),
      dy = shifts[, 1],
      dx = shifts[, 2]
    )

    if (type == "shifts") {
      p <- ggplot2::ggplot(df) +
        ggplot2::geom_line(ggplot2::aes(x = frame, y = dy, color = "Y shift")) +
        ggplot2::geom_line(ggplot2::aes(x = frame, y = dx, color = "X shift")) +
        ggplot2::labs(
          title = "Motion Correction Shifts",
          x = "Frame",
          y = "Shift (pixels)",
          color = "Direction"
        ) +
        ggplot2::theme_minimal()

    } else if (type == "trajectory") {
      p <- ggplot2::ggplot(df, ggplot2::aes(x = dx, y = dy, color = frame)) +
        ggplot2::geom_path() +
        ggplot2::geom_point(size = 0.5) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::labs(
          title = "Motion Trajectory",
          x = "X shift (pixels)",
          y = "Y shift (pixels)",
          color = "Frame"
        ) +
        ggplot2::coord_equal() +
        ggplot2::theme_minimal()
    }

    return(p)
  }

  message("Piecewise shifts not yet supported for plotting")
  return(NULL)
}
