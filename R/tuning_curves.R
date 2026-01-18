#' Tuning Curve Analysis
#'
#' Fit and analyze neural tuning curves for various stimulus dimensions.
#' Supports orientation/direction tuning (von Mises), contrast tuning (sigmoid),
#' spatial tuning (place fields), and general parametric tuning.
#'
#' @name tuning_curves
NULL

#' Fit Orientation Tuning Curves
#'
#' Fit von Mises (circular Gaussian) tuning curves to orientation-selective responses.
#'
#' @param responses Matrix of responses (cells x orientations) or (cells x trials x orientations)
#' @param orientations Vector of stimulus orientations in degrees (0-360 or 0-180)
#' @param method Fitting method: "vonmises", "gaussian", "halfwave_rectified"
#' @param bootstrap_n Number of bootstrap iterations for confidence intervals (0 = none)
#'
#' @return Data frame with tuning parameters per cell:
#' \itemize{
#'   \item{preferred_orientation}{Peak orientation}
#'   \item{tuning_width}{Width (FWHM or kappa)}
#'   \item{amplitude}{Response amplitude}
#'   \item{baseline}{Baseline response}
#'   \item{OSI}{Orientation selectivity index}
#'   \item{DSI}{Direction selectivity index (if applicable)}
#'   \item{r_squared}{Goodness of fit}
#' }
#'
#' @keywords internal
#' @note Use \code{\link{RunTuning}} for CaExperiment workflow
#'
#' @examples
#' \dontrun{
#' # responses: cells x orientations
#' orientations <- seq(0, 330, by = 30)  # 12 orientations
#' tuning <- fit_orientation_tuning(responses, orientations)
#'
#' # Plot tuning curve for cell 1
#' plot_tuning_curve(tuning, cell = 1)
#' }
fit_orientation_tuning <- function(responses, orientations,
                                   method = c("vonmises", "gaussian", "halfwave_rectified"),
                                   bootstrap_n = 0) {
  method <- match.arg(method)

  # Handle 3D input (cells x trials x orientations)
  if (length(dim(responses)) == 3) {
    # Average over trials
    responses <- apply(responses, c(1, 3), mean)
  }

  n_cells <- nrow(responses)
  n_ori <- length(orientations)

  # Convert to radians
  ori_rad <- orientations * pi / 180

  results <- data.frame(
    cell_id = seq_len(n_cells),
    preferred_orientation = NA_real_,
    tuning_width = NA_real_,
    amplitude = NA_real_,
    baseline = NA_real_,
    OSI = NA_real_,
    DSI = NA_real_,
    gOSI = NA_real_,
    r_squared = NA_real_,
    p_value = NA_real_
  )

  for (i in seq_len(n_cells)) {
    resp <- responses[i, ]

    # Compute selectivity indices
    results$OSI[i] <- compute_osi(resp, orientations)
    results$DSI[i] <- compute_dsi(resp, orientations)
    results$gOSI[i] <- compute_gosi(resp, orientations)

    # Fit tuning curve
    fit <- tryCatch({
      switch(method,
             vonmises = .fit_vonmises(resp, ori_rad),
             gaussian = .fit_gaussian_tuning(resp, orientations),
             halfwave_rectified = .fit_halfwave(resp, ori_rad))
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      results$preferred_orientation[i] <- fit$preferred * 180 / pi
      results$tuning_width[i] <- fit$width
      results$amplitude[i] <- fit$amplitude
      results$baseline[i] <- fit$baseline
      results$r_squared[i] <- fit$r_squared

      # Permutation test for significance
      if (bootstrap_n > 0) {
        null_r2 <- replicate(bootstrap_n, {
          shuffled <- sample(resp)
          fit_null <- tryCatch(.fit_vonmises(shuffled, ori_rad), error = function(e) list(r_squared = 0))
          fit_null$r_squared
        })
        results$p_value[i] <- mean(null_r2 >= fit$r_squared)
      }
    }
  }

  # Wrap preferred orientation to valid range
  results$preferred_orientation <- results$preferred_orientation %% 360

  class(results) <- c("tuning_result", "data.frame")
  attr(results, "orientations") <- orientations
  attr(results, "responses") <- responses
  attr(results, "method") <- method

  results
}

#' Fit von Mises function
#' @keywords internal
.fit_vonmises <- function(resp, ori_rad) {
  # Initial estimates
  baseline <- min(resp)
  amplitude <- max(resp) - baseline
  preferred <- ori_rad[which.max(resp)]

  # Grid search for kappa
  kappa_grid <- c(0.5, 1, 2, 4, 8)
  best_r2 <- -Inf
  best_params <- NULL

  for (kappa in kappa_grid) {
    # Optimize
    opt <- tryCatch({
      optim(
        par = c(preferred, kappa, amplitude, baseline),
        fn = function(p) {
          pred <- p[4] + p[3] * exp(p[2] * (cos(2 * (ori_rad - p[1])) - 1))
          sum((resp - pred)^2)
        },
        method = "L-BFGS-B",
        lower = c(-pi, 0.1, 0, -Inf),
        upper = c(pi, 50, Inf, Inf)
      )
    }, error = function(e) NULL)

    if (!is.null(opt)) {
      pred <- opt$par[4] + opt$par[3] * exp(opt$par[2] * (cos(2 * (ori_rad - opt$par[1])) - 1))
      ss_res <- sum((resp - pred)^2)
      ss_tot <- sum((resp - mean(resp))^2)
      r2 <- 1 - ss_res / ss_tot

      if (r2 > best_r2) {
        best_r2 <- r2
        best_params <- opt$par
      }
    }
  }

  if (is.null(best_params)) {
    return(list(preferred = preferred, width = NA, amplitude = amplitude,
                baseline = baseline, r_squared = 0))
  }

  # Convert kappa to width (FWHM in degrees)
  kappa <- best_params[2]
  fwhm <- 2 * acos(1 - log(2) / kappa) * 180 / pi / 2  # Half-width at half-max

  list(
    preferred = best_params[1],
    width = fwhm,
    amplitude = best_params[3],
    baseline = best_params[4],
    kappa = kappa,
    r_squared = best_r2
  )
}

#' Fit Gaussian tuning
#' @keywords internal
.fit_gaussian_tuning <- function(resp, orientations) {
  baseline <- min(resp)
  amplitude <- max(resp) - baseline
  preferred <- orientations[which.max(resp)]
  sigma <- 30

  opt <- tryCatch({
    optim(
      par = c(preferred, sigma, amplitude, baseline),
      fn = function(p) {
        pred <- p[4] + p[3] * exp(-0.5 * ((orientations - p[1]) / p[2])^2)
        sum((resp - pred)^2)
      },
      method = "L-BFGS-B",
      lower = c(min(orientations), 1, 0, -Inf),
      upper = c(max(orientations), 180, Inf, Inf)
    )
  }, error = function(e) NULL)

  if (is.null(opt)) {
    return(list(preferred = preferred * pi / 180, width = NA, amplitude = amplitude,
                baseline = baseline, r_squared = 0))
  }

  pred <- opt$par[4] + opt$par[3] * exp(-0.5 * ((orientations - opt$par[1]) / opt$par[2])^2)
  r2 <- 1 - sum((resp - pred)^2) / sum((resp - mean(resp))^2)

  list(
    preferred = opt$par[1] * pi / 180,
    width = opt$par[2] * 2.355,  # FWHM
    amplitude = opt$par[3],
    baseline = opt$par[4],
    r_squared = r2
  )
}

#' Fit half-wave rectified sinusoid
#' @keywords internal
.fit_halfwave <- function(resp, ori_rad) {
  baseline <- min(resp)
  amplitude <- max(resp) - baseline
  preferred <- ori_rad[which.max(resp)]

  opt <- tryCatch({
    optim(
      par = c(preferred, amplitude, baseline),
      fn = function(p) {
        pred <- p[3] + p[2] * pmax(0, cos(2 * (ori_rad - p[1])))
        sum((resp - pred)^2)
      },
      method = "L-BFGS-B",
      lower = c(-pi, 0, -Inf),
      upper = c(pi, Inf, Inf)
    )
  }, error = function(e) NULL)

  if (is.null(opt)) {
    return(list(preferred = preferred, width = 90, amplitude = amplitude,
                baseline = baseline, r_squared = 0))
  }

  pred <- opt$par[3] + opt$par[2] * pmax(0, cos(2 * (ori_rad - opt$par[1])))
  r2 <- 1 - sum((resp - pred)^2) / sum((resp - mean(resp))^2)

  list(
    preferred = opt$par[1],
    width = 90,  # Fixed for half-wave
    amplitude = opt$par[2],
    baseline = opt$par[3],
    r_squared = r2
  )
}

#' Compute Orientation Selectivity Index (OSI)
#'
#' Classical OSI: (R_pref - R_orth) / (R_pref + R_orth)
#'
#' @param responses Vector of responses at each orientation
#' @param orientations Vector of orientations in degrees
#'
#' @return OSI value (0-1)
#' @keywords internal
compute_osi <- function(responses, orientations) {
  pref_idx <- which.max(responses)
  pref_ori <- orientations[pref_idx]

  # Find orthogonal orientation
  orth_ori <- (pref_ori + 90) %% 180
  if (max(orientations) > 180) {
    orth_ori <- (pref_ori + 90) %% 360
  }

  # Find closest orientation to orthogonal
  orth_idx <- which.min(abs((orientations %% 180) - (orth_ori %% 180)))

  R_pref <- responses[pref_idx]
  R_orth <- responses[orth_idx]

  (R_pref - R_orth) / (R_pref + R_orth + 1e-10)
}

#' Compute Direction Selectivity Index (DSI)
#'
#' DSI: (R_pref - R_null) / (R_pref + R_null)
#'
#' @param responses Vector of responses at each direction
#' @param directions Vector of directions in degrees (0-360)
#'
#' @return DSI value (0-1)
#' @keywords internal
compute_dsi <- function(responses, directions) {
  if (max(directions) <= 180) {
    return(NA_real_)  # Cannot compute DSI without full 360 degrees
  }

  pref_idx <- which.max(responses)
  pref_dir <- directions[pref_idx]

  # Null direction is opposite
  null_dir <- (pref_dir + 180) %% 360
  null_idx <- which.min(abs(directions - null_dir))

  R_pref <- responses[pref_idx]
  R_null <- responses[null_idx]

  (R_pref - R_null) / (R_pref + R_null + 1e-10)
}

#' Compute Global OSI (gOSI)
#'
#' Vector-averaged orientation selectivity (1 - circular variance).
#'
#' @param responses Vector of responses at each orientation
#' @param orientations Vector of orientations in degrees
#'
#' @return gOSI value (0-1)
#' @keywords internal
compute_gosi <- function(responses, orientations) {
  # Convert to radians (double for orientation)
  theta <- orientations * 2 * pi / 180

  # Normalize responses
  R <- responses - min(responses)
  R <- R / (sum(R) + 1e-10)

  # Vector sum
  x <- sum(R * cos(theta))
  y <- sum(R * sin(theta))

  sqrt(x^2 + y^2)
}

#' Fit Contrast Response Function
#'
#' Fit sigmoidal (Naka-Rushton) contrast response function.
#'
#' @param responses Matrix (cells x contrasts) or vector
#' @param contrasts Vector of contrast values (0-1 or 0-100)
#' @param method "naka_rushton" or "hyperbolic_ratio"
#'
#' @return Data frame with parameters: Rmax, C50, n, baseline, r_squared
#' @keywords internal
fit_contrast_response <- function(responses, contrasts,
                                  method = c("naka_rushton", "hyperbolic_ratio")) {
  method <- match.arg(method)

  if (is.vector(responses)) {
    responses <- matrix(responses, nrow = 1)
  }

  # Normalize contrasts to 0-100 if needed
  if (max(contrasts) <= 1) {
    contrasts <- contrasts * 100
  }

  n_cells <- nrow(responses)

  results <- data.frame(
    cell_id = seq_len(n_cells),
    Rmax = NA_real_,
    C50 = NA_real_,
    n = NA_real_,
    baseline = NA_real_,
    r_squared = NA_real_
  )

  for (i in seq_len(n_cells)) {
    resp <- responses[i, ]

    # Initial estimates
    baseline <- resp[which.min(contrasts)]
    Rmax <- max(resp) - baseline
    C50 <- 50
    n <- 2

    opt <- tryCatch({
      optim(
        par = c(Rmax, C50, n, baseline),
        fn = function(p) {
          # Naka-Rushton: R = baseline + Rmax * C^n / (C^n + C50^n)
          pred <- p[4] + p[1] * (contrasts^p[3]) / (contrasts^p[3] + p[2]^p[3])
          sum((resp - pred)^2)
        },
        method = "L-BFGS-B",
        lower = c(0, 1, 0.1, -Inf),
        upper = c(Inf, 100, 10, Inf)
      )
    }, error = function(e) NULL)

    if (!is.null(opt)) {
      pred <- opt$par[4] + opt$par[1] * (contrasts^opt$par[3]) / (contrasts^opt$par[3] + opt$par[2]^opt$par[3])
      r2 <- 1 - sum((resp - pred)^2) / sum((resp - mean(resp))^2)

      results$Rmax[i] <- opt$par[1]
      results$C50[i] <- opt$par[2]
      results$n[i] <- opt$par[3]
      results$baseline[i] <- opt$par[4]
      results$r_squared[i] <- r2
    }
  }

  class(results) <- c("contrast_tuning", "data.frame")
  attr(results, "contrasts") <- contrasts
  attr(results, "responses") <- responses

  results
}

#' Fit Place Field
#'
#' Fit 2D Gaussian place field to spatial firing rate maps.
#'
#' @param rate_map Matrix of firing rates (x_bins x y_bins)
#' @param x_coords Vector of x bin centers
#' @param y_coords Vector of y bin centers
#' @param min_peak Minimum peak rate to consider as place cell
#'
#' @return List with place field parameters
#' @keywords internal
fit_place_field <- function(rate_map, x_coords = NULL, y_coords = NULL,
                            min_peak = 1) {
  n_x <- nrow(rate_map)
  n_y <- ncol(rate_map)

  if (is.null(x_coords)) x_coords <- seq_len(n_x)
  if (is.null(y_coords)) y_coords <- seq_len(n_y)

  # Check for valid place field
  peak_rate <- max(rate_map, na.rm = TRUE)
  if (peak_rate < min_peak) {
    return(list(
      is_place_cell = FALSE,
      peak_rate = peak_rate,
      center_x = NA, center_y = NA,
      sigma_x = NA, sigma_y = NA,
      r_squared = NA
    ))
  }

  # Find peak location
  peak_idx <- which(rate_map == peak_rate, arr.ind = TRUE)[1, ]
  center_x <- x_coords[peak_idx[1]]
  center_y <- y_coords[peak_idx[2]]

  # Create coordinate matrices
  X <- matrix(rep(x_coords, n_y), n_x, n_y)
  Y <- matrix(rep(y_coords, each = n_x), n_x, n_y)

  # Initial sigma estimate
  sigma_init <- (max(x_coords) - min(x_coords)) / 4

  # Fit 2D Gaussian
  opt <- tryCatch({
    optim(
      par = c(center_x, center_y, sigma_init, sigma_init, peak_rate, 0),
      fn = function(p) {
        pred <- p[6] + p[5] * exp(-0.5 * (((X - p[1]) / p[3])^2 + ((Y - p[2]) / p[4])^2))
        sum((rate_map - pred)^2, na.rm = TRUE)
      },
      method = "L-BFGS-B",
      lower = c(min(x_coords), min(y_coords), 0.1, 0.1, 0, 0),
      upper = c(max(x_coords), max(y_coords), diff(range(x_coords)), diff(range(y_coords)), Inf, Inf)
    )
  }, error = function(e) NULL)

  if (is.null(opt)) {
    return(list(
      is_place_cell = FALSE,
      peak_rate = peak_rate,
      center_x = center_x, center_y = center_y,
      sigma_x = NA, sigma_y = NA,
      r_squared = NA
    ))
  }

  # Compute fit quality
  pred <- opt$par[6] + opt$par[5] * exp(-0.5 * (((X - opt$par[1]) / opt$par[3])^2 +
                                                  ((Y - opt$par[2]) / opt$par[4])^2))
  ss_res <- sum((rate_map - pred)^2, na.rm = TRUE)
  ss_tot <- sum((rate_map - mean(rate_map, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- 1 - ss_res / ss_tot

  list(
    is_place_cell = r2 > 0.3,
    peak_rate = opt$par[5],
    center_x = opt$par[1],
    center_y = opt$par[2],
    sigma_x = opt$par[3],
    sigma_y = opt$par[4],
    baseline = opt$par[6],
    r_squared = r2,
    field_size = pi * opt$par[3] * opt$par[4]  # Approximate area
  )
}

#' Compute Spatial Information
#'
#' Skaggs spatial information in bits per spike.
#'
#' @param rate_map Firing rate map (x_bins x y_bins)
#' @param occupancy Occupancy map (time spent in each bin)
#'
#' @return Spatial information in bits/spike
#' @keywords internal
compute_spatial_information <- function(rate_map, occupancy) {
  # Normalize occupancy to probability
  p_x <- occupancy / sum(occupancy, na.rm = TRUE)

  # Mean firing rate
  mean_rate <- sum(rate_map * p_x, na.rm = TRUE)

  if (mean_rate <= 0) return(0)

  # Spatial information: sum(p_x * lambda_x * log2(lambda_x / mean_lambda))
  valid <- rate_map > 0 & p_x > 0
  SI <- sum(p_x[valid] * rate_map[valid] * log2(rate_map[valid] / mean_rate))

  SI / mean_rate  # bits/spike
}

#' Fit Generic Tuning Curve
#'
#' Fit parametric tuning curve to any stimulus dimension.
#'
#' @param responses Response matrix (cells x stimulus_values)
#' @param stimulus_values Vector of stimulus parameter values
#' @param model Tuning model: "gaussian", "sigmoid", "polynomial", "linear"
#' @param circular Whether stimulus dimension is circular (e.g., orientation)
#'
#' @return Data frame with fitted parameters
#' @keywords internal
fit_tuning_curve <- function(responses, stimulus_values,
                             model = c("gaussian", "sigmoid", "polynomial", "linear"),
                             circular = FALSE) {
  model <- match.arg(model)

  if (is.vector(responses)) {
    responses <- matrix(responses, nrow = 1)
  }

  n_cells <- nrow(responses)

  results <- list()

  for (i in seq_len(n_cells)) {
    resp <- responses[i, ]

    fit <- switch(model,
                  gaussian = .fit_generic_gaussian(resp, stimulus_values, circular),
                  sigmoid = .fit_generic_sigmoid(resp, stimulus_values),
                  polynomial = .fit_generic_poly(resp, stimulus_values),
                  linear = .fit_generic_linear(resp, stimulus_values))

    fit$cell_id <- i
    results[[i]] <- fit
  }

  do.call(rbind, lapply(results, as.data.frame))
}

#' Fit generic Gaussian
#' @keywords internal
.fit_generic_gaussian <- function(resp, stim, circular = FALSE) {
  baseline <- min(resp)
  amplitude <- max(resp) - baseline
  preferred <- stim[which.max(resp)]
  sigma <- diff(range(stim)) / 4

  opt <- tryCatch({
    optim(
      par = c(preferred, sigma, amplitude, baseline),
      fn = function(p) {
        if (circular) {
          # Circular distance
          d <- abs(stim - p[1])
          d <- pmin(d, 360 - d)
          pred <- p[4] + p[3] * exp(-0.5 * (d / p[2])^2)
        } else {
          pred <- p[4] + p[3] * exp(-0.5 * ((stim - p[1]) / p[2])^2)
        }
        sum((resp - pred)^2)
      },
      method = "L-BFGS-B",
      lower = c(min(stim), 0.1, 0, -Inf),
      upper = c(max(stim), diff(range(stim)), Inf, Inf)
    )
  }, error = function(e) NULL)

  if (is.null(opt)) {
    return(list(preferred = preferred, width = NA, amplitude = amplitude,
                baseline = baseline, r_squared = 0))
  }

  if (circular) {
    d <- abs(stim - opt$par[1])
    d <- pmin(d, 360 - d)
    pred <- opt$par[4] + opt$par[3] * exp(-0.5 * (d / opt$par[2])^2)
  } else {
    pred <- opt$par[4] + opt$par[3] * exp(-0.5 * ((stim - opt$par[1]) / opt$par[2])^2)
  }
  r2 <- 1 - sum((resp - pred)^2) / sum((resp - mean(resp))^2)

  list(
    model = "gaussian",
    preferred = opt$par[1],
    width = opt$par[2] * 2.355,
    amplitude = opt$par[3],
    baseline = opt$par[4],
    r_squared = r2
  )
}

#' Fit generic sigmoid
#' @keywords internal
.fit_generic_sigmoid <- function(resp, stim) {
  R_min <- min(resp)
  R_max <- max(resp)
  x50 <- median(stim)
  slope <- 1

  opt <- tryCatch({
    optim(
      par = c(R_min, R_max, x50, slope),
      fn = function(p) {
        pred <- p[1] + (p[2] - p[1]) / (1 + exp(-p[4] * (stim - p[3])))
        sum((resp - pred)^2)
      },
      method = "L-BFGS-B",
      lower = c(-Inf, -Inf, min(stim), 0.001),
      upper = c(Inf, Inf, max(stim), 100)
    )
  }, error = function(e) NULL)

  if (is.null(opt)) {
    return(list(model = "sigmoid", R_min = R_min, R_max = R_max,
                x50 = x50, slope = slope, r_squared = 0))
  }

  pred <- opt$par[1] + (opt$par[2] - opt$par[1]) / (1 + exp(-opt$par[4] * (stim - opt$par[3])))
  r2 <- 1 - sum((resp - pred)^2) / sum((resp - mean(resp))^2)

  list(
    model = "sigmoid",
    R_min = opt$par[1],
    R_max = opt$par[2],
    x50 = opt$par[3],
    slope = opt$par[4],
    r_squared = r2
  )
}

#' Fit polynomial
#' @keywords internal
.fit_generic_poly <- function(resp, stim, degree = 3) {
  fit <- lm(resp ~ poly(stim, degree))
  pred <- predict(fit)
  r2 <- summary(fit)$r.squared

  list(
    model = "polynomial",
    degree = degree,
    coefficients = coef(fit),
    r_squared = r2
  )
}

#' Fit linear
#' @keywords internal
.fit_generic_linear <- function(resp, stim) {
  fit <- lm(resp ~ stim)
  coefs <- coef(fit)
  r2 <- summary(fit)$r.squared

  list(
    model = "linear",
    intercept = coefs[1],
    slope = coefs[2],
    r_squared = r2
  )
}

#' Plot Tuning Curve
#'
#' Visualize fitted tuning curve with data points.
#'
#' @param tuning_result Result from fit_orientation_tuning or fit_contrast_response
#' @param cell Which cell to plot (index or "all")
#' @param show_fit Whether to show fitted curve
#' @param show_ci Whether to show confidence intervals (requires bootstrap)
#'
#' @return ggplot object
#' @keywords internal
plot_tuning_curve <- function(tuning_result, cell = 1, show_fit = TRUE, show_ci = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  responses <- attr(tuning_result, "responses")

  if (inherits(tuning_result, "tuning_result")) {
    # Orientation tuning
    orientations <- attr(tuning_result, "orientations")

    if (cell == "all") {
      cells <- seq_len(nrow(responses))
    } else {
      cells <- cell
    }

    plot_data <- do.call(rbind, lapply(cells, function(i) {
      data.frame(
        cell = i,
        orientation = orientations,
        response = responses[i, ]
      )
    }))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = orientation, y = response)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "Orientation (degrees)", y = "Response") +
      ggplot2::theme_minimal()

    if (show_fit && length(cells) == 1) {
      params <- tuning_result[cell, ]
      if (!is.na(params$preferred_orientation)) {
        ori_fine <- seq(min(orientations), max(orientations), length.out = 100)
        # Reconstruct fitted curve (von Mises)
        pref_rad <- params$preferred_orientation * pi / 180
        ori_rad <- ori_fine * pi / 180
        kappa <- 2  # Approximate from width
        pred <- params$baseline + params$amplitude *
                exp(kappa * (cos(2 * (ori_rad - pref_rad)) - 1))

        fit_data <- data.frame(orientation = ori_fine, response = pred)
        p <- p + ggplot2::geom_line(data = fit_data, color = "red", linewidth = 1)
      }
    }

    if (length(cells) > 1) {
      p <- p + ggplot2::facet_wrap(~cell)
    }

  } else if (inherits(tuning_result, "contrast_tuning")) {
    # Contrast response
    contrasts <- attr(tuning_result, "contrasts")

    plot_data <- data.frame(
      contrast = contrasts,
      response = responses[cell, ]
    )

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = contrast, y = response)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "Contrast (%)", y = "Response") +
      ggplot2::theme_minimal()

    if (show_fit) {
      params <- tuning_result[cell, ]
      if (!is.na(params$Rmax)) {
        c_fine <- seq(min(contrasts), max(contrasts), length.out = 100)
        pred <- params$baseline + params$Rmax * (c_fine^params$n) /
                (c_fine^params$n + params$C50^params$n)

        fit_data <- data.frame(contrast = c_fine, response = pred)
        p <- p + ggplot2::geom_line(data = fit_data, color = "red", linewidth = 1)
      }
    }
  }

  p
}

#' Population Tuning Summary
#'
#' Summarize tuning properties across a population of neurons.
#'
#' @param tuning_result Result from tuning curve fitting
#' @param min_r2 Minimum R-squared to include cell
#'
#' @return Summary statistics
#' @keywords internal
tuning_summary <- function(tuning_result, min_r2 = 0.3) {
  # Filter by fit quality
  good_fits <- tuning_result$r_squared >= min_r2

  if (inherits(tuning_result, "tuning_result")) {
    list(
      n_cells = nrow(tuning_result),
      n_tuned = sum(good_fits, na.rm = TRUE),
      fraction_tuned = mean(good_fits, na.rm = TRUE),
      mean_OSI = mean(tuning_result$OSI[good_fits], na.rm = TRUE),
      sd_OSI = sd(tuning_result$OSI[good_fits], na.rm = TRUE),
      mean_DSI = mean(tuning_result$DSI[good_fits], na.rm = TRUE),
      sd_DSI = sd(tuning_result$DSI[good_fits], na.rm = TRUE),
      mean_width = mean(tuning_result$tuning_width[good_fits], na.rm = TRUE),
      sd_width = sd(tuning_result$tuning_width[good_fits], na.rm = TRUE),
      preferred_distribution = hist(tuning_result$preferred_orientation[good_fits],
                                    breaks = seq(0, 360, 30), plot = FALSE)
    )
  } else if (inherits(tuning_result, "contrast_tuning")) {
    list(
      n_cells = nrow(tuning_result),
      n_tuned = sum(good_fits, na.rm = TRUE),
      fraction_tuned = mean(good_fits, na.rm = TRUE),
      mean_C50 = mean(tuning_result$C50[good_fits], na.rm = TRUE),
      sd_C50 = sd(tuning_result$C50[good_fits], na.rm = TRUE),
      mean_n = mean(tuning_result$n[good_fits], na.rm = TRUE),
      sd_n = sd(tuning_result$n[good_fits], na.rm = TRUE)
    )
  }
}
