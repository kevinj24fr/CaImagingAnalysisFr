#' Oscillation and Spectral Analysis
#'
#' Frequency domain analysis of calcium imaging data.
#' Power spectra, coherence, and oscillation detection.
#'
#' @name oscillation_analysis
#' @keywords internal
NULL

#' Compute power spectrum
#'
#' Estimate power spectral density of calcium traces.
#'
#' @param traces Matrix of traces (cells x time) or single trace
#' @param frame_rate Sampling rate in Hz
#' @param method Spectrum method: "periodogram", "welch", "multitaper"
#' @param window_size Window size for Welch method (seconds)
#' @param overlap Overlap fraction for Welch method
#'
#' @return List with frequencies and power values
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' psd <- compute_power_spectrum(traces, frame_rate = 10)
#' plot_spectrum(psd)
#' }
compute_power_spectrum <- function(traces, frame_rate = 10,
                                    method = "welch",
                                    window_size = 10,
                                    overlap = 0.5) {
  if (is.vector(traces)) {
    traces <- matrix(traces, nrow = 1)
  }

  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  if (method == "periodogram") {
    # Simple periodogram
    spectra <- apply(traces, 1, function(x) {
      x_centered <- x - mean(x)
      fft_result <- fft(x_centered)
      power <- Mod(fft_result)^2 / n_time
      power[1:(n_time %/% 2)]
    })

    freqs <- seq(0, frame_rate / 2, length.out = n_time %/% 2)

  } else if (method == "welch") {
    # Welch's method
    window_frames <- round(window_size * frame_rate)
    step <- round(window_frames * (1 - overlap))

    spectra <- apply(traces, 1, function(x) {
      .welch_spectrum(x, window_frames, step, frame_rate)
    })

    freqs <- spectra[[1]]$freq
    spectra <- sapply(spectra, `[[`, "power")

  } else if (method == "multitaper") {
    if (!requireNamespace("multitaper", quietly = TRUE)) {
      warning("multitaper package not available, using Welch method")
      return(compute_power_spectrum(traces, frame_rate, method = "welch"))
    }

    spectra <- apply(traces, 1, function(x) {
      spec <- multitaper::spec.mtm(x, Fs = frame_rate, plot = FALSE)
      list(freq = spec$freq, power = spec$spec)
    })

    freqs <- spectra[[1]]$freq
    spectra <- sapply(spectra, `[[`, "power")
  }

  if (!is.matrix(spectra)) {
    spectra <- matrix(spectra, ncol = 1)
  }

  structure(
    list(
      frequencies = freqs,
      power = t(spectra),
      mean_power = colMeans(spectra),
      method = method,
      frame_rate = frame_rate,
      n_cells = n_cells
    ),
    class = "power_spectrum"
  )
}

.welch_spectrum <- function(x, window_size, step, fs) {
  n <- length(x)
  n_windows <- floor((n - window_size) / step) + 1

  # Hanning window
  window <- 0.5 * (1 - cos(2 * pi * seq(0, window_size - 1) / (window_size - 1)))

  power_sum <- rep(0, window_size %/% 2)

  for (i in seq_len(n_windows)) {
    start <- (i - 1) * step + 1
    segment <- x[start:(start + window_size - 1)]
    segment <- (segment - mean(segment)) * window

    fft_result <- fft(segment)
    power <- Mod(fft_result)^2 / window_size
    power_sum <- power_sum + power[1:(window_size %/% 2)]
  }

  power_avg <- power_sum / n_windows
  freqs <- seq(0, fs / 2, length.out = window_size %/% 2)

  list(freq = freqs, power = power_avg)
}

#' Detect oscillation frequency bands
#'
#' Find dominant oscillation frequencies.
#'
#' @param psd Power spectrum object
#' @param min_prominence Minimum peak prominence
#' @param freq_range Frequency range to search (Hz)
#'
#' @return Data frame of detected oscillations
#' @keywords internal
detect_oscillations <- function(psd, min_prominence = 0.1,
                                 freq_range = c(0.01, 5)) {
  freqs <- psd$frequencies
  power <- psd$mean_power

  # Restrict to frequency range
  valid <- freqs >= freq_range[1] & freqs <= freq_range[2]
  freqs <- freqs[valid]
  power <- power[valid]

  # Find peaks
  peaks <- .find_spectrum_peaks(power, min_prominence)

  if (length(peaks) == 0) {
    return(data.frame(
      frequency = numeric(),
      power = numeric(),
      prominence = numeric()
    ))
  }

  # Compute peak properties
  df <- data.frame(
    frequency = freqs[peaks],
    power = power[peaks],
    prominence = sapply(peaks, function(p) {
      .peak_prominence(power, p)
    })
  )

  # Sort by power
  df <- df[order(df$power, decreasing = TRUE), ]
  df
}

.find_spectrum_peaks <- function(power, min_prominence) {
  n <- length(power)
  peaks <- c()

  for (i in 2:(n - 1)) {
    if (power[i] > power[i - 1] && power[i] > power[i + 1]) {
      prom <- .peak_prominence(power, i)
      if (prom > min_prominence * max(power)) {
        peaks <- c(peaks, i)
      }
    }
  }

  peaks
}

.peak_prominence <- function(power, peak_idx) {
  n <- length(power)
  peak_val <- power[peak_idx]

  # Find lowest point between peak and higher peaks on each side
  left_min <- peak_val
  for (i in (peak_idx - 1):1) {
    if (power[i] > peak_val) break
    left_min <- min(left_min, power[i])
  }

  right_min <- peak_val
  for (i in (peak_idx + 1):n) {
    if (power[i] > peak_val) break
    right_min <- min(right_min, power[i])
  }

  peak_val - max(left_min, right_min)
}

#' Compute band power
#'
#' Get power in specific frequency bands.
#'
#' @param psd Power spectrum object
#' @param bands Named list of frequency bands (e.g., list(slow = c(0.01, 0.1)))
#'
#' @return Data frame with band powers per cell
#' @keywords internal
compute_band_power <- function(psd, bands = NULL) {
  if (is.null(bands)) {
    bands <- list(
      ultra_slow = c(0.001, 0.01),
      slow = c(0.01, 0.1),
      medium = c(0.1, 0.5),
      fast = c(0.5, 2)
    )
  }

  freqs <- psd$frequencies
  power <- psd$power

  band_power <- sapply(bands, function(band) {
    idx <- freqs >= band[1] & freqs <= band[2]
    if (sum(idx) == 0) return(rep(NA, nrow(power)))
    rowSums(power[, idx, drop = FALSE])
  })

  as.data.frame(band_power)
}

#' Compute coherence between cells
#'
#' Cross-spectral coherence analysis.
#'
#' @param trace1 First trace
#' @param trace2 Second trace
#' @param frame_rate Sampling rate
#' @param window_size Window for spectral estimation
#'
#' @return List with frequencies and coherence values
#' @keywords internal
compute_coherence <- function(trace1, trace2, frame_rate = 10,
                               window_size = 256) {
  n <- length(trace1)

  # Detrend
  trace1 <- trace1 - mean(trace1)
  trace2 <- trace2 - mean(trace2)

  # Compute cross-spectrum using Welch
  step <- window_size %/% 2
  n_windows <- floor((n - window_size) / step) + 1

  window <- 0.5 * (1 - cos(2 * pi * seq(0, window_size - 1) / (window_size - 1)))

  pxx <- pxy <- pyy <- rep(0 + 0i, window_size %/% 2)

  for (i in seq_len(n_windows)) {
    start <- (i - 1) * step + 1

    seg1 <- (trace1[start:(start + window_size - 1)] - mean(trace1)) * window
    seg2 <- (trace2[start:(start + window_size - 1)] - mean(trace2)) * window

    fft1 <- fft(seg1)[1:(window_size %/% 2)]
    fft2 <- fft(seg2)[1:(window_size %/% 2)]

    pxx <- pxx + fft1 * Conj(fft1)
    pyy <- pyy + fft2 * Conj(fft2)
    pxy <- pxy + fft1 * Conj(fft2)
  }

  pxx <- Re(pxx / n_windows)
  pyy <- Re(pyy / n_windows)

  # Coherence = |Pxy|^2 / (Pxx * Pyy)
  coherence <- (Mod(pxy)^2) / (pxx * pyy + 1e-10)

  freqs <- seq(0, frame_rate / 2, length.out = window_size %/% 2)

  list(
    frequencies = freqs,
    coherence = Re(coherence),
    cross_spectrum = pxy / n_windows
  )
}

#' Compute coherence matrix
#'
#' Pairwise coherence at specified frequency band.
#'
#' @param traces Matrix of traces (cells x time)
#' @param frame_rate Sampling rate
#' @param freq_band Frequency band for coherence (Hz)
#'
#' @return Coherence matrix
#' @keywords internal
coherence_matrix <- function(traces, frame_rate = 10,
                              freq_band = c(0.1, 0.5)) {
  n_cells <- nrow(traces)

  coh_mat <- matrix(0, n_cells, n_cells)

  for (i in seq_len(n_cells)) {
    for (j in i:n_cells) {
      if (i == j) {
        coh_mat[i, j] <- 1
      } else {
        coh <- compute_coherence(traces[i, ], traces[j, ], frame_rate)
        # Average coherence in band
        idx <- coh$frequencies >= freq_band[1] & coh$frequencies <= freq_band[2]
        coh_mat[i, j] <- mean(coh$coherence[idx])
        coh_mat[j, i] <- coh_mat[i, j]
      }
    }
  }

  coh_mat
}

#' Compute phase locking value
#'
#' Measure phase synchronization between cells.
#'
#' @param trace1 First trace
#' @param trace2 Second trace
#' @param freq_band Frequency band for filtering
#' @param frame_rate Sampling rate
#'
#' @return Phase locking value (0-1)
#' @keywords internal
phase_locking_value <- function(trace1, trace2, freq_band = c(0.1, 0.5),
                                 frame_rate = 10) {
  # Bandpass filter
  trace1_filt <- .bandpass_filter(trace1, freq_band, frame_rate)
  trace2_filt <- .bandpass_filter(trace2, freq_band, frame_rate)

  # Extract instantaneous phase via Hilbert transform
  phase1 <- .instantaneous_phase(trace1_filt)
  phase2 <- .instantaneous_phase(trace2_filt)

  # PLV = |mean(exp(i * (phase1 - phase2)))|
  phase_diff <- phase1 - phase2
  plv <- abs(mean(exp(1i * phase_diff)))

  plv
}

.bandpass_filter <- function(x, freq_band, fs) {
  # Simple FIR bandpass
  n <- length(x)
  fft_x <- fft(x)

  freqs <- seq(0, fs, length.out = n)
  mask <- (freqs >= freq_band[1] & freqs <= freq_band[2]) |
    (freqs >= fs - freq_band[2] & freqs <= fs - freq_band[1])

  fft_x[!mask] <- 0

  Re(fft(fft_x, inverse = TRUE)) / n
}

.instantaneous_phase <- function(x) {
  # Hilbert transform approximation
  n <- length(x)
  fft_x <- fft(x)

  h <- rep(0, n)
  h[1] <- 1
  if (n %% 2 == 0) {
    h[2:(n/2)] <- 2
    h[n/2 + 1] <- 1
  } else {
    h[2:((n + 1)/2)] <- 2
  }

  analytic <- fft(fft_x * h, inverse = TRUE) / n
  Arg(analytic)
}

#' Compare spectra between conditions
#'
#' Statistical comparison of power spectra.
#'
#' @param psd1 Power spectrum for condition 1
#' @param psd2 Power spectrum for condition 2
#' @param method Test method: "permutation", "bootstrap"
#' @param n_perm Number of permutations
#'
#' @return List with frequency-specific p-values
#' @keywords internal
compare_spectra <- function(psd1, psd2, method = "permutation", n_perm = 1000) {
  power1 <- psd1$power
  power2 <- psd2$power
  freqs <- psd1$frequencies

  # Observed difference in mean power
  diff_observed <- colMeans(power2) - colMeans(power1)

  # Permutation test
  n1 <- nrow(power1)
  n2 <- nrow(power2)
  combined <- rbind(power1, power2)

  null_diffs <- replicate(n_perm, {
    perm_idx <- sample(n1 + n2)
    perm1 <- combined[perm_idx[1:n1], , drop = FALSE]
    perm2 <- combined[perm_idx[(n1 + 1):(n1 + n2)], , drop = FALSE]
    colMeans(perm2) - colMeans(perm1)
  })

  # Two-tailed p-values
  p_values <- sapply(seq_along(freqs), function(i) {
    2 * min(
      mean(null_diffs[i, ] >= diff_observed[i]),
      mean(null_diffs[i, ] <= diff_observed[i])
    )
  })

  # FDR correction
  p_adj <- p.adjust(p_values, method = "fdr")

  list(
    frequencies = freqs,
    difference = diff_observed,
    p_values = p_values,
    p_adjusted = p_adj,
    significant = p_adj < 0.05
  )
}

#' Plot power spectrum
#'
#' @param psd Power spectrum object
#' @param log_scale Use log scale for power
#' @param show_cells Show individual cells or just mean
#' @param freq_range Frequency range to display
#'
#' @return ggplot object
#' @keywords internal
plot_spectrum <- function(psd, log_scale = TRUE, show_cells = FALSE,
                          freq_range = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  freqs <- psd$frequencies

  if (!is.null(freq_range)) {
    idx <- freqs >= freq_range[1] & freqs <= freq_range[2]
  } else {
    idx <- rep(TRUE, length(freqs))
  }

  if (show_cells && nrow(psd$power) > 1) {
    df <- data.frame(
      frequency = rep(freqs[idx], nrow(psd$power)),
      power = as.vector(t(psd$power[, idx])),
      cell = factor(rep(seq_len(nrow(psd$power)), each = sum(idx)))
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = frequency, y = power, group = cell)) +
      ggplot2::geom_line(alpha = 0.3)

  } else {
    df <- data.frame(
      frequency = freqs[idx],
      power = psd$mean_power[idx]
    )

    # Add SEM if multiple cells
    if (nrow(psd$power) > 1) {
      df$sem <- apply(psd$power[, idx, drop = FALSE], 2, function(x) sd(x) / sqrt(length(x)))
    }

    p <- ggplot2::ggplot(df, ggplot2::aes(x = frequency, y = power))

    if ("sem" %in% names(df)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = power - sem, ymax = power + sem),
        alpha = 0.3
      )
    }

    p <- p + ggplot2::geom_line()
  }

  if (log_scale) {
    p <- p + ggplot2::scale_y_log10()
  }

  p +
    ggplot2::labs(
      title = "Power Spectrum",
      x = "Frequency (Hz)",
      y = "Power"
    ) +
    ggplot2::theme_minimal()
}
