#' data.table Operations for Calcium Imaging
#'
#' High-performance operations using data.table for large trace matrices.
#' These functions provide significant speedups for common operations on
#' fluorescence traces.
#'
#' @name data_table_ops
#' @keywords internal
NULL

#' Convert trace matrix to data.table format
#'
#' Converts a cells x time matrix to long-format data.table for
#' efficient grouped operations.
#'
#' @param traces Matrix of traces (cells x time)
#' @param time_vector Optional numeric vector of timestamps
#' @param cell_ids Optional character vector of cell identifiers
#'
#' @return data.table in long format with columns: cell_id, time, frame, value
#' @export
#'
#' @examples
#' traces <- matrix(rnorm(10 * 100), nrow = 10)
#' dt <- traces_to_dt(traces)
traces_to_dt <- function(traces, time_vector = NULL, cell_ids = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required. Install with: install.packages('data.table')")
  }

  n_cells <- nrow(traces)
  n_frames <- ncol(traces)

  if (is.null(time_vector)) {
    time_vector <- seq_len(n_frames)
  }

  if (is.null(cell_ids)) {
    cell_ids <- paste0("cell_", seq_len(n_cells))
  }

  # Efficient conversion using data.table
  dt <- data.table::data.table(
    cell_id = rep(cell_ids, each = n_frames),
    frame = rep(seq_len(n_frames), times = n_cells),
    time = rep(time_vector, times = n_cells),
    value = as.vector(t(traces))
  )

  data.table::setkey(dt, cell_id, frame)
  dt
}

#' Convert data.table back to trace matrix
#'
#' @param dt data.table with columns cell_id, frame, value
#'
#' @return Matrix of traces (cells x time)
#' @export
dt_to_traces <- function(dt) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required")
  }

  # Use data.table's efficient dcast
  wide <- data.table::dcast(dt, cell_id ~ frame, value.var = "value")
  cell_ids <- wide$cell_id
  mat <- as.matrix(wide[, -1, with = FALSE])
  rownames(mat) <- cell_ids
  mat
}

#' Compute delta F/F using data.table
#'
#' Fast baseline-normalized fluorescence calculation.
#'
#' @param dt data.table from traces_to_dt
#' @param baseline_frames Integer vector of frames to use as baseline
#' @param method Baseline method: "mean", "median", or "rolling"
#' @param rolling_window Window size for rolling baseline
#'
#' @return data.table with additional dff column
#' @export
dt_compute_dff <- function(dt, baseline_frames = NULL, method = "mean",
                           rolling_window = 100) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required")
  }

  # Make a copy to avoid modifying input
  dt <- data.table::copy(dt)

  if (method == "rolling") {
    # Rolling baseline using data.table's frollmean
    dt[, baseline := data.table::frollmean(value, n = rolling_window, align = "right", fill = NA),
       by = cell_id]
    dt[is.na(baseline), baseline := mean(value, na.rm = TRUE), by = cell_id]
  } else if (is.null(baseline_frames)) {
    # Use first 10% as baseline
    dt[, baseline := {
      n <- .N
      bl_end <- max(1, floor(n * 0.1))
      if (method == "median") {
        median(value[1:bl_end])
      } else {
        mean(value[1:bl_end])
      }
    }, by = cell_id]
  } else {
    # Use specified frames
    dt[, baseline := {
      if (method == "median") {
        median(value[frame %in% baseline_frames])
      } else {
        mean(value[frame %in% baseline_frames])
      }
    }, by = cell_id]
  }

  dt[, dff := (value - baseline) / baseline]
  dt[, baseline := NULL]  # Clean up

  dt
}

#' Z-score normalization using data.table
#'
#' @param dt data.table from traces_to_dt
#' @param by_cell Logical, normalize each cell separately
#'
#' @return data.table with additional zscore column
#' @export
dt_zscore <- function(dt, by_cell = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required")
  }

  dt <- data.table::copy(dt)

  if (by_cell) {
    dt[, zscore := (value - mean(value)) / sd(value), by = cell_id]
  } else {
    mu <- mean(dt$value)
    sigma <- sd(dt$value)
    dt[, zscore := (value - mu) / sigma]
  }

  dt
}

#' Detect events/spikes using data.table
#'
#' Fast threshold-based event detection.
#'
#' @param dt data.table with value or zscore column
#' @param threshold Threshold for detection (in SD if use_zscore=TRUE)
#' @param use_zscore Use zscore column if available
#' @param min_gap Minimum frames between events
#'
#' @return data.table with event column (0/1)
#' @export
dt_detect_events <- function(dt, threshold = 2.5, use_zscore = TRUE,
                             min_gap = 3) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required")
  }

  dt <- data.table::copy(dt)

  # Compute zscore if needed
  if (use_zscore && !"zscore" %in% names(dt)) {
    dt <- dt_zscore(dt, by_cell = TRUE)
  }

  col <- if (use_zscore && "zscore" %in% names(dt)) "zscore" else "value"

  # Threshold detection with minimum gap
  dt[, event := {
    vals <- get(col)
    above <- vals > threshold
    events <- rep(0L, length(vals))

    # Find peaks above threshold
    if (any(above)) {
      idx <- which(above)
      events[idx[1]] <- 1L

      if (length(idx) > 1) {
        for (i in 2:length(idx)) {
          if (idx[i] - idx[i-1] >= min_gap || events[idx[i-1]] == 0L) {
            events[idx[i]] <- 1L
          }
        }
      }
    }
    events
  }, by = cell_id]

  dt
}

#' Compute summary statistics per cell
#'
#' @param dt data.table from traces_to_dt
#'
#' @return data.table with summary statistics per cell
#' @export
dt_cell_summary <- function(dt) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required")
  }

  dt[, .(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE),
    range_value = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
    skewness = {
      x <- value - mean(value)
      n <- length(x)
      m3 <- sum(x^3) / n
      s3 <- (sum(x^2) / n)^(3/2)
      m3 / s3
    },
    n_frames = .N
  ), by = cell_id]
}

#' Compute pairwise correlations efficiently
#'
#' Uses data.table for memory-efficient correlation computation.
#'
#' @param dt data.table from traces_to_dt
#' @param method Correlation method: "pearson", "spearman"
#'
#' @return Matrix of pairwise correlations
#' @export
dt_pairwise_cor <- function(dt, method = "pearson") {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required")
  }

  # Convert to wide matrix for correlation
  traces <- dt_to_traces(dt)
  cor(t(traces), method = method)
}

#' Bin traces by time windows
#'
#' Aggregate traces into time bins for reduced temporal resolution.
#'
#' @param dt data.table from traces_to_dt
#' @param bin_size Number of frames per bin
#' @param fun Aggregation function
#'
#' @return data.table with binned values
#' @export
dt_bin_traces <- function(dt, bin_size = 10, fun = mean) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required")
  }

  dt <- data.table::copy(dt)
  dt[, bin := (frame - 1) %/% bin_size + 1]

  dt[, .(
    value = fun(value, na.rm = TRUE),
    time = mean(time),
    frame = min(frame)
  ), by = .(cell_id, bin)]
}

#' Apply function to traces using data.table
#'
#' Efficiently apply a function to each cell's trace.
#'
#' @param dt data.table from traces_to_dt
#' @param fun Function to apply (receives vector, returns vector of same length)
#' @param ... Additional arguments to fun
#'
#' @return data.table with transformed values
#' @export
dt_transform <- function(dt, fun, ...) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required")
  }

  dt <- data.table::copy(dt)
  dt[, value := fun(value, ...), by = cell_id]
  dt
}

#' Smooth traces using data.table
#'
#' @param dt data.table from traces_to_dt
#' @param window_size Smoothing window size
#' @param method Smoothing method: "mean", "median", "gaussian"
#'
#' @return data.table with smoothed values
#' @export
dt_smooth <- function(dt, window_size = 5, method = "mean") {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package required")
  }

  dt <- data.table::copy(dt)

  if (method == "mean") {
    dt[, value := data.table::frollmean(value, n = window_size, align = "center", fill = NA),
       by = cell_id]
  } else if (method == "median") {
    dt[, value := data.table::frollapply(value, n = window_size, FUN = median, align = "center", fill = NA),
       by = cell_id]
  } else if (method == "gaussian") {
    # Gaussian weights
    half_w <- window_size %/% 2
    sigma <- window_size / 4
    weights <- dnorm(seq(-half_w, half_w), sd = sigma)
    weights <- weights / sum(weights)

    dt[, value := {
      x <- value
      n <- length(x)
      result <- rep(NA_real_, n)
      for (i in (half_w + 1):(n - half_w)) {
        result[i] <- sum(x[(i - half_w):(i + half_w)] * weights)
      }
      result
    }, by = cell_id]
  }

  # Fill NAs at edges with original values
  dt[is.na(value), value := dt[!is.na(value), value[1]], by = cell_id]

  dt
}
