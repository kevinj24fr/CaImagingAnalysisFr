#' Recurrence Quantification Analysis for Neural Dynamics
#'
#' Characterize dynamics complexity, determinism, and regime transitions
#' using recurrence plots and derived statistics.
#'
#' @name recurrence
NULL

# ============================================================================
# Recurrence Plot Construction
# ============================================================================

#' Compute Recurrence Plot
#'
#' Construct a recurrence plot from neural trajectories to analyze dynamics.
#'
#' @param x Time series or trajectory matrix (time x dimensions)
#' @param embed_dim Embedding dimension for time series
#' @param embed_delay Embedding delay (tau)
#' @param threshold Distance threshold for recurrence
#' @param threshold_method How to set threshold ("fixed", "fan", "quantile")
#' @param threshold_value Value for threshold method
#' @param norm Distance norm ("euclidean", "max", "manhattan")
#'
#' @return List with:
#'   - rp: Recurrence matrix (binary)
#'   - distance: Distance matrix
#'   - threshold: Used threshold value
#'   - embedding: Embedded trajectory (if applicable)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze neural trajectory dynamics
#' pca <- prcomp(t(traces))
#' rp <- recurrence_plot(pca$x[, 1:3], threshold_method = "fan")
#' image(rp$rp)
#' }
recurrence_plot <- function(x,
                             embed_dim = NULL,
                             embed_delay = 1,
                             threshold = NULL,
                             threshold_method = c("fixed", "fan", "quantile"),
                             threshold_value = 0.1,
                             norm = c("euclidean", "max", "manhattan")) {

  threshold_method <- match.arg(threshold_method)
  norm <- match.arg(norm)

  # Handle time series vs trajectory
  if (is.vector(x) || (is.matrix(x) && ncol(x) == 1)) {
    x <- as.numeric(x)

    # Time-delay embedding
    if (is.null(embed_dim)) {
      embed_dim <- estimate_embedding_dim(x)
    }

    embedding <- embed_time_series(x, embed_dim, embed_delay)
  } else {
    embedding <- x
  }

  n <- nrow(embedding)

  # Compute distance matrix
  D <- compute_distance_matrix(embedding, norm)

  # Determine threshold
  if (is.null(threshold)) {
    threshold <- switch(threshold_method,
      "fixed" = threshold_value,
      "fan" = quantile(D[upper.tri(D)], threshold_value),
      "quantile" = quantile(D[upper.tri(D)], threshold_value)
    )
  }

  # Construct recurrence matrix
  R <- D <= threshold

  structure(
    list(
      rp = R,
      distance = D,
      threshold = threshold,
      embedding = embedding,
      embed_dim = embed_dim,
      embed_delay = embed_delay,
      n = n
    ),
    class = "recurrence_plot"
  )
}

#' Time-Delay Embedding
#' @keywords internal
embed_time_series <- function(x, dim, delay) {
  n <- length(x)
  n_embedded <- n - (dim - 1) * delay

  if (n_embedded < 10) {
    stop("Time series too short for given embedding parameters")
  }

  embedding <- matrix(0, n_embedded, dim)
  for (d in 1:dim) {
    start <- (d - 1) * delay + 1
    embedding[, d] <- x[start:(start + n_embedded - 1)]
  }

  embedding
}

#' Estimate Optimal Embedding Dimension
#' @keywords internal
estimate_embedding_dim <- function(x, max_dim = 10, threshold = 0.1) {
  # False nearest neighbors method (simplified)
  n <- length(x)

  for (dim in 1:max_dim) {
    E <- embed_time_series(x, dim, delay = 1)
    E_next <- embed_time_series(x, dim + 1, delay = 1)

    n_e <- nrow(E)
    fnn_ratio <- 0

    for (i in 1:min(100, n_e - 1)) {
      # Find nearest neighbor in dim
      dists <- sqrt(rowSums((E[-i, ] - E[i, ])^2))
      nn_idx <- which.min(dists)
      nn_dist <- dists[nn_idx]

      if (nn_dist < 1e-10) next

      # Check if still neighbors in dim+1
      dist_next <- sqrt(sum((E_next[i, ] - E_next[nn_idx, ])^2))

      if ((dist_next - nn_dist) / nn_dist > threshold) {
        fnn_ratio <- fnn_ratio + 1
      }
    }

    if (fnn_ratio / min(100, n_e - 1) < 0.01) {
      return(dim)
    }
  }

  max_dim
}

#' Compute Distance Matrix
#' @keywords internal
compute_distance_matrix <- function(X, norm) {
  n <- nrow(X)

  if (norm == "euclidean") {
    as.matrix(dist(X, method = "euclidean"))
  } else if (norm == "max") {
    as.matrix(dist(X, method = "maximum"))
  } else {
    as.matrix(dist(X, method = "manhattan"))
  }
}

# ============================================================================
# Recurrence Quantification Analysis (RQA)
# ============================================================================

#' Compute RQA Metrics
#'
#' Extract quantitative measures from a recurrence plot.
#'
#' @param rp Recurrence plot object or binary matrix
#' @param min_line_length Minimum length for diagonal/vertical lines
#'
#' @return List of RQA metrics:
#'   - RR: Recurrence rate
#'   - DET: Determinism
#'   - LMAX: Maximum diagonal line length
#'   - L: Average diagonal line length
#'   - ENTR: Entropy of diagonal line distribution
#'   - LAM: Laminarity
#'   - TT: Trapping time
#'   - VMAX: Maximum vertical line length
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rp <- recurrence_plot(trajectory)
#' metrics <- rqa_metrics(rp)
#' print(metrics$DET)  # Determinism
#' }
rqa_metrics <- function(rp, min_line_length = 2) {

  if (inherits(rp, "recurrence_plot")) {
    R <- rp$rp
  } else {
    R <- as.matrix(rp)
  }

  n <- nrow(R)

  # Basic recurrence rate
  RR <- sum(R) / n^2

  # Diagonal line structures (for determinism)
  diag_lines <- extract_diagonal_lines(R, min_line_length)

  if (length(diag_lines) > 0) {
    DET <- sum(diag_lines) / sum(R[upper.tri(R) | lower.tri(R)])
    LMAX <- max(diag_lines)
    L <- mean(diag_lines)

    # Entropy of line length distribution
    line_hist <- table(diag_lines)
    p <- as.numeric(line_hist) / sum(line_hist)
    ENTR <- -sum(p * log(p + 1e-10))
  } else {
    DET <- 0
    LMAX <- 0
    L <- 0
    ENTR <- 0
  }

  # Vertical line structures (for laminarity)
  vert_lines <- extract_vertical_lines(R, min_line_length)

  if (length(vert_lines) > 0) {
    LAM <- sum(vert_lines) / sum(R)
    TT <- mean(vert_lines)
    VMAX <- max(vert_lines)
  } else {
    LAM <- 0
    TT <- 0
    VMAX <- 0
  }

  # Divergence (inverse of LMAX)
  DIV <- if (LMAX > 0) 1 / LMAX else Inf

  # Ratio (DET/RR)
  RATIO <- if (RR > 0) DET / RR else 0

  structure(
    list(
      RR = RR,
      DET = DET,
      L = L,
      LMAX = LMAX,
      ENTR = ENTR,
      LAM = LAM,
      TT = TT,
      VMAX = VMAX,
      DIV = DIV,
      RATIO = RATIO,
      n = n
    ),
    class = "rqa_metrics"
  )
}

#' Extract Diagonal Line Lengths
#' @keywords internal
extract_diagonal_lines <- function(R, min_length) {
  n <- nrow(R)
  lines <- integer(0)

  # Upper diagonals
  for (k in 1:(n-min_length)) {
    diag_k <- diag(R[1:(n-k), (k+1):n])
    lines <- c(lines, extract_runs(diag_k, min_length))
  }

  # Lower diagonals
  for (k in 1:(n-min_length)) {
    diag_k <- diag(R[(k+1):n, 1:(n-k)])
    lines <- c(lines, extract_runs(diag_k, min_length))
  }

  lines
}

#' Extract Vertical Line Lengths
#' @keywords internal
extract_vertical_lines <- function(R, min_length) {
  n <- nrow(R)
  lines <- integer(0)

  for (j in 1:n) {
    lines <- c(lines, extract_runs(R[, j], min_length))
  }

  lines
}

#' Extract Run Lengths from Binary Vector
#' @keywords internal
extract_runs <- function(x, min_length) {
  if (length(x) == 0) return(integer(0))

  rle_result <- rle(x)
  lengths <- rle_result$lengths[rle_result$values == TRUE]
  lengths[lengths >= min_length]
}

#' @export
print.rqa_metrics <- function(x, ...) {
  cat("Recurrence Quantification Analysis\n")
  cat("==================================\n")
  cat(sprintf("Recurrence Rate (RR): %.4f\n", x$RR))
  cat(sprintf("Determinism (DET): %.4f\n", x$DET))
  cat(sprintf("Avg Diagonal Length (L): %.2f\n", x$L))
  cat(sprintf("Max Diagonal Length (LMAX): %d\n", x$LMAX))
  cat(sprintf("Entropy (ENTR): %.4f\n", x$ENTR))
  cat(sprintf("Laminarity (LAM): %.4f\n", x$LAM))
  cat(sprintf("Trapping Time (TT): %.2f\n", x$TT))
  cat(sprintf("Max Vertical Length (VMAX): %d\n", x$VMAX))
  invisible(x)
}

# ============================================================================
# Cross-Recurrence Analysis
# ============================================================================

#' Cross-Recurrence Plot
#'
#' Analyze coupling between two neural signals.
#'
#' @param x First time series or trajectory
#' @param y Second time series or trajectory
#' @param embed_dim Embedding dimension
#' @param embed_delay Embedding delay
#' @param threshold Distance threshold
#' @param norm Distance norm
#'
#' @return Cross-recurrence plot object
#' @export
cross_recurrence_plot <- function(x, y,
                                   embed_dim = 3,
                                   embed_delay = 1,
                                   threshold = NULL,
                                   norm = "euclidean") {

  # Embed both series
  if (is.vector(x)) {
    Ex <- embed_time_series(x, embed_dim, embed_delay)
  } else {
    Ex <- x
  }

  if (is.vector(y)) {
    Ey <- embed_time_series(y, embed_dim, embed_delay)
  } else {
    Ey <- y
  }

  # Ensure same length
  n <- min(nrow(Ex), nrow(Ey))
  Ex <- Ex[1:n, ]
  Ey <- Ey[1:n, ]

  # Cross-distance matrix
  D <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      D[i, j] <- switch(norm,
        "euclidean" = sqrt(sum((Ex[i, ] - Ey[j, ])^2)),
        "max" = max(abs(Ex[i, ] - Ey[j, ])),
        "manhattan" = sum(abs(Ex[i, ] - Ey[j, ]))
      )
    }
  }

  if (is.null(threshold)) {
    threshold <- quantile(D, 0.1)
  }

  R <- D <= threshold

  structure(
    list(
      crp = R,
      distance = D,
      threshold = threshold,
      n = n
    ),
    class = "cross_recurrence_plot"
  )
}

#' Cross-RQA Metrics
#'
#' @param crp Cross-recurrence plot object
#' @param min_line_length Minimum diagonal line length
#'
#' @return Cross-RQA metrics
#' @export
cross_rqa_metrics <- function(crp, min_line_length = 2) {
  R <- crp$crp
  n <- nrow(R)

  # Cross-recurrence rate
  RR <- sum(R) / n^2

  # Diagonal lines (measure of coupling)
  diag_lines <- integer(0)

  # Main diagonal and off-diagonals
  for (k in (-(n-min_line_length)):(n-min_line_length)) {
    if (k >= 0) {
      diag_k <- diag(R[1:(n-k), (k+1):n, drop = FALSE])
    } else {
      diag_k <- diag(R[(-k+1):n, 1:(n+k), drop = FALSE])
    }
    diag_lines <- c(diag_lines, extract_runs(diag_k, min_line_length))
  }

  if (length(diag_lines) > 0) {
    DET <- sum(diag_lines) / sum(R)
    LMAX <- max(diag_lines)
    L <- mean(diag_lines)
  } else {
    DET <- 0
    LMAX <- 0
    L <- 0
  }

  list(
    RR = RR,
    DET = DET,
    L = L,
    LMAX = LMAX,
    n = n
  )
}

# ============================================================================
# Windowed RQA
# ============================================================================

#' Windowed RQA for Time-Varying Dynamics
#'
#' Compute RQA metrics in sliding windows to detect regime changes.
#'
#' @param x Time series
#' @param window_size Window size
#' @param step Step size between windows
#' @param embed_dim Embedding dimension
#' @param embed_delay Embedding delay
#' @param threshold Recurrence threshold
#'
#' @return Data frame of time-varying RQA metrics
#' @export
windowed_rqa <- function(x, window_size = 100, step = 10,
                          embed_dim = 3, embed_delay = 1,
                          threshold = NULL) {

  x <- as.numeric(x)
  n <- length(x)

  if (is.null(threshold)) {
    # Global threshold
    rp_full <- recurrence_plot(x, embed_dim, embed_delay,
                                threshold_method = "quantile",
                                threshold_value = 0.1)
    threshold <- rp_full$threshold
  }

  # Window positions
  starts <- seq(1, n - window_size + 1, by = step)
  n_windows <- length(starts)

  results <- data.frame(
    window = seq_len(n_windows),
    center = starts + window_size / 2,
    RR = numeric(n_windows),
    DET = numeric(n_windows),
    L = numeric(n_windows),
    LMAX = numeric(n_windows),
    ENTR = numeric(n_windows),
    LAM = numeric(n_windows)
  )

  for (i in seq_len(n_windows)) {
    x_win <- x[starts[i]:(starts[i] + window_size - 1)]

    rp <- recurrence_plot(x_win, embed_dim, embed_delay,
                          threshold = threshold,
                          threshold_method = "fixed")
    metrics <- rqa_metrics(rp)

    results$RR[i] <- metrics$RR
    results$DET[i] <- metrics$DET
    results$L[i] <- metrics$L
    results$LMAX[i] <- metrics$LMAX
    results$ENTR[i] <- metrics$ENTR
    results$LAM[i] <- metrics$LAM
  }

  results
}

#' Plot Recurrence Plot
#'
#' @param x Recurrence plot object
#' @param ... Additional arguments
#'
#' @export
plot.recurrence_plot <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  n <- nrow(x$rp)
  df <- expand.grid(i = 1:n, j = 1:n)
  df$recur <- as.vector(x$rp)

  ggplot2::ggplot(df, ggplot2::aes(x = i, y = j, fill = recur)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_manual(values = c("white", "black"), guide = "none") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "Recurrence Plot", x = "Time", y = "Time") +
    ggplot2::theme_minimal()
}

#' Detect Regime Transitions via RQA
#'
#' @param x Time series
#' @param window_size Window size
#' @param threshold_quantile Quantile for detecting transitions
#'
#' @return List with transition times and windowed metrics
#' @export
detect_regime_transitions <- function(x, window_size = 100,
                                       threshold_quantile = 0.9) {

  wrqa <- windowed_rqa(x, window_size = window_size)

  # Detect transitions as drops in determinism
  det_diff <- abs(diff(wrqa$DET))
  threshold <- quantile(det_diff, threshold_quantile)

  transitions <- wrqa$center[which(det_diff > threshold) + 1]

  list(
    transitions = transitions,
    windowed_metrics = wrqa,
    det_changes = det_diff
  )
}
