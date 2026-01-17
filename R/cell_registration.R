#' Cell Registration Functions
#'
#' Functions for registering and tracking cells across multiple imaging sessions.
#'
#' @name cell_registration
NULL

#' Register Cells Across Sessions
#'
#' Match cells detected in different imaging sessions based on spatial location
#' and other features.
#'
#' @param sessions List of session data, each containing ROI information
#' @param method Registration method ("centroid", "overlap", "correlation", "hungarian")
#' @param max_distance Maximum allowed distance between matched cells (pixels)
#' @param min_overlap Minimum required overlap ratio (for "overlap" method)
#' @param reference_session Which session to use as reference (default: 1)
#' @param alignment Optional spatial alignment transforms between sessions
#' @param verbose Show progress messages
#' @return List with registration results
#'
#' @details
#' Each session should be a list with at least one of:
#' \describe{
#'   \item{rois}{List of ROI masks or labeled matrix}
#'   \item{centroids}{Matrix of [x, y] centroid coordinates}
#'   \item{mean_image}{Mean projection image (for correlation method)}
#' }
#'
#' Methods:
#' \describe{
#'   \item{centroid}{Match based on centroid distance}
#'   \item{overlap}{Match based on spatial overlap of ROIs}
#'   \item{correlation}{Match based on spatial correlation of ROI footprints}
#'   \item{hungarian}{Optimal assignment using Hungarian algorithm}
#' }
#'
#' @examples
#' \dontrun{
#' # Register cells from two sessions
#' session1 <- list(centroids = centroids_day1, rois = rois_day1)
#' session2 <- list(centroids = centroids_day2, rois = rois_day2)
#'
#' registration <- register_cells(
#'   sessions = list(day1 = session1, day2 = session2),
#'   method = "centroid",
#'   max_distance = 10
#' )
#'
#' # Get matched cells
#' matches <- registration$matches
#' }
#'
#' @export
register_cells <- function(sessions,
                           method = c("centroid", "overlap", "correlation", "hungarian"),
                           max_distance = 10,
                           min_overlap = 0.3,
                           reference_session = 1,
                           alignment = NULL,
                           verbose = TRUE) {

  method <- match.arg(method)
  n_sessions <- length(sessions)

  if (n_sessions < 2) {
    stop("At least 2 sessions required for registration")
  }

  if (verbose) {
    message("Registering cells across ", n_sessions, " sessions")
    message("Method: ", method)
    message("Reference session: ", reference_session)
  }

  # Extract or compute centroids for each session
  centroids_list <- lapply(seq_along(sessions), function(i) {
    session <- sessions[[i]]

    if (!is.null(session$centroids)) {
      return(session$centroids)
    }

    if (!is.null(session$rois)) {
      return(compute_roi_centroids(session$rois))
    }

    stop("Session ", i, " must have 'centroids' or 'rois'")
  })

  # Apply alignment transforms if provided
  if (!is.null(alignment)) {
    centroids_list <- apply_alignment(centroids_list, alignment, reference_session)
  }

  # Get reference centroids
  ref_centroids <- centroids_list[[reference_session]]
  n_ref <- nrow(ref_centroids)

  if (verbose) message("Reference session has ", n_ref, " cells")

  # Register each session to reference
  registrations <- list()

  for (i in seq_len(n_sessions)) {
    if (i == reference_session) {
      # Reference maps to itself
      registrations[[i]] <- data.frame(
        ref_cell = 1:n_ref,
        session_cell = 1:n_ref,
        distance = rep(0, n_ref),
        confidence = rep(1, n_ref)
      )
      next
    }

    target_centroids <- centroids_list[[i]]
    n_target <- nrow(target_centroids)

    if (verbose) message("Registering session ", i, " (", n_target, " cells) to reference")

    # Compute registration based on method
    if (method == "centroid" || method == "hungarian") {
      reg <- register_by_centroid(ref_centroids, target_centroids, max_distance,
                                  use_hungarian = (method == "hungarian"))
    } else if (method == "overlap") {
      if (is.null(sessions[[reference_session]]$rois) || is.null(sessions[[i]]$rois)) {
        stop("Overlap method requires ROI masks")
      }
      reg <- register_by_overlap(sessions[[reference_session]]$rois,
                                 sessions[[i]]$rois,
                                 min_overlap, max_distance)
    } else if (method == "correlation") {
      if (is.null(sessions[[reference_session]]$rois) || is.null(sessions[[i]]$rois)) {
        stop("Correlation method requires ROI masks")
      }
      reg <- register_by_correlation(sessions[[reference_session]]$rois,
                                     sessions[[i]]$rois,
                                     max_distance)
    }

    registrations[[i]] <- reg

    if (verbose) {
      n_matched <- sum(!is.na(reg$session_cell))
      message("  Matched ", n_matched, "/", n_ref, " reference cells")
    }
  }

  # Build unified registration table
  unified <- build_unified_registration(registrations, n_sessions)

  # Compute statistics
  stats <- compute_registration_stats(registrations, n_sessions)

  result <- list(
    registrations = registrations,
    unified = unified,
    centroids = centroids_list,
    stats = stats,
    method = method,
    parameters = list(
      max_distance = max_distance,
      min_overlap = min_overlap,
      reference_session = reference_session
    )
  )

  class(result) <- c("cell_registration", "list")
  return(result)
}

#' Compute ROI Centroids
#'
#' Compute centroid coordinates for a set of ROIs.
#'
#' @param rois List of ROI masks or labeled matrix
#' @return Matrix with columns [x, y]
#'
#' @export
compute_roi_centroids <- function(rois) {
  # Convert labeled matrix to list if needed
  if (is.matrix(rois) && !is.logical(rois)) {
    labels <- unique(as.vector(rois))
    labels <- labels[labels > 0]
    roi_list <- lapply(labels, function(l) rois == l)
  } else if (is.list(rois)) {
    roi_list <- rois
  } else {
    stop("ROIs must be a list of masks or a labeled matrix")
  }

  centroids <- matrix(0, nrow = length(roi_list), ncol = 2)
  colnames(centroids) <- c("x", "y")

  for (i in seq_along(roi_list)) {
    mask <- roi_list[[i]]
    coords <- which(mask, arr.ind = TRUE)

    if (nrow(coords) > 0) {
      centroids[i, 1] <- mean(coords[, 2])  # x = column
      centroids[i, 2] <- mean(coords[, 1])  # y = row
    }
  }

  return(centroids)
}

#' Register by Centroid Distance
#' @keywords internal
register_by_centroid <- function(ref_centroids, target_centroids, max_distance,
                                 use_hungarian = FALSE) {
  n_ref <- nrow(ref_centroids)
  n_target <- nrow(target_centroids)

  # Compute distance matrix
  dist_matrix <- matrix(Inf, nrow = n_ref, ncol = n_target)

  for (i in seq_len(n_ref)) {
    for (j in seq_len(n_target)) {
      dist_matrix[i, j] <- sqrt(sum((ref_centroids[i, ] - target_centroids[j, ])^2))
    }
  }

  # Apply distance threshold
  dist_matrix[dist_matrix > max_distance] <- Inf

  if (use_hungarian) {
    # Use Hungarian algorithm for optimal assignment
    assignment <- hungarian_assignment(dist_matrix)
    matches <- assignment$matches
    distances <- assignment$costs
  } else {
    # Greedy nearest neighbor
    matches <- rep(NA, n_ref)
    distances <- rep(NA, n_ref)

    for (i in seq_len(n_ref)) {
      min_dist <- min(dist_matrix[i, ])
      if (is.finite(min_dist)) {
        j <- which.min(dist_matrix[i, ])
        matches[i] <- j
        distances[i] <- min_dist

        # Remove matched target from consideration
        dist_matrix[, j] <- Inf
      }
    }
  }

  # Compute confidence based on distance
  confidence <- rep(NA, n_ref)
  matched <- !is.na(matches)
  if (any(matched)) {
    confidence[matched] <- 1 - distances[matched] / max_distance
  }

  data.frame(
    ref_cell = 1:n_ref,
    session_cell = matches,
    distance = distances,
    confidence = confidence
  )
}

#' Hungarian Algorithm for Optimal Assignment
#' @keywords internal
hungarian_assignment <- function(cost_matrix) {
  n_rows <- nrow(cost_matrix)
  n_cols <- ncol(cost_matrix)
  n <- max(n_rows, n_cols)

  # Pad matrix to square
  padded <- matrix(Inf, n, n)
  padded[1:n_rows, 1:n_cols] <- cost_matrix

  # Replace Inf with large value
  max_val <- max(padded[is.finite(padded)], na.rm = TRUE) * n + 1
  padded[!is.finite(padded)] <- max_val

  # Step 1: Subtract row minimum
  for (i in 1:n) {
    padded[i, ] <- padded[i, ] - min(padded[i, ])
  }

  # Step 2: Subtract column minimum
  for (j in 1:n) {
    padded[, j] <- padded[, j] - min(padded[, j])
  }

  # Simplified assignment: greedy on zero-cost entries
  matches <- rep(NA, n_rows)
  costs <- rep(NA, n_rows)
  used_cols <- logical(n_cols)

  # Sort rows by number of zeros (least first for better assignment)
  zero_counts <- apply(padded[1:n_rows, 1:n_cols, drop = FALSE], 1,
                       function(x) sum(x == 0))
  row_order <- order(zero_counts)

  for (i in row_order) {
    for (j in which(padded[i, 1:n_cols] == 0)) {
      if (!used_cols[j]) {
        matches[i] <- j
        costs[i] <- cost_matrix[i, j]
        used_cols[j] <- TRUE
        break
      }
    }
  }

  # For unmatched, find minimum
  for (i in which(is.na(matches))) {
    available <- which(!used_cols)
    if (length(available) > 0) {
      min_idx <- available[which.min(cost_matrix[i, available])]
      if (cost_matrix[i, min_idx] < max_val) {
        matches[i] <- min_idx
        costs[i] <- cost_matrix[i, min_idx]
        used_cols[min_idx] <- TRUE
      }
    }
  }

  # Set Inf costs to NA
  matches[costs >= max_val] <- NA
  costs[costs >= max_val] <- NA

  list(matches = matches, costs = costs)
}

#' Register by ROI Overlap
#' @keywords internal
register_by_overlap <- function(ref_rois, target_rois, min_overlap, max_distance) {
  # Convert to lists if needed
  if (is.matrix(ref_rois) && !is.logical(ref_rois)) {
    labels <- unique(as.vector(ref_rois))[unique(as.vector(ref_rois)) > 0]
    ref_list <- lapply(labels, function(l) ref_rois == l)
  } else {
    ref_list <- ref_rois
  }

  if (is.matrix(target_rois) && !is.logical(target_rois)) {
    labels <- unique(as.vector(target_rois))[unique(as.vector(target_rois)) > 0]
    target_list <- lapply(labels, function(l) target_rois == l)
  } else {
    target_list <- target_rois
  }

  n_ref <- length(ref_list)
  n_target <- length(target_list)

  # Compute overlap matrix
  overlap_matrix <- matrix(0, nrow = n_ref, ncol = n_target)

  for (i in seq_len(n_ref)) {
    ref_mask <- ref_list[[i]]
    ref_size <- sum(ref_mask)

    for (j in seq_len(n_target)) {
      target_mask <- target_list[[j]]
      target_size <- sum(target_mask)

      intersection <- sum(ref_mask & target_mask)

      # Jaccard index
      union_size <- ref_size + target_size - intersection
      if (union_size > 0) {
        overlap_matrix[i, j] <- intersection / union_size
      }
    }
  }

  # Find best matches above threshold
  matches <- rep(NA, n_ref)
  overlaps <- rep(NA, n_ref)

  for (i in seq_len(n_ref)) {
    max_overlap <- max(overlap_matrix[i, ])
    if (max_overlap >= min_overlap) {
      j <- which.max(overlap_matrix[i, ])
      matches[i] <- j
      overlaps[i] <- max_overlap

      # Remove from consideration
      overlap_matrix[, j] <- 0
    }
  }

  data.frame(
    ref_cell = 1:n_ref,
    session_cell = matches,
    distance = NA,  # Not applicable for overlap
    confidence = overlaps
  )
}

#' Register by Correlation
#' @keywords internal
register_by_correlation <- function(ref_rois, target_rois, max_distance) {
  # Get centroids first
  ref_centroids <- compute_roi_centroids(ref_rois)
  target_centroids <- compute_roi_centroids(target_rois)

  # Convert to lists
  if (is.matrix(ref_rois) && !is.logical(ref_rois)) {
    labels <- unique(as.vector(ref_rois))[unique(as.vector(ref_rois)) > 0]
    ref_list <- lapply(labels, function(l) ref_rois == l * 1)  # Convert to numeric
  } else {
    ref_list <- lapply(ref_rois, function(m) m * 1)
  }

  if (is.matrix(target_rois) && !is.logical(target_rois)) {
    labels <- unique(as.vector(target_rois))[unique(as.vector(target_rois)) > 0]
    target_list <- lapply(labels, function(l) target_rois == l * 1)
  } else {
    target_list <- lapply(target_rois, function(m) m * 1)
  }

  n_ref <- length(ref_list)
  n_target <- length(target_list)

  matches <- rep(NA, n_ref)
  correlations <- rep(NA, n_ref)
  distances <- rep(NA, n_ref)

  for (i in seq_len(n_ref)) {
    ref_mask <- ref_list[[i]]
    ref_centroid <- ref_centroids[i, ]

    best_cor <- -1
    best_j <- NA
    best_dist <- NA

    for (j in seq_len(n_target)) {
      # Check distance first
      dist <- sqrt(sum((ref_centroid - target_centroids[j, ])^2))
      if (dist > max_distance) next

      target_mask <- target_list[[j]]

      # Correlation of spatial footprints
      cor_val <- cor(as.vector(ref_mask), as.vector(target_mask))

      if (!is.na(cor_val) && cor_val > best_cor) {
        best_cor <- cor_val
        best_j <- j
        best_dist <- dist
      }
    }

    if (!is.na(best_j) && best_cor > 0.3) {  # Minimum correlation threshold
      matches[i] <- best_j
      correlations[i] <- best_cor
      distances[i] <- best_dist
    }
  }

  data.frame(
    ref_cell = 1:n_ref,
    session_cell = matches,
    distance = distances,
    confidence = correlations
  )
}

#' Build Unified Registration Table
#' @keywords internal
build_unified_registration <- function(registrations, n_sessions) {
  n_ref <- nrow(registrations[[1]])

  # Create unified table
  unified <- data.frame(matrix(NA, nrow = n_ref, ncol = n_sessions + 1))
  colnames(unified) <- c("unified_id", paste0("session_", 1:n_sessions))
  unified$unified_id <- 1:n_ref

  for (s in seq_len(n_sessions)) {
    reg <- registrations[[s]]
    unified[[paste0("session_", s)]] <- reg$session_cell
  }

  # Count how many sessions each cell is found in
  unified$n_sessions <- rowSums(!is.na(unified[, 2:(n_sessions + 1)]))

  return(unified)
}

#' Compute Registration Statistics
#' @keywords internal
compute_registration_stats <- function(registrations, n_sessions) {
  n_ref <- nrow(registrations[[1]])

  # Number of matches per session
  matches_per_session <- sapply(registrations, function(r) sum(!is.na(r$session_cell)))

  # Cells found in all sessions
  unified <- build_unified_registration(registrations, n_sessions)
  cells_all_sessions <- sum(unified$n_sessions == n_sessions)

  # Mean distance for matched cells
  mean_distances <- sapply(registrations, function(r) {
    mean(r$distance, na.rm = TRUE)
  })

  list(
    n_reference_cells = n_ref,
    matches_per_session = matches_per_session,
    cells_in_all_sessions = cells_all_sessions,
    fraction_in_all = cells_all_sessions / n_ref,
    mean_match_distance = mean_distances
  )
}

#' Apply Alignment Transforms
#' @keywords internal
apply_alignment <- function(centroids_list, alignment, reference_session) {
  for (i in seq_along(centroids_list)) {
    if (i == reference_session) next

    if (i <= length(alignment) && !is.null(alignment[[i]])) {
      transform <- alignment[[i]]

      # Apply affine transform
      centroids <- centroids_list[[i]]
      if (!is.null(transform$rotation)) {
        theta <- transform$rotation
        R <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
        centroids <- t(R %*% t(centroids))
      }
      if (!is.null(transform$translation)) {
        centroids[, 1] <- centroids[, 1] + transform$translation[1]
        centroids[, 2] <- centroids[, 2] + transform$translation[2]
      }
      if (!is.null(transform$scale)) {
        centroids <- centroids * transform$scale
      }

      centroids_list[[i]] <- centroids
    }
  }

  return(centroids_list)
}

#' Plot Cell Registration Results
#'
#' Visualize cell registration across sessions.
#'
#' @param registration Registration result from register_cells()
#' @param type Plot type ("scatter", "matching", "overlap")
#' @param sessions Which sessions to plot (default: all)
#' @return ggplot object
#'
#' @export
plot_cell_registration <- function(registration,
                                   type = c("scatter", "matching", "overlap"),
                                   sessions = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }

  type <- match.arg(type)
  centroids <- registration$centroids

  if (is.null(sessions)) {
    sessions <- seq_along(centroids)
  }

  if (type == "scatter") {
    # Scatter plot of cell positions across sessions
    df_list <- lapply(sessions, function(s) {
      data.frame(
        x = centroids[[s]][, 1],
        y = centroids[[s]][, 2],
        session = as.factor(s)
      )
    })
    df <- do.call(rbind, df_list)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = session)) +
      ggplot2::geom_point(alpha = 0.6, size = 2) +
      ggplot2::labs(
        title = "Cell Positions Across Sessions",
        x = "X (pixels)",
        y = "Y (pixels)",
        color = "Session"
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal()

  } else if (type == "matching") {
    # Show matched cells with connecting lines
    ref_session <- registration$parameters$reference_session
    ref_centroids <- centroids[[ref_session]]

    # Build data for matching lines
    df_lines <- data.frame()

    for (s in sessions[sessions != ref_session]) {
      reg <- registration$registrations[[s]]
      matched <- !is.na(reg$session_cell)

      if (any(matched)) {
        target_centroids <- centroids[[s]]

        lines_df <- data.frame(
          x_start = ref_centroids[matched, 1],
          y_start = ref_centroids[matched, 2],
          x_end = target_centroids[reg$session_cell[matched], 1],
          y_end = target_centroids[reg$session_cell[matched], 2],
          session = as.factor(s)
        )
        df_lines <- rbind(df_lines, lines_df)
      }
    }

    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = df_lines,
        ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = session),
        alpha = 0.5
      ) +
      ggplot2::geom_point(
        data = data.frame(x = ref_centroids[, 1], y = ref_centroids[, 2]),
        ggplot2::aes(x = x, y = y),
        color = "black", size = 2
      ) +
      ggplot2::labs(
        title = "Cell Matching Across Sessions",
        x = "X (pixels)",
        y = "Y (pixels)",
        color = "Target Session"
      ) +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal()
  }

  return(p)
}

#' Extract Registered Traces
#'
#' Extract traces for cells that are matched across sessions.
#'
#' @param registration Registration result from register_cells()
#' @param traces_list List of trace data frames (one per session)
#' @param sessions Which sessions to include (default: all)
#' @param require_all_sessions Only include cells found in all requested sessions
#' @return List with registered traces
#'
#' @export
extract_registered_traces <- function(registration,
                                      traces_list,
                                      sessions = NULL,
                                      require_all_sessions = TRUE) {
  unified <- registration$unified
  n_sessions <- length(registration$centroids)

  if (is.null(sessions)) {
    sessions <- 1:n_sessions
  }

  # Filter to cells present in requested sessions
  present_cols <- paste0("session_", sessions)

  if (require_all_sessions) {
    keep_rows <- rowSums(!is.na(unified[, present_cols, drop = FALSE])) == length(sessions)
  } else {
    keep_rows <- rowSums(!is.na(unified[, present_cols, drop = FALSE])) >= 1
  }

  unified_filtered <- unified[keep_rows, ]
  n_cells <- nrow(unified_filtered)

  if (n_cells == 0) {
    warning("No cells found in the requested sessions")
    return(NULL)
  }

  # Extract traces for each session
  result <- list()

  for (s in sessions) {
    session_cells <- unified_filtered[[paste0("session_", s)]]
    traces <- traces_list[[s]]

    if (is.data.frame(traces)) {
      # Get cell columns
      cell_cols <- grep("^Cell_", names(traces), value = TRUE)

      if (length(cell_cols) > 0) {
        # Map session cell indices to column names
        matched_cols <- cell_cols[session_cells]
        matched_cols <- matched_cols[!is.na(matched_cols)]

        if (length(matched_cols) > 0) {
          result[[paste0("session_", s)]] <- traces[, matched_cols, drop = FALSE]
        }
      }
    } else if (is.matrix(traces)) {
      valid_cells <- session_cells[!is.na(session_cells)]
      result[[paste0("session_", s)]] <- traces[, valid_cells, drop = FALSE]
    }
  }

  result$unified_ids <- unified_filtered$unified_id
  result$n_cells <- n_cells

  return(result)
}
