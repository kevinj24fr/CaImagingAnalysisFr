#' Publication-Ready Calcium Imaging Visualizations
#'
#' Cutting-edge visualization functions for calcium imaging data,
#' designed for Cell, Nature, and similar high-impact journals.
#'
#' @name publication_plots
NULL

# =============================================================================
# NEURAL ACTIVITY RASTERS
# =============================================================================

#' Neural Population Raster with State Annotations
#'
#' Creates a publication-quality raster plot showing population activity
#' with optional behavioral state annotations and event markers.
#'
#' @param traces Trace matrix (cells x time), CalciumTraces object, or CaExperiment object
#' @param frame_rate Frame rate in Hz (auto-extracted from CaExperiment)
#' @param states Optional data frame with columns: start, end, state
#' @param events Optional vector of event times (in seconds)
#' @param sort_by How to sort cells: "none", "activity", "peak_time", "correlation"
#' @param normalize Normalize each cell: "zscore", "minmax", "none"
#' @param show_colorbar Show color scale bar
#' @param time_range Time range to plot (c(start, end) in seconds)
#' @param assay Assay to use when traces is a CaExperiment (default: "dff")
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic raster from CaExperiment
#' plot_neural_raster(ca)
#'
#' # Basic raster from matrix
#' plot_neural_raster(traces, frame_rate = 10)
#'
#' # With behavioral states
#' states <- data.frame(start = c(0, 30, 60), end = c(30, 60, 90),
#'                      state = c("baseline", "stimulus", "recovery"))
#' plot_neural_raster(traces, frame_rate = 10, states = states)
#' }
plot_neural_raster <- function(traces, frame_rate = NULL, states = NULL,
                               events = NULL, sort_by = "activity",
                               normalize = "zscore", show_colorbar = TRUE,
                               time_range = NULL, assay = "dff") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Handle CaExperiment objects
  if (inherits(traces, "CaExperiment")) {
    ca <- traces
    frame_rate <- ca$frame_rate %||% 10
    traces <- GetTraces(ca, assay = assay)
  }

  # Handle S7 objects
  if (inherits(traces, "CalciumTraces")) {
    frame_rate <- traces@frame_rate
    traces <- traces@data
  }

  # Default frame_rate if still NULL
  if (is.null(frame_rate)) frame_rate <- 10

  traces <- as.matrix(traces)
  n_cells <- nrow(traces)
  n_time <- ncol(traces)
  time_vec <- seq(0, (n_time - 1) / frame_rate, length.out = n_time)

  # Apply time range
  if (!is.null(time_range)) {
    idx <- time_vec >= time_range[1] & time_vec <= time_range[2]
    traces <- traces[, idx, drop = FALSE]
    time_vec <- time_vec[idx]
    n_time <- ncol(traces)
  }

  # Normalize
  if (normalize == "zscore") {
    traces <- t(scale(t(traces)))
  } else if (normalize == "minmax") {
    traces <- t(apply(traces, 1, function(x) (x - min(x)) / (max(x) - min(x) + 1e-10)))
  }

  # Sort cells
  if (sort_by == "activity") {
    cell_order <- order(rowMeans(traces), decreasing = TRUE)
  } else if (sort_by == "peak_time") {
    cell_order <- order(apply(traces, 1, which.max))
  } else if (sort_by == "correlation") {
    mean_trace <- colMeans(traces)
    cors <- apply(traces, 1, function(x) cor(x, mean_trace))
    cell_order <- order(cors, decreasing = TRUE)
  } else {
    cell_order <- seq_len(n_cells)
  }
  traces <- traces[cell_order, , drop = FALSE]

  # Build data frame
  plot_data <- expand.grid(cell = seq_len(n_cells), time_idx = seq_len(n_time))
  plot_data$time <- time_vec[plot_data$time_idx]
  plot_data$activity <- as.vector(t(traces))

  # Clamp extreme values for visualization
  zlim <- quantile(plot_data$activity, c(0.01, 0.99), na.rm = TRUE)
  plot_data$activity <- pmax(pmin(plot_data$activity, zlim[2]), zlim[1])

  # Base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = cell, fill = activity)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    scale_fill_publication("diverging", midpoint = 0) +
    ggplot2::labs(
      x = "Time (s)",
      y = "Cell",
      fill = if (normalize == "zscore") "Z-score" else "Activity"
    ) +
    theme_publication(grid = FALSE) +
    ggplot2::theme(
      legend.position = if (show_colorbar) "right" else "none"
    )

  # Add behavioral states
  if (!is.null(states)) {
    state_colors <- pub_colors(length(unique(states$state)), "categorical")
    names(state_colors) <- unique(states$state)

    for (i in seq_len(nrow(states))) {
      p <- p + ggplot2::annotate(
        "rect",
        xmin = states$start[i], xmax = states$end[i],
        ymin = n_cells + 0.5, ymax = n_cells + 2,
        fill = state_colors[states$state[i]],
        alpha = 0.8
      )
    }
    p <- p + ggplot2::coord_cartesian(ylim = c(0.5, n_cells + 2.5), expand = FALSE)
  }

  # Add event markers
  if (!is.null(events)) {
    events <- events[events >= min(time_vec) & events <= max(time_vec)]
    for (ev in events) {
      p <- p + ggplot2::geom_vline(xintercept = ev, linetype = "dashed",
                                    color = "#FF7F00", linewidth = 0.5, alpha = 0.7)
    }
  }

  p
}

#' Spike Raster Plot
#'
#' Publication-quality spike raster with trial structure.
#'
#' @param spikes Spike matrix (cells x time), CaExperiment object, or list of spike times
#' @param frame_rate Frame rate (auto-extracted from CaExperiment)
#' @param trial_starts Vector of trial start times (seconds)
#' @param trial_duration Duration of each trial (seconds)
#' @param color_by Color points by: "cell", "trial", "uniform"
#' @param point_size Size of spike markers
#' @param method Spike method to use when spikes is a CaExperiment (default: first available)
#'
#' @return A ggplot object
#' @export
plot_spike_raster_pub <- function(spikes, frame_rate = NULL, trial_starts = NULL,
                                  trial_duration = NULL, color_by = "cell",
                                  point_size = 0.5, method = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Handle CaExperiment objects
  if (inherits(spikes, "CaExperiment")) {
    ca <- spikes
    frame_rate <- ca$frame_rate %||% 10
    spikes <- GetSpikes(ca, method = method)
  }

  # Default frame_rate if still NULL
  if (is.null(frame_rate)) frame_rate <- 10

  # Convert matrix to spike times
  if (is.matrix(spikes)) {
    n_cells <- nrow(spikes)
    n_time <- ncol(spikes)
    time_vec <- seq(0, (n_time - 1) / frame_rate, length.out = n_time)

    spike_data <- do.call(rbind, lapply(seq_len(n_cells), function(i) {
      spike_idx <- which(spikes[i, ] > 0)
      if (length(spike_idx) > 0) {
        data.frame(cell = i, time = time_vec[spike_idx])
      } else {
        NULL
      }
    }))
  } else {
    stop("spikes must be a matrix or CaExperiment object")
  }

  if (is.null(spike_data) || nrow(spike_data) == 0) {
    warning("No spikes detected")
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  # Assign to trials if provided
  if (!is.null(trial_starts) && !is.null(trial_duration)) {
    spike_data$trial <- NA
    for (i in seq_along(trial_starts)) {
      in_trial <- spike_data$time >= trial_starts[i] &
                  spike_data$time < trial_starts[i] + trial_duration
      spike_data$trial[in_trial] <- i
      spike_data$trial_time[in_trial] <- spike_data$time[in_trial] - trial_starts[i]
    }
    spike_data <- spike_data[!is.na(spike_data$trial), ]
    use_trial_time <- TRUE
    x_lab <- "Time from trial onset (s)"
  } else {
    use_trial_time <- FALSE
    x_lab <- "Time (s)"
  }

  # Build plot using modern aes() with .data pronoun
  if (use_trial_time) {
    p <- ggplot2::ggplot(spike_data, ggplot2::aes(x = .data$trial_time, y = .data$cell))
  } else {
    p <- ggplot2::ggplot(spike_data, ggplot2::aes(x = .data$time, y = .data$cell))
  }

  if (color_by == "cell") {
    n_unique_cells <- length(unique(spike_data$cell))
    p <- p + ggplot2::geom_point(ggplot2::aes(color = factor(.data$cell)),
                                  size = point_size, shape = "|") +
      scale_color_publication(n = n_unique_cells) +
      ggplot2::theme(legend.position = "none")
  } else if (color_by == "trial" && !is.null(trial_starts)) {
    n_unique_trials <- length(unique(spike_data$trial))
    p <- p + ggplot2::geom_point(ggplot2::aes(color = factor(.data$trial)),
                                  size = point_size, shape = "|") +
      scale_color_publication(n = n_unique_trials)
  } else {
    p <- p + ggplot2::geom_point(size = point_size, shape = "|", color = "#377EB8")
  }

  p +
    ggplot2::labs(x = x_lab, y = "Cell") +
    theme_publication(grid = FALSE)
}

# =============================================================================
# TRAJECTORY AND STATE SPACE VISUALIZATIONS
# =============================================================================

#' GPFA Neural Trajectory Plot
#'
#' Publication-quality visualization of GPFA latent trajectories.
#'
#' @param gpfa_result Result from fit_gpfa() or CaExperiment object
#' @param dims Which dimensions to plot (2 or 3)
#' @param trials Which trials to include (NULL for all)
#' @param conditions Condition labels for coloring
#' @param condition_colors Named vector of colors
#' @param show_variance Show explained variance in axis labels
#' @param arrow_interval Add directional arrows every N points (0 = none)
#' @param start_marker Shape for trajectory start
#' @param end_marker Shape for trajectory end
#' @param reduction Name of reduction to use when gpfa_result is a CaExperiment (default: "gpfa")
#'
#' @return A ggplot object
#' @export
plot_gpfa_trajectory <- function(gpfa_result, dims = c(1, 2), trials = NULL,
                                  conditions = NULL, condition_colors = NULL,
                                  show_variance = TRUE, arrow_interval = 0,
                                  start_marker = 16, end_marker = 17,
                                  reduction = "gpfa") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Handle CaExperiment objects
  if (inherits(gpfa_result, "CaExperiment")) {
    ca <- gpfa_result
    gpfa_result <- GetReduction(ca, reduction)
    if (is.null(gpfa_result)) {
      stop("No GPFA reduction found. Run RunGPFA() first.")
    }
  }

  trajectories <- gpfa_result$trajectories
  if (is.null(trials)) trials <- seq_along(trajectories)

  # Build plot data
  plot_data <- do.call(rbind, lapply(trials, function(i) {
    traj <- trajectories[[i]]
    n_time <- ncol(traj)
    data.frame(
      trial = i,
      time = seq_len(n_time),
      dim1 = traj[dims[1], ],
      dim2 = traj[dims[2], ],
      dim3 = if (length(dims) == 3) traj[dims[3], ] else NA
    )
  }))

  # Add conditions
  if (!is.null(conditions)) {
    if (length(conditions) == length(trials)) {
      plot_data$condition <- factor(conditions[match(plot_data$trial, trials)])
    } else {
      plot_data$condition <- factor(conditions[plot_data$trial])
    }
  } else {
    plot_data$condition <- factor(plot_data$trial)
  }

  # Colors
  if (is.null(condition_colors)) {
    n_cond <- length(unique(plot_data$condition))
    condition_colors <- pub_colors(n_cond, "categorical")
    names(condition_colors) <- levels(plot_data$condition)
  }

  # Labels
  dim_labels <- paste0("Latent ", dims)
  if (show_variance && !is.null(gpfa_result$explained_var)) {
    # Approximate variance explained
    dim_labels <- paste0(dim_labels, " (", round(gpfa_result$explained_var * 100), "%)")
  }

  # 2D Plot
  if (length(dims) == 2) {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = dim1, y = dim2,
                                                   color = condition, group = trial)) +
      ggplot2::geom_path(linewidth = 0.8, alpha = 0.7) +
      ggplot2::scale_color_manual(values = condition_colors) +
      theme_publication() +
      ggplot2::labs(x = dim_labels[1], y = dim_labels[2], color = "Condition")

    # Add start/end markers
    starts <- plot_data[plot_data$time == 1, ]
    ends <- do.call(rbind, lapply(unique(plot_data$trial), function(tr) {
      subset_trial <- plot_data[plot_data$trial == tr, ]
      subset_trial[nrow(subset_trial), ]
    }))

    p <- p +
      ggplot2::geom_point(data = starts, shape = start_marker, size = 3) +
      ggplot2::geom_point(data = ends, shape = end_marker, size = 3)

    # Add arrows
    if (arrow_interval > 0) {
      arrow_data <- plot_data[plot_data$time %% arrow_interval == 0, ]
      p <- p + ggplot2::geom_segment(
        data = arrow_data,
        ggplot2::aes(xend = dim1 + 0.1, yend = dim2),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")),
        linewidth = 0.3
      )
    }

    return(p)
  }

  # 3D would need plotly - return 2D projection instead
  warning("3D plots require plotly; returning 2D projection")
  plot_gpfa_trajectory(gpfa_result, dims = dims[1:2], trials = trials,
                       conditions = conditions, condition_colors = condition_colors)
}

#' dPCA Component Trajectories
#'
#' Visualize demixed PCA components across conditions.
#'
#' @param dpca_result Result from fit_dpca()
#' @param marginalization Which marginalization to plot
#' @param components Which components (1:3 by default)
#' @param condition_labels Labels for conditions
#' @param show_variance Show variance explained
#'
#' @return A ggplot object
#' @export
plot_dpca_components <- function(dpca_result, marginalization = NULL,
                                  components = 1:3, condition_labels = NULL,
                                  show_variance = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  if (is.null(marginalization)) {
    marginalization <- dpca_result$marg_names[1]
  }

  proj <- dpca_result$projections[[marginalization]]
  var_exp <- dpca_result$explained_var[[marginalization]]

  n_comp <- min(length(components), dim(proj)[1])
  n_time <- dim(proj)[2]
  n_cond <- dim(proj)[3]

  # Build data
  plot_data <- do.call(rbind, lapply(1:n_comp, function(c) {
    do.call(rbind, lapply(1:n_cond, function(cond) {
      data.frame(
        component = c,
        time = seq_len(n_time),
        condition = cond,
        value = proj[components[c], , cond]
      )
    }))
  }))

  if (!is.null(condition_labels) && length(condition_labels) == n_cond) {
    plot_data$condition <- factor(plot_data$condition, labels = condition_labels)
  } else {
    plot_data$condition <- factor(plot_data$condition)
  }

  # Component labels with variance
  comp_labels <- paste0("dPC", components)
  if (show_variance && length(var_exp) >= max(components)) {
    comp_labels <- paste0(comp_labels, " (", round(var_exp[components] * 100, 1), "%)")
  }
  plot_data$component_label <- factor(plot_data$component,
                                       labels = comp_labels[1:n_comp])

  # Plot
  ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = value,
                                           color = condition, group = condition)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_wrap(~component_label, scales = "free_y", ncol = 1) +
    scale_color_publication() +
    theme_publication() +
    ggplot2::labs(
      x = "Time bin",
      y = "Projection",
      color = "Condition",
      title = paste(marginalization, "components")
    )
}

#' State Space Flow Field
#'
#' Visualize population dynamics as a flow field in 2D state space.
#'
#' @param trajectories Result from compute_trajectories() or GPFA
#' @param dims Which dimensions to use
#' @param n_grid Number of grid points per dimension
#' @param arrow_scale Scale factor for arrows
#' @param show_trajectories Overlay actual trajectories
#'
#' @return A ggplot object
#' @export
plot_flow_field <- function(trajectories, dims = c(1, 2), n_grid = 15,
                             arrow_scale = 1, show_trajectories = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Extract trajectory data
  if (is.list(trajectories) && !is.null(trajectories$trajectories)) {
    # GPFA result
    traj_list <- trajectories$trajectories
    all_data <- do.call(rbind, lapply(seq_along(traj_list), function(i) {
      traj <- traj_list[[i]]
      data.frame(
        trial = i,
        time = seq_len(ncol(traj)),
        x = traj[dims[1], ],
        y = traj[dims[2], ]
      )
    }))
  } else if (is.list(trajectories) && !is.null(trajectories$trajectories_2d)) {
    # compute_trajectories result
    traj <- trajectories$trajectories_2d
    all_data <- data.frame(
      trial = 1,
      time = seq_len(ncol(traj)),
      x = traj[1, ],
      y = traj[2, ]
    )
  } else {
    stop("Unrecognized trajectory format")
  }

  # Compute velocity field on grid
  x_range <- range(all_data$x) + c(-0.1, 0.1) * diff(range(all_data$x))
  y_range <- range(all_data$y) + c(-0.1, 0.1) * diff(range(all_data$y))

  grid_x <- seq(x_range[1], x_range[2], length.out = n_grid)
  grid_y <- seq(y_range[1], y_range[2], length.out = n_grid)
  grid_df <- expand.grid(x = grid_x, y = grid_y)

  # For each grid point, find nearby trajectory points and average velocity
  grid_df$dx <- 0
  grid_df$dy <- 0
  bandwidth_x <- diff(x_range) / (n_grid - 1) * 2
  bandwidth_y <- diff(y_range) / (n_grid - 1) * 2

  # Compute velocities from trajectories
  all_data$dx <- c(diff(all_data$x), NA)
  all_data$dy <- c(diff(all_data$y), NA)
  all_data <- all_data[!is.na(all_data$dx), ]

  for (i in seq_len(nrow(grid_df))) {
    weights <- exp(-0.5 * (((all_data$x - grid_df$x[i]) / bandwidth_x)^2 +
                           ((all_data$y - grid_df$y[i]) / bandwidth_y)^2))
    if (sum(weights) > 0.01) {
      grid_df$dx[i] <- sum(weights * all_data$dx) / sum(weights)
      grid_df$dy[i] <- sum(weights * all_data$dy) / sum(weights)
    }
  }

  # Normalize arrow lengths
  arrow_mag <- sqrt(grid_df$dx^2 + grid_df$dy^2)
  max_mag <- quantile(arrow_mag, 0.95)
  grid_df$dx <- grid_df$dx / max_mag * arrow_scale * diff(x_range) / n_grid
  grid_df$dy <- grid_df$dy / max_mag * arrow_scale * diff(y_range) / n_grid

  # Plot
  p <- ggplot2::ggplot()

  # Flow field
  p <- p + ggplot2::geom_segment(
    data = grid_df,
    ggplot2::aes(x = x, y = y, xend = x + dx, yend = y + dy),
    arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm")),
    color = "grey40",
    linewidth = 0.3,
    alpha = 0.6
  )

  # Overlay trajectories
  if (show_trajectories) {
    p <- p + ggplot2::geom_path(
      data = all_data,
      ggplot2::aes(x = x, y = y, group = trial),
      color = "#E41A1C",
      linewidth = 0.6,
      alpha = 0.5
    )
  }

  p +
    theme_publication() +
    ggplot2::labs(x = paste("Dim", dims[1]), y = paste("Dim", dims[2])) +
    ggplot2::coord_fixed()
}

# =============================================================================
# NETWORK AND CONNECTIVITY VISUALIZATIONS
# =============================================================================

#' Connectivity Matrix Heatmap
#'
#' Publication-quality heatmap for functional connectivity matrices.
#'
#' @param connectivity Connectivity matrix, CaExperiment object, or functional_connectivity result
#' @param order_cells How to order cells: "none", "hierarchical", "community"
#' @param communities Optional community assignments for annotation
#' @param show_diagonal Show diagonal values
#' @param symmetric For asymmetric matrices, show full or just triangle
#' @param title Plot title
#' @param graph Name of graph to use when connectivity is a CaExperiment (default: first available)
#'
#' @return A ggplot object
#' @export
plot_connectivity_matrix <- function(connectivity, order_cells = "hierarchical",
                                      communities = NULL, show_diagonal = FALSE,
                                      symmetric = TRUE, title = NULL, graph = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Handle CaExperiment objects
  if (inherits(connectivity, "CaExperiment")) {
    ca <- connectivity
    connectivity <- GetGraph(ca, graph)
    if (is.null(connectivity)) {
      stop("No connectivity graph found. Run RunConnectivity() first.")
    }
  }

  # Extract matrix if it's a list result
  if (is.list(connectivity) && !is.null(connectivity$connectivity_matrix)) {
    conn_mat <- connectivity$connectivity_matrix
  } else {
    conn_mat <- as.matrix(connectivity)
  }

  n_cells <- nrow(conn_mat)

  # Handle NAs
  conn_mat[is.na(conn_mat)] <- 0

  # Remove diagonal if requested
  if (!show_diagonal) {
    diag(conn_mat) <- NA
  }

  # Order cells
  if (order_cells == "hierarchical") {
    d <- as.dist(1 - abs(conn_mat))
    d[is.na(d)] <- 1
    hc <- hclust(d, method = "ward.D2")
    cell_order <- hc$order
  } else if (order_cells == "community" && !is.null(communities)) {
    cell_order <- order(communities)
  } else {
    cell_order <- seq_len(n_cells)
  }

  conn_mat <- conn_mat[cell_order, cell_order]
  if (!is.null(communities)) {
    communities <- communities[cell_order]
  }

  # Build plot data
  plot_data <- expand.grid(row = seq_len(n_cells), col = seq_len(n_cells))
  plot_data$value <- as.vector(conn_mat)

  if (symmetric) {
    # Show only lower triangle
    plot_data$value[plot_data$row < plot_data$col] <- NA
  }

  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = col, y = row, fill = value)) +
    ggplot2::geom_tile() +
    scale_fill_publication("diverging", midpoint = 0, na.value = "white") +
    ggplot2::scale_y_reverse() +
    theme_publication(grid = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "Cell", y = "Cell", fill = "Correlation", title = title) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  # Add community boundaries
  if (!is.null(communities)) {
    boundaries <- which(diff(communities) != 0) + 0.5
    for (b in boundaries) {
      p <- p +
        ggplot2::geom_hline(yintercept = b, color = "black", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = b, color = "black", linewidth = 0.5)
    }
  }

  p
}

#' Circular Network Graph (Chord Diagram Style)
#'
#' Visualize network connectivity in a circular layout.
#'
#' @param connectivity Connectivity matrix, CaExperiment object, or functional_connectivity result
#' @param threshold Only show connections above this threshold
#' @param communities Community assignments for grouping
#' @param show_labels Show cell labels
#' @param edge_alpha Transparency of edges
#' @param graph Name of graph to use when connectivity is a CaExperiment (default: first available)
#'
#' @return A ggplot object
#' @export
plot_circular_network <- function(connectivity, threshold = 0.3,
                                   communities = NULL, show_labels = FALSE,
                                   edge_alpha = 0.5, graph = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Handle CaExperiment objects
  if (inherits(connectivity, "CaExperiment")) {
    ca <- connectivity
    connectivity <- GetGraph(ca, graph)
    if (is.null(connectivity)) {
      stop("No connectivity graph found. Run RunConnectivity() first.")
    }
  }

  # Extract matrix
  if (is.list(connectivity)) {
    conn_mat <- connectivity$connectivity_matrix
  } else {
    conn_mat <- as.matrix(connectivity)
  }

  n_cells <- nrow(conn_mat)
  conn_mat[is.na(conn_mat)] <- 0
  diag(conn_mat) <- 0

  # Order by community
  if (!is.null(communities)) {
    cell_order <- order(communities)
    conn_mat <- conn_mat[cell_order, cell_order]
    communities <- communities[cell_order]
  } else {
    cell_order <- seq_len(n_cells)
    communities <- rep(1, n_cells)
  }

  # Node positions (circular layout)
  angles <- seq(0, 2 * pi * (n_cells - 1) / n_cells, length.out = n_cells)
  node_data <- data.frame(
    cell = seq_len(n_cells),
    x = cos(angles),
    y = sin(angles),
    community = factor(communities)
  )

  # Edge data
  edges <- which(abs(conn_mat) > threshold & upper.tri(conn_mat), arr.ind = TRUE)
  if (nrow(edges) == 0) {
    warning("No edges above threshold")
    edge_data <- data.frame()
  } else {
    edge_data <- data.frame(
      from = edges[, 1],
      to = edges[, 2],
      weight = conn_mat[edges]
    )
    edge_data$x <- node_data$x[edge_data$from]
    edge_data$y <- node_data$y[edge_data$from]
    edge_data$xend <- node_data$x[edge_data$to]
    edge_data$yend <- node_data$y[edge_data$to]
    edge_data$sign <- ifelse(edge_data$weight > 0, "positive", "negative")
  }

  # Plot
  p <- ggplot2::ggplot() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()

  # Add edges
  if (nrow(edge_data) > 0) {
    p <- p + ggplot2::geom_segment(
      data = edge_data,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = sign,
                   linewidth = abs(weight)),
      alpha = edge_alpha
    ) +
      ggplot2::scale_linewidth_continuous(range = c(0.2, 1.5), guide = "none") +
      ggplot2::scale_color_manual(
        values = c("positive" = "#4DAF4A", "negative" = "#E41A1C"),
        name = "Connection"
      )
  }

  # Add nodes
  p <- p + ggplot2::geom_point(
    data = node_data,
    ggplot2::aes(x = x, y = y, fill = community),
    shape = 21, size = 3, color = "black", stroke = 0.5
  ) +
    scale_fill_publication_d()

  # Add labels
  if (show_labels) {
    p <- p + ggplot2::geom_text(
      data = node_data,
      ggplot2::aes(x = x * 1.15, y = y * 1.15, label = cell),
      size = 2.5
    )
  }

  p + ggplot2::labs(fill = "Community")
}

#' Assembly Activation Timeline
#'
#' Visualize assembly activations over time.
#'
#' @param assembly_result Result from detect_assemblies() or CaExperiment object
#' @param frame_rate Frame rate (auto-extracted from CaExperiment)
#' @param time_range Optional time range to plot
#' @param show_events Mark individual activation events
#' @param threshold Activation threshold for event detection
#' @param name Name of assembly result when assembly_result is a CaExperiment (default: "default")
#'
#' @return A ggplot object
#' @export
plot_assembly_timeline <- function(assembly_result, frame_rate = NULL,
                                    time_range = NULL, show_events = TRUE,
                                    threshold = 2, name = "default") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Handle CaExperiment objects
  if (inherits(assembly_result, "CaExperiment")) {
    ca <- assembly_result
    frame_rate <- ca$frame_rate %||% 10
    assembly_result <- GetAssemblies(ca, name)
    if (is.null(assembly_result)) {
      stop("No assembly results found. Run RunAssemblies() first.")
    }
  }

  # Default frame_rate if still NULL
  if (is.null(frame_rate)) frame_rate <- 10

  # Extract activation time series
  if (!is.null(assembly_result$activations)) {
    activations <- assembly_result$activations
  } else if (!is.null(assembly_result$activation_strength)) {
    activations <- assembly_result$activation_strength
  } else {
    stop("No activation data found in assembly result")
  }

  activations <- as.matrix(activations)
  n_assemblies <- nrow(activations)
  n_time <- ncol(activations)
  time_vec <- seq(0, (n_time - 1) / frame_rate, length.out = n_time)

  # Apply time range
  if (!is.null(time_range)) {
    idx <- time_vec >= time_range[1] & time_vec <= time_range[2]
    activations <- activations[, idx, drop = FALSE]
    time_vec <- time_vec[idx]
  }

  # Build data
  plot_data <- do.call(rbind, lapply(seq_len(n_assemblies), function(a) {
    data.frame(
      assembly = factor(a),
      time = time_vec,
      activation = activations[a, ]
    )
  }))

  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = activation,
                                                 color = assembly)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(~assembly, ncol = 1, scales = "free_y",
                        labeller = ggplot2::labeller(assembly = function(x) paste("Assembly", x))) +
    scale_color_publication() +
    theme_publication() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "Time (s)", y = "Activation")

  # Add event markers
  if (show_events) {
    events <- plot_data[plot_data$activation > threshold, ]
    if (nrow(events) > 0) {
      p <- p + ggplot2::geom_point(data = events, shape = 21,
                                    fill = "#FF7F00", color = "black", size = 2)
    }
  }

  # Add threshold line
  p + ggplot2::geom_hline(yintercept = threshold, linetype = "dashed",
                           color = "grey50", linewidth = 0.5)
}

# =============================================================================
# TUNING CURVE VISUALIZATIONS
# =============================================================================

#' Polar Tuning Curve Plot
#'
#' Visualize orientation/direction tuning on a polar plot.
#'
#' @param tuning_result Result from fit_orientation_tuning()
#' @param cells Which cells to plot
#' @param normalize Normalize responses to [0,1]
#' @param show_fit Show fitted curve
#' @param ncol Number of columns for multi-cell plots
#'
#' @return A ggplot object
#' @export
plot_tuning_polar <- function(tuning_result, cells = 1:4, normalize = TRUE,
                               show_fit = TRUE, ncol = 2) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  orientations <- attr(tuning_result, "orientations")
  responses <- attr(tuning_result, "responses")

  cells <- cells[cells <= nrow(responses)]

  # Build data
  plot_data <- do.call(rbind, lapply(cells, function(cell) {
    resp <- responses[cell, ]
    if (normalize) {
      resp <- (resp - min(resp)) / (max(resp) - min(resp) + 1e-10)
    }
    data.frame(
      cell = factor(cell),
      orientation = orientations,
      response = resp
    )
  }))

  # Close the circle by adding first point at end
  plot_data_closed <- do.call(rbind, lapply(cells, function(cell) {
    sub_data <- plot_data[plot_data$cell == cell, ]
    first_row <- sub_data[1, ]
    first_row$orientation <- max(orientations) + diff(orientations)[1]
    rbind(sub_data, first_row)
  }))

  # Polar plot
  p <- ggplot2::ggplot(plot_data_closed, ggplot2::aes(x = orientation, y = response,
                                                        color = cell, fill = cell)) +
    ggplot2::geom_polygon(alpha = 0.2, linewidth = 0.8) +
    ggplot2::geom_point(data = plot_data, size = 2) +
    ggplot2::coord_polar(start = -pi/2, direction = -1) +
    scale_color_publication() +
    scale_fill_publication_d() +
    theme_publication() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey80", linewidth = 0.3)
    ) +
    ggplot2::labs(x = "Orientation (deg)", y = NULL, color = "Cell", fill = "Cell")

  if (length(cells) > 1) {
    p <- p + ggplot2::facet_wrap(~cell, ncol = ncol)
  }

  p
}

#' Population Tuning Curve Summary
#'
#' Heatmap showing tuning curves for entire population.
#'
#' @param tuning_result Result from fit_orientation_tuning()
#' @param sort_by How to sort: "none", "preferred", "selectivity"
#' @param show_pref_markers Show preferred orientation markers
#'
#' @return A ggplot object
#' @export
plot_population_tuning <- function(tuning_result, sort_by = "preferred",
                                    show_pref_markers = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  orientations <- attr(tuning_result, "orientations")
  responses <- attr(tuning_result, "responses")

  n_cells <- nrow(responses)

  # Normalize each cell
  responses_norm <- t(apply(responses, 1, function(x) {
    (x - min(x)) / (max(x) - min(x) + 1e-10)
  }))

  # Sort cells
  if (sort_by == "preferred") {
    cell_order <- order(tuning_result$preferred_orientation)
  } else if (sort_by == "selectivity") {
    cell_order <- order(tuning_result$OSI, decreasing = TRUE)
  } else {
    cell_order <- seq_len(n_cells)
  }

  responses_norm <- responses_norm[cell_order, ]
  prefs <- tuning_result$preferred_orientation[cell_order]

  # Build data
  plot_data <- expand.grid(cell = seq_len(n_cells), ori_idx = seq_along(orientations))
  plot_data$orientation <- orientations[plot_data$ori_idx]
  plot_data$response <- as.vector(responses_norm)

  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = orientation, y = cell, fill = response)) +
    ggplot2::geom_tile() +
    scale_fill_publication("sequential") +
    theme_publication(grid = FALSE) +
    ggplot2::labs(x = "Orientation (deg)", y = "Cell", fill = "Norm.\nresponse")

  # Add preferred orientation markers
  if (show_pref_markers) {
    pref_data <- data.frame(
      cell = seq_len(n_cells),
      orientation = prefs
    )
    pref_data <- pref_data[!is.na(pref_data$orientation), ]

    if (nrow(pref_data) > 0) {
      p <- p + ggplot2::geom_point(
        data = pref_data,
        ggplot2::aes(x = orientation, y = cell),
        shape = 4, color = "#E41A1C", size = 1, stroke = 0.8,
        inherit.aes = FALSE
      )
    }
  }

  p
}

#' Selectivity Distribution Plot
#'
#' Histogram of selectivity indices across the population.
#'
#' @param tuning_result Result from fit_orientation_tuning()
#' @param index Which index: "OSI", "DSI", "gOSI"
#' @param bins Number of histogram bins
#' @param show_median Show median line
#'
#' @return A ggplot object
#' @export
plot_selectivity_distribution <- function(tuning_result, index = "OSI",
                                           bins = 20, show_median = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  values <- tuning_result[[index]]
  values <- values[!is.na(values)]

  plot_data <- data.frame(value = values)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = bins, fill = "#377EB8", color = "black",
                             alpha = 0.7, linewidth = 0.4) +
    theme_publication() +
    ggplot2::labs(x = index, y = "Count",
                  subtitle = sprintf("n = %d cells", length(values)))

  if (show_median) {
    med <- median(values)
    p <- p +
      ggplot2::geom_vline(xintercept = med, linetype = "dashed",
                           color = "#E41A1C", linewidth = 0.8) +
      ggplot2::annotate("text", x = med, y = Inf, vjust = 2,
                         label = sprintf("Median = %.2f", med),
                         size = 3.5, fontface = "italic")
  }

  p
}

# =============================================================================
# PHARMACOLOGY AND RESPONSE VISUALIZATIONS
# =============================================================================

#' Drug Response Waterfall Plot
#'
#' Waterfall plot showing individual cell responses sorted by magnitude.
#'
#' @param baseline_activity Baseline activity values (one per cell)
#' @param treatment_activity Treatment activity values (one per cell)
#' @param responder_class Optional classification ("activated", "inhibited", "non-responder")
#' @param threshold Threshold for significance (fold-change or z-score)
#' @param log_scale Use log2 fold-change
#' @param title Plot title
#'
#' @return A ggplot object
#' @export
plot_response_waterfall <- function(baseline_activity, treatment_activity,
                                     responder_class = NULL, threshold = 1.5,
                                     log_scale = TRUE, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Calculate response
  if (log_scale) {
    response <- log2((treatment_activity + 0.01) / (baseline_activity + 0.01))
    y_lab <- expression(log[2]~"(Treatment / Baseline)")
    thresh_vals <- c(-log2(threshold), log2(threshold))
  } else {
    response <- treatment_activity - baseline_activity
    y_lab <- "Response (Treatment - Baseline)"
    thresh_vals <- c(-threshold, threshold)
  }

  # Sort by response
  cell_order <- order(response, decreasing = TRUE)
  response <- response[cell_order]

  # Build data
  plot_data <- data.frame(
    rank = seq_along(response),
    response = response
  )

  # Add responder class
  if (!is.null(responder_class)) {
    plot_data$class <- factor(responder_class[cell_order])
  } else {
    plot_data$class <- factor(ifelse(response > thresh_vals[2], "activated",
                                      ifelse(response < thresh_vals[1], "inhibited",
                                             "non-responder")))
  }

  # Colors
  class_colors <- c("activated" = "#E41A1C", "inhibited" = "#377EB8",
                    "non-responder" = "#BDBDBD")

  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = rank, y = response, fill = class)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = NA) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = thresh_vals, linetype = "dashed",
                         color = "#757575", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = class_colors, name = "Response") +
    theme_publication(grid = FALSE) +
    ggplot2::labs(x = "Cell rank", y = y_lab, title = title)

  # Add count annotations
  counts <- table(plot_data$class)
  subtitle <- sprintf("Activated: %d | Inhibited: %d | Non-responder: %d",
                      counts["activated"], counts["inhibited"], counts["non-responder"])
  p + ggplot2::labs(subtitle = subtitle)
}

#' Response Magnitude Comparison
#'
#' Violin plot comparing response magnitudes across conditions or groups.
#'
#' @param response_data Data frame with columns: cell, condition, response
#' @param x Column name for x-axis grouping
#' @param y Column name for response values
#' @param fill Column for fill color (optional)
#' @param add_points Show individual data points
#' @param add_boxplot Overlay boxplot
#' @param comparisons List of pairs to compare (for significance)
#'
#' @return A ggplot object
#' @export
plot_response_comparison <- function(response_data, x = "condition", y = "response",
                                      fill = NULL, add_points = TRUE,
                                      add_boxplot = TRUE, comparisons = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  if (is.null(fill)) fill <- x

  p <- ggplot2::ggplot(response_data, ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    ggplot2::geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
    scale_fill_publication_d() +
    theme_publication() +
    ggplot2::labs(x = NULL, y = "Response")

  if (add_boxplot) {
    p <- p + ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5)
  }

  if (add_points) {
    p <- p + ggplot2::geom_jitter(width = 0.1, size = 0.8, alpha = 0.5, color = "black")
  }

  # Add significance comparisons
  if (!is.null(comparisons) && requireNamespace("ggsignif", quietly = TRUE)) {
    p <- p + ggsignif::geom_signif(comparisons = comparisons, map_signif_level = TRUE,
                                    textsize = 3.5, tip_length = 0.01)
  }

  p
}

#' Temporal Response Profile
#'
#' Visualize response dynamics over time with condition overlays.
#'
#' @param traces Trace matrix (cells x time)
#' @param frame_rate Frame rate in Hz
#' @param conditions Condition labels per cell
#' @param time_windows Named list of time windows for shading
#' @param show_individual Show individual cell traces
#' @param show_mean Show condition means with error bands
#' @param error_type "se", "sd", or "ci"
#'
#' @return A ggplot object
#' @export
plot_temporal_response <- function(traces, frame_rate = 10, conditions = NULL,
                                    time_windows = NULL, show_individual = TRUE,
                                    show_mean = TRUE, error_type = "se") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  traces <- as.matrix(traces)
  n_cells <- nrow(traces)
  n_time <- ncol(traces)
  time_vec <- seq(0, (n_time - 1) / frame_rate, length.out = n_time)

  if (is.null(conditions)) {
    conditions <- rep("All", n_cells)
  }
  conditions <- factor(conditions)

  # Build individual trace data
  trace_data <- do.call(rbind, lapply(seq_len(n_cells), function(i) {
    data.frame(
      cell = i,
      condition = conditions[i],
      time = time_vec,
      value = traces[i, ]
    )
  }))

  # Compute condition means and errors
  mean_data <- aggregate(value ~ condition + time, data = trace_data, FUN = mean)
  names(mean_data)[3] <- "mean"

  error_data <- aggregate(value ~ condition + time, data = trace_data,
                          FUN = function(x) {
                            if (error_type == "se") sd(x) / sqrt(length(x))
                            else if (error_type == "sd") sd(x)
                            else qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
                          })
  names(error_data)[3] <- "error"

  mean_data <- merge(mean_data, error_data)
  mean_data$lower <- mean_data$mean - mean_data$error
  mean_data$upper <- mean_data$mean + mean_data$error

  # Plot
  p <- ggplot2::ggplot()

  # Add time window shading
  if (!is.null(time_windows)) {
    window_colors <- c("#DEEBF7", "#FEE0D2", "#E5F5E0")  # Light blue, red, green
    for (i in seq_along(time_windows)) {
      w <- time_windows[[i]]
      p <- p + ggplot2::annotate(
        "rect",
        xmin = w[1], xmax = w[2],
        ymin = -Inf, ymax = Inf,
        fill = window_colors[((i - 1) %% 3) + 1],
        alpha = 0.3
      )
    }
  }

  # Individual traces
  if (show_individual) {
    p <- p + ggplot2::geom_line(
      data = trace_data,
      ggplot2::aes(x = time, y = value, group = cell, color = condition),
      alpha = 0.2, linewidth = 0.3
    )
  }

  # Mean with error band
  if (show_mean) {
    p <- p +
      ggplot2::geom_ribbon(
        data = mean_data,
        ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = condition),
        alpha = 0.3
      ) +
      ggplot2::geom_line(
        data = mean_data,
        ggplot2::aes(x = time, y = mean, color = condition),
        linewidth = 1
      )
  }

  p +
    scale_color_publication() +
    scale_fill_publication_d() +
    theme_publication() +
    ggplot2::labs(x = "Time (s)", y = "Response", color = "Condition", fill = "Condition")
}

#' Event-Triggered Average Plot
#'
#' Publication-quality peri-event time histogram.
#'
#' @param aligned_data Result from align_to_events() or matrix (cells x time)
#' @param window Window around event (c(pre, post) in seconds)
#' @param frame_rate Frame rate
#' @param conditions Optional condition labels for grouping
#' @param show_individual Show individual traces
#' @param baseline_window Frames for baseline normalization (NULL = no normalization)
#'
#' @return A ggplot object
#' @export
plot_event_triggered <- function(aligned_data, window = c(-1, 2), frame_rate = 10,
                                  conditions = NULL, show_individual = FALSE,
                                  baseline_window = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Extract data
  if (is.list(aligned_data) && !is.null(aligned_data$aligned)) {
    traces <- aligned_data$aligned
  } else {
    traces <- as.matrix(aligned_data)
  }

  n_cells <- nrow(traces)
  n_time <- ncol(traces)
  time_vec <- seq(window[1], window[2], length.out = n_time)

  # Baseline normalize
  if (!is.null(baseline_window)) {
    baseline_idx <- which(time_vec >= baseline_window[1] & time_vec <= baseline_window[2])
    if (length(baseline_idx) > 0) {
      baseline_mean <- rowMeans(traces[, baseline_idx, drop = FALSE])
      traces <- traces - baseline_mean
    }
  }

  if (is.null(conditions)) {
    conditions <- rep("All", n_cells)
  }

  # Build data
  trace_data <- do.call(rbind, lapply(seq_len(n_cells), function(i) {
    data.frame(
      cell = i,
      condition = conditions[i],
      time = time_vec,
      value = traces[i, ]
    )
  }))

  # Summary stats
  mean_data <- aggregate(value ~ condition + time, data = trace_data, FUN = mean)
  se_data <- aggregate(value ~ condition + time, data = trace_data,
                       FUN = function(x) sd(x) / sqrt(length(x)))
  mean_data$se <- se_data$value
  mean_data$lower <- mean_data$value - mean_data$se
  mean_data$upper <- mean_data$value + mean_data$se

  # Plot
  p <- ggplot2::ggplot()

  # Individual traces
  if (show_individual) {
    p <- p + ggplot2::geom_line(
      data = trace_data,
      ggplot2::aes(x = time, y = value, group = cell, color = condition),
      alpha = 0.15, linewidth = 0.3
    )
  }

  # Mean with SEM ribbon
  p <- p +
    ggplot2::geom_ribbon(
      data = mean_data,
      ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = condition),
      alpha = 0.3
    ) +
    ggplot2::geom_line(
      data = mean_data,
      ggplot2::aes(x = time, y = value, color = condition),
      linewidth = 1
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#757575", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "grey80", linewidth = 0.5) +
    scale_color_publication() +
    scale_fill_publication_d() +
    theme_publication() +
    ggplot2::labs(
      x = "Time from event (s)",
      y = expression(Delta*"F/F"),
      color = "Condition",
      fill = "Condition"
    )

  p
}

#' Dose-Response Curve
#'
#' Publication-quality dose-response visualization.
#'
#' @param dose_response Result from fit_dose_response()
#' @param show_individual Show individual cell responses
#' @param show_hill Show fitted Hill curve
#' @param log_dose Plot dose on log scale
#' @param ec50_marker Show EC50 marker
#'
#' @return A ggplot object
#' @export
plot_dose_response_pub <- function(dose_response, show_individual = TRUE,
                                    show_hill = TRUE, log_dose = TRUE,
                                    ec50_marker = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  # Extract data from result
  doses <- dose_response$doses
  responses <- dose_response$responses

  if (is.matrix(responses)) {
    mean_resp <- colMeans(responses)
    se_resp <- apply(responses, 2, function(x) sd(x) / sqrt(length(x)))

    # Individual data
    if (show_individual) {
      ind_data <- do.call(rbind, lapply(seq_len(nrow(responses)), function(i) {
        data.frame(dose = doses, response = responses[i, ], cell = i)
      }))
    }
  } else {
    mean_resp <- responses
    se_resp <- rep(0, length(responses))
    ind_data <- NULL
  }

  mean_data <- data.frame(
    dose = doses,
    response = mean_resp,
    lower = mean_resp - se_resp,
    upper = mean_resp + se_resp
  )

  # Plot
  p <- ggplot2::ggplot(mean_data, ggplot2::aes(x = dose, y = response))

  # Individual points
  if (show_individual && !is.null(ind_data)) {
    p <- p + ggplot2::geom_point(
      data = ind_data,
      ggplot2::aes(x = dose, y = response),
      alpha = 0.3, size = 1, color = "grey50"
    )
  }

  # Mean with error bars
  p <- p +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.1, linewidth = 0.6, color = "black"
    ) +
    ggplot2::geom_point(size = 3, color = "#377EB8")

  # Hill curve
  if (show_hill && !is.null(dose_response$fit_params)) {
    params <- dose_response$fit_params
    dose_fine <- seq(min(doses), max(doses), length.out = 100)
    hill_curve <- params$baseline + params$emax *
                  (dose_fine^params$n) / (dose_fine^params$n + params$ec50^params$n)

    fit_data <- data.frame(dose = dose_fine, response = hill_curve)
    p <- p + ggplot2::geom_line(data = fit_data, color = "#E41A1C", linewidth = 1)

    # EC50 marker
    if (ec50_marker) {
      ec50 <- params$ec50
      ec50_resp <- params$baseline + params$emax / 2
      p <- p +
        ggplot2::geom_segment(x = ec50, xend = ec50, y = -Inf, yend = ec50_resp,
                               linetype = "dotted", color = "#757575") +
        ggplot2::geom_segment(x = -Inf, xend = ec50, y = ec50_resp, yend = ec50_resp,
                               linetype = "dotted", color = "#757575") +
        ggplot2::annotate("text", x = ec50, y = ec50_resp, label = sprintf("EC50 = %.2f", ec50),
                           hjust = -0.1, vjust = -0.5, size = 3.5, fontface = "italic")
    }
  }

  p <- p + theme_publication() +
    ggplot2::labs(x = "Dose", y = "Response")

  if (log_dose) {
    p <- p + ggplot2::scale_x_log10()
  }

  p
}

# =============================================================================
# SUMMARY AND OVERVIEW VISUALIZATIONS
# =============================================================================

#' Create Analysis Summary Figure
#'
#' Multi-panel figure summarizing key analysis results.
#'
#' @param traces Original trace data, CaExperiment object, or trace matrix
#' @param spikes Spike detection results (auto-extracted from CaExperiment)
#' @param connectivity Connectivity results (auto-extracted from CaExperiment)
#' @param assemblies Assembly detection results (auto-extracted from CaExperiment)
#' @param frame_rate Frame rate (auto-extracted from CaExperiment)
#' @param assay Assay to use when traces is a CaExperiment (default: "dff")
#'
#' @return A patchwork object
#' @export
plot_analysis_summary <- function(traces, spikes = NULL, connectivity = NULL,
                                   assemblies = NULL, frame_rate = NULL, assay = "dff") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package required")
  }

  # Handle CaExperiment objects
  if (inherits(traces, "CaExperiment")) {
    ca <- traces
    frame_rate <- ca$frame_rate %||% 10
    traces <- GetTraces(ca, assay = assay)
    # Auto-extract available results
    if (is.null(spikes) && length(ca$spikes) > 0) {
      spikes <- GetSpikes(ca)
    }
    if (is.null(connectivity) && length(ca$graphs) > 0) {
      connectivity <- GetGraph(ca)
    }
    if (is.null(assemblies) && length(ca$assemblies) > 0) {
      assemblies <- GetAssemblies(ca)
    }
  }

  # Default frame_rate if still NULL
  if (is.null(frame_rate)) frame_rate <- 10

  plots <- list()

  # Panel A: Activity raster
  plots$A <- plot_neural_raster(traces, frame_rate = frame_rate,
                                 sort_by = "activity", show_colorbar = TRUE) +
    ggplot2::labs(title = "Population Activity")

  # Panel B: Mean activity trace
  time_vec <- seq(0, (ncol(traces) - 1) / frame_rate, length.out = ncol(traces))
  mean_trace <- data.frame(time = time_vec, activity = colMeans(traces))
  plots$B <- ggplot2::ggplot(mean_trace, ggplot2::aes(x = .data$time, y = .data$activity)) +
    ggplot2::geom_line(color = "#377EB8", linewidth = 0.6) +
    theme_publication() +
    ggplot2::labs(x = "Time (s)", y = "Mean activity", title = "Population Average")

  # Panel C: Connectivity (if provided)
  if (!is.null(connectivity)) {
    plots$C <- plot_connectivity_matrix(connectivity, order_cells = "hierarchical",
                                         title = "Functional Connectivity")
  }

  # Panel D: Assembly structure (if provided)
  if (!is.null(assemblies)) {
    plots$D <- plot_assembly_timeline(assemblies, frame_rate = frame_rate,
                                       time_range = c(0, min(60, max(time_vec))))
  }

  # Combine panels
  if (length(plots) == 4) {
    combined <- (plots$A | plots$C) / (plots$B | plots$D)
  } else if (length(plots) == 3) {
    combined <- (plots$A | plots$C) / plots$B
  } else {
    combined <- plots$A / plots$B
  }

  combined +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))
}
