#' Enhanced Visualization Functions
#'
#' Advanced visualization and plotting utilities for calcium imaging data.
#'
#' @name enhanced_visualization
NULL

#' Interactive Calcium Trace Plot
#'
#' Create an interactive plot of calcium traces using plotly.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param cell_ids Vector of cell identifiers (default: NULL, use column names)
#' @param time_points Vector of time points (default: NULL, use 1:ncol(traces))
#' @param selected_cells Vector of cells to highlight (default: NULL)
#' @param color_palette Color palette for cells (default: "viridis")
#' @param ... Additional arguments
#' @return plotly object
#' @export
interactive_calcium_plot <- function(traces, cell_ids = NULL, time_points = NULL, selected_cells = NULL, color_palette = "viridis", ...) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for interactive plots. Install with: install.packages('plotly')")
  }
  
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for color palettes. Install with: install.packages('viridis')")
  }
  
  # Prepare data
  if (is.null(cell_ids)) {
    cell_ids <- if (!is.null(colnames(traces))) colnames(traces) else paste0("Cell_", 1:ncol(traces))
  }
  
  if (is.null(time_points)) {
    time_points <- 1:nrow(traces)
  }
  
  # Create plot
  p <- plotly::plot_ly()
  
  # Add traces for each cell
  colors <- viridis::viridis(length(cell_ids))
  
  for (i in 1:length(cell_ids)) {
    cell_id <- cell_ids[i]
    trace_data <- traces[, i]
    
    # Determine line width and opacity based on selection
    line_width <- if (cell_id %in% selected_cells) 3 else 1
    opacity <- if (cell_id %in% selected_cells) 1 else 0.7
    
    p <- plotly::add_trace(
      p,
      x = time_points,
      y = trace_data,
      type = "scatter",
      mode = "lines",
      name = cell_id,
      line = list(color = colors[i], width = line_width),
      opacity = opacity,
      hovertemplate = paste(
        "<b>%{fullData.name}</b><br>",
        "Time: %{x}<br>",
        "Signal: %{y:.3f}<br>",
        "<extra></extra>"
      )
    )
  }
  
  # Update layout
  p <- plotly::layout(
    p,
    title = "Interactive Calcium Traces",
    xaxis = list(title = "Time", showgrid = TRUE),
    yaxis = list(title = "Calcium Signal (Delta F/F)", showgrid = TRUE),
    hovermode = "closest",
    showlegend = TRUE,
    legend = list(
      orientation = "h",
      x = 0.5,
      y = -0.1
    )
  )
  
  return(p)
}

#' 3D Calcium Imaging Visualization
#'
#' Create a 3D visualization of calcium imaging data with spatial information.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param spatial_coords Matrix of spatial coordinates (cells x 2 or 3)
#' @param time_points Vector of time points (default: NULL)
#' @param cell_ids Vector of cell identifiers (default: NULL)
#' @param color_by Variable to color by ("time", "cell", "activity")
#' @param ... Additional arguments
#' @return plotly 3D object
#' @export
visualize_3d_calcium <- function(traces, spatial_coords, time_points = NULL, cell_ids = NULL, color_by = "activity", ...) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required for 3D visualization")
  }
  
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("Package 'viridis' is required for color palettes")
  }
  
  # Prepare data
  if (is.null(time_points)) {
    time_points <- 1:ncol(traces)
  }
  
  if (is.null(cell_ids)) {
    cell_ids <- paste0("Cell_", 1:nrow(traces))
  }
  
  # Ensure spatial coordinates are 3D
  if (ncol(spatial_coords) == 2) {
    spatial_coords <- cbind(spatial_coords, 0)
  }
  
  # Create 3D scatter plot
  if (color_by == "activity") {
    # Color by mean activity
    activity_levels <- rowMeans(traces, na.rm = TRUE)
    colors <- viridis::viridis(100)[cut(activity_levels, 100)]
  } else if (color_by == "cell") {
    colors <- viridis::viridis(nrow(traces))
  } else {
    colors <- rep("blue", nrow(traces))
  }
  
  p <- plotly::plot_ly(
    x = spatial_coords[, 1],
    y = spatial_coords[, 2],
    z = spatial_coords[, 3],
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 8,
      color = colors,
      opacity = 0.8
    ),
    text = cell_ids,
    hovertemplate = paste(
      "<b>%{text}</b><br>",
      "X: %{x:.2f}<br>",
      "Y: %{y:.2f}<br>",
      "Z: %{z:.2f}<br>",
      "<extra></extra>"
    )
  )
  
  # Update layout
  p <- plotly::layout(
    p,
    title = "3D Calcium Imaging Visualization",
    scene = list(
      xaxis = list(title = "X Position"),
      yaxis = list(title = "Y Position"),
      zaxis = list(title = "Z Position")
    ),
    showlegend = FALSE
  )
  
  return(p)
}

#' Publication-Ready Calcium Trace Plot
#'
#' Create a publication-ready plot of calcium traces using ggplot2.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param cell_ids Vector of cell identifiers (default: NULL)
#' @param time_points Vector of time points (default: NULL)
#' @param selected_cells Vector of cells to highlight (default: NULL)
#' @param color_palette Color palette (default: "viridis")
#' @param line_size Line size (default: 0.5)
#' @param alpha Transparency (default: 0.8)
#' @param theme Theme to use (default: "minimal")
#' @param ... Additional arguments
#' @return ggplot object
#' @export
publication_calcium_plot <- function(traces, cell_ids = NULL, time_points = NULL, selected_cells = NULL, color_palette = "viridis", line_size = 0.5, alpha = 0.8, theme = "minimal", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for publication plots")
  }
  
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("viridis is required for color palettes")
  }
  
  # Prepare data
  if (is.null(cell_ids)) {
    cell_ids <- if (!is.null(colnames(traces))) colnames(traces) else paste0("Cell_", 1:ncol(traces))
  }
  
  if (is.null(time_points)) {
    time_points <- 1:nrow(traces)
  }
  
  # Create long format data
  plot_data <- data.frame()
  for (i in 1:length(cell_ids)) {
    cell_data <- data.frame(
      Time = time_points,
      Signal = traces[, i],
      Cell = cell_ids[i],
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, cell_data)
  }
  
  # Add selection indicator
  plot_data$Selected <- plot_data$Cell %in% selected_cells
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Time, y = Signal, color = Cell, alpha = Selected)) +
    ggplot2::geom_line(size = line_size) +
    viridis::scale_color_viridis(discrete = TRUE) +
    ggplot2::scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = alpha)) +
    ggplot2::labs(
      title = "Calcium Imaging Traces",
      x = "Time",
      y = "Calcium Signal (Delta F/F)",
      color = "Cell"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 16, hjust = 0.5)
    )
  
  return(p)
}

#' Heatmap of Calcium Activity
#'
#' Create a heatmap visualization of calcium activity across cells and time.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param cell_ids Vector of cell identifiers (default: NULL)
#' @param time_points Vector of time points (default: NULL)
#' @param color_scale Color scale (default: "viridis")
#' @param cluster_cells Whether to cluster cells (default: TRUE)
#' @param cluster_time Whether to cluster time points (default: FALSE)
#' @param ... Additional arguments
#' @return ggplot object
#' @export
calcium_activity_heatmap <- function(traces, cell_ids = NULL, time_points = NULL, color_scale = "viridis", cluster_cells = TRUE, cluster_time = FALSE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for heatmap")
  }
  
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("viridis is required for color scales")
  }
  
  # Prepare data
  if (is.null(cell_ids)) {
    cell_ids <- if (!is.null(colnames(traces))) colnames(traces) else paste0("Cell_", 1:ncol(traces))
  }
  
  if (is.null(time_points)) {
    time_points <- 1:nrow(traces)
  }
  
  # Cluster if requested
  if (cluster_cells) {
    cell_order <- hclust(dist(t(traces)))$order
    cell_ids <- cell_ids[cell_order]
    traces <- traces[, cell_order]
  }
  
  if (cluster_time) {
    time_order <- hclust(dist(traces))$order
    time_points <- time_points[time_order]
    traces <- traces[time_order, ]
  }
  
  # Create long format data
  plot_data <- expand.grid(
    Time = time_points,
    Cell = cell_ids
  )
  plot_data$Signal <- as.vector(traces)
  
  # Create heatmap
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Time, y = Cell, fill = Signal)) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis() +
    ggplot2::labs(
      title = "Calcium Activity Heatmap",
      x = "Time",
      y = "Cell",
      fill = "Signal"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      axis.text.x = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
      legend.position = "right"
    )
  
  return(p)
}

#' Interactive Network Visualization
#'
#' Create an interactive network visualization using plotly.
#'
#' @param adjacency_matrix Adjacency matrix of the network
#' @param node_attributes Data frame with node attributes (optional)
#' @param edge_attributes Data frame with edge attributes (optional)
#' @param layout Network layout ("spring", "circular", "random")
#' @param node_size Variable for node size (default: "degree")
#' @param node_color Variable for node color (default: "community")
#' @param ... Additional arguments
#' @return plotly object
#' @export
interactive_network_plot <- function(adjacency_matrix, node_attributes = NULL, edge_attributes = NULL, layout = "spring", node_size = "degree", node_color = "community", ...) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("plotly is required for interactive network plots")
  }
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required for network analysis")
  }
  
  # Create igraph object
  g <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE)
  
  # Calculate layout
  if (layout == "spring") {
    coords <- igraph::layout_with_fr(g)
  } else if (layout == "circular") {
    coords <- igraph::layout_in_circle(g)
  } else {
    coords <- igraph::layout_randomly(g)
  }
  
  # Prepare node data
  node_data <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    id = 1:igraph::vcount(g),
    degree = igraph::degree(g)
  )
  
  # Add node attributes if provided
  if (!is.null(node_attributes)) {
    node_data <- cbind(node_data, node_attributes)
  }
  
  # Prepare edge data
  edges <- igraph::get.edgelist(g)
  edge_data <- data.frame(
    x0 = coords[edges[, 1], 1],
    y0 = coords[edges[, 1], 2],
    x1 = coords[edges[, 2], 1],
    y1 = coords[edges[, 2], 2],
    weight = igraph::E(g)$weight
  )
  
  # Create plot
  p <- plotly::plot_ly()
  
  # Add edges
  p <- plotly::add_segments(
    p,
    data = edge_data,
    x = ~x0, y = ~y0,
    xend = ~x1, yend = ~y1,
    line = list(color = "gray", width = 1),
    showlegend = FALSE,
    hoverinfo = "skip"
  )
  
  # Add nodes
  p <- plotly::add_markers(
    p,
    data = node_data,
    x = ~x, y = ~y,
    size = ~degree,
    color = ~degree,
    colors = "viridis",
    text = ~id,
    hovertemplate = paste(
      "<b>Node %{text}</b><br>",
      "Degree: %{marker.size}<br>",
      "<extra></extra>"
    )
  )
  
  # Update layout
  p <- plotly::layout(
    p,
    title = "Interactive Network Visualization",
    xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
    showlegend = FALSE
  )
  
  return(p)
}

#' Multi-Panel Calcium Analysis Plot
#'
#' Create a multi-panel figure showing various aspects of calcium imaging analysis.
#'
#' @param traces Matrix of calcium traces (cells x time)
#' @param spike_results List containing spike detection results (optional)
#' @param network_results List containing network analysis results (optional)
#' @param cell_ids Vector of cell identifiers (default: NULL)
#' @param time_points Vector of time points (default: NULL)
#' @param ... Additional arguments
#' @return ggplot object with multiple panels
#' @export
multi_panel_calcium_plot <- function(traces, spike_results = NULL, network_results = NULL, cell_ids = NULL, time_points = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for multi-panel plots")
  }
  
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra is required for multi-panel layouts")
  }
  
  # Prepare data
  if (is.null(cell_ids)) {
    cell_ids <- if (!is.null(colnames(traces))) colnames(traces) else paste0("Cell_", 1:ncol(traces))
  }
  
  if (is.null(time_points)) {
    time_points <- 1:nrow(traces)
  }
  
  # Create individual plots
  plots <- list()
  
  # Panel 1: Calcium traces
  plots[[1]] <- publication_calcium_plot(traces, cell_ids, time_points, ...) +
    ggplot2::labs(title = "A) Calcium Traces")
  
  # Panel 2: Activity heatmap
  plots[[2]] <- calcium_activity_heatmap(traces, cell_ids, time_points, ...) +
    ggplot2::labs(title = "B) Activity Heatmap")
  
  # Panel 3: Spike raster (if available)
  if (!is.null(spike_results)) {
    spikes <- spike_results$spikes
    # Ensure spikes is a matrix with correct dimensions (cells x time)
    if (is.null(dim(spikes))) {
      # If vector, convert to matrix with one row
      spikes <- matrix(spikes, nrow = 1)
    }
    if (nrow(spikes) != length(cell_ids) && ncol(spikes) == length(cell_ids)) {
      # Transpose if needed
      spikes <- t(spikes)
    }
    if (nrow(spikes) != length(cell_ids)) {
      stop("spike_results$spikes must have one row per cell")
    }
    if (ncol(spikes) != length(time_points)) {
      stop("spike_results$spikes must have one column per time point")
    }
    # Build spike data frame
    spike_data <- data.frame(
      Time = rep(time_points, each = length(cell_ids)),
      Cell = rep(cell_ids, times = length(time_points)),
      Spike = as.vector(spikes)
    )
    plots[[3]] <- ggplot2::ggplot(spike_data[spike_data$Spike == 1, ], 
                                 ggplot2::aes(x = Time, y = Cell)) +
      ggplot2::geom_point(size = 0.5) +
      ggplot2::labs(title = "C) Spike Raster") +
      ggplot2::theme_minimal() +
      ggplot2::scale_y_discrete(limits = rev(cell_ids))
  }
  
  # Panel 4: Network visualization (if available)
  if (!is.null(network_results)) {
    # Create network plot
    plots[[4]] <- ggplot2::ggplot() +
      ggplot2::geom_blank() +
      ggplot2::labs(title = "D) Network Connectivity") +
      ggplot2::theme_minimal()
  }
  
  # Combine plots
  n_plots <- length(plots)
  if (n_plots == 2) {
    combined_plot <- gridExtra::grid.arrange(plots[[1]], plots[[2]], ncol = 2)
  } else if (n_plots == 3) {
    combined_plot <- gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol = 3)
  } else if (n_plots == 4) {
    combined_plot <- gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2)
  } else {
    combined_plot <- plots[[1]]
  }
  
  return(combined_plot)
}

#' Save Publication-Ready Figure
#'
#' Save a ggplot object as a publication-ready figure with proper formatting.
#'
#' @param plot ggplot object to save
#' @param filename Output filename
#' @param width Figure width in inches (default: 10)
#' @param height Figure height in inches (default: 8)
#' @param dpi Resolution in DPI (default: 300)
#' @param format Output format ("pdf", "png", "tiff", "svg")
#' @param ... Additional arguments
#' @return Path to saved file
#' @export
save_publication_figure <- function(plot, filename, width = 10, height = 8, dpi = 300, format = "pdf", ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for saving figures")
  }
  
  # Ensure filename has correct extension
  if (!grepl(paste0("\\.", format, "$"), filename)) {
    filename <- paste0(filename, ".", format)
  }
  
  # Save figure
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    ...
  )
  
  message("Figure saved: ", filename)
  return(filename)
} 