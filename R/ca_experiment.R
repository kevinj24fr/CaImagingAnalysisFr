#' CaExperiment: Unified Calcium Imaging Analysis Object
#'
#' A Seurat-style unified object for calcium imaging analysis workflows.
#' Stores all data, analysis results, and maintains full provenance tracking.
#'
#' @name CaExperiment
#' @aliases CaExperiment-class
NULL

# ============================================================================
# CaExperiment S7 Class Definition
# ============================================================================

#' Create a CaExperiment Object
#'
#' A unified container for calcium imaging analysis that stores raw data,
#' processed traces, spike inference results, dimensionality reductions,
#' network graphs, and metadata in a single object with full provenance tracking.
#'
#' @param traces A numeric matrix of calcium traces (cells x timepoints) or CalciumTraces object
#' @param movie Optional CalciumMovie object
#' @param rois Optional ROISet object with cell spatial information
#' @param meta.data Optional data frame with cell-level metadata
#' @param frame_rate Imaging frame rate in Hz
#' @param experiment_id Unique identifier for the experiment
#' @param ... Additional metadata to store
#'
#' @return A CaExperiment object
#' @export
#'
#' @examples
#' \dontrun
#' # Create from trace matrix
#' traces <- matrix(rnorm(100 * 1000), nrow = 100)
#' ca <- CaExperiment(traces, frame_rate = 30)
#'
#' # Analysis pipeline with method chaining
#' ca <- ca |>
#'   RunCorrection(method = "neuropil") |>
#'   RunDFF() |>
#'   RunSpikes(method = "oasis") |>
#'   RunPCA(n_components = 20) |>
#'   RunConnectivity() |>
#'   RunAssemblies()
#'
#' # Access results
#' GetTraces(ca, assay = "dff")
#' GetSpikes(ca)
#' GetReduction(ca, "pca")
#' }
CaExperiment <- function(traces,
                         movie = NULL,
                         rois = NULL,
                         meta.data = NULL,
                         frame_rate = 30,
                         experiment_id = NULL,
                         ...) {

  # Handle input types
  if (inherits(traces, "CalciumTraces")) {
    trace_matrix <- traces@traces
    if (is.null(frame_rate) || frame_rate == 30) {
      frame_rate <- traces@frame_rate
    }
  } else if (is.matrix(traces)) {
    trace_matrix <- traces
  } else if (is.data.frame(traces)) {
    trace_matrix <- as.matrix(traces)
  } else {
    stop("traces must be a matrix, data.frame, or CalciumTraces object")
  }

  # Validate dimensions
  if (!is.numeric(trace_matrix)) {
    stop("traces must contain numeric values")
  }

  n_cells <- nrow(trace_matrix)
  n_frames <- ncol(trace_matrix)

  # Generate cell IDs if not present
  if (is.null(rownames(trace_matrix))) {
    rownames(trace_matrix) <- paste0("cell_", seq_len(n_cells))
  }
  cell_ids <- rownames(trace_matrix)

  # Generate time vector
  time_vec <- seq(0, (n_frames - 1) / frame_rate, length.out = n_frames)

  # Initialize metadata
  if (is.null(meta.data)) {
    meta.data <- data.frame(
      cell_id = cell_ids,
      row.names = cell_ids,
      stringsAsFactors = FALSE
    )
  } else {
    if (nrow(meta.data) != n_cells) {
      stop("meta.data must have same number of rows as cells in traces")
    }
    rownames(meta.data) <- cell_ids
  }

  # Generate experiment ID
  if (is.null(experiment_id)) {
    experiment_id <- paste0("ca_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }

  # Create assays list with raw data
  assays <- list(
    raw = trace_matrix
  )

  # Initialize empty containers
  spikes <- list()
  reductions <- list()
  graphs <- list()
  clusters <- list()
  transients <- list()
  assemblies <- list()
  sequences <- list()
  trajectories <- list()
  tuning <- list()

  # Command log for provenance
  commands <- list(
    list(
      command = "CaExperiment",
      timestamp = Sys.time(),
      args = list(
        n_cells = n_cells,
        n_frames = n_frames,
        frame_rate = frame_rate
      )
    )
  )

  # Additional info from ...
  misc <- list(...)

  # Create the object
  obj <- structure(
    list(
      assays = assays,
      spikes = spikes,
      reductions = reductions,
      graphs = graphs,
      clusters = clusters,
      transients = transients,
      assemblies = assemblies,
      sequences = sequences,
      trajectories = trajectories,
      tuning = tuning,
      meta.data = meta.data,
      commands = commands,
      misc = misc,
      movie = movie,
      rois = rois,
      frame_rate = frame_rate,
      time = time_vec,
      experiment_id = experiment_id,
      active_assay = "raw",
      version = packageVersion("CaImagingAnalysisFr")
    ),
    class = "CaExperiment"
  )

  obj
}

# ============================================================================
# Print and Summary Methods
# ============================================================================

#' @export
print.CaExperiment <- function(x, ...) {
  n_cells <- nrow(x$assays[[x$active_assay]])
  n_frames <- ncol(x$assays[[x$active_assay]])
  duration <- n_frames / x$frame_rate


  cat("CaExperiment object\n")
  cat("-------------------\n")
  cat(sprintf("  Cells: %d\n", n_cells))
  cat(sprintf("  Frames: %d (%.1f sec @ %.1f Hz)\n", n_frames, duration, x$frame_rate))
  cat(sprintf("  Experiment ID: %s\n", x$experiment_id))
  cat("\n")

  # Assays
  cat(sprintf("Assays (%d): %s\n", length(x$assays), paste(names(x$assays), collapse = ", ")))
  cat(sprintf("  Active: %s\n", x$active_assay))

  # Spikes
  if (length(x$spikes) > 0) {
    cat(sprintf("Spikes (%d): %s\n", length(x$spikes), paste(names(x$spikes), collapse = ", ")))
  }

  # Reductions
  if (length(x$reductions) > 0) {
    cat(sprintf("Reductions (%d): %s\n", length(x$reductions), paste(names(x$reductions), collapse = ", ")))
  }

  # Graphs
  if (length(x$graphs) > 0) {
    cat(sprintf("Graphs (%d): %s\n", length(x$graphs), paste(names(x$graphs), collapse = ", ")))
  }

  # Clusters
  if (length(x$clusters) > 0) {
    cat(sprintf("Clusters (%d): %s\n", length(x$clusters), paste(names(x$clusters), collapse = ", ")))
  }

  # Assemblies

  if (length(x$assemblies) > 0) {
    cat(sprintf("Assemblies (%d): %s\n", length(x$assemblies), paste(names(x$assemblies), collapse = ", ")))
  }

  # Metadata columns
  meta_cols <- setdiff(names(x$meta.data), "cell_id")
  if (length(meta_cols) > 0) {
    cat(sprintf("Metadata (%d): %s\n", length(meta_cols),
                paste(head(meta_cols, 5), collapse = ", ")))
    if (length(meta_cols) > 5) cat("  ...\n")
  }

  # Commands
  cat(sprintf("\nCommands logged: %d\n", length(x$commands)))

  invisible(x)
}

#' @export
summary.CaExperiment <- function(object, ...) {
  cat("=== CaExperiment Summary ===\n\n")

  print(object)

  cat("\n--- Trace Statistics ---\n")
  traces <- object$assays[[object$active_assay]]
  cat(sprintf("  Mean: %.3f\n", mean(traces, na.rm = TRUE)))
  cat(sprintf("  SD: %.3f\n", sd(traces, na.rm = TRUE)))
  cat(sprintf("  Range: [%.3f, %.3f]\n", min(traces, na.rm = TRUE), max(traces, na.rm = TRUE)))

  if (length(object$spikes) > 0) {
    cat("\n--- Spike Summary ---\n")
    for (name in names(object$spikes)) {
      spk <- object$spikes[[name]]
      if (is.matrix(spk)) {
        total_spikes <- sum(spk > 0, na.rm = TRUE)
        mean_rate <- total_spikes / (ncol(spk) / object$frame_rate) / nrow(spk)
        cat(sprintf("  %s: %.2f Hz mean rate\n", name, mean_rate))
      }
    }
  }

  invisible(object)
}

# ============================================================================
# Accessor Functions (Getters)
# ============================================================================

#' Get Traces from CaExperiment
#'
#' @param object CaExperiment object
#' @param assay Which assay to retrieve (default: active assay)
#' @param cells Cell IDs or indices to subset (default: all)
#' @param frames Frame indices to subset (default: all)
#'
#' @return Numeric matrix of traces
#' @export
GetTraces <- function(object, assay = NULL, cells = NULL, frames = NULL) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.null(assay)) {
    assay <- object$active_assay
  }

  if (!assay %in% names(object$assays)) {
    stop(sprintf("Assay '%s' not found. Available: %s",
                 assay, paste(names(object$assays), collapse = ", ")))
  }

  traces <- object$assays[[assay]]

  # Subset cells
  if (!is.null(cells)) {
    if (is.character(cells)) {
      cells <- match(cells, rownames(traces))
    }
    traces <- traces[cells, , drop = FALSE]
  }

  # Subset frames
  if (!is.null(frames)) {
    traces <- traces[, frames, drop = FALSE]
  }

  traces
}

#' Get Spike Data from CaExperiment
#'
#' @param object CaExperiment object
#' @param method Which spike inference method (default: first available)
#' @param cells Cell IDs or indices to subset
#'
#' @return Numeric matrix of spike data
#' @export
GetSpikes <- function(object, method = NULL, cells = NULL) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (length(object$spikes) == 0) {
    stop("No spike data available. Run RunSpikes() first.")
  }

  if (is.null(method)) {
    method <- names(object$spikes)[1]
  }

  if (!method %in% names(object$spikes)) {
    stop(sprintf("Spike method '%s' not found. Available: %s",
                 method, paste(names(object$spikes), collapse = ", ")))
  }

  spikes <- object$spikes[[method]]

  if (!is.null(cells)) {
    if (is.character(cells)) {
      cells <- match(cells, rownames(spikes))
    }
    spikes <- spikes[cells, , drop = FALSE]
  }

  spikes
}

#' Get Dimensionality Reduction from CaExperiment
#'
#' @param object CaExperiment object
#' @param reduction Name of reduction (e.g., "pca", "umap", "gpfa")
#'
#' @return Reduction result object
#' @export
GetReduction <- function(object, reduction) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (!reduction %in% names(object$reductions)) {
    stop(sprintf("Reduction '%s' not found. Available: %s",
                 reduction, paste(names(object$reductions), collapse = ", ")))
  }

  object$reductions[[reduction]]
}

#' Get Embeddings from Reduction
#'
#' @param object CaExperiment object
#' @param reduction Name of reduction
#' @param dims Which dimensions to return
#'
#' @return Matrix of embeddings (cells x dimensions)
#' @export
GetEmbeddings <- function(object, reduction, dims = NULL) {
  red <- GetReduction(object, reduction)

  embeddings <- if ("embeddings" %in% names(red)) {
    red$embeddings
  } else if ("x" %in% names(red)) {
    red$x
  } else if ("cell_scores" %in% names(red)) {
    red$cell_scores
  } else {
    stop("Could not find embeddings in reduction object")
  }

  if (!is.null(dims)) {
    embeddings <- embeddings[, dims, drop = FALSE]
  }

  embeddings
}

#' Get Graph from CaExperiment
#'
#' @param object CaExperiment object
#' @param graph Name of graph
#'
#' @return Graph object (matrix or igraph)
#' @export
GetGraph <- function(object, graph) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (!graph %in% names(object$graphs)) {
    stop(sprintf("Graph '%s' not found. Available: %s",
                 graph, paste(names(object$graphs), collapse = ", ")))
  }

  object$graphs[[graph]]
}

#' Get Cell Metadata from CaExperiment
#'
#' @param object CaExperiment object
#' @param vars Variables to retrieve (default: all)
#'
#' @return Data frame of metadata
#' @export
GetMetaData <- function(object, vars = NULL) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  meta <- object$meta.data

  if (!is.null(vars)) {
    missing <- setdiff(vars, names(meta))
    if (length(missing) > 0) {
      warning(sprintf("Variables not found: %s", paste(missing, collapse = ", ")))
    }
    vars <- intersect(vars, names(meta))
    meta <- meta[, vars, drop = FALSE]
  }

  meta
}

#' Get Assemblies from CaExperiment
#'
#' @param object CaExperiment object
#' @param name Name of assembly detection result
#'
#' @return Assembly detection result
#' @export
GetAssemblies <- function(object, name = NULL) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (length(object$assemblies) == 0) {
    stop("No assembly data available. Run RunAssemblies() first.")
  }

  if (is.null(name)) {
    name <- names(object$assemblies)[1]
  }

  object$assemblies[[name]]
}

#' Get Transient Events from CaExperiment
#'
#' @param object CaExperiment object
#' @param name Name of transient detection result
#'
#' @return Transient detection result
#' @export
GetTransients <- function(object, name = NULL) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (length(object$transients) == 0) {
    stop("No transient data available. Run RunTransients() first.")
  }

  if (is.null(name)) {
    name <- names(object$transients)[1]
  }

  object$transients[[name]]
}

#' Get Tuning Data from CaExperiment
#'
#' @param object CaExperiment object
#' @param type Type of tuning ("orientation", "contrast", "place", etc.)
#'
#' @return Tuning analysis result
#' @export
GetTuning <- function(object, type = NULL) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (length(object$tuning) == 0) {
    stop("No tuning data available. Run RunTuning() first.")
  }

  if (is.null(type)) {
    type <- names(object$tuning)[1]
  }

  object$tuning[[type]]
}

#' Get Command History
#'
#' @param object CaExperiment object
#' @param n Number of recent commands to show (default: all)
#'
#' @return List of commands
#' @export
GetCommands <- function(object, n = NULL) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  commands <- object$commands

  if (!is.null(n)) {
    n <- min(n, length(commands))
    commands <- commands[(length(commands) - n + 1):length(commands)]
  }

  commands
}

# ============================================================================
# Setter Functions
# ============================================================================

#' Set Active Assay
#'
#' @param object CaExperiment object
#' @param assay Name of assay to set as active
#'
#' @return Modified CaExperiment object
#' @export
SetActiveAssay <- function(object, assay) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (!assay %in% names(object$assays)) {
    stop(sprintf("Assay '%s' not found", assay))
  }

  object$active_assay <- assay
  object
}

#' Add Metadata to CaExperiment
#'
#' @param object CaExperiment object
#' @param metadata Data frame or named vector of metadata
#' @param col_name Column name if metadata is a vector
#'
#' @return Modified CaExperiment object
#' @export
AddMetaData <- function(object, metadata, col_name = NULL) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.vector(metadata) && !is.list(metadata)) {
    if (is.null(col_name)) {
      stop("col_name required when metadata is a vector")
    }
    if (length(metadata) != nrow(object$meta.data)) {
      stop("metadata vector length must match number of cells")
    }
    object$meta.data[[col_name]] <- metadata
  } else if (is.data.frame(metadata)) {
    if (nrow(metadata) != nrow(object$meta.data)) {
      stop("metadata rows must match number of cells")
    }
    for (col in names(metadata)) {
      object$meta.data[[col]] <- metadata[[col]]
    }
  } else {
    stop("metadata must be a vector or data.frame")
  }

  object
}

# ============================================================================
# Subsetting
# ============================================================================

#' Subset CaExperiment Object
#'
#' @param x CaExperiment object
#' @param i Cell indices or names
#' @param j Frame indices
#' @param ... Additional arguments (ignored)
#' @param drop Ignored
#'
#' @return Subsetted CaExperiment object
#' @export
`[.CaExperiment` <- function(x, i = NULL, j = NULL, ..., drop = FALSE) {

  # Get current dimensions
  n_cells <- nrow(x$assays[[x$active_assay]])
  n_frames <- ncol(x$assays[[x$active_assay]])
  cell_ids <- rownames(x$assays[[x$active_assay]])

  # Handle cell selection
  if (is.null(i)) {
    i <- seq_len(n_cells)
  } else if (is.character(i)) {
    i <- match(i, cell_ids)
    i <- i[!is.na(i)]
  } else if (is.logical(i)) {
    i <- which(i)
  }

  # Handle frame selection
  if (is.null(j)) {
    j <- seq_len(n_frames)
  } else if (is.logical(j)) {
    j <- which(j)
  }

  # Subset assays
  new_assays <- lapply(x$assays, function(a) {
    a[i, j, drop = FALSE]
  })

  # Subset spikes
  new_spikes <- lapply(x$spikes, function(s) {
    if (is.matrix(s)) {
      s[i, j, drop = FALSE]
    } else {
      s  # Keep as-is if not a matrix
    }
  })

  # Subset metadata
  new_meta <- x$meta.data[i, , drop = FALSE]

  # Subset time vector
  new_time <- x$time[j]

  # Subset reductions (cell dimension only)
  new_reductions <- lapply(x$reductions, function(r) {
    if (is.list(r)) {
      if ("embeddings" %in% names(r)) {
        r$embeddings <- r$embeddings[i, , drop = FALSE]
      }
      if ("x" %in% names(r)) {
        r$x <- r$x[i, , drop = FALSE]
      }
      if ("cell_scores" %in% names(r)) {
        r$cell_scores <- r$cell_scores[i, , drop = FALSE]
      }
    }
    r
  })

  # Subset graphs
  new_graphs <- lapply(x$graphs, function(g) {
    if (is.matrix(g)) {
      g[i, i, drop = FALSE]
    } else {
      g  # Keep igraph objects as-is (would need special handling)
    }
  })

  # Log the subset operation
  new_commands <- c(x$commands, list(list(
    command = "subset",
    timestamp = Sys.time(),
    args = list(
      cells = length(i),
      frames = length(j)
    )
  )))

  # Create new object
  structure(
    list(
      assays = new_assays,
      spikes = new_spikes,
      reductions = new_reductions,
      graphs = new_graphs,
      clusters = x$clusters,
      transients = x$transients,
      assemblies = x$assemblies,
      sequences = x$sequences,
      trajectories = x$trajectories,
      tuning = x$tuning,
      meta.data = new_meta,
      commands = new_commands,
      misc = x$misc,
      movie = x$movie,
      rois = x$rois,
      frame_rate = x$frame_rate,
      time = new_time,
      experiment_id = x$experiment_id,
      active_assay = x$active_assay,
      version = x$version
    ),
    class = "CaExperiment"
  )
}

#' Subset Cells by Metadata
#'
#' @param object CaExperiment object
#' @param ... Logical expressions for filtering (e.g., cell_type == "excitatory")
#'
#' @return Subsetted CaExperiment object
#' @export
SubsetCells <- function(object, ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  # Evaluate filter in context of metadata
  filter_expr <- substitute(...)
  keep <- eval(filter_expr, object$meta.data, parent.frame())

  if (!is.logical(keep)) {
    stop("Filter expression must evaluate to logical")
  }

  object[keep, ]
}

# ============================================================================
# Dimension Functions
# ============================================================================

#' Get Number of Cells
#'
#' @param object CaExperiment object
#'
#' @return Integer number of cells
#' @export
ncells <- function(object) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }
  nrow(object$assays[[object$active_assay]])
}

#' Get Number of Frames
#'
#' @param object CaExperiment object
#'
#' @return Integer number of frames
#' @export
nframes <- function(object) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }
  ncol(object$assays[[object$active_assay]])
}

#' Get Cell IDs
#'
#' @param object CaExperiment object
#'
#' @return Character vector of cell IDs
#' @export
CellIDs <- function(object) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }
  rownames(object$assays[[object$active_assay]])
}

#' Get Recording Duration
#'
#' @param object CaExperiment object
#'
#' @return Duration in seconds
#' @export
Duration <- function(object) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }
  nframes(object) / object$frame_rate
}

# ============================================================================
# Run* Analysis Wrapper Functions
# ============================================================================

#' Log Command to CaExperiment
#'
#' @param object CaExperiment object
#' @param command Command name
#' @param args List of arguments
#'
#' @return Modified object with logged command
#' @keywords internal
log_command <- function(object, command, args = list()) {
  object$commands <- c(object$commands, list(list(
    command = command,
    timestamp = Sys.time(),
    args = args
  )))
  object
}

#' Run Neuropil Correction
#'
#' @param object CaExperiment object
#' @param neuropil_traces Matrix of neuropil traces (optional)
#' @param coefficient Neuropil coefficient (default: 0.7)
#' @param method Correction method
#' @param assay_name Name for the corrected assay
#'
#' @return Modified CaExperiment object
#' @export
RunCorrection <- function(object,
                          neuropil_traces = NULL,
                          coefficient = 0.7,
                          method = c("subtract", "regression"),
                          assay_name = "corrected") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  method <- match.arg(method)
  traces <- GetTraces(object, assay = "raw")

  if (is.null(neuropil_traces)) {
    # If no neuropil provided, just copy raw
    corrected <- traces
    message("No neuropil traces provided; using raw traces as corrected")
  } else {
    corrected <- neuropil_correct(traces, neuropil_traces,
                                  coefficient = coefficient,
                                  method = method)
  }

  object$assays[[assay_name]] <- corrected
  object$active_assay <- assay_name

  log_command(object, "RunCorrection", list(
    method = method,
    coefficient = coefficient,
    assay_name = assay_name
  ))
}

#' Compute dF/F
#'
#' @param object CaExperiment object
#' @param method Baseline method ("rolling", "percentile", "mode")
#' @param window_size Window size for rolling baseline
#' @param percentile Percentile for baseline
#' @param assay Source assay (default: active)
#' @param assay_name Name for dF/F assay
#'
#' @return Modified CaExperiment object
#' @export
RunDFF <- function(object,
                   method = c("rolling", "percentile", "mode"),
                   window_size = 300,
                   percentile = 0.08,
                   assay = NULL,
                   assay_name = "dff") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  method <- match.arg(method)

  if (is.null(assay)) {
    assay <- object$active_assay
  }

  traces <- GetTraces(object, assay = assay)

 # Compute dF/F based on method
  n_cells <- nrow(traces)
  n_time <- ncol(traces)
  dff <- matrix(0, nrow = n_cells, ncol = n_time)
  rownames(dff) <- rownames(traces)
  colnames(dff) <- colnames(traces)

  for (i in seq_len(n_cells)) {
    trace <- traces[i, ]

    # Compute baseline based on method
    baseline <- switch(method,
      rolling = {
        # Rolling minimum baseline
        half_win <- window_size %/% 2
        bl <- numeric(n_time)
        for (t in seq_len(n_time)) {
          start <- max(1, t - half_win)
          end <- min(n_time, t + half_win)
          bl[t] <- quantile(trace[start:end], percentile, na.rm = TRUE)
        }
        bl
      },
      percentile = {
        # Global percentile baseline
        rep(quantile(trace, percentile, na.rm = TRUE), n_time)
      },
      mode = {
        # Mode-based baseline (histogram peak)
        h <- hist(trace, breaks = 50, plot = FALSE)
        mode_val <- h$mids[which.max(h$counts)]
        rep(mode_val, n_time)
      }
    )

    # Compute dF/F, avoiding division by zero
    baseline[baseline < 1e-10] <- 1e-10
    dff[i, ] <- (trace - baseline) / baseline
  }

  object$assays[[assay_name]] <- dff
  object$active_assay <- assay_name

  log_command(object, "RunDFF", list(
    method = method,
    window_size = window_size,
    percentile = percentile,
    source_assay = assay,
    assay_name = assay_name
  ))
}

#' Run Spike Inference
#'
#' @param object CaExperiment object
#' @param method Spike inference method
#' @param assay Source assay (default: active)
#' @param name Name for spike result
#' @param ... Additional arguments passed to inference function
#'
#' @return Modified CaExperiment object
#' @export
RunSpikes <- function(object,
                      method = c("oasis", "threshold", "bayesian", "statistical"),
                      assay = NULL,
                      name = NULL,
                      threshold_sd = 2.5,
                      ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  method <- match.arg(method)

  if (is.null(assay)) {
    assay <- object$active_assay
  }

  if (is.null(name)) {
    name <- method
  }

  traces <- GetTraces(object, assay = assay)

  # Run spike inference based on method
  result <- switch(method,
    threshold = {
      # Use threshold_spike_detection
      threshold_spike_detection(traces, threshold_sd = threshold_sd, ...)
    },
    statistical = {
      # Use statistical_spike_inference
      statistical_spike_inference(traces, ...)
    },
    oasis = {
      # Apply infer_spikes to each cell
      n_cells <- nrow(traces)
      spike_matrix <- matrix(0, nrow = n_cells, ncol = ncol(traces))
      rownames(spike_matrix) <- rownames(traces)
      for (i in seq_len(n_cells)) {
        res <- infer_spikes(traces[i, ], method = "oasis", verbose = FALSE, ...)
        if (is.data.frame(res) && "spike" %in% names(res)) {
          spike_matrix[i, ] <- res$spike
        }
      }
      list(spike_predictions = spike_matrix, method = "oasis")
    },
    bayesian = {
      # Use bayesian_spike_inference for each cell
      n_cells <- nrow(traces)
      spike_matrix <- matrix(0, nrow = n_cells, ncol = ncol(traces))
      rownames(spike_matrix) <- rownames(traces)
      for (i in seq_len(n_cells)) {
        res <- bayesian_spike_inference(traces[i, ], ...)
        if (is.list(res) && "spikes" %in% names(res)) {
          spike_matrix[i, ] <- res$spikes
        }
      }
      list(spike_predictions = spike_matrix, method = "bayesian")
    }
  )

  # Store spikes - handle different result formats
  if (is.list(result)) {
    if ("spike_predictions" %in% names(result)) {
      object$spikes[[name]] <- result$spike_predictions
    } else if ("spikes" %in% names(result)) {
      object$spikes[[name]] <- result$spikes
    } else {
      object$spikes[[name]] <- result
    }
    object$misc[[paste0("spikes_", name)]] <- result
  } else {
    object$spikes[[name]] <- result
  }

  # Add spike rate to metadata
  spike_matrix <- object$spikes[[name]]
  if (is.matrix(spike_matrix)) {
    spike_rates <- rowSums(spike_matrix > 0, na.rm = TRUE) / Duration(object)
    object <- AddMetaData(object, spike_rates, paste0("spike_rate_", name))
  }

  log_command(object, "RunSpikes", list(
    method = method,
    threshold_sd = threshold_sd,
    source_assay = assay,
    name = name
  ))
}

#' Run PCA Dimensionality Reduction
#'
#' @param object CaExperiment object
#' @param n_components Number of principal components
#' @param assay Source assay
#' @param name Name for reduction
#' @param scale Scale data before PCA
#'
#' @return Modified CaExperiment object
#' @export
RunPCA <- function(object,
                   n_components = 20,
                   assay = NULL,
                   name = "pca",
                   scale = TRUE) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.null(assay)) {
    assay <- object$active_assay
  }

  traces <- GetTraces(object, assay = assay)
  n_cells <- nrow(traces)
  n_frames <- ncol(traces)

  # Limit components to valid range
  max_components <- min(n_cells, n_frames) - 1
  n_components <- min(n_components, max_components)

  # Run PCA on cells (cells as observations, time points as features)
  # This gives us cell embeddings in PC space
  pca_result <- prcomp(traces, center = TRUE, scale. = scale, rank. = n_components)

  # pca_result$x contains cell embeddings (n_cells x n_components)
  # pca_result$rotation contains loadings (n_frames x n_components)
  cell_embeddings <- pca_result$x

  object$reductions[[name]] <- list(
    embeddings = cell_embeddings,
    rotation = pca_result$rotation,
    sdev = pca_result$sdev,
    center = pca_result$center,
    scale = if (scale) pca_result$scale else NULL,
    variance_explained = pca_result$sdev^2 / sum(pca_result$sdev^2),
    method = "pca"
  )

  log_command(object, "RunPCA", list(
    n_components = n_components,
    source_assay = assay,
    scale = scale,
    name = name
  ))
}

#' Run GPFA (Gaussian Process Factor Analysis)
#'
#' @param object CaExperiment object
#' @param n_latents Number of latent dimensions
#' @param trial_ids Vector of trial identifiers
#' @param assay Source assay
#' @param name Name for reduction
#' @param ... Additional arguments to fit_gpfa
#'
#' @return Modified CaExperiment object
#' @export
RunGPFA <- function(object,
                    n_latents = 3,
                    trial_ids = NULL,
                    assay = NULL,
                    name = "gpfa",
                    ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.null(assay)) {
    assay <- object$active_assay
  }

  traces <- GetTraces(object, assay = assay)

  # Run GPFA
  gpfa_result <- fit_gpfa(traces, n_latents = n_latents,
                          trial_ids = trial_ids, ...)

  object$reductions[[name]] <- gpfa_result

  log_command(object, "RunGPFA", list(
    n_latents = n_latents,
    source_assay = assay,
    name = name
  ))
}

#' Run dPCA (Demixed PCA)
#'
#' @param object CaExperiment object
#' @param trial_data 3D array (neurons x time x trials)
#' @param condition_labels Condition labels for each trial
#' @param n_components Number of components per marginalization
#' @param name Name for reduction
#' @param ... Additional arguments to fit_dpca
#'
#' @return Modified CaExperiment object
#' @export
RunDPCA <- function(object,
                    trial_data,
                    condition_labels,
                    n_components = 10,
                    name = "dpca",
                    ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  # Run dPCA
  dpca_result <- fit_dpca(trial_data, condition_labels,
                          n_components = n_components, ...)

  object$reductions[[name]] <- dpca_result

  log_command(object, "RunDPCA", list(
    n_components = n_components,
    n_conditions = length(unique(condition_labels)),
    name = name
  ))
}

#' Run Functional Connectivity Analysis
#'
#' @param object CaExperiment object
#' @param method Connectivity method
#' @param assay Source assay
#' @param name Name for graph
#' @param ... Additional arguments to functional_connectivity
#'
#' @return Modified CaExperiment object
#' @export
RunConnectivity <- function(object,
                            method = c("correlation", "partial", "mutual_information"),
                            assay = NULL,
                            name = NULL,
                            ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  method <- match.arg(method)

  if (is.null(assay)) {
    assay <- object$active_assay
  }

  if (is.null(name)) {
    name <- paste0("connectivity_", method)
  }

  traces <- GetTraces(object, assay = assay)

  # Run connectivity analysis
  conn_result <- functional_connectivity(t(traces), method = method, ...)

  # Extract connectivity matrix
  if (is.list(conn_result)) {
    conn_matrix <- conn_result[[1]]  # First element is usually the matrix
  } else {
    conn_matrix <- conn_result
  }

  rownames(conn_matrix) <- CellIDs(object)
  colnames(conn_matrix) <- CellIDs(object)

  object$graphs[[name]] <- conn_matrix

  log_command(object, "RunConnectivity", list(
    method = method,
    source_assay = assay,
    name = name
  ))
}

#' Run Graph Metrics
#'
#' @param object CaExperiment object
#' @param graph Name of connectivity graph to use
#' @param threshold Threshold for binarizing connectivity
#'
#' @return Modified CaExperiment object with metrics added to metadata
#' @export
RunGraphMetrics <- function(object, graph = NULL, threshold = 0.3) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (length(object$graphs) == 0) {
    stop("No graphs available. Run RunConnectivity() first.")
  }

  if (is.null(graph)) {
    graph <- names(object$graphs)[1]
  }

  conn_matrix <- GetGraph(object, graph)

  # Compute graph metrics
  metrics <- graph_metrics(conn_matrix, threshold = threshold)

  # Add to metadata
  object <- AddMetaData(object, metrics$degree, "graph_degree")
  object <- AddMetaData(object, metrics$strength, "graph_strength")
  object <- AddMetaData(object, metrics$clustering, "graph_clustering")
  if (!is.null(metrics$betweenness)) {
    object <- AddMetaData(object, metrics$betweenness, "graph_betweenness")
  }

  # Store global metrics
  object$misc$graph_metrics <- list(
    density = metrics$density,
    transitivity = metrics$transitivity,
    modularity = metrics$modularity
  )

  log_command(object, "RunGraphMetrics", list(
    graph = graph,
    threshold = threshold
  ))
}

#' Run Neural Assembly Detection
#'
#' @param object CaExperiment object
#' @param method Detection method
#' @param assay Source assay
#' @param name Name for assembly result
#' @param ... Additional arguments to detect_assemblies
#'
#' @return Modified CaExperiment object
#' @export
RunAssemblies <- function(object,
                          method = c("ica", "pca", "nmf", "correlation"),
                          assay = NULL,
                          name = "default",
                          ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  method <- match.arg(method)

  if (is.null(assay)) {
    assay <- object$active_assay
  }

  traces <- GetTraces(object, assay = assay)

  # Run assembly detection
  if (method == "correlation") {
    assembly_result <- detect_assemblies_correlation(traces, ...)
  } else {
    assembly_result <- detect_assemblies(traces, method = method, ...)
  }

  object$assemblies[[name]] <- assembly_result

  # Add assembly membership to metadata
  if (!is.null(assembly_result$membership)) {
    object <- AddMetaData(object, assembly_result$membership,
                          paste0("assembly_", name))
  }

  log_command(object, "RunAssemblies", list(
    method = method,
    source_assay = assay,
    name = name
  ))
}

#' Run Transient Detection
#'
#' @param object CaExperiment object
#' @param assay Source assay
#' @param name Name for transient result
#' @param ... Additional arguments to detect_transients
#'
#' @return Modified CaExperiment object
#' @export
RunTransients <- function(object,
                          assay = NULL,
                          name = "default",
                          ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.null(assay)) {
    assay <- object$active_assay
  }

  traces <- GetTraces(object, assay = assay)

  # Run transient detection
  transient_result <- detect_transients(traces, frame_rate = object$frame_rate, ...)

  object$transients[[name]] <- transient_result

  # Add transient counts to metadata
  if (!is.null(transient_result$events)) {
    event_counts <- table(transient_result$events$cell_id)
    counts <- rep(0, ncells(object))
    names(counts) <- CellIDs(object)
    counts[names(event_counts)] <- as.numeric(event_counts)
    object <- AddMetaData(object, counts, paste0("n_transients_", name))
  }

  log_command(object, "RunTransients", list(
    source_assay = assay,
    name = name
  ))
}

#' Run Tuning Curve Analysis
#'
#' @param object CaExperiment object
#' @param stimulus Stimulus values for each frame
#' @param type Type of tuning ("orientation", "contrast", "place")
#' @param assay Source assay
#' @param ... Additional arguments to tuning functions
#'
#' @return Modified CaExperiment object
#' @export
RunTuning <- function(object,
                      stimulus,
                      type = c("orientation", "contrast", "place", "generic"),
                      assay = NULL,
                      ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  type <- match.arg(type)

  if (is.null(assay)) {
    assay <- object$active_assay
  }

  traces <- GetTraces(object, assay = assay)

  # Run appropriate tuning analysis
  tuning_result <- switch(type,
    orientation = fit_orientation_tuning(traces, orientations = stimulus, ...),
    contrast = fit_contrast_response(traces, contrasts = stimulus, ...),
    place = fit_place_field(traces, positions = stimulus, ...),
    generic = fit_tuning_curve(traces, stimulus = stimulus, ...)
  )

  object$tuning[[type]] <- tuning_result

  # Add selectivity indices to metadata
  if (type == "orientation") {
    object <- AddMetaData(object, tuning_result$osi, "osi")
    object <- AddMetaData(object, tuning_result$dsi, "dsi")
    object <- AddMetaData(object, tuning_result$preferred, "pref_orientation")
  } else if (type == "place") {
    if (!is.null(tuning_result$spatial_info)) {
      object <- AddMetaData(object, tuning_result$spatial_info, "spatial_info")
    }
  }

  log_command(object, "RunTuning", list(
    type = type,
    source_assay = assay
  ))
}

#' Run State Space Trajectories
#'
#' @param object CaExperiment object
#' @param reduction Name of reduction to use for trajectories
#' @param dims Dimensions to use
#' @param smooth Smoothing parameter
#' @param name Name for trajectory result
#'
#' @return Modified CaExperiment object
#' @export
RunTrajectories <- function(object,
                            reduction = "pca",
                            dims = 1:3,
                            smooth = 5,
                            name = "default") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  # Get embeddings
  embeddings <- GetEmbeddings(object, reduction, dims = dims)

  # Compute trajectories
  traj_result <- compute_trajectories(embeddings, smooth = smooth)

  object$trajectories[[name]] <- traj_result

  log_command(object, "RunTrajectories", list(
    reduction = reduction,
    dims = paste(dims, collapse = ","),
    smooth = smooth,
    name = name
  ))
}

#' Run Sequence Detection
#'
#' @param object CaExperiment object
#' @param assay Source assay
#' @param name Name for sequence result
#' @param ... Additional arguments to detect_sequences
#'
#' @return Modified CaExperiment object
#' @export
RunSequences <- function(object,
                         assay = NULL,
                         name = "default",
                         ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.null(assay)) {
    assay <- object$active_assay
  }

  traces <- GetTraces(object, assay = assay)

  # Run sequence detection
  seq_result <- detect_sequences(traces, ...)

  object$sequences[[name]] <- seq_result

  log_command(object, "RunSequences", list(
    source_assay = assay,
    name = name
  ))
}

# ============================================================================
# Advanced Analysis Run* Methods
# ============================================================================

#' Run Hidden Markov Model
#'
#' @param object CaExperiment object
#' @param n_states Number of hidden states
#' @param assay Source assay
#' @param name Name for result
#' @param ... Additional arguments to fit_hmm
#'
#' @return Modified CaExperiment object
#' @export
RunHMM <- function(object, n_states = 3, assay = NULL, name = "hmm", ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.null(assay)) assay <- object$active_assay
  traces <- GetTraces(object, assay = assay)

  hmm_result <- fit_hmm(traces, n_states = n_states, ...)

  object$reductions[[name]] <- hmm_result
  object <- AddMetaData(object, rep(NA, ncells(object)), paste0("state_", name))

  # Add state sequence to misc
  object$misc[[paste0("hmm_", name)]] <- hmm_result

  log_command(object, "RunHMM", list(n_states = n_states, source_assay = assay, name = name))
}

#' Run Switching Linear Dynamical System
#'
#' @param object CaExperiment object
#' @param n_states Number of discrete states
#' @param latent_dim Latent dimension
#' @param assay Source assay
#' @param name Name for result
#' @param ... Additional arguments to fit_slds
#'
#' @return Modified CaExperiment object
#' @export
RunSLDS <- function(object, n_states = 2, latent_dim = 3, assay = NULL, name = "slds", ...) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.null(assay)) assay <- object$active_assay
  traces <- GetTraces(object, assay = assay)

  slds_result <- fit_slds(traces, n_states = n_states, latent_dim = latent_dim, ...)

  object$reductions[[name]] <- list(
    embeddings = t(slds_result$latent_states),
    discrete_states = slds_result$discrete_states,
    method = "slds"
  )
  object$misc[[paste0("slds_", name)]] <- slds_result

  log_command(object, "RunSLDS", list(n_states = n_states, latent_dim = latent_dim, name = name))
}

#' Run Tensor Decomposition
#'
#' @param object CaExperiment object
#' @param trial_starts Trial start indices
#' @param trial_length Length of each trial
#' @param rank Decomposition rank
#' @param method "cp" or "tucker"
#' @param name Name for result
#'
#' @return Modified CaExperiment object
#' @export
RunTensorDecomp <- function(object, trial_starts, trial_length, rank = 5,
                             method = c("cp", "tucker"), name = "tensor") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  method <- match.arg(method)
  traces <- GetTraces(object)

  # Create trial tensor
  X <- create_trial_tensor(traces, trial_starts, trial_length)

  # Decompose
  result <- if (method == "cp") {
    cp_decompose(X, rank = rank)
  } else {
    tucker_decompose(X, ranks = rank)
  }

  object$reductions[[name]] <- result
  object <- AddMetaData(object, result$factors$neuron[, 1], paste0("tensor_loading_", name))

  log_command(object, "RunTensorDecomp", list(method = method, rank = rank, name = name))
}

#' Run Changepoint Detection
#'
#' @param object CaExperiment object
#' @param method Detection method
#' @param type Type of change
#' @param assay Source assay
#' @param name Name for result
#'
#' @return Modified CaExperiment object
#' @export
RunChangepoints <- function(object, method = "pelt", type = "mean",
                             assay = NULL, name = "changepoints") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.null(assay)) assay <- object$active_assay
  traces <- GetTraces(object, assay = assay)

  cp_result <- detect_population_changepoints(traces, method = method, aggregate = "pc1")

  object$misc[[name]] <- cp_result

  log_command(object, "RunChangepoints", list(method = method, type = type, name = name))
}

#' Run Topological Data Analysis
#'
#' @param object CaExperiment object
#' @param reduction Name of reduction to use
#' @param max_dim Maximum homology dimension
#' @param name Name for result
#'
#' @return Modified CaExperiment object
#' @export
RunTDA <- function(object, reduction = "pca", max_dim = 1, name = "tda") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  embeddings <- GetEmbeddings(object, reduction)
  ph_result <- compute_persistent_homology(embeddings, max_dim = max_dim)

  object$misc[[name]] <- ph_result

  log_command(object, "RunTDA", list(reduction = reduction, max_dim = max_dim, name = name))
}

#' Run Recurrence Quantification Analysis
#'
#' @param object CaExperiment object
#' @param reduction Name of reduction to use
#' @param embed_dim Embedding dimension
#' @param name Name for result
#'
#' @return Modified CaExperiment object
#' @export
RunRQA <- function(object, reduction = "pca", embed_dim = 3, name = "rqa") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  embeddings <- GetEmbeddings(object, reduction)
  rp <- recurrence_plot(embeddings, embed_dim = embed_dim)
  rqa <- rqa_metrics(rp)

  object$misc[[name]] <- list(recurrence_plot = rp, metrics = rqa)

  log_command(object, "RunRQA", list(reduction = reduction, embed_dim = embed_dim, name = name))
}

#' Run Causal Discovery
#'
#' @param object CaExperiment object
#' @param method "pc", "fci", or "ges"
#' @param alpha Significance level
#' @param assay Source assay
#' @param name Name for result
#'
#' @return Modified CaExperiment object
#' @export
RunCausalDiscovery <- function(object, method = "pc", alpha = 0.01,
                                assay = NULL, name = "causal") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.null(assay)) assay <- object$active_assay
  traces <- GetTraces(object, assay = assay)

  causal_result <- discover_neural_causality(traces, method = method, alpha = alpha)

  object$graphs[[name]] <- causal_result$adjacency %||% causal_result$skeleton
  object$misc[[paste0("causal_", name)]] <- causal_result

  log_command(object, "RunCausalDiscovery", list(method = method, alpha = alpha, name = name))
}

#' Run Neural ODE
#'
#' @param object CaExperiment object
#' @param latent_dim Latent dimension
#' @param n_epochs Training epochs
#' @param reduction Use existing reduction as input
#' @param name Name for result
#'
#' @return Modified CaExperiment object
#' @export
RunNeuralODE <- function(object, latent_dim = 3, n_epochs = 50,
                          reduction = NULL, name = "node") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  traces <- GetTraces(object)
  time <- object$time

  node_result <- fit_latent_node(traces, time, latent_dim = latent_dim, n_epochs = n_epochs)

  object$reductions[[name]] <- list(
    embeddings = node_result$pca$x[, 1:latent_dim],
    method = "neural_ode",
    model = node_result
  )
  object$misc[[paste0("node_", name)]] <- node_result

  log_command(object, "RunNeuralODE", list(latent_dim = latent_dim, n_epochs = n_epochs, name = name))
}

#' Run Gaussian Process Smoothing
#'
#' @param object CaExperiment object
#' @param kernel GP kernel type
#' @param assay Source assay
#' @param assay_name Name for smoothed assay
#'
#' @return Modified CaExperiment object
#' @export
RunGPSmooth <- function(object, kernel = "rbf", assay = NULL, assay_name = "gp_smooth") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  if (is.null(assay)) assay <- object$active_assay
  traces <- GetTraces(object, assay = assay)
  time <- object$time

  gp_fit <- fit_gp_traces(time, traces, kernel = kernel)
  pred <- predict_gp_traces(gp_fit, time, return_var = FALSE)

  object$assays[[assay_name]] <- pred$mean
  object$misc$gp_fit <- gp_fit

  log_command(object, "RunGPSmooth", list(kernel = kernel, assay_name = assay_name))
}

#' Run Optimal Transport Comparison
#'
#' @param object CaExperiment object
#' @param condition_var Metadata variable for conditions
#' @param method "wasserstein" or "gromov_wasserstein"
#' @param name Name for result
#'
#' @return Modified CaExperiment object
#' @export
RunOTCompare <- function(object, condition_var, method = "wasserstein", name = "ot") {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  meta <- GetMetaData(object)
  if (!condition_var %in% names(meta)) {
    stop(sprintf("Condition variable '%s' not found in metadata", condition_var))
  }

  conditions <- unique(meta[[condition_var]])
  traces <- GetTraces(object)

  # Compute pairwise distances
  distances <- matrix(0, length(conditions), length(conditions))
  rownames(distances) <- colnames(distances) <- conditions

  for (i in seq_along(conditions)) {
    for (j in seq_along(conditions)) {
      if (i < j) {
        idx_i <- which(meta[[condition_var]] == conditions[i])
        idx_j <- which(meta[[condition_var]] == conditions[j])

        d <- compare_representations_ot(
          traces[idx_i, , drop = FALSE],
          traces[idx_j, , drop = FALSE],
          method = method,
          by = "population"
        )
        distances[i, j] <- d
        distances[j, i] <- d
      }
    }
  }

  object$misc[[name]] <- list(distances = distances, conditions = conditions, method = method)

  log_command(object, "RunOTCompare", list(condition_var = condition_var, method = method, name = name))
}

# ============================================================================
# Utility Functions
# ============================================================================

#' Merge CaExperiment Objects
#'
#' @param x First CaExperiment object
#' @param y Second CaExperiment object
#' @param add_batch_id Add batch identifier to metadata
#'
#' @return Merged CaExperiment object
#' @export
MergeCaExperiments <- function(x, y, add_batch_id = TRUE) {
  if (!inherits(x, "CaExperiment") || !inherits(y, "CaExperiment")) {
    stop("Both objects must be CaExperiment")
  }

  if (x$frame_rate != y$frame_rate) {
    warning("Frame rates differ; using frame rate from first object")
  }

  # Merge assays (cells dimension)
  shared_assays <- intersect(names(x$assays), names(y$assays))
  if (length(shared_assays) == 0) {
    stop("No shared assays to merge")
  }

  new_assays <- lapply(shared_assays, function(a) {
    # Get minimum frame count
    n_frames <- min(ncol(x$assays[[a]]), ncol(y$assays[[a]]))
    rbind(
      x$assays[[a]][, 1:n_frames],
      y$assays[[a]][, 1:n_frames]
    )
  })
  names(new_assays) <- shared_assays

  # Merge metadata
  new_meta <- rbind(
    x$meta.data,
    y$meta.data
  )

  if (add_batch_id) {
    new_meta$batch <- c(
      rep(x$experiment_id, nrow(x$meta.data)),
      rep(y$experiment_id, nrow(y$meta.data))
    )
  }

  # Create merged object
  merged <- CaExperiment(
    traces = new_assays[[1]],
    meta.data = new_meta,
    frame_rate = x$frame_rate,
    experiment_id = paste(x$experiment_id, y$experiment_id, sep = "_")
  )

  # Add other assays
  for (a in shared_assays[-1]) {
    merged$assays[[a]] <- new_assays[[a]]
  }

  # Log merge
  merged <- log_command(merged, "MergeCaExperiments", list(
    experiment_1 = x$experiment_id,
    experiment_2 = y$experiment_id
  ))

  merged
}

#' Export Command History as Script
#'
#' @param object CaExperiment object
#' @param file Output file path
#'
#' @export
ExportCommands <- function(object, file = NULL) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  lines <- c(
    "# CaExperiment Analysis Pipeline",
    sprintf("# Generated: %s", Sys.time()),
    sprintf("# Experiment ID: %s", object$experiment_id),
    "",
    "library(CaImagingAnalysisFr)",
    ""
  )

  for (cmd in object$commands) {
    args_str <- paste(names(cmd$args), "=", sapply(cmd$args, deparse), collapse = ", ")
    lines <- c(lines, sprintf("# %s", cmd$timestamp))
    lines <- c(lines, sprintf("ca <- %s(ca, %s)", cmd$command, args_str))
    lines <- c(lines, "")
  }

  script <- paste(lines, collapse = "\n")

  if (!is.null(file)) {
    writeLines(script, file)
    message(sprintf("Commands exported to: %s", file))
  }

  invisible(script)
}

#' Check CaExperiment Object
#'
#' Validate internal consistency of a CaExperiment object.
#'
#' @param object CaExperiment object
#'
#' @return TRUE if valid, otherwise throws error
#' @export
ValidateCaExperiment <- function(object) {
  if (!inherits(object, "CaExperiment")) {
    stop("object must be a CaExperiment")
  }

  # Check assays consistency
  n_cells <- nrow(object$assays[[1]])
  n_frames <- ncol(object$assays[[1]])

  for (name in names(object$assays)) {
    if (nrow(object$assays[[name]]) != n_cells) {
      stop(sprintf("Assay '%s' has inconsistent number of cells", name))
    }
    if (ncol(object$assays[[name]]) != n_frames) {
      stop(sprintf("Assay '%s' has inconsistent number of frames", name))
    }
  }

  # Check spikes consistency
  for (name in names(object$spikes)) {
    spk <- object$spikes[[name]]
    if (is.matrix(spk)) {
      if (nrow(spk) != n_cells) {
        stop(sprintf("Spike '%s' has inconsistent number of cells", name))
      }
    }
  }

  # Check metadata
  if (nrow(object$meta.data) != n_cells) {
    stop("Metadata has inconsistent number of rows")
  }

  # Check graphs
  for (name in names(object$graphs)) {
    g <- object$graphs[[name]]
    if (is.matrix(g)) {
      if (nrow(g) != n_cells || ncol(g) != n_cells) {
        stop(sprintf("Graph '%s' has inconsistent dimensions", name))
      }
    }
  }

  message("CaExperiment object is valid")
  TRUE
}

#' Convert Legacy Data to CaExperiment
#'
#' @param traces Trace matrix
#' @param spikes Optional spike matrix
#' @param metadata Optional metadata data frame
#' @param ... Additional arguments to CaExperiment
#'
#' @return CaExperiment object
#' @export
as_ca_experiment <- function(traces, spikes = NULL, metadata = NULL, ...) {
  ca <- CaExperiment(traces, meta.data = metadata, ...)

  if (!is.null(spikes)) {
    ca$spikes[["imported"]] <- spikes
  }

  ca
}

#' @export
is_ca_experiment <- function(x) {
  inherits(x, "CaExperiment")
}
