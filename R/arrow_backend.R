#' Arrow/DuckDB Backend for Large Datasets
#'
#' Out-of-core processing for calcium imaging datasets that don't fit
#' in memory. Uses Apache Arrow for efficient columnar storage and
#' DuckDB for SQL-like queries on large trace data.
#'
#' @name arrow_backend
#' @keywords internal
NULL

#' Check Arrow availability
#'
#' @return Logical
#' @keywords internal
.arrow_available <- function() {
  requireNamespace("arrow", quietly = TRUE)
}

#' Check DuckDB availability
#'
#' @return Logical
#' @keywords internal
.duckdb_available <- function() {
  requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("DBI", quietly = TRUE)
}

#' Save traces to Arrow/Parquet format
#'
#' Saves trace data in efficient columnar format for large datasets.
#'
#' @param traces Matrix of traces (cells x time) or data.table
#' @param path Output file path (.parquet extension)
#' @param time_vector Optional time vector
#' @param cell_ids Optional cell identifiers
#' @param metadata Optional named list of metadata
#' @param compression Compression codec: "snappy", "gzip", "zstd", "lz4"
#'
#' @return Path to saved file
#' @export
#'
#' @examples
#' \dontrun{
#' traces <- matrix(rnorm(1000 * 10000), nrow = 1000)
#' save_traces_parquet(traces, "large_traces.parquet")
#' }
save_traces_parquet <- function(traces, path, time_vector = NULL,
                                 cell_ids = NULL, metadata = NULL,
                                 compression = "snappy") {
  if (!.arrow_available()) {
    stop("arrow package required. Install with: install.packages('arrow')")
  }

  # Convert to long format for efficient storage
  if (is.matrix(traces)) {
    n_cells <- nrow(traces)
    n_frames <- ncol(traces)

    if (is.null(time_vector)) time_vector <- seq_len(n_frames)
    if (is.null(cell_ids)) cell_ids <- paste0("cell_", seq_len(n_cells))

    df <- data.frame(
      cell_id = rep(cell_ids, each = n_frames),
      frame = rep(seq_len(n_frames), times = n_cells),
      time = rep(time_vector, times = n_cells),
      value = as.vector(t(traces))
    )
  } else {
    df <- as.data.frame(traces)
  }

  # Convert to Arrow table
  tbl <- arrow::arrow_table(df)

  # Add metadata
  if (!is.null(metadata)) {
    tbl <- arrow::arrow_table(
      df,
      metadata = list(custom = jsonlite::toJSON(metadata))
    )
  }

  # Write to Parquet
  arrow::write_parquet(
    tbl,
    path,
    compression = compression
  )

  message(sprintf("Saved traces to %s (%.2f MB)",
                  path, file.size(path) / 1024^2))
  invisible(path)
}

#' Load traces from Parquet file
#'
#' Efficiently load trace data with optional filtering.
#'
#' @param path Path to Parquet file
#' @param cells Character vector of cell IDs to load (NULL for all)
#' @param frames Integer vector of frames to load (NULL for all)
#' @param as_matrix Return as matrix (TRUE) or data.frame (FALSE)
#'
#' @return Matrix or data.frame of traces
#' @export
load_traces_parquet <- function(path, cells = NULL, frames = NULL,
                                 as_matrix = TRUE) {
  if (!.arrow_available()) {
    stop("arrow package required")
  }

  # Open dataset (lazy)
  ds <- arrow::open_dataset(path)

  # Apply filters
  if (!is.null(cells) || !is.null(frames)) {
    if (!is.null(cells)) {
      ds <- dplyr::filter(ds, cell_id %in% cells)
    }
    if (!is.null(frames)) {
      ds <- dplyr::filter(ds, frame %in% frames)
    }
  }

  # Collect data
  df <- dplyr::collect(ds)

  if (as_matrix) {
    # Convert back to wide matrix
    wide <- tidyr::pivot_wider(
      df,
      names_from = frame,
      values_from = value,
      id_cols = cell_id
    )
    mat <- as.matrix(wide[, -1])
    rownames(mat) <- wide$cell_id
    return(mat)
  }

  as.data.frame(df)
}

#' Create DuckDB connection for trace data
#'
#' Opens a DuckDB database for SQL queries on trace data.
#'
#' @param path Path to Parquet file or directory
#' @param memory In-memory database (faster but uses more RAM)
#'
#' @return DuckDB connection object
#' @export
#'
#' @examples
#' \dontrun{
#' conn <- duckdb_traces("traces.parquet")
#' result <- DBI::dbGetQuery(conn, "
#'   SELECT cell_id, AVG(value) as mean_value
#'   FROM traces
#'   GROUP BY cell_id
#' ")
#' close_duckdb_traces(conn)
#' }
duckdb_traces <- function(path, memory = TRUE) {
  if (!.duckdb_available()) {
    stop("duckdb and DBI packages required")
  }

  # Create connection
  if (memory) {
    conn <- DBI::dbConnect(duckdb::duckdb(), ":memory:")
  } else {
    conn <- DBI::dbConnect(duckdb::duckdb(), paste0(path, ".duckdb"))
  }

  # Register Parquet file as table
  query <- sprintf("CREATE VIEW traces AS SELECT * FROM '%s'", path)
  DBI::dbExecute(conn, query)

  message("DuckDB connection opened. Query with: DBI::dbGetQuery(conn, 'SELECT ...')")
  conn
}

#' Close DuckDB connection
#'
#' @param conn DuckDB connection
#'
#' @export
close_duckdb_traces <- function(conn) {
  if (!.duckdb_available()) return(invisible(NULL))
  DBI::dbDisconnect(conn, shutdown = TRUE)
  invisible(NULL)
}

#' Query traces with SQL
#'
#' Execute SQL query on trace data stored in Parquet.
#'
#' @param path Path to Parquet file
#' @param query SQL query string (use 'traces' as table name)
#'
#' @return data.frame with query results
#' @export
#'
#' @examples
#' \dontrun{
#' # Get mean fluorescence per cell
#' result <- query_traces("traces.parquet", "
#'   SELECT cell_id, AVG(value) as mean_f, MAX(value) as max_f
#'   FROM traces
#'   GROUP BY cell_id
#' ")
#'
#' # Get time points with high activity
#' result <- query_traces("traces.parquet", "
#'   SELECT frame, AVG(value) as population_mean
#'   FROM traces
#'   GROUP BY frame
#'   HAVING AVG(value) > 2
#'   ORDER BY frame
#' ")
#' }
query_traces <- function(path, query) {
  conn <- duckdb_traces(path)
  on.exit(close_duckdb_traces(conn))

  result <- DBI::dbGetQuery(conn, query)
  as.data.frame(result)
}

#' Stream-process large trace files
#'
#' Process trace data in batches for memory efficiency.
#'
#' @param path Path to Parquet file
#' @param fun Function to apply to each batch
#' @param batch_size Number of cells per batch
#' @param ... Additional arguments to fun
#'
#' @return List of results from each batch
#' @export
#'
#' @examples
#' \dontrun{
#' # Compute per-cell statistics in batches
#' results <- stream_process_traces(
#'   "large_traces.parquet",
#'   function(batch) {
#'     data.frame(
#'       cell_id = unique(batch$cell_id),
#'       mean = tapply(batch$value, batch$cell_id, mean)
#'     )
#'   },
#'   batch_size = 100
#' )
#' stats <- do.call(rbind, results)
#' }
stream_process_traces <- function(path, fun, batch_size = 100, ...) {
  if (!.arrow_available()) {
    stop("arrow package required")
  }

  # Get cell IDs
  ds <- arrow::open_dataset(path)
  all_cells <- dplyr::collect(
    dplyr::distinct(dplyr::select(ds, cell_id))
  )$cell_id

  # Process in batches
  n_cells <- length(all_cells)
  n_batches <- ceiling(n_cells / batch_size)

  results <- vector("list", n_batches)

  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_cells)
    batch_cells <- all_cells[start_idx:end_idx]

    # Load batch
    batch_data <- load_traces_parquet(path, cells = batch_cells, as_matrix = FALSE)

    # Apply function
    results[[i]] <- fun(batch_data, ...)

    if (i %% 10 == 0) {
      message(sprintf("Processed %d/%d batches", i, n_batches))
    }
  }

  results
}

#' Create partitioned dataset for multi-session data
#'
#' Organize trace data from multiple sessions into partitioned storage.
#'
#' @param session_list Named list of trace matrices
#' @param path Output directory
#' @param partition_by Partitioning column(s)
#'
#' @return Path to dataset directory
#' @export
create_partitioned_dataset <- function(session_list, path,
                                        partition_by = "session_id") {
  if (!.arrow_available()) {
    stop("arrow package required")
  }

  # Create output directory
  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  # Process each session
  all_data <- lapply(names(session_list), function(session_id) {
    traces <- session_list[[session_id]]

    if (is.matrix(traces)) {
      n_cells <- nrow(traces)
      n_frames <- ncol(traces)

      data.frame(
        session_id = session_id,
        cell_id = rep(paste0(session_id, "_cell_", seq_len(n_cells)), each = n_frames),
        frame = rep(seq_len(n_frames), times = n_cells),
        value = as.vector(t(traces))
      )
    } else {
      traces$session_id <- session_id
      traces
    }
  })

  combined <- do.call(rbind, all_data)

  # Write partitioned dataset
  arrow::write_dataset(
    combined,
    path,
    format = "parquet",
    partitioning = partition_by
  )

  message(sprintf("Created partitioned dataset at %s", path))
  invisible(path)
}

#' Open partitioned dataset
#'
#' @param path Path to partitioned dataset directory
#'
#' @return Arrow dataset object
#' @export
open_partitioned_dataset <- function(path) {
  if (!.arrow_available()) {
    stop("arrow package required")
  }

  arrow::open_dataset(path)
}

#' Compute summary statistics on Arrow dataset
#'
#' Efficient grouped statistics on out-of-core data.
#'
#' @param dataset Arrow dataset or path to Parquet
#' @param group_by Grouping column(s)
#'
#' @return data.frame with summary statistics
#' @export
arrow_cell_summary <- function(dataset, group_by = "cell_id") {
  if (!.arrow_available()) {
    stop("arrow package required")
  }

  if (is.character(dataset)) {
    dataset <- arrow::open_dataset(dataset)
  }

  result <- dataset |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) |>
    dplyr::summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      min_value = min(value, na.rm = TRUE),
      max_value = max(value, na.rm = TRUE),
      n_frames = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::collect()

  as.data.frame(result)
}

#' Convert between formats
#'
#' Convert trace data between matrix, data.frame, and Parquet.
#'
#' @param data Input data (matrix, data.frame, or path)
#' @param to Target format: "matrix", "data.frame", "parquet"
#' @param path Output path (required if to = "parquet")
#'
#' @return Converted data or path
#' @export
convert_trace_format <- function(data, to, path = NULL) {
  # Determine input format
  if (is.character(data) && file.exists(data)) {
    input_format <- "parquet"
    df <- load_traces_parquet(data, as_matrix = FALSE)
  } else if (is.matrix(data)) {
    input_format <- "matrix"
    n_cells <- nrow(data)
    n_frames <- ncol(data)
    df <- data.frame(
      cell_id = rep(seq_len(n_cells), each = n_frames),
      frame = rep(seq_len(n_frames), times = n_cells),
      value = as.vector(t(data))
    )
  } else {
    input_format <- "data.frame"
    df <- as.data.frame(data)
  }

  # Convert to target format
  if (to == "matrix") {
    wide <- tidyr::pivot_wider(
      df,
      names_from = frame,
      values_from = value,
      id_cols = cell_id
    )
    mat <- as.matrix(wide[, -1])
    rownames(mat) <- wide$cell_id
    return(mat)
  } else if (to == "data.frame") {
    return(df)
  } else if (to == "parquet") {
    if (is.null(path)) {
      stop("path required for Parquet output")
    }
    save_traces_parquet(df, path)
    return(path)
  }

  stop("Unknown target format: ", to)
}
