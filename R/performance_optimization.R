#' Process Large Datasets in Chunks
#' 
#' Processes large calcium imaging datasets in chunks to manage memory usage.
#' 
#' @param raw_df Input data frame
#' @param chunk_size Size of each chunk (default: from config)
#' @param process_function Function to apply to each chunk
#' @param combine_function Function to combine chunk results
#' @param verbose Whether to show progress messages
#' @return Combined results from all chunks
#' @export
process_in_chunks <- function(raw_df,
                             chunk_size = NULL,
                             process_function,
                             combine_function = NULL,
                             verbose = TRUE) {
  
  config <- get_config()
  
  if (is.null(chunk_size)) {
    chunk_size <- config$chunk_size
  }
  
  n_rows <- nrow(raw_df)
  n_chunks <- ceiling(n_rows / chunk_size)
  
  if (verbose) {
    message("Processing ", n_rows, " rows in ", n_chunks, " chunks of size ", chunk_size)
  }
  
  # Check memory usage
  estimated_size_mb <- (n_rows * ncol(raw_df) * 8) / (1024 * 1024)  # 8 bytes per numeric
  check_memory_usage(estimated_size_mb)
  
  chunk_results <- list()
  
  for (i in seq_len(n_chunks)) {
    if (verbose) {
      message("  Processing chunk ", i, "/", n_chunks)
    }
    
    # Define chunk indices
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_rows)
    
    # Extract chunk
    chunk_data <- raw_df[start_idx:end_idx, ]
    
    # Process chunk
    tryCatch({
      chunk_result <- process_function(chunk_data)
      chunk_results[[i]] <- chunk_result
    }, error = function(e) {
      warning("Error processing chunk ", i, ": ", e$message)
      chunk_results[[i]] <<- NULL
    })
  }
  
  # Combine results
  if (!is.null(combine_function)) {
    final_result <- combine_function(chunk_results)
  } else {
    final_result <- chunk_results
  }
  
  if (verbose) {
    message("Chunk processing completed")
  }
  
  return(final_result)
}

#' Memory-efficient calcium correction
#' 
#' @param raw_df Input data frame
#' @param method Correction method
#' @param span Smoothing span
#' @param normalize Whether to normalize
#' @param chunk_size Chunk size for processing
#' @param verbose Whether to show progress messages
#' @return Corrected data frame
#' @export
calcium_correction_chunked <- function(raw_df,
                                      method = "modern",
                                      span = NULL,
                                      normalize = TRUE,
                                      chunk_size = NULL,
                                      verbose = TRUE) {
  
  if (is.null(span)) {
    span <- get_config()$default_span
  }
  
  # Define chunk processing function
  process_chunk <- function(chunk_data) {
    calcium_correction(chunk_data, method = method, span = span, 
                      normalize = normalize, verbose = FALSE)
  }
  
  # Define combine function
  combine_chunks <- function(chunk_results) {
    # Remove NULL results
    valid_results <- chunk_results[!sapply(chunk_results, is.null)]
    
    if (length(valid_results) == 0) {
      stop("No valid chunk results to combine")
    }
    
    # Combine data frames
    do.call(rbind, valid_results)
  }
  
  # Process in chunks
  result <- process_in_chunks(
    raw_df = raw_df,
    chunk_size = chunk_size,
    process_function = process_chunk,
    combine_function = combine_chunks,
    verbose = verbose
  )
  
  return(result)
}

#' Memory-efficient spike inference
#' 
#' @param corrected_df Corrected data frame
#' @param method Spike inference method
#' @param chunk_size Chunk size for processing
#' @param verbose Whether to show progress messages
#' @return List of spike results
#' @export
infer_spikes_chunked <- function(corrected_df,
                                method = "oasis",
                                chunk_size = NULL,
                                verbose = TRUE) {
  
  config <- get_config()
  cell_cols <- names(corrected_df)[grepl(config$cell_pattern, names(corrected_df))]
  
  if (verbose) {
    message("Performing chunked spike inference for ", length(cell_cols), " cells")
  }
  
  results <- list()
  
  for (cell in cell_cols) {
    if (verbose) {
      message("  Processing cell: ", cell)
    }
    
    trace <- corrected_df[[cell]]
    
    # Process trace in chunks if it's very long
    if (length(trace) > config$chunk_size) {
      chunk_results <- process_trace_in_chunks(trace, method, chunk_size, verbose = FALSE)
      results[[cell]] <- combine_trace_chunks(chunk_results)
    } else {
      results[[cell]] <- infer_spikes(trace, method = method, fallback = TRUE, verbose = FALSE)
    }
  }
  
  return(results)
}

#' Process a single trace in chunks
#' 
#' @param trace Numeric vector
#' @param method Spike inference method
#' @param chunk_size Chunk size
#' @param verbose Whether to show progress messages
#' @return List of chunk results
#' @keywords internal
process_trace_in_chunks <- function(trace, method, chunk_size, verbose = TRUE) {
  
  n_samples <- length(trace)
  n_chunks <- ceiling(n_samples / chunk_size)
  
  chunk_results <- list()
  
  for (i in seq_len(n_chunks)) {
    if (verbose && i %% 10 == 0) {
      message("    Trace chunk ", i, "/", n_chunks)
    }
    
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_samples)
    
    chunk_trace <- trace[start_idx:end_idx]
    
    tryCatch({
      chunk_result <- infer_spikes(chunk_trace, method = method, fallback = TRUE, verbose = FALSE)
      chunk_results[[i]] <- list(
        start_idx = start_idx,
        end_idx = end_idx,
        result = chunk_result
      )
    }, error = function(e) {
      warning("Error processing trace chunk ", i, ": ", e$message)
      chunk_results[[i]] <<- NULL
    })
  }
  
  return(chunk_results)
}

#' Combine trace chunk results
#' 
#' @param chunk_results List of chunk results
#' @return Combined spike results
#' @keywords internal
combine_trace_chunks <- function(chunk_results) {
  
  # Remove NULL results
  valid_results <- chunk_results[!sapply(chunk_results, is.null)]
  
  if (length(valid_results) == 0) {
    stop("No valid trace chunk results to combine")
  }
  
  # Combine all chunks
  all_fit <- numeric(0)
  all_spike <- numeric(0)
  
  for (chunk in valid_results) {
    all_fit <- c(all_fit, chunk$result$fit)
    all_spike <- c(all_spike, chunk$result$spike)
  }
  
  data.frame(
    fit = all_fit,
    spike = all_spike,
    stringsAsFactors = FALSE
  )
}

#' Monitor memory usage
#' 
#' @param operation_name Name of the operation for logging
#' @param verbose Whether to show memory usage
#' @return Memory usage information
#' @export
monitor_memory_usage <- function(operation_name = "operation", verbose = TRUE) {
  
  mem_used <- 0  # Default value
  
  if (verbose) {
    # Get memory usage (platform dependent)
    tryCatch({
      if (.Platform$OS.type == "windows") {
        mem_info <- gc()
        mem_used <- sum(mem_info[, "used"]) / 1024  # MB
      } else {
        # Unix-like systems
        mem_info <- system("ps -o rss= -p $PPID", intern = TRUE)
        mem_used <- as.numeric(mem_info) / 1024  # MB
      }
      
      message("Memory usage for ", operation_name, ": ", round(mem_used, 1), " MB")
    }, error = function(e) {
      message("Could not determine memory usage for ", operation_name)
      mem_used <<- 0
    })
  }
  
  invisible(list(operation = operation_name, memory_mb = mem_used))
}

#' Optimize data structure for memory efficiency
#' 
#' @param raw_df Input data frame
#' @param optimize_types Whether to optimize data types
#' @param remove_na Whether to remove NA values
#' @param verbose Whether to show progress messages
#' @return Optimized data frame
#' @export
optimize_data_structure <- function(raw_df,
                                   optimize_types = TRUE,
                                   remove_na = FALSE,
                                   verbose = TRUE) {
  
  if (verbose) {
    message("Optimizing data structure for memory efficiency...")
  }
  
  optimized_df <- raw_df
  
  if (optimize_types) {
    # Optimize numeric columns
    numeric_cols <- sapply(optimized_df, is.numeric)
    
    for (col in names(optimized_df)[numeric_cols]) {
      values <- optimized_df[[col]]
      
      # Check if we can use integer instead of numeric
      if (all(values == round(values), na.rm = TRUE)) {
        optimized_df[[col]] <- as.integer(values)
      }
      
      # Check if we can use single precision
      if (max(abs(values), na.rm = TRUE) < 1e6) {
        optimized_df[[col]] <- as.numeric(values)  # Keep as double for precision
      }
    }
  }
  
  if (remove_na) {
    # Remove rows with too many NA values
    config <- get_config()
    max_na_fraction <- config$quality_control$max_missing_fraction
    
    na_counts <- rowSums(is.na(optimized_df))
    na_fraction <- na_counts / ncol(optimized_df)
    
    rows_to_keep <- na_fraction <= max_na_fraction
    optimized_df <- optimized_df[rows_to_keep, ]
    
    if (verbose) {
      removed_rows <- sum(!rows_to_keep)
      message("Removed ", removed_rows, " rows with too many NA values")
    }
  }
  
  if (verbose) {
    original_size <- object.size(raw_df)
    optimized_size <- object.size(optimized_df)
    reduction <- (original_size - optimized_size) / original_size * 100
    
    message("Memory reduction: ", round(reduction, 1), "%")
  }
  
  return(optimized_df)
}

#' Stream processing for very large datasets
#' 
#' @param file_path Path to data file
#' @param chunk_size Number of rows to read at once
#' @param process_function Function to apply to each chunk
#' @param output_file Path to output file (optional)
#' @param verbose Whether to show progress messages
#' @return Processing results
#' @export
stream_process_data <- function(file_path,
                               chunk_size = 1000,
                               process_function,
                               output_file = NULL,
                               verbose = TRUE) {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  if (verbose) {
    message("Starting stream processing of: ", file_path)
  }
  
  # Count total rows
  total_rows <- length(count.fields(file_path, sep = ",")) - 1  # Subtract header
  
  # Read header
  header <- read.csv(file_path, nrows = 1, header = TRUE)
  
  # Initialize output
  if (!is.null(output_file)) {
    write.csv(header, output_file, row.names = FALSE)
  }
  
  results <- list()
  chunk_count <- 0
  
  # Process in chunks
  for (start_row in seq(1, total_rows, by = chunk_size)) {
    chunk_count <- chunk_count + 1
    end_row <- min(start_row + chunk_size - 1, total_rows)
    
    if (verbose) {
      message("  Processing chunk ", chunk_count, " (rows ", start_row, "-", end_row, ")")
    }
    
    # Read chunk
    chunk_data <- read.csv(file_path, skip = start_row, nrows = chunk_size, 
                          header = FALSE, col.names = names(header))
    
    # Process chunk
    tryCatch({
      chunk_result <- process_function(chunk_data)
      
      # Write to output file if specified
      if (!is.null(output_file)) {
        write.csv(chunk_result, output_file, append = TRUE, row.names = FALSE)
      }
      
      results[[chunk_count]] <- chunk_result
      
    }, error = function(e) {
      warning("Error processing chunk ", chunk_count, ": ", e$message)
      results[[chunk_count]] <<- NULL
    })
  }
  
  if (verbose) {
    message("Stream processing completed. Processed ", chunk_count, " chunks.")
  }
  
  return(results)
}

#' Performance benchmarking
#' 
#' @param raw_df Input data frame
#' @param operations List of operations to benchmark
#' @param n_repeats Number of times to repeat each operation
#' @param verbose Whether to show progress messages
#' @return Benchmark results
#' @export
benchmark_performance <- function(raw_df,
                                 operations = list(
                                   correction = function(x) calcium_correction(x, verbose = FALSE),
                                   spikes = function(x) {
                                     corrected <- calcium_correction(x, verbose = FALSE)
                                     cell_cols <- names(corrected)[grepl("^Cell_", names(corrected))]
                                     lapply(cell_cols, function(col) {
                                       infer_spikes(corrected[[col]], verbose = FALSE)
                                     })
                                   }
                                 ),
                                 n_repeats = 3,
                                 verbose = TRUE) {
  
  if (verbose) {
    message("Starting performance benchmark...")
  }
  
  results <- list()
  
  for (op_name in names(operations)) {
    if (verbose) {
      message("Benchmarking: ", op_name)
    }
    
    times <- numeric(n_repeats)
    memory_usage <- numeric(n_repeats)
    
    for (i in seq_len(n_repeats)) {
      if (verbose) {
        message("  Repeat ", i, "/", n_repeats)
      }
      
      # Measure time and memory
      start_time <- Sys.time()
      start_mem <- monitor_memory_usage(paste0(op_name, "_start"), verbose = FALSE)
      
      # Run operation
      tryCatch({
        operation_result <- operations[[op_name]](raw_df)
        
        end_time <- Sys.time()
        end_mem <- monitor_memory_usage(paste0(op_name, "_end"), verbose = FALSE)
        
        times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
        memory_usage[i] <- end_mem$memory_mb - start_mem$memory_mb
        
      }, error = function(e) {
        warning("Operation ", op_name, " failed on repeat ", i, ": ", e$message)
        times[i] <<- NA
        memory_usage[i] <<- NA
      })
    }
    
    results[[op_name]] <- list(
      mean_time = mean(times, na.rm = TRUE),
      sd_time = sd(times, na.rm = TRUE),
      mean_memory = mean(memory_usage, na.rm = TRUE),
      sd_memory = sd(memory_usage, na.rm = TRUE),
      times = times,
      memory_usage = memory_usage
    )
  }
  
  if (verbose) {
    message("Benchmark completed")
    for (op_name in names(results)) {
      result <- results[[op_name]]
      message(op_name, ": ", round(result$mean_time, 3), " ± ", 
              round(result$sd_time, 3), " seconds, ",
              round(result$mean_memory, 1), " ± ", 
              round(result$sd_memory, 1), " MB")
    }
  }
  
  return(results)
} 