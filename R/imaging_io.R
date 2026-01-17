#' Imaging I/O Functions
#'
#' Functions for reading and writing calcium imaging data from various file formats.
#'
#' @name imaging_io
NULL

#' Read TIFF Stack
#'
#' Read a multi-page TIFF file (image stack) into a 3D array.
#'
#' @param path Path to TIFF file
#' @param frames Frame indices to read (NULL = all frames)
#' @param as_array Return as 3D array (TRUE) or list of matrices (FALSE)
#' @param verbose Show progress messages
#' @return 3D array [height x width x frames] or list of matrices
#'
#' @details
#' This function uses the `tiff` package if available, otherwise falls back to
#' a basic implementation. For large files, consider using `frames` parameter
#' to load subsets.
#'
#' @examples
#' \dontrun{
#' # Read entire stack
#' movie <- read_tiff_stack("recording.tif")
#'
#' # Read specific frames
#' movie <- read_tiff_stack("recording.tif", frames = 1:100)
#' }
#'
#' @export
read_tiff_stack <- function(path, frames = NULL, as_array = TRUE, verbose = TRUE) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  # Check for tiff package
 if (requireNamespace("tiff", quietly = TRUE)) {
    if (verbose) message("Reading TIFF with tiff package...")
    return(read_tiff_with_package(path, frames, as_array, verbose))
  }

  # Check for EBImage (Bioconductor)
  if (requireNamespace("EBImage", quietly = TRUE)) {
    if (verbose) message("Reading TIFF with EBImage...")
    return(read_tiff_with_ebimage(path, frames, as_array, verbose))
  }

  # Fallback to basic binary reading (limited support)
  stop("No TIFF reading package available. Install with:\n",
       "  install.packages('tiff')           # CRAN\n",
       "  BiocManager::install('EBImage')    # Bioconductor")
}

#' Read TIFF using tiff package
#' @keywords internal
read_tiff_with_package <- function(path, frames, as_array, verbose) {
  # Read all images from TIFF
  imgs <- tiff::readTIFF(path, all = TRUE, as.is = TRUE)

  if (!is.list(imgs)) {
    imgs <- list(imgs)
  }

  n_frames <- length(imgs)
  if (verbose) message("  Found ", n_frames, " frames")

  # Subset frames if requested
  if (!is.null(frames)) {
    frames <- frames[frames >= 1 & frames <= n_frames]
    imgs <- imgs[frames]
    if (verbose) message("  Selected ", length(imgs), " frames")
  }

  if (as_array && length(imgs) > 0) {
    # Convert to 3D array
    dims <- dim(imgs[[1]])
    arr <- array(0, dim = c(dims[1], dims[2], length(imgs)))
    for (i in seq_along(imgs)) {
      arr[, , i] <- imgs[[i]]
    }
    return(arr)
  }

  return(imgs)
}

#' Read TIFF using EBImage
#' @keywords internal
read_tiff_with_ebimage <- function(path, frames, as_array, verbose) {
  img <- EBImage::readImage(path, all = TRUE)

  # EBImage returns [width x height x frames]
  n_frames <- dim(img)[3]
  if (verbose) message("  Found ", n_frames, " frames")

  # Subset frames if requested
  if (!is.null(frames)) {
    frames <- frames[frames >= 1 & frames <= n_frames]
    img <- img[, , frames]
    if (verbose) message("  Selected ", length(frames), " frames")
  }

  if (as_array) {
    # Transpose to [height x width x frames] for consistency
    return(aperm(img, c(2, 1, 3)))
  }

  # Return as list of matrices
  lapply(seq_len(dim(img)[3]), function(i) t(img[, , i]))
}

#' Write TIFF Stack
#'
#' Write a 3D array to a multi-page TIFF file.
#'
#' @param data 3D array [height x width x frames] or list of matrices
#' @param path Output file path
#' @param bits Bit depth (8 or 16)
#' @param compression Compression type ("none", "LZW", "deflate")
#' @param verbose Show progress messages
#' @return Path to written file (invisibly)
#'
#' @export
write_tiff_stack <- function(data, path, bits = 16, compression = "none", verbose = TRUE) {
  if (!requireNamespace("tiff", quietly = TRUE)) {
    stop("Package 'tiff' required. Install with: install.packages('tiff')")
  }

  # Convert list to array if needed
  if (is.list(data)) {
    dims <- dim(data[[1]])
    arr <- array(0, dim = c(dims[1], dims[2], length(data)))
    for (i in seq_along(data)) {
      arr[, , i] <- data[[i]]
    }
    data <- arr
  }

  if (verbose) {
    message("Writing TIFF: ", dim(data)[1], "x", dim(data)[2], " x ", dim(data)[3], " frames")
  }

  # Normalize to 0-1 range for tiff package
  data_range <- range(data, na.rm = TRUE)
  if (data_range[2] > 1 || data_range[1] < 0) {
    data <- (data - data_range[1]) / (data_range[2] - data_range[1])
  }

  # Write frame by frame
  for (i in seq_len(dim(data)[3])) {
    tiff::writeTIFF(
      data[, , i],
      path,
      bits.per.sample = bits,
      compression = compression,
      reduce = FALSE,
      append = (i > 1)
    )
  }

  if (verbose) message("Written to: ", path)
  invisible(path)
}

#' Read Imaging Metadata
#'
#' Extract metadata from imaging files (ScanImage, Prairie View, etc.)
#'
#' @param path Path to imaging file or metadata file
#' @param format Format hint ("auto", "scanimage", "prairie", "micromanager")
#' @return List of metadata fields
#'
#' @export
read_imaging_metadata <- function(path, format = "auto") {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  # Auto-detect format
  if (format == "auto") {
    ext <- tolower(tools::file_ext(path))
    if (ext == "xml") {
      format <- "prairie"
    } else if (ext %in% c("tif", "tiff")) {
      format <- "scanimage"
    } else if (ext == "json") {
      format <- "micromanager"
    }
  }

  metadata <- switch(format,
    "scanimage" = read_scanimage_metadata(path),
    "prairie" = read_prairie_metadata(path),
    "micromanager" = read_micromanager_metadata(path),
    list(format = "unknown", warning = "Could not detect metadata format")
  )

  metadata$source_file <- path
  metadata$format <- format
  return(metadata)
}

#' Read ScanImage TIFF metadata
#' @keywords internal
read_scanimage_metadata <- function(path) {
  metadata <- list()

  # Try to read TIFF tags
  if (requireNamespace("tiff", quietly = TRUE)) {
    tryCatch({
      # Read first frame to get tags
      info <- tiff::readTIFF(path, payload = FALSE)

      if (!is.null(attr(info, "description"))) {
        desc <- attr(info, "description")

        # Parse ScanImage header
        if (grepl("SI\\.", desc)) {
          lines <- strsplit(desc, "\n")[[1]]

          for (line in lines) {
            if (grepl("=", line)) {
              parts <- strsplit(line, "=")[[1]]
              key <- trimws(parts[1])
              value <- trimws(parts[2])
              metadata[[key]] <- tryCatch(eval(parse(text = value)), error = function(e) value)
            }
          }
        }
      }
    }, error = function(e) NULL)
  }

  # Add basic file info
  metadata$file_size_mb <- file.info(path)$size / 1024^2

  return(metadata)
}

#' Read Prairie View XML metadata
#' @keywords internal
read_prairie_metadata <- function(path) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    return(list(error = "Package 'xml2' required for Prairie metadata"))
  }

  metadata <- list()

  tryCatch({
    doc <- xml2::read_xml(path)

    # Extract common fields
    metadata$date <- xml2::xml_text(xml2::xml_find_first(doc, "//Date"))
    metadata$time <- xml2::xml_text(xml2::xml_find_first(doc, "//Time"))

    # Frame rate
    frame_period <- xml2::xml_find_first(doc, "//PVStateValue[@key='framePeriod']")
    if (!is.na(frame_period)) {
      metadata$frame_period <- as.numeric(xml2::xml_attr(frame_period, "value"))
      metadata$frame_rate <- 1 / metadata$frame_period
    }

    # Image dimensions
    pixels_x <- xml2::xml_find_first(doc, "//PVStateValue[@key='pixelsPerLine']")
    pixels_y <- xml2::xml_find_first(doc, "//PVStateValue[@key='linesPerFrame']")
    if (!is.na(pixels_x)) metadata$width <- as.integer(xml2::xml_attr(pixels_x, "value"))
    if (!is.na(pixels_y)) metadata$height <- as.integer(xml2::xml_attr(pixels_y, "value"))

    # Zoom/magnification
    zoom <- xml2::xml_find_first(doc, "//PVStateValue[@key='opticalZoom']")
    if (!is.na(zoom)) metadata$zoom <- as.numeric(xml2::xml_attr(zoom, "value"))

  }, error = function(e) {
    metadata$error <- e$message
  })

  return(metadata)
}

#' Read Micro-Manager metadata
#' @keywords internal
read_micromanager_metadata <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(list(error = "Package 'jsonlite' required for Micro-Manager metadata"))
  }

  tryCatch({
    jsonlite::fromJSON(path)
  }, error = function(e) {
    list(error = e$message)
  })
}

#' Read Multiple Imaging Formats
#'
#' Unified interface for reading calcium imaging data from various formats.
#'
#' @param path Path to file or directory
#' @param format Format ("auto", "tiff", "hdf5", "nwb", "mat", "csv")
#' @param dataset For HDF5/NWB, the dataset path
#' @param verbose Show progress messages
#' @return List with 'data' (array or data frame) and 'metadata'
#'
#' @details
#' Supported formats:
#' - TIFF stacks (.tif, .tiff)
#' - HDF5 files (.h5, .hdf5) - requires rhdf5 package
#' - NWB files (.nwb) - basic HDF5 support
#' - MATLAB files (.mat) - requires R.matlab package
#' - CSV/TSV files (.csv, .tsv)
#'
#' @export
read_calcium_imaging <- function(path, format = "auto", dataset = NULL, verbose = TRUE) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  # Auto-detect format
  if (format == "auto") {
    ext <- tolower(tools::file_ext(path))
    format <- switch(ext,
      "tif" = "tiff", "tiff" = "tiff",
      "h5" = "hdf5", "hdf5" = "hdf5",
      "nwb" = "nwb",
      "mat" = "mat",
      "csv" = "csv", "tsv" = "csv",
      "unknown"
    )
  }

  result <- switch(format,
    "tiff" = {
      list(
        data = read_tiff_stack(path, verbose = verbose),
        metadata = read_imaging_metadata(path, format = "scanimage")
      )
    },
    "hdf5" = read_hdf5_calcium(path, dataset, verbose),
    "nwb" = read_nwb_calcium(path, verbose),
    "mat" = read_mat_calcium(path, verbose),
    "csv" = {
      list(
        data = utils::read.csv(path, stringsAsFactors = FALSE),
        metadata = list(format = "csv", file = path)
      )
    },
    stop("Unsupported format: ", format)
  )

  result$source <- path
  result$format <- format
  return(result)
}

#' Read HDF5 calcium imaging data
#' @keywords internal
read_hdf5_calcium <- function(path, dataset, verbose) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' required. Install with: BiocManager::install('rhdf5')")
  }

  # List datasets if not specified
  if (is.null(dataset)) {
    datasets <- rhdf5::h5ls(path)
    if (verbose) {
      message("Available datasets:")
      print(datasets[, c("group", "name", "dim")])
    }

    # Try to find common dataset names
    common_names <- c("data", "mov", "movie", "imaging", "fluorescence", "dff", "F")
    for (name in common_names) {
      matches <- datasets[grepl(name, datasets$name, ignore.case = TRUE), ]
      if (nrow(matches) > 0) {
        dataset <- paste0(matches$group[1], "/", matches$name[1])
        dataset <- gsub("^/+", "/", dataset)
        if (verbose) message("Auto-selected dataset: ", dataset)
        break
      }
    }

    if (is.null(dataset)) {
      stop("Specify 'dataset' parameter. Use rhdf5::h5ls('", path, "') to list datasets.")
    }
  }

  data <- rhdf5::h5read(path, dataset)
  rhdf5::h5closeAll()

  list(
    data = data,
    metadata = list(dataset = dataset)
  )
}

#' Read NWB calcium imaging data
#' @keywords internal
read_nwb_calcium <- function(path, verbose) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("Package 'rhdf5' required. Install with: BiocManager::install('rhdf5')")
  }

  # NWB files store data in specific locations
  result <- list(metadata = list())

  # Try common NWB paths for calcium imaging
  nwb_paths <- c(
    "/processing/ophys/Fluorescence/RoiResponseSeries/data",
    "/processing/ophys/DfOverF/RoiResponseSeries/data",
    "/acquisition/TwoPhotonSeries/data"
  )

  data <- NULL
  for (nwb_path in nwb_paths) {
    tryCatch({
      data <- rhdf5::h5read(path, nwb_path)
      result$metadata$nwb_path <- nwb_path
      if (verbose) message("Found data at: ", nwb_path)
      break
    }, error = function(e) NULL)
  }

  if (is.null(data)) {
    if (verbose) {
      message("Could not find standard NWB paths. Listing available datasets:")
      print(rhdf5::h5ls(path))
    }
    stop("Could not find calcium imaging data in NWB file")
  }

  # Try to read ROI info
  tryCatch({
    rois <- rhdf5::h5read(path, "/processing/ophys/ImageSegmentation/PlaneSegmentation")
    result$metadata$n_rois <- nrow(rois)
  }, error = function(e) NULL)

  rhdf5::h5closeAll()

  result$data <- data
  return(result)
}

#' Read MATLAB calcium imaging data
#' @keywords internal
read_mat_calcium <- function(path, verbose) {
  if (!requireNamespace("R.matlab", quietly = TRUE)) {
    stop("Package 'R.matlab' required. Install with: install.packages('R.matlab')")
  }

  mat_data <- R.matlab::readMat(path)

  # Look for common variable names
  common_names <- c("data", "F", "dff", "C", "traces", "fluorescence", "mov", "Y")
  data <- NULL

  for (name in names(mat_data)) {
    if (tolower(name) %in% tolower(common_names)) {
      data <- mat_data[[name]]
      if (verbose) message("Using variable: ", name)
      break
    }
  }

  if (is.null(data)) {
    # Use first numeric array
    for (name in names(mat_data)) {
      if (is.numeric(mat_data[[name]]) && length(dim(mat_data[[name]])) >= 2) {
        data <- mat_data[[name]]
        if (verbose) message("Using first numeric array: ", name)
        break
      }
    }
  }

  list(
    data = data,
    metadata = list(variables = names(mat_data)),
    all_data = mat_data
  )
}

#' Extract Traces from ROIs
#'
#' Extract fluorescence traces from a movie given ROI masks.
#'
#' @param movie 3D array [height x width x frames]
#' @param rois List of ROI masks (each a logical matrix) or labeled matrix
#' @param method Extraction method ("mean", "median", "weighted", "neuropil_corrected")
#' @param neuropil_radius For neuropil_corrected, radius of neuropil ring
#' @param neuropil_coefficient For neuropil_corrected, subtraction coefficient
#' @param verbose Show progress
#' @return Data frame with Time and Cell_* columns
#'
#' @export
extract_traces_from_rois <- function(movie,
                                      rois,
                                      method = "mean",
                                      neuropil_radius = 15,
                                      neuropil_coefficient = 0.7,
                                      verbose = TRUE) {
  n_frames <- dim(movie)[3]
  height <- dim(movie)[1]
  width <- dim(movie)[2]

  # Convert labeled matrix to list of masks
  if (is.matrix(rois) && !is.logical(rois)) {
    labels <- unique(as.vector(rois))
    labels <- labels[labels > 0]
    roi_list <- lapply(labels, function(l) rois == l)
    names(roi_list) <- paste0("Cell_", seq_along(roi_list))
    rois <- roi_list
  }

  n_rois <- length(rois)
  if (verbose) message("Extracting traces from ", n_rois, " ROIs across ", n_frames, " frames")

  # Initialize output
  traces <- matrix(0, nrow = n_frames, ncol = n_rois)
  colnames(traces) <- if (!is.null(names(rois))) names(rois) else paste0("Cell_", 1:n_rois)

  # Extract traces
  for (i in seq_along(rois)) {
    mask <- rois[[i]]

    if (method == "neuropil_corrected") {
      # Create neuropil mask (ring around cell)
      dilated <- dilate_mask(mask, neuropil_radius)
      eroded <- erode_mask(mask, 2)
      neuropil_mask <- dilated & !mask

      for (t in seq_len(n_frames)) {
        frame <- movie[, , t]
        cell_signal <- mean(frame[mask], na.rm = TRUE)
        neuropil_signal <- mean(frame[neuropil_mask], na.rm = TRUE)
        traces[t, i] <- cell_signal - neuropil_coefficient * neuropil_signal
      }
    } else {
      for (t in seq_len(n_frames)) {
        frame <- movie[, , t]
        pixels <- frame[mask]

        traces[t, i] <- switch(method,
          "mean" = mean(pixels, na.rm = TRUE),
          "median" = median(pixels, na.rm = TRUE),
          "weighted" = {
            weights <- pixels / sum(pixels, na.rm = TRUE)
            sum(pixels * weights, na.rm = TRUE)
          },
          mean(pixels, na.rm = TRUE)
        )
      }
    }

    if (verbose && i %% 10 == 0) message("  Processed ", i, "/", n_rois, " ROIs")
  }

  # Convert to data frame
  result <- as.data.frame(traces)
  result$Time <- seq_len(n_frames)

  return(result)
}

#' Helper: Dilate mask
#' @keywords internal
dilate_mask <- function(mask, radius) {
  if (requireNamespace("EBImage", quietly = TRUE)) {
    kernel <- EBImage::makeBrush(2 * radius + 1, shape = "disc")
    return(EBImage::dilate(mask, kernel) > 0)
  }

  # Simple dilation without EBImage
  result <- mask
  for (dx in -radius:radius) {
    for (dy in -radius:radius) {
      if (dx^2 + dy^2 <= radius^2) {
        shifted <- matrix(FALSE, nrow(mask), ncol(mask))
        rows <- max(1, 1 - dy):min(nrow(mask), nrow(mask) - dy)
        cols <- max(1, 1 - dx):min(ncol(mask), ncol(mask) - dx)
        shifted[rows + dy, cols + dx] <- mask[rows, cols]
        result <- result | shifted
      }
    }
  }
  return(result)
}

#' Helper: Erode mask
#' @keywords internal
erode_mask <- function(mask, radius) {
  if (requireNamespace("EBImage", quietly = TRUE)) {
    kernel <- EBImage::makeBrush(2 * radius + 1, shape = "disc")
    return(EBImage::erode(mask, kernel) > 0)
  }

  # Simple erosion: a pixel stays TRUE only if all neighbors are TRUE
  result <- mask
  for (dx in -radius:radius) {
    for (dy in -radius:radius) {
      if (dx^2 + dy^2 <= radius^2) {
        shifted <- matrix(TRUE, nrow(mask), ncol(mask))
        rows <- max(1, 1 - dy):min(nrow(mask), nrow(mask) - dy)
        cols <- max(1, 1 - dx):min(ncol(mask), ncol(mask) - dx)
        shifted[rows + dy, cols + dx] <- mask[rows, cols]
        result <- result & shifted
      }
    }
  }
  return(result)
}
