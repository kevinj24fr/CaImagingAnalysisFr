#' S7 Class Definitions for Calcium Imaging Analysis
#'
#' Type-safe class system using R's modern S7 OOP framework.
#' These classes provide runtime validation, method dispatch, and
#' better integration with modern R tooling.
#'
#' @name s7_classes
#' @keywords internal
NULL

# Check if S7 is available
.s7_available <- function() {

requireNamespace("S7", quietly = TRUE)
}

#' CalciumMovie Class
#'
#' Represents a calcium imaging movie with metadata.
#'
#' @param data 3D array (height x width x frames) or file path
#' @param frame_rate Acquisition frame rate in Hz
#' @param pixel_size Pixel size in micrometers (optional)
#' @param metadata Named list of additional metadata
#'
#' @return A CalciumMovie S7 object
#' @export
#'
#' @examples
#' \dontrun{
#' movie <- CalciumMovie(
#'   data = array(rnorm(100*100*500), c(100, 100, 500)),
#'   frame_rate = 30
#' )
#' }
CalciumMovie <- NULL

#' SpikeTrains Class
#'
#' Represents inferred spike trains from calcium traces.
#'
#' @param spikes Binary matrix (cells x time) or sparse representation
#' @param times Numeric vector of time points
#' @param cell_ids Character vector of cell identifiers
#' @param method Method used for spike inference
#'
#' @return A SpikeTrains S7 object
#' @export
SpikeTrains <- NULL
#' ROISet Class
#'
#' Represents a collection of regions of interest (ROIs).
#'
#' @param masks List of binary masks or coordinate matrices
#' @param labels Character vector of ROI labels
#' @param properties data.frame of ROI properties (area, centroid, etc.)
#' @param image_dims Integer vector of source image dimensions
#'
#' @return An ROISet S7 object
#' @export
ROISet <- NULL

#' CalciumTraces Class
#'
#' Represents extracted fluorescence traces.
#'
#' @param traces Matrix of traces (cells x time)
#' @param time_vector Numeric vector of timestamps
#' @param cell_ids Character vector of cell identifiers
#' @param trace_type Type of trace ("raw", "dff", "denoised", "deconvolved")
#'
#' @return A CalciumTraces S7 object
#' @export
CalciumTraces <- NULL

#' ExperimentSession Class
#'
#' Container for a complete imaging session.
#'
#' @param movie CalciumMovie object (optional, can be lazy-loaded)
#' @param rois ROISet object
#' @param traces CalciumTraces object
#' @param spikes SpikeTrains object (optional)
#' @param metadata Named list of session metadata
#'
#' @return An ExperimentSession S7 object
#' @export
ExperimentSession <- NULL

# Initialize S7 classes when package loads
.init_s7_classes <- function() {
  if (!.s7_available()) {
    return(invisible(NULL))
  }

  S7 <- asNamespace("S7")

  # CalciumMovie class
  CalciumMovie <<- S7$new_class(
    name = "CalciumMovie",
    properties = list(
      data = S7$class_any,
      frame_rate = S7$new_property(
        class = S7$class_numeric,
        validator = function(value) {
          if (length(value) != 1 || value <= 0) {
            return("frame_rate must be a positive scalar")
          }
          NULL
        }
      ),
      pixel_size = S7$new_property(
        class = S7$class_numeric,
        default = NA_real_
      ),
      n_frames = S7$new_property(
        class = S7$class_integer,
        getter = function(self) {
          if (is.array(self@data)) {
            as.integer(dim(self@data)[3])
          } else {
            NA_integer_
          }
        }
      ),
      dims = S7$new_property(
        class = S7$class_integer,
        getter = function(self) {
          if (is.array(self@data)) {
            as.integer(dim(self@data)[1:2])
          } else {
            c(NA_integer_, NA_integer_)
          }
        }
      ),
      duration = S7$new_property(
        class = S7$class_numeric,
        getter = function(self) {
          if (!is.na(self@n_frames)) {
            self@n_frames / self@frame_rate
          } else {
            NA_real_
          }
        }
      ),
      metadata = S7$new_property(
        class = S7$class_list,
        default = list()
      )
    ),
    validator = function(self) {
      if (is.array(self@data) && length(dim(self@data)) != 3) {
        return("data must be a 3D array (height x width x frames)")
      }
      NULL
    }
  )

  # CalciumTraces class
  CalciumTraces <<- S7$new_class(
    name = "CalciumTraces",
    properties = list(
      traces = S7$new_property(
        class = S7$class_any,
        validator = function(value) {
          if (!is.matrix(value) && !inherits(value, "data.table")) {
            return("traces must be a matrix or data.table")
          }
          NULL
        }
      ),
      time_vector = S7$new_property(
        class = S7$class_numeric,
        default = NULL
      ),
      cell_ids = S7$new_property(
        class = S7$class_character,
        default = NULL
      ),
      trace_type = S7$new_property(
        class = S7$class_character,
        default = "raw",
        validator = function(value) {
          valid <- c("raw", "dff", "denoised", "deconvolved", "corrected")
          if (!value %in% valid) {
            return(paste("trace_type must be one of:", paste(valid, collapse = ", ")))
          }
          NULL
        }
      ),
      n_cells = S7$new_property(
        class = S7$class_integer,
        getter = function(self) {
          as.integer(nrow(self@traces))
        }
      ),
      n_frames = S7$new_property(
        class = S7$class_integer,
        getter = function(self) {
          as.integer(ncol(self@traces))
        }
      )
    )
  )

  # SpikeTrains class
  SpikeTrains <<- S7$new_class(
    name = "SpikeTrains",
    properties = list(
      spikes = S7$new_property(
        class = S7$class_any,
        validator = function(value) {
          if (!is.matrix(value) && !inherits(value, "sparseMatrix")) {
            return("spikes must be a matrix or sparse matrix")
          }
          NULL
        }
      ),
      times = S7$new_property(
        class = S7$class_numeric,
        default = NULL
      ),
      cell_ids = S7$new_property(
        class = S7$class_character,
        default = NULL
      ),
      method = S7$new_property(
        class = S7$class_character,
        default = "unknown"
      ),
      n_cells = S7$new_property(
        class = S7$class_integer,
        getter = function(self) {
          as.integer(nrow(self@spikes))
        }
      ),
      spike_counts = S7$new_property(
        class = S7$class_numeric,
        getter = function(self) {
          rowSums(self@spikes > 0)
        }
      )
    )
  )

  # ROISet class
  ROISet <<- S7$new_class(
    name = "ROISet",
    properties = list(
      masks = S7$new_property(
        class = S7$class_list,
        validator = function(value) {
          if (length(value) == 0) {
            return("masks cannot be empty")
          }
          NULL
        }
      ),
      labels = S7$new_property(
        class = S7$class_character,
        default = NULL
      ),
      properties = S7$new_property(
        class = S7$class_any,
        default = NULL
      ),
      image_dims = S7$new_property(
        class = S7$class_integer,
        default = c(NA_integer_, NA_integer_)
      ),
      n_rois = S7$new_property(
        class = S7$class_integer,
        getter = function(self) {
          as.integer(length(self@masks))
        }
      )
    )
  )

  # ExperimentSession class
  ExperimentSession <<- S7$new_class(
    name = "ExperimentSession",
    properties = list(
      movie = S7$new_property(
        class = S7$class_any,
        default = NULL
      ),
      rois = S7$new_property(
        class = S7$class_any,
        default = NULL
      ),
      traces = S7$new_property(
        class = S7$class_any,
        default = NULL
      ),
      spikes = S7$new_property(
        class = S7$class_any,
        default = NULL
      ),
      metadata = S7$new_property(
        class = S7$class_list,
        default = list()
      ),
      session_id = S7$new_property(
        class = S7$class_character,
        default = NULL
      )
    )
  )

  invisible(NULL)
}

#' Convert matrix to CalciumTraces object
#'
#' @param traces Matrix of traces (cells x time)
#' @param frame_rate Frame rate for generating time vector
#' @param trace_type Type of trace
#'
#' @return CalciumTraces object
#' @export
as_calcium_traces <- function(traces, frame_rate = 30, trace_type = "raw") {
  if (!.s7_available()) {
    stop("S7 package required for class-based interface. Install with: install.packages('S7')")
  }

  if (is.null(CalciumTraces)) {
    .init_s7_classes()
  }

  time_vector <- seq(0, ncol(traces) - 1) / frame_rate
  cell_ids <- paste0("cell_", seq_len(nrow(traces)))

  CalciumTraces(
    traces = traces,
    time_vector = time_vector,
    cell_ids = cell_ids,
    trace_type = trace_type
  )
}

#' Convert array to CalciumMovie object
#'
#' @param data 3D array (height x width x frames)
#' @param frame_rate Acquisition frame rate
#' @param pixel_size Pixel size in micrometers
#'
#' @return CalciumMovie object
#' @export
as_calcium_movie <- function(data, frame_rate = 30, pixel_size = NA_real_) {
  if (!.s7_available()) {
    stop("S7 package required for class-based interface. Install with: install.packages('S7')")
  }

  if (is.null(CalciumMovie)) {
    .init_s7_classes()
  }

  CalciumMovie(
    data = data,
    frame_rate = frame_rate,
    pixel_size = pixel_size
  )
}

#' Convert list of masks to ROISet object
#'
#' @param masks List of binary masks
#' @param image_dims Dimensions of source image
#'
#' @return ROISet object
#' @export
as_roi_set <- function(masks, image_dims = NULL) {
  if (!.s7_available()) {
    stop("S7 package required for class-based interface. Install with: install.packages('S7')")
  }

  if (is.null(ROISet)) {
    .init_s7_classes()
  }

  labels <- paste0("roi_", seq_along(masks))

  if (is.null(image_dims) && length(masks) > 0) {
    if (is.matrix(masks[[1]])) {
      image_dims <- as.integer(dim(masks[[1]]))
    }
  }

  ROISet(
    masks = masks,
    labels = labels,
    image_dims = image_dims %||% c(NA_integer_, NA_integer_)
  )
}

#' Check if object is a CaImaging S7 class
#'
#' @param x Object to check
#'
#' @return Logical
#' @export
is_calcium_object <- function(x) {
  if (!.s7_available()) {
    return(FALSE)
  }
  inherits(x, c("CalciumMovie", "CalciumTraces", "SpikeTrains", "ROISet", "ExperimentSession"))
}

# Note: %||% operator defined in aaa_utils.R
