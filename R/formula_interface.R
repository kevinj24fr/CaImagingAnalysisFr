#' Formula Interface for Calcium Imaging Analysis
#'
#' R-native formula interface for specifying analysis models and
#' relationships. Enables intuitive specification of statistical models.
#'
#' @name formula_interface
#' @keywords internal
NULL

#' Analyze traces using formula interface
#'
#' Specify trace analysis using R formulas for intuitive model specification.
#'
#' @param formula Formula specifying the analysis. Left side is the response,
#'   right side specifies predictors and grouping.
#' @param data data.frame or data.table with trace data
#' @param method Analysis method: "lm", "glm", "mixed", "correlation"
#' @param ... Additional arguments passed to underlying model function
#'
#' @return Analysis result object
#' @export
#'
#' @examples
#' \dontrun{
#' # Linear model of fluorescence by time
#' result <- analyze_traces(dff ~ time | cell_id, data = trace_data)
#'
#' # Correlation analysis
#' result <- analyze_traces(cell_1 ~ cell_2, data = wide_data, method = "correlation")
#'
#' # Mixed effects model
#' result <- analyze_traces(response ~ condition + (1 | cell_id), data = data, method = "mixed")
#' }
analyze_traces <- function(formula, data, method = "lm", ...) {
  # Parse the formula
  parsed <- parse_analysis_formula(formula)

  if (method == "lm") {
    # Standard linear model per group
    if (!is.null(parsed$group)) {
      results <- by(data, data[[parsed$group]], function(subset) {
        lm(parsed$base_formula, data = subset, ...)
      })
      class(results) <- c("trace_lm_list", class(results))
      return(results)
    } else {
      return(lm(formula, data = data, ...))
    }
  } else if (method == "glm") {
    if (!is.null(parsed$group)) {
      results <- by(data, data[[parsed$group]], function(subset) {
        glm(parsed$base_formula, data = subset, ...)
      })
      class(results) <- c("trace_glm_list", class(results))
      return(results)
    } else {
      return(glm(formula, data = data, ...))
    }
  } else if (method == "mixed") {
    if (!requireNamespace("lme4", quietly = TRUE)) {
      stop("lme4 package required for mixed effects models")
    }
    return(lme4::lmer(formula, data = data, ...))
  } else if (method == "correlation") {
    return(compute_formula_correlation(formula, data, ...))
  }

  stop("Unknown method: ", method)
}

#' Parse analysis formula
#'
#' Internal function to parse formula with grouping.
#'
#' @param formula Formula to parse
#'
#' @return List with base_formula, response, predictors, group
#' @keywords internal
parse_analysis_formula <- function(formula) {
  # Get formula parts
  formula_str <- deparse(formula)

  # Check for grouping operator |
  group <- NULL
  if (grepl("\\|", formula_str)) {
    parts <- strsplit(formula_str, "\\|")[[1]]
    base_str <- trimws(parts[1])
    group <- trimws(parts[2])
    formula <- as.formula(base_str)
  }

  # Extract response and predictors
  terms_obj <- terms(formula)
  response <- all.vars(formula)[1]
  predictors <- attr(terms_obj, "term.labels")

  list(
    base_formula = formula,
    response = response,
    predictors = predictors,
    group = group,
    original = formula_str
  )
}

#' Compute correlation from formula
#'
#' @param formula Formula specifying variables
#' @param data data.frame with data
#' @param method Correlation method
#'
#' @return Correlation result
#' @keywords internal
compute_formula_correlation <- function(formula, data, method = "pearson", ...) {
  vars <- all.vars(formula)
  if (length(vars) < 2) {
    stop("Correlation requires at least two variables")
  }

  x <- data[[vars[1]]]
  y <- data[[vars[2]]]

  cor.test(x, y, method = method, ...)
}

#' Fit response model using formula
#'
#' Model neural responses to stimuli using formula interface.
#'
#' @param formula Formula: response ~ stimulus + covariates | cell_id
#' @param data data.frame with aligned response data
#' @param baseline_window Numeric vector c(start, end) for baseline
#' @param response_window Numeric vector c(start, end) for response
#'
#' @return Model fit results
#' @export
#'
#' @examples
#' \dontrun{
#' # Model response amplitude by stimulus type
#' result <- fit_response_model(
#'   amplitude ~ stimulus_type + trial_number | cell_id,
#'   data = aligned_data,
#'   baseline_window = c(-0.5, 0),
#'   response_window = c(0, 1)
#' )
#' }
fit_response_model <- function(formula, data, baseline_window = c(-1, 0),
                               response_window = c(0, 2)) {
  parsed <- parse_analysis_formula(formula)

  # Compute response metrics if needed
  if ("amplitude" %in% parsed$response || "response" %in% parsed$response) {
    data <- compute_response_from_data(data, baseline_window, response_window)
  }

  # Fit model
  if (!is.null(parsed$group)) {
    # Fit per cell
    cells <- unique(data[[parsed$group]])

    results <- lapply(cells, function(cell) {
      subset <- data[data[[parsed$group]] == cell, ]
      tryCatch(
        lm(parsed$base_formula, data = subset),
        error = function(e) NULL
      )
    })
    names(results) <- cells

    # Summarize
    coefs <- do.call(rbind, lapply(results, function(m) {
      if (is.null(m)) return(NULL)
      coef(summary(m))[, c("Estimate", "Std. Error", "Pr(>|t|)")]
    }))

    structure(
      list(
        models = results,
        coefficients = coefs,
        formula = formula
      ),
      class = "response_model_fit"
    )
  } else {
    lm(parsed$base_formula, data = data)
  }
}

#' Compute response from aligned data
#'
#' @param data data.frame with time and value columns
#' @param baseline_window Baseline time window
#' @param response_window Response time window
#'
#' @return data.frame with response metrics added
#' @keywords internal
compute_response_from_data <- function(data, baseline_window, response_window) {
  if (!all(c("time", "value") %in% names(data))) {
    return(data)
  }

  # Compute baseline and response per trial
  data$baseline <- NA
  data$response <- NA
  data$amplitude <- NA

  trials <- unique(data$trial_id)

  for (trial in trials) {
    idx <- data$trial_id == trial
    times <- data$time[idx]
    vals <- data$value[idx]

    bl_idx <- times >= baseline_window[1] & times <= baseline_window[2]
    resp_idx <- times >= response_window[1] & times <= response_window[2]

    baseline <- mean(vals[bl_idx], na.rm = TRUE)
    response <- mean(vals[resp_idx], na.rm = TRUE)

    data$baseline[idx] <- baseline
    data$response[idx] <- response
    data$amplitude[idx] <- response - baseline
  }

  data
}

#' Specify spike detection parameters using formula
#'
#' @param formula Formula: spikes ~ trace + threshold(2.5) + method("adaptive")
#' @param traces Matrix or data.frame of traces
#'
#' @return Spike detection result
#' @export
detect_spikes_formula <- function(formula, traces) {
  # Parse formula for parameters
  formula_str <- deparse(formula)

  # Extract parameters using regex
  threshold <- 2.5
  method <- "threshold"

  if (grepl("threshold\\(", formula_str)) {
    threshold <- as.numeric(gsub(".*threshold\\(([0-9.]+)\\).*", "\\1", formula_str))
  }

  if (grepl("method\\(", formula_str)) {
    method <- gsub(".*method\\([\"']([^\"']+)[\"']\\).*", "\\1", formula_str)
  }

  # Run detection
  if (is.data.frame(traces) || inherits(traces, "data.table")) {
    traces <- dt_to_traces(traces)
  }

  if (method == "adaptive") {
    adaptive_threshold_detection(traces, threshold_sd = threshold)
  } else {
    threshold_spike_detection(traces, threshold_sd = threshold)
  }
}

#' Define analysis pipeline using formula
#'
#' Create reusable analysis pipeline from formula specification.
#'
#' @param formula Pipeline formula: output ~ step1 + step2 + step3
#' @param ... Named arguments for pipeline steps
#'
#' @return Pipeline object
#' @export
#'
#' @examples
#' \dontrun{
#' # Define preprocessing pipeline
#' pipeline <- define_pipeline(
#'   processed ~ correct_motion + extract_traces + compute_dff,
#'   correct_motion = list(max_shift = 20),
#'   extract_traces = list(method = "mean"),
#'   compute_dff = list(baseline_frames = 1:100)
#' )
#'
#' # Apply pipeline
#' result <- run_pipeline(pipeline, movie)
#' }
define_pipeline <- function(formula, ...) {
  parsed <- parse_analysis_formula(formula)
  steps <- parsed$predictors
  params <- list(...)

  structure(
    list(
      steps = steps,
      params = params,
      formula = formula,
      output_name = parsed$response
    ),
    class = "analysis_pipeline"
  )
}

#' Run defined pipeline
#'
#' @param pipeline Pipeline object from define_pipeline
#' @param data Input data (movie, traces, etc.)
#'
#' @return Processed data
#' @export
run_pipeline <- function(pipeline, data) {
  if (!inherits(pipeline, "analysis_pipeline")) {
    stop("pipeline must be an analysis_pipeline object")
  }

  result <- data

  for (step in pipeline$steps) {
    step_params <- pipeline$params[[step]] %||% list()

    # Map step names to functions
    step_fun <- switch(step,
      correct_motion = motion_correct,
      extract_traces = extract_traces,
      compute_dff = compute_delta_f,
      detect_spikes = threshold_spike_detection,
      denoise = wavelet_denoise,
      neuropil_correct = neuropil_correct,
      zscore = function(x, ...) scale(t(x)),
      NULL
    )

    if (is.null(step_fun)) {
      # Try to find function by name
      if (exists(step, mode = "function")) {
        step_fun <- get(step, mode = "function")
      } else {
        warning("Unknown pipeline step: ", step)
        next
      }
    }

    result <- do.call(step_fun, c(list(result), step_params))
  }

  result
}

#' Print method for analysis pipeline
#'
#' @param x Pipeline object
#' @param ... Ignored
#'
#' @export
print.analysis_pipeline <- function(x, ...) {
  cat("Analysis Pipeline\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Steps:\n")
  for (i in seq_along(x$steps)) {
    step <- x$steps[i]
    params <- x$params[[step]]
    param_str <- if (length(params) > 0) {
      paste(names(params), "=", params, collapse = ", ")
    } else {
      "(default)"
    }
    cat(sprintf("  %d. %s %s\n", i, step, param_str))
  }
  invisible(x)
}

# Note: %||% operator defined in aaa_utils.R
