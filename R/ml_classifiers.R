#' Classical Machine Learning for Calcium Imaging
#'
#' Proven ML methods that R excels at: random forests, gradient boosting,
#' SVMs, and regularized regression for cell classification, neural decoding,
#' and event detection.
#'
#' @name ml_classifiers
#' @keywords internal
NULL

# ============================================================================
# Cell Type Classification
# ============================================================================

#' Classify cell types from activity patterns
#'
#' Uses classical ML to classify neurons into functional types based on
#' their calcium activity features.
#'
#' @param traces Matrix of traces (cells x time) or features (cells x features)
#' @param labels Factor or character vector of known cell type labels (for training)
#' @param method Classification method: "ranger", "xgboost", "svm", "glmnet"
#' @param features Character vector of features to extract, or NULL if traces
#'   already contains features. Options: "temporal", "spectral", "statistical", "all"
#' @param train_fraction Fraction of data for training (rest for validation)
#' @param ... Additional arguments passed to the classifier
#'
#' @return List with model, predictions, accuracy, and feature importance
#' @export
#'
#' @examples
#' \dontrun{
#' # With known labels for some cells
#' labels <- c("excitatory", "inhibitory", "excitatory", NA, NA, "inhibitory", NA)
#' result <- classify_cell_types(traces, labels, method = "ranger")
#'
#' # Predict unlabeled cells
#' predicted <- result$predictions
#' }
classify_cell_types <- function(traces, labels = NULL, method = "ranger",
                                 features = "all", train_fraction = 0.7, ...) {
  # Extract features if needed
  if (!is.null(features)) {
    feature_matrix <- extract_trace_features(traces, features)
  } else {
    feature_matrix <- traces
  }

  # Handle labels
  if (is.null(labels)) {
    stop("labels required for supervised classification. Use cluster_cell_types() for unsupervised.")
  }

  labels <- as.factor(labels)
  labeled_idx <- which(!is.na(labels))

  if (length(labeled_idx) < 10) {
    stop("Need at least 10 labeled examples for training")
  }

  # Split data
  n_train <- floor(length(labeled_idx) * train_fraction)
  train_idx <- sample(labeled_idx, n_train)
  test_idx <- setdiff(labeled_idx, train_idx)

  train_X <- feature_matrix[train_idx, , drop = FALSE]
  train_y <- labels[train_idx]
  test_X <- feature_matrix[test_idx, , drop = FALSE]
  test_y <- labels[test_idx]

  # Train model
  model <- switch(method,
    ranger = .train_ranger_classifier(train_X, train_y, ...),
    xgboost = .train_xgboost_classifier(train_X, train_y, ...),
    svm = .train_svm_classifier(train_X, train_y, ...),
    glmnet = .train_glmnet_classifier(train_X, train_y, ...),
    stop("Unknown method: ", method)
  )

  # Predictions on test set
  test_pred <- predict_cell_types(model, test_X, method = method)

  # Accuracy
  accuracy <- mean(test_pred == test_y)

  # Predict all cells (including unlabeled)
  all_pred <- predict_cell_types(model, feature_matrix, method = method)

  # Feature importance
  importance <- .get_feature_importance(model, method, colnames(feature_matrix))

  list(
    model = model,
    method = method,
    predictions = all_pred,
    test_accuracy = accuracy,
    confusion_matrix = table(Predicted = test_pred, Actual = test_y),
    feature_importance = importance,
    features_used = colnames(feature_matrix),
    labeled_indices = labeled_idx,
    train_indices = train_idx,
    test_indices = test_idx
  )
}

#' Predict cell types using trained model
#'
#' @param model Trained model from classify_cell_types
#' @param new_data Feature matrix for new cells
#' @param method Method used (must match training)
#'
#' @return Factor of predicted cell types
#' @export
predict_cell_types <- function(model, new_data, method = "ranger") {
  switch(method,
    ranger = {
      if (!requireNamespace("ranger", quietly = TRUE)) {
        stop("ranger package required")
      }
      predict(model, data = as.data.frame(new_data))$predictions
    },
    xgboost = {
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("xgboost package required")
      }
      probs <- predict(model, as.matrix(new_data))
      if (is.matrix(probs)) {
        factor(model$classes[apply(probs, 1, which.max)], levels = model$classes)
      } else {
        factor(ifelse(probs > 0.5, model$classes[2], model$classes[1]), levels = model$classes)
      }
    },
    svm = {
      if (!requireNamespace("e1071", quietly = TRUE)) {
        stop("e1071 package required")
      }
      predict(model, as.data.frame(new_data))
    },
    glmnet = {
      if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("glmnet package required")
      }
      as.factor(predict(model$model, as.matrix(new_data), s = model$lambda, type = "class"))
    }
  )
}

# Training helpers
.train_ranger_classifier <- function(X, y, num.trees = 500, importance = "impurity", ...) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("ranger package required. Install with: install.packages('ranger')")
  }

  df <- as.data.frame(X)
  df$.label <- y

  ranger::ranger(
    .label ~ .,
    data = df,
    num.trees = num.trees,
    importance = importance,
    ...
  )
}

.train_xgboost_classifier <- function(X, y, nrounds = 100, ...) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("xgboost package required. Install with: install.packages('xgboost')")
  }

  classes <- levels(y)
  y_numeric <- as.integer(y) - 1

  if (length(classes) == 2) {
    model <- xgboost::xgboost(
      data = as.matrix(X),
      label = y_numeric,
      nrounds = nrounds,
      objective = "binary:logistic",
      verbose = 0,
      ...
    )
  } else {
    model <- xgboost::xgboost(
      data = as.matrix(X),
      label = y_numeric,
      nrounds = nrounds,
      objective = "multi:softprob",
      num_class = length(classes),
      verbose = 0,
      ...
    )
  }

  model$classes <- classes
  model
}

.train_svm_classifier <- function(X, y, kernel = "radial", ...) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("e1071 package required. Install with: install.packages('e1071')")
  }

  e1071::svm(
    x = as.data.frame(X),
    y = y,
    kernel = kernel,
    probability = TRUE,
    ...
  )
}

.train_glmnet_classifier <- function(X, y, alpha = 0.5, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet package required. Install with: install.packages('glmnet')")
  }

  cv_fit <- glmnet::cv.glmnet(
    x = as.matrix(X),
    y = y,
    family = "multinomial",
    alpha = alpha,
    ...
  )

  list(
    model = cv_fit,
    lambda = cv_fit$lambda.min
  )
}

.get_feature_importance <- function(model, method, feature_names) {
  imp <- switch(method,
    ranger = {
      if (!is.null(model$variable.importance)) {
        sort(model$variable.importance, decreasing = TRUE)
      } else {
        NULL
      }
    },
    xgboost = {
      if (requireNamespace("xgboost", quietly = TRUE)) {
        imp <- xgboost::xgb.importance(model = model)
        setNames(imp$Gain, imp$Feature)
      } else {
        NULL
      }
    },
    svm = NULL,  # SVM doesn't have built-in importance
    glmnet = {
      coefs <- coef(model$model, s = model$lambda)
      if (is.list(coefs)) {
        # Multinomial: average absolute coefficients across classes
        avg_coef <- Reduce(`+`, lapply(coefs, function(x) abs(as.matrix(x)))) / length(coefs)
        imp <- avg_coef[-1, 1]  # Remove intercept
        sort(imp, decreasing = TRUE)
      } else {
        imp <- abs(as.matrix(coefs))[-1, 1]
        sort(imp, decreasing = TRUE)
      }
    }
  )

  imp
}

# ============================================================================
# Neural Decoding / Behavior Prediction
# ============================================================================

#' Decode behavior from neural activity
#'
#' Predict behavioral variables from calcium traces using classical ML.
#'
#' @param traces Matrix of traces (cells x time) or (time x cells)
#' @param behavior Vector or matrix of behavioral variables to predict
#' @param method Decoding method: "ridge", "lasso", "elastic_net", "ranger", "xgboost"
#' @param cv_folds Number of cross-validation folds
#' @param time_lags Integer vector of time lags to include (e.g., -2:2)
#' @param transpose If TRUE, traces is (time x cells); if FALSE, (cells x time)
#' @param ... Additional arguments to the method
#'
#' @return List with model, predictions, R-squared, and weights
#' @export
#'
#' @examples
#' \dontrun{
#' # Decode velocity from neural activity
#' result <- decode_behavior(traces, velocity, method = "ridge")
#'
#' # With time lags
#' result <- decode_behavior(traces, position, method = "elastic_net", time_lags = -3:3)
#' }
decode_behavior <- function(traces, behavior, method = "ridge", cv_folds = 5,
                            time_lags = 0, transpose = FALSE, ...) {
  # Ensure traces is time x cells
  if (!transpose) {
    traces <- t(traces)
  }

  n_time <- nrow(traces)
  n_cells <- ncol(traces)

  # Handle behavior
  if (is.vector(behavior)) {
    behavior <- matrix(behavior, ncol = 1)
  }

  if (nrow(behavior) != n_time) {
    stop("behavior must have same number of time points as traces")
  }

  # Create lagged features
  if (length(time_lags) > 1 || time_lags[1] != 0) {
    X <- create_lagged_features(traces, time_lags)
    # Trim behavior to match
    max_lag <- max(abs(time_lags))
    valid_idx <- (max_lag + 1):(n_time - max_lag)
    y <- behavior[valid_idx, , drop = FALSE]
  } else {
    X <- traces
    y <- behavior
  }

  # Train-test split (80-20)
  n <- nrow(X)
  train_idx <- seq_len(floor(n * 0.8))
  test_idx <- setdiff(seq_len(n), train_idx)

  train_X <- X[train_idx, , drop = FALSE]
  train_y <- y[train_idx, , drop = FALSE]
  test_X <- X[test_idx, , drop = FALSE]
  test_y <- y[test_idx, , drop = FALSE]

  # Fit model
  model <- switch(method,
    ridge = .fit_glmnet_decoder(train_X, train_y, alpha = 0, cv_folds = cv_folds, ...),
    lasso = .fit_glmnet_decoder(train_X, train_y, alpha = 1, cv_folds = cv_folds, ...),
    elastic_net = .fit_glmnet_decoder(train_X, train_y, alpha = 0.5, cv_folds = cv_folds, ...),
    ranger = .fit_ranger_decoder(train_X, train_y, ...),
    xgboost = .fit_xgboost_decoder(train_X, train_y, ...),
    stop("Unknown method: ", method)
  )

  # Predictions
  train_pred <- predict_decoder(model, train_X, method)
  test_pred <- predict_decoder(model, test_X, method)

  # Metrics
  train_r2 <- .compute_r2(train_y, train_pred)
  test_r2 <- .compute_r2(test_y, test_pred)
  test_corr <- cor(as.vector(test_y), as.vector(test_pred))

  # Get weights for linear models
  weights <- .get_decoder_weights(model, method, n_cells, time_lags)

  list(
    model = model,
    method = method,
    train_predictions = train_pred,
    test_predictions = test_pred,
    train_r2 = train_r2,
    test_r2 = test_r2,
    test_correlation = test_corr,
    weights = weights,
    time_lags = time_lags,
    n_cells = n_cells,
    train_indices = train_idx,
    test_indices = test_idx
  )
}

#' Predict behavior using trained decoder
#'
#' @param model Trained decoder model
#' @param new_traces New trace data (time x cells)
#' @param method Method used in training
#'
#' @return Predicted behavioral values
#' @export
predict_decoder <- function(model, new_traces, method = "ridge") {
  switch(method,
    ridge = , lasso = , elastic_net = {
      predict(model$model, as.matrix(new_traces), s = model$lambda)
    },
    ranger = {
      predict(model, data = as.data.frame(new_traces))$predictions
    },
    xgboost = {
      predict(model, as.matrix(new_traces))
    }
  )
}

.fit_glmnet_decoder <- function(X, y, alpha, cv_folds, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet package required")
  }

  cv_fit <- glmnet::cv.glmnet(
    x = as.matrix(X),
    y = as.matrix(y),
    alpha = alpha,
    nfolds = cv_folds,
    ...
  )

  list(model = cv_fit, lambda = cv_fit$lambda.min, alpha = alpha)
}

.fit_ranger_decoder <- function(X, y, num.trees = 500, ...) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("ranger package required")
  }

  df <- as.data.frame(X)
  df$.target <- y[, 1]

  ranger::ranger(
    .target ~ .,
    data = df,
    num.trees = num.trees,
    ...
  )
}

.fit_xgboost_decoder <- function(X, y, nrounds = 100, ...) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("xgboost package required")
  }

  xgboost::xgboost(
    data = as.matrix(X),
    label = y[, 1],
    nrounds = nrounds,
    objective = "reg:squarederror",
    verbose = 0,
    ...
  )
}

.compute_r2 <- function(actual, predicted) {
  ss_res <- sum((actual - predicted)^2)
  ss_tot <- sum((actual - mean(actual))^2)
  1 - ss_res / ss_tot
}

.get_decoder_weights <- function(model, method, n_cells, time_lags) {
  if (method %in% c("ridge", "lasso", "elastic_net")) {
    coefs <- as.matrix(coef(model$model, s = model$lambda))[-1, , drop = FALSE]

    if (length(time_lags) > 1 || time_lags[1] != 0) {
      # Reshape to cells x lags
      matrix(coefs, nrow = n_cells, ncol = length(time_lags))
    } else {
      coefs
    }
  } else {
    NULL
  }
}

#' Create lagged features for time series
#'
#' @param X Matrix (time x features)
#' @param lags Integer vector of lags
#'
#' @return Matrix with lagged features
#' @export
create_lagged_features <- function(X, lags) {
  n_time <- nrow(X)
  n_features <- ncol(X)
  n_lags <- length(lags)
  max_lag <- max(abs(lags))

  valid_range <- (max_lag + 1):(n_time - max_lag)
  n_valid <- length(valid_range)

  result <- matrix(0, nrow = n_valid, ncol = n_features * n_lags)

  for (i in seq_along(lags)) {
    lag <- lags[i]
    col_idx <- ((i - 1) * n_features + 1):(i * n_features)
    result[, col_idx] <- X[valid_range - lag, ]
  }

  colnames(result) <- paste0(
    rep(colnames(X) %||% paste0("V", seq_len(n_features)), n_lags),
    "_lag",
    rep(lags, each = n_features)
  )

  result
}

# ============================================================================
# Event/Spike Classification
# ============================================================================

#' Train event classifier
#'
#' Train a classifier to detect calcium events or spikes.
#'
#' @param traces Matrix of traces (cells x time)
#' @param events Binary matrix of known events (same dims as traces)
#' @param method Classification method: "ranger", "xgboost", "logistic"
#' @param window_size Window size for feature extraction around each time point
#' @param balance_classes Balance positive/negative examples
#' @param ... Additional arguments to classifier
#'
#' @return Trained event classifier
#' @export
#'
#' @examples
#' \dontrun{
#' # Train on manually annotated events
#' classifier <- train_event_classifier(traces, manual_events, method = "xgboost")
#'
#' # Detect events in new data
#' detected <- detect_events_ml(new_traces, classifier)
#' }
train_event_classifier <- function(traces, events, method = "ranger",
                                   window_size = 5, balance_classes = TRUE, ...) {
  # Extract training examples
  examples <- .extract_event_examples(traces, events, window_size)

  X <- examples$features
  y <- factor(examples$labels, levels = c(0, 1), labels = c("no_event", "event"))

  # Balance classes if requested
  if (balance_classes) {
    balanced <- .balance_classes(X, y)
    X <- balanced$X
    y <- balanced$y
  }

  # Train model
  model <- switch(method,
    ranger = .train_ranger_classifier(X, y, ...),
    xgboost = .train_xgboost_classifier(X, y, ...),
    logistic = {
      df <- as.data.frame(X)
      df$.event <- y
      glm(.event ~ ., data = df, family = binomial())
    },
    stop("Unknown method: ", method)
  )

  structure(
    list(
      model = model,
      method = method,
      window_size = window_size,
      feature_names = colnames(X)
    ),
    class = "event_classifier"
  )
}

#' Detect events using trained classifier
#'
#' @param traces Matrix of traces (cells x time)
#' @param classifier Trained event classifier
#' @param threshold Probability threshold for event detection
#'
#' @return Binary matrix of detected events
#' @export
detect_events_ml <- function(traces, classifier, threshold = 0.5) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)
  window_size <- classifier$window_size
  half_w <- window_size %/% 2

  events <- matrix(0, nrow = n_cells, ncol = n_time)

  for (cell in seq_len(n_cells)) {
    trace <- traces[cell, ]

    for (t in (half_w + 1):(n_time - half_w)) {
      features <- .extract_window_features(trace, t, window_size)
      features_df <- as.data.frame(matrix(features, nrow = 1))
      names(features_df) <- classifier$feature_names

      prob <- switch(classifier$method,
        ranger = predict(classifier$model, data = features_df)$predictions,
        xgboost = {
          pred <- predict(classifier$model, as.matrix(features_df))
          if (length(classifier$model$classes) == 2) pred else pred
        },
        logistic = predict(classifier$model, features_df, type = "response")
      )

      if (classifier$method == "ranger") {
        events[cell, t] <- as.integer(prob == "event")
      } else {
        events[cell, t] <- as.integer(prob > threshold)
      }
    }
  }

  events
}

.extract_event_examples <- function(traces, events, window_size) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)
  half_w <- window_size %/% 2

  features_list <- list()
  labels_list <- list()

  for (cell in seq_len(n_cells)) {
    trace <- traces[cell, ]
    event_vec <- events[cell, ]

    for (t in (half_w + 1):(n_time - half_w)) {
      features <- .extract_window_features(trace, t, window_size)
      features_list <- c(features_list, list(features))
      labels_list <- c(labels_list, event_vec[t])
    }
  }

  list(
    features = do.call(rbind, features_list),
    labels = unlist(labels_list)
  )
}

.extract_window_features <- function(trace, center, window_size) {
  half_w <- window_size %/% 2
  window <- trace[(center - half_w):(center + half_w)]

  c(
    window,  # Raw values
    mean(window),
    sd(window),
    max(window) - min(window),
    window[half_w + 1] - mean(window[1:half_w]),  # Rise from baseline
    diff(window)  # First differences
  )
}

.balance_classes <- function(X, y) {
  pos_idx <- which(y == levels(y)[2])
  neg_idx <- which(y == levels(y)[1])

  n_pos <- length(pos_idx)
  n_neg <- length(neg_idx)

  if (n_pos < n_neg) {
    neg_idx <- sample(neg_idx, n_pos)
  } else {
    pos_idx <- sample(pos_idx, n_neg, replace = TRUE)
  }

  idx <- c(pos_idx, neg_idx)
  list(X = X[idx, , drop = FALSE], y = y[idx])
}

# ============================================================================
# Feature Extraction
# ============================================================================

#' Extract features from calcium traces
#'
#' Compute features for machine learning from raw traces.
#'
#' @param traces Matrix of traces (cells x time)
#' @param feature_types Types of features: "temporal", "spectral", "statistical", "all"
#' @param normalize Normalize features
#'
#' @return Feature matrix (cells x features)
#' @export
#'
#' @examples
#' \dontrun{
#' features <- extract_trace_features(traces, "all")
#' }
extract_trace_features <- function(traces, feature_types = "all",
                                   normalize = TRUE) {
  n_cells <- nrow(traces)

  if ("all" %in% feature_types) {
    feature_types <- c("temporal", "spectral", "statistical")
  }

  features_list <- list()

  if ("statistical" %in% feature_types) {
    features_list$statistical <- .extract_statistical_features(traces)
  }

  if ("temporal" %in% feature_types) {
    features_list$temporal <- .extract_temporal_features(traces)
  }

  if ("spectral" %in% feature_types) {
    features_list$spectral <- .extract_spectral_features(traces)
  }

  features <- do.call(cbind, features_list)

  if (normalize) {
    features <- scale(features)
    features[is.na(features)] <- 0
  }

  features
}

.extract_statistical_features <- function(traces) {
  n_cells <- nrow(traces)

  data.frame(
    mean = rowMeans(traces),
    sd = apply(traces, 1, sd),
    skewness = apply(traces, 1, function(x) {
      n <- length(x)
      m <- mean(x)
      s <- sd(x)
      sum((x - m)^3) / (n * s^3)
    }),
    kurtosis = apply(traces, 1, function(x) {
      n <- length(x)
      m <- mean(x)
      s <- sd(x)
      sum((x - m)^4) / (n * s^4) - 3
    }),
    range = apply(traces, 1, function(x) max(x) - min(x)),
    iqr = apply(traces, 1, IQR),
    median = apply(traces, 1, median),
    mad = apply(traces, 1, mad),
    cv = apply(traces, 1, function(x) sd(x) / mean(x)),
    above_mean_frac = apply(traces, 1, function(x) mean(x > mean(x)))
  )
}

.extract_temporal_features <- function(traces) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Autocorrelation at different lags
  acf_lags <- c(1, 5, 10, 20)
  acf_features <- sapply(acf_lags, function(lag) {
    apply(traces, 1, function(x) {
      if (lag >= length(x)) return(0)
      acf(x, lag.max = lag, plot = FALSE)$acf[lag + 1]
    })
  })
  colnames(acf_features) <- paste0("acf_lag", acf_lags)

  # Derivative statistics
  diffs <- t(apply(traces, 1, diff))
  diff_features <- data.frame(
    diff_mean = rowMeans(diffs),
    diff_sd = apply(diffs, 1, sd),
    diff_max = apply(diffs, 1, max),
    diff_min = apply(diffs, 1, min),
    zero_crossings = apply(diffs, 1, function(x) sum(diff(sign(x)) != 0))
  )

  # Peak statistics
  peak_features <- t(apply(traces, 1, function(x) {
    threshold <- mean(x) + 2 * sd(x)
    peaks <- which(x > threshold & c(FALSE, diff(x) > 0) & c(diff(x) < 0, FALSE))
    c(
      n_peaks = length(peaks),
      peak_rate = length(peaks) / length(x),
      mean_peak_height = if (length(peaks) > 0) mean(x[peaks]) else 0
    )
  }))

  cbind(acf_features, diff_features, peak_features)
}

.extract_spectral_features <- function(traces) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Simple spectral features using FFT
  t(apply(traces, 1, function(x) {
    # FFT
    fft_result <- Mod(fft(x - mean(x)))
    freqs <- fft_result[1:(length(x) %/% 2)]
    freqs <- freqs / sum(freqs)

    # Spectral features
    n_freq <- length(freqs)
    freq_idx <- seq_len(n_freq)

    c(
      spectral_centroid = sum(freq_idx * freqs) / sum(freqs),
      spectral_spread = sqrt(sum((freq_idx - sum(freq_idx * freqs) / sum(freqs))^2 * freqs) / sum(freqs)),
      spectral_entropy = -sum(freqs * log(freqs + 1e-10)),
      low_freq_power = sum(freqs[1:min(10, n_freq)]),
      mid_freq_power = sum(freqs[min(11, n_freq):min(50, n_freq)]),
      high_freq_power = sum(freqs[min(51, n_freq):n_freq]),
      dominant_freq = which.max(freqs)
    )
  }))
}

# ============================================================================
# Ensemble Methods
# ============================================================================

#' Ensemble spike detection
#'
#' Combine multiple spike detection methods using ML.
#'
#' @param traces Matrix of traces (cells x time)
#' @param methods Vector of detection methods to combine
#' @param ground_truth Optional binary matrix of known spikes for training
#' @param ensemble_method How to combine: "vote", "stack_ranger", "stack_xgboost"
#' @param ... Additional arguments
#'
#' @return List with ensemble predictions and individual method results
#' @export
#'
#' @examples
#' \dontrun{
#' result <- ensemble_spike_detection(
#'   traces,
#'   methods = c("threshold", "adaptive", "derivative"),
#'   ensemble_method = "vote"
#' )
#' }
ensemble_spike_detection <- function(traces, methods = c("threshold", "adaptive"),
                                      ground_truth = NULL,
                                      ensemble_method = "vote", ...) {
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Run individual methods
  method_results <- lapply(methods, function(m) {
    if (m == "threshold") {
      threshold_spike_detection(traces, ...)$spike_predictions
    } else if (m == "adaptive") {
      adaptive_threshold_detection(traces, ...)$spike_predictions
    } else if (m == "derivative") {
      .derivative_spike_detection(traces)
    } else {
      stop("Unknown method: ", m)
    }
  })
  names(method_results) <- methods

  # Combine results
  if (ensemble_method == "vote") {
    # Simple majority voting
    combined <- Reduce(`+`, method_results)
    ensemble_pred <- (combined > length(methods) / 2) * 1
  } else if (ensemble_method %in% c("stack_ranger", "stack_xgboost")) {
    if (is.null(ground_truth)) {
      stop("ground_truth required for stacking ensemble")
    }

    # Create stacking features
    stack_features <- do.call(cbind, lapply(method_results, as.vector))
    colnames(stack_features) <- methods
    stack_labels <- as.vector(ground_truth)

    # Train stacking model
    stack_method <- sub("stack_", "", ensemble_method)
    stack_model <- switch(stack_method,
      ranger = .train_ranger_classifier(stack_features, factor(stack_labels), ...),
      xgboost = .train_xgboost_classifier(stack_features, factor(stack_labels), ...)
    )

    # Predict
    ensemble_pred <- matrix(
      as.integer(predict_cell_types(stack_model, stack_features, stack_method) == "1"),
      nrow = n_cells,
      ncol = n_time
    )
  }

  list(
    ensemble_predictions = ensemble_pred,
    individual_predictions = method_results,
    method = ensemble_method,
    n_methods = length(methods)
  )
}

.derivative_spike_detection <- function(traces, threshold_sd = 2) {
  # Simple derivative-based detection
  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  spikes <- matrix(0, n_cells, n_time)

  for (i in seq_len(n_cells)) {
    trace <- traces[i, ]
    deriv <- c(0, diff(trace))
    threshold <- mean(deriv) + threshold_sd * sd(deriv)
    spikes[i, ] <- as.integer(deriv > threshold)
  }

  spikes
}

# ============================================================================
# Cross-Validation Utilities
# ============================================================================

#' Cross-validate decoder performance
#'
#' @param traces Matrix of traces
#' @param behavior Behavioral variable to predict
#' @param method Decoding method
#' @param n_folds Number of CV folds
#' @param ... Additional arguments
#'
#' @return Cross-validation results
#' @export
cv_decode_behavior <- function(traces, behavior, method = "ridge",
                               n_folds = 5, ...) {
  n <- ncol(traces)
  fold_size <- floor(n / n_folds)
  fold_ids <- rep(1:n_folds, each = fold_size)[1:n]

  results <- lapply(1:n_folds, function(fold) {
    test_idx <- which(fold_ids == fold)
    train_idx <- which(fold_ids != fold)

    train_traces <- traces[, train_idx]
    test_traces <- traces[, test_idx]
    train_behavior <- behavior[train_idx]
    test_behavior <- behavior[test_idx]

    model <- decode_behavior(
      train_traces,
      matrix(train_behavior, ncol = 1),
      method = method,
      transpose = FALSE,
      ...
    )

    test_pred <- predict_decoder(model$model, t(test_traces), method)

    list(
      fold = fold,
      r2 = .compute_r2(test_behavior, test_pred),
      correlation = cor(test_behavior, test_pred)
    )
  })

  r2_values <- sapply(results, `[[`, "r2")
  corr_values <- sapply(results, `[[`, "correlation")

  list(
    mean_r2 = mean(r2_values),
    sd_r2 = sd(r2_values),
    mean_correlation = mean(corr_values),
    sd_correlation = sd(corr_values),
    fold_results = results
  )
}

# Null coalesce
`%||%` <- function(x, y) if (is.null(x)) y else x
