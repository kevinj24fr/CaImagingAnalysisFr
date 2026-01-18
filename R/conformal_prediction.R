#' Conformal Prediction for Neural Data
#'
#' Distribution-free uncertainty quantification for neural predictions.
#' Provides prediction intervals with guaranteed coverage.
#'
#' @name conformal_prediction
NULL

# ============================================================================
# Conformal Prediction Framework
# ============================================================================

#' Fit Conformal Predictor
#'
#' Create a conformal predictor with distribution-free coverage guarantees.
#'
#' @param model Fitted model (any model with predict method)
#' @param X_calib Calibration features (n_calib x features)
#' @param y_calib Calibration responses
#' @param method Conformal method ("split", "full", "cv")
#' @param alpha Significance level (1 - alpha = coverage)
#' @param score_fn Custom score function (optional)
#'
#' @return Conformal predictor object
#' @export
#'
#' @examples
#' \dontrun{
#' # Train a decoder, then calibrate for uncertainty
#' decoder <- decode_behavior(neural_features, behavior)
#' conf <- conformal_predictor(decoder, X_calib, y_calib, alpha = 0.1)
#'
#' # Get prediction intervals
#' pred <- conformal_predict(conf, X_test)
#' # pred$lower and pred$upper have 90% coverage guarantee
#' }
conformal_predictor <- function(model,
                                 X_calib,
                                 y_calib,
                                 method = c("split", "full", "cv"),
                                 alpha = 0.1,
                                 score_fn = NULL) {

  method <- match.arg(method)

  if (is.vector(X_calib)) {
    X_calib <- matrix(X_calib, ncol = 1)
  }

  n_calib <- nrow(X_calib)

  # Default score function: absolute residual
  if (is.null(score_fn)) {
    score_fn <- function(y_true, y_pred) {
      abs(y_true - y_pred)
    }
  }

  if (method == "split") {
    # Split conformal: compute scores on calibration set
    y_pred <- predict(model, newdata = as.data.frame(X_calib))
    if (is.list(y_pred)) y_pred <- y_pred$predictions %||% y_pred[[1]]

    scores <- score_fn(y_calib, y_pred)

    # Quantile for coverage
    q_level <- ceiling((n_calib + 1) * (1 - alpha)) / n_calib
    q_hat <- quantile(scores, probs = min(q_level, 1))

    result <- list(
      model = model,
      scores = scores,
      q_hat = q_hat,
      alpha = alpha,
      n_calib = n_calib,
      method = "split",
      score_fn = score_fn
    )

  } else if (method == "cv") {
    # Cross-validation conformal
    n_folds <- 5
    folds <- sample(rep(1:n_folds, length.out = n_calib))
    scores <- numeric(n_calib)

    for (k in 1:n_folds) {
      train_idx <- folds != k
      test_idx <- folds == k

      # Refit model on fold
      model_k <- refit_model(model, X_calib[train_idx, ], y_calib[train_idx])

      # Predict on held-out
      y_pred_k <- predict(model_k, newdata = as.data.frame(X_calib[test_idx, , drop = FALSE]))
      if (is.list(y_pred_k)) y_pred_k <- y_pred_k$predictions %||% y_pred_k[[1]]

      scores[test_idx] <- score_fn(y_calib[test_idx], y_pred_k)
    }

    q_level <- ceiling((n_calib + 1) * (1 - alpha)) / n_calib
    q_hat <- quantile(scores, probs = min(q_level, 1))

    result <- list(
      model = model,
      scores = scores,
      q_hat = q_hat,
      alpha = alpha,
      n_calib = n_calib,
      method = "cv",
      score_fn = score_fn
    )

  } else {
    # Full conformal (computationally expensive)
    result <- list(
      model = model,
      X_calib = X_calib,
      y_calib = y_calib,
      alpha = alpha,
      n_calib = n_calib,
      method = "full",
      score_fn = score_fn
    )
  }

  structure(result, class = "conformal_predictor")
}

#' Refit Model on Subset
#' @keywords internal
refit_model <- function(model, X, y) {
  # Try common model types
  if (inherits(model, "lm")) {
    data <- as.data.frame(X)
    data$y <- y
    lm(y ~ ., data = data)
  } else if (inherits(model, "glm")) {
    data <- as.data.frame(X)
    data$y <- y
    glm(y ~ ., data = data, family = model$family)
  } else if (inherits(model, "ranger")) {
    if (requireNamespace("ranger", quietly = TRUE)) {
      data <- as.data.frame(X)
      data$y <- y
      ranger::ranger(y ~ ., data = data)
    } else {
      model
    }
  } else {
    # Return original model if can't refit
    model
  }
}

#' Conformal Prediction Intervals
#'
#' @param object Conformal predictor
#' @param newdata New data for prediction
#' @param ... Additional arguments
#'
#' @return List with predictions, lower bounds, and upper bounds
#' @export
conformal_predict <- function(object, newdata, ...) {
  if (!inherits(object, "conformal_predictor")) {
    stop("object must be a conformal_predictor")
  }

  if (is.vector(newdata)) {
    newdata <- matrix(newdata, ncol = 1)
  }

  if (object$method %in% c("split", "cv")) {
    # Point predictions
    y_pred <- predict(object$model, newdata = as.data.frame(newdata))
    if (is.list(y_pred)) y_pred <- y_pred$predictions %||% y_pred[[1]]

    # Prediction intervals using calibrated quantile
    lower <- y_pred - object$q_hat
    upper <- y_pred + object$q_hat

    list(
      pred = y_pred,
      lower = lower,
      upper = upper,
      coverage = 1 - object$alpha,
      width = upper - lower
    )

  } else {
    # Full conformal (much slower)
    n_new <- nrow(newdata)
    y_pred <- numeric(n_new)
    lower <- numeric(n_new)
    upper <- numeric(n_new)

    for (i in seq_len(n_new)) {
      # Grid search for interval
      y_range <- range(object$y_calib)
      y_grid <- seq(y_range[1] - diff(y_range), y_range[2] + diff(y_range), length.out = 100)

      in_interval <- sapply(y_grid, function(y_trial) {
        # Augment calibration set
        X_aug <- rbind(object$X_calib, newdata[i, ])
        y_aug <- c(object$y_calib, y_trial)

        # Compute all scores (would need to refit model for full conformal)
        # Simplified: use original model
        y_pred_aug <- predict(object$model, newdata = as.data.frame(X_aug))
        if (is.list(y_pred_aug)) y_pred_aug <- y_pred_aug$predictions %||% y_pred_aug[[1]]

        scores <- object$score_fn(y_aug, y_pred_aug)

        # Compute p-value
        p_value <- mean(scores >= scores[length(scores)])
        p_value >= object$alpha
      })

      valid_y <- y_grid[in_interval]
      if (length(valid_y) > 0) {
        lower[i] <- min(valid_y)
        upper[i] <- max(valid_y)
        y_pred[i] <- median(valid_y)
      } else {
        y_pred[i] <- predict(object$model, newdata = as.data.frame(newdata[i, , drop = FALSE]))[[1]]
        lower[i] <- y_pred[i] - sd(object$y_calib)
        upper[i] <- y_pred[i] + sd(object$y_calib)
      }
    }

    list(
      pred = y_pred,
      lower = lower,
      upper = upper,
      coverage = 1 - object$alpha
    )
  }
}

#' @export
print.conformal_predictor <- function(x, ...) {
  cat("Conformal Predictor\n")
  cat("===================\n")
  cat(sprintf("Method: %s\n", x$method))
  cat(sprintf("Calibration samples: %d\n", x$n_calib))
  cat(sprintf("Target coverage: %.1f%%\n", 100 * (1 - x$alpha)))
  if (!is.null(x$q_hat)) {
    cat(sprintf("Calibrated quantile: %.4f\n", x$q_hat))
  }
  invisible(x)
}

# ============================================================================
# Conformal Methods for Specific Tasks
# ============================================================================

#' Conformal Regression for Neural Decoding
#'
#' Add prediction intervals to behavioral decoding.
#'
#' @param traces Neural traces (cells x time)
#' @param behavior Behavior variable to decode
#' @param train_idx Training indices
#' @param calib_idx Calibration indices
#' @param test_idx Test indices
#' @param alpha Significance level
#' @param method Decoder method
#'
#' @return Decoding results with conformal intervals
#' @export
conformal_decode <- function(traces, behavior, train_idx, calib_idx, test_idx,
                              alpha = 0.1, method = "ridge") {

  X <- t(traces)
  y <- behavior

  # Train decoder
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]

  decoder <- train_decoder(X_train, y_train, method)

  # Calibrate
  X_calib <- X[calib_idx, ]
  y_calib <- y[calib_idx]

  conf <- conformal_predictor(decoder, X_calib, y_calib, alpha = alpha)

  # Predict with intervals
  X_test <- X[test_idx, ]
  y_test <- y[test_idx]

  pred <- conformal_predict(conf, X_test)

  # Evaluate coverage
  covered <- y_test >= pred$lower & y_test <= pred$upper
  actual_coverage <- mean(covered)

  list(
    predictions = pred$pred,
    lower = pred$lower,
    upper = pred$upper,
    true = y_test,
    target_coverage = 1 - alpha,
    actual_coverage = actual_coverage,
    mean_width = mean(pred$width),
    mse = mean((pred$pred - y_test)^2)
  )
}

#' Train Simple Decoder
#' @keywords internal
train_decoder <- function(X, y, method) {
  data <- as.data.frame(X)
  data$y <- y

  if (method == "ridge") {
    # Ridge regression
    if (requireNamespace("glmnet", quietly = TRUE)) {
      cv_fit <- glmnet::cv.glmnet(as.matrix(X), y, alpha = 0)
      return(cv_fit)
    }
  }

  # Fallback to OLS
  lm(y ~ ., data = data)
}

#' Conformal Classification
#'
#' Prediction sets for cell type classification with coverage guarantees.
#'
#' @param model Classification model
#' @param X_calib Calibration features
#' @param y_calib Calibration labels
#' @param alpha Significance level
#'
#' @return Conformal classifier
#' @export
conformal_classifier <- function(model, X_calib, y_calib, alpha = 0.1) {

  if (is.vector(X_calib)) {
    X_calib <- matrix(X_calib, ncol = 1)
  }

  n_calib <- nrow(X_calib)
  classes <- unique(y_calib)

  # Get predicted probabilities
  probs <- predict(model, newdata = as.data.frame(X_calib), type = "prob")
  if (!is.matrix(probs) && !is.data.frame(probs)) {
    # Model doesn't give probabilities, use predicted class
    probs <- predict(model, newdata = as.data.frame(X_calib))
    probs <- sapply(classes, function(c) as.numeric(probs == c))
  }

  # Compute conformity scores (1 - probability of true class)
  scores <- numeric(n_calib)
  for (i in seq_len(n_calib)) {
    true_class <- as.character(y_calib[i])
    if (true_class %in% colnames(probs)) {
      scores[i] <- 1 - probs[i, true_class]
    } else {
      scores[i] <- 1
    }
  }

  # Calibrated quantile
  q_level <- ceiling((n_calib + 1) * (1 - alpha)) / n_calib
  q_hat <- quantile(scores, probs = min(q_level, 1))

  structure(
    list(
      model = model,
      scores = scores,
      q_hat = q_hat,
      alpha = alpha,
      classes = classes,
      n_calib = n_calib
    ),
    class = "conformal_classifier"
  )
}

#' Conformal Classification Prediction Sets
#'
#' @param object Conformal classifier
#' @param newdata New data
#' @param ... Additional arguments
#'
#' @return List of prediction sets (one per observation)
#' @export
conformal_predict_set <- function(object, newdata, ...) {
  if (!inherits(object, "conformal_classifier")) {
    stop("object must be a conformal_classifier")
  }

  if (is.vector(newdata)) {
    newdata <- matrix(newdata, ncol = 1)
  }

  n_new <- nrow(newdata)

  # Get probabilities
  probs <- predict(object$model, newdata = as.data.frame(newdata), type = "prob")
  if (!is.matrix(probs) && !is.data.frame(probs)) {
    probs <- predict(object$model, newdata = as.data.frame(newdata))
    probs <- sapply(object$classes, function(c) as.numeric(probs == c))
  }

  # Build prediction sets
  prediction_sets <- lapply(seq_len(n_new), function(i) {
    # Include classes where 1 - prob <= q_hat
    included <- colnames(probs)[1 - probs[i, ] <= object$q_hat]
    if (length(included) == 0) {
      # Return most likely class if empty
      included <- colnames(probs)[which.max(probs[i, ])]
    }
    included
  })

  # Summary statistics
  set_sizes <- sapply(prediction_sets, length)

  list(
    sets = prediction_sets,
    set_sizes = set_sizes,
    mean_size = mean(set_sizes),
    coverage = 1 - object$alpha
  )
}

# ============================================================================
# Adaptive Conformal Inference
# ============================================================================

#' Adaptive Conformal Inference for Time Series
#'
#' Conformal prediction that adapts to distribution shift over time.
#'
#' @param model Predictive model
#' @param y_stream Stream of observations
#' @param X_stream Stream of features (optional)
#' @param alpha Target miscoverage rate
#' @param gamma Learning rate for adaptation
#'
#' @return Adaptive conformal results
#' @export
adaptive_conformal <- function(model, y_stream, X_stream = NULL,
                                alpha = 0.1, gamma = 0.01) {

  n <- length(y_stream)

  # Initialize
  alphas <- numeric(n)
  alphas[1] <- alpha

  q_hats <- numeric(n)
  errors <- logical(n)

  lower <- numeric(n)
  upper <- numeric(n)

  for (t in 2:n) {
    # Current prediction
    if (!is.null(X_stream)) {
      y_pred <- predict(model, newdata = as.data.frame(X_stream[t, , drop = FALSE]))[[1]]
    } else {
      y_pred <- predict(model, n.ahead = 1)$pred[1]
    }

    # Current quantile (from past errors)
    past_resid <- abs(y_stream[1:(t-1)] -
                       if (!is.null(X_stream)) {
                         predict(model, newdata = as.data.frame(X_stream[1:(t-1), , drop = FALSE]))
                       } else {
                         fitted(model)[1:(t-1)]
                       })

    q_hats[t] <- quantile(past_resid, 1 - alphas[t-1], na.rm = TRUE)

    lower[t] <- y_pred - q_hats[t]
    upper[t] <- y_pred + q_hats[t]

    # Check coverage
    errors[t] <- y_stream[t] < lower[t] | y_stream[t] > upper[t]

    # Adapt alpha
    alphas[t] <- alphas[t-1] + gamma * (errors[t] - alpha)
    alphas[t] <- max(0.01, min(0.5, alphas[t]))  # Clip to valid range
  }

  list(
    lower = lower,
    upper = upper,
    alphas = alphas,
    errors = errors,
    actual_coverage = 1 - mean(errors[-1]),
    target_coverage = 1 - alpha
  )
}

#' Plot Conformal Predictions
#'
#' @param x Conformal prediction result from conformal_predict
#' @param y_true True values (optional)
#' @param ... Additional arguments
#'
#' @export
plot_conformal <- function(x, y_true = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 required")
  }

  n <- length(x$pred)
  df <- data.frame(
    idx = 1:n,
    pred = x$pred,
    lower = x$lower,
    upper = x$upper
  )

  if (!is.null(y_true)) {
    df$true <- y_true
    df$covered <- y_true >= x$lower & y_true <= x$upper
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = idx)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "steelblue", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = pred), color = "steelblue", linewidth = 1)

  if (!is.null(y_true)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(y = true, color = covered), size = 2) +
      ggplot2::scale_color_manual(values = c("FALSE" = "red", "TRUE" = "black"),
                                   name = "Covered")
  }

  p + ggplot2::labs(title = sprintf("Conformal Prediction (%.0f%% coverage)",
                                     100 * x$coverage),
                    x = "Index", y = "Value") +
    ggplot2::theme_minimal()
}

#' Evaluate Conformal Predictor Performance
#'
#' @param pred Predictions with intervals
#' @param y_true True values
#'
#' @return Coverage and efficiency metrics
#' @export
evaluate_conformal <- function(pred, y_true) {
  covered <- y_true >= pred$lower & y_true <= pred$upper
  width <- pred$upper - pred$lower

  list(
    coverage = mean(covered),
    mean_width = mean(width),
    median_width = median(width),
    mse = mean((pred$pred - y_true)^2),
    mae = mean(abs(pred$pred - y_true)),
    n_uncovered = sum(!covered),
    conditional_coverage = tapply(covered, cut(y_true, 5), mean)
  )
}
