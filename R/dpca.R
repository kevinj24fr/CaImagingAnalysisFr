#' Demixed Principal Component Analysis (dPCA)
#'
#' Decompose neural population activity into components associated with
#' different task parameters. dPCA finds linear combinations of neurons
#' that capture variance explained by specific task variables.
#'
#' @name dpca
#' @references
#' Kobak, D., et al. (2016). Demixed principal component analysis of neural
#' population data. eLife.
NULL

#' Fit dPCA Model
#'
#' Perform demixed PCA to separate neural variance by task parameters.
#'
#' @param data 4D array of neural data: (neurons x time x condition1 x condition2) or
#'        list with $data (neurons x time x trials) and $labels (data frame of trial labels)
#' @param labels If data is 3D array, data frame with trial labels (columns are conditions)
#' @param n_components Number of components per marginalization
#' @param regularization Regularization parameter for decoder (0 = no regularization)
#' @param center Whether to center data
#' @param verbose Print progress
#'
#' @return List containing:
#' \itemize{
#'   \item{components}{List of components per marginalization}
#'   \item{explained_var}{Variance explained by each marginalization}
#'   \item{projections}{Data projected onto each set of components}
#'   \item{marginalizations}{Names of marginalizations}
#'   \item{encoders}{Encoding matrices (neurons -> components)}
#'   \item{decoders}{Decoding matrices (components -> neurons)}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun
#' # Create trial data with stimulus and time conditions
#' # data: neurons x time x trials
#' # labels: data frame with 'stimulus' column
#'
#' dpca_result <- fit_dpca(data, labels, n_components = 5)
#'
#' # Plot stimulus-related components
#' plot_dpca(dpca_result, marginalization = "stimulus")
#' }
fit_dpca <- function(data, labels = NULL, n_components = 10,
                     regularization = 0, center = TRUE, verbose = TRUE) {

  # Parse input format
  parsed <- .dpca_parse_input(data, labels)
  X <- parsed$X  # neurons x time x conditions (averaged)
  trial_data <- parsed$trial_data  # neurons x time x trials (if available)
  condition_labels <- parsed$condition_labels
  marginalizations <- parsed$marginalizations

  n_neurons <- dim(X)[1]
  n_time <- dim(X)[2]
  n_conditions <- prod(dim(X)[-(1:2)])

  if (verbose) {
    cat(sprintf("Fitting dPCA: %d neurons, %d time points, %d conditions\n",
                n_neurons, n_time, n_conditions))
    cat(sprintf("Marginalizations: %s\n", paste(names(marginalizations), collapse = ", ")))
  }

  # Center data
  if (center) {
    X_mean <- apply(X, 1, mean)
    X_centered <- sweep(X, 1, X_mean, "-")
  } else {
    X_mean <- rep(0, n_neurons)
    X_centered <- X
  }

  # Compute total covariance
  X_flat <- matrix(X_centered, n_neurons, n_time * n_conditions)
  C_total <- tcrossprod(X_flat) / (ncol(X_flat) - 1)

  # Compute marginalized covariances
  marg_covs <- lapply(names(marginalizations), function(marg_name) {
    .dpca_marginalize_cov(X_centered, marginalizations[[marg_name]], marg_name)
  })
  names(marg_covs) <- names(marginalizations)

  # Compute dPCA axes for each marginalization
  results <- list()

  for (marg_name in names(marginalizations)) {
    if (verbose) cat(sprintf("  Computing %s components...\n", marg_name))

    C_marg <- marg_covs[[marg_name]]

    # Regularized decoder: D = (C_total + lambda * I)^-1 * C_marg
    if (regularization > 0) {
      C_reg <- C_total + regularization * diag(n_neurons)
    } else {
      C_reg <- C_total + 1e-6 * diag(n_neurons)  # Small ridge for stability
    }

    # Solve for decoder
    M <- solve(C_reg, C_marg)

    # SVD of the product to get components
    svd_result <- svd(M, nu = n_components, nv = n_components)

    # Decoder (F): maps from neural space to component space
    decoder <- svd_result$u[, 1:min(n_components, ncol(svd_result$u)), drop = FALSE]

    # Encoder (D): maps from component space to neural space
    encoder <- svd_result$v[, 1:min(n_components, ncol(svd_result$v)), drop = FALSE]

    # Variance explained
    var_explained <- svd_result$d[1:min(n_components, length(svd_result$d))]^2
    var_explained <- var_explained / sum(diag(C_marg))

    # Project data
    projections <- array(0, dim = c(n_components, dim(X)[2:length(dim(X))]))
    for (t in seq_len(n_time)) {
      X_t <- matrix(X_centered[, t, ], n_neurons)
      projections[, t, ] <- t(decoder) %*% X_t
    }

    results[[marg_name]] <- list(
      decoder = decoder,
      encoder = encoder,
      var_explained = var_explained,
      projections = projections
    )
  }

  # Compile output
  list(
    components = lapply(results, function(r) r$decoder),
    encoders = lapply(results, function(r) r$encoder),
    decoders = lapply(results, function(r) r$decoder),
    explained_var = lapply(results, function(r) r$var_explained),
    projections = lapply(results, function(r) r$projections),
    marginalizations = marginalizations,
    marg_names = names(marginalizations),
    X_mean = X_mean,
    n_components = n_components,
    condition_labels = condition_labels
  )
}

#' Parse dPCA input
#' @keywords internal
.dpca_parse_input <- function(data, labels) {
  if (is.list(data) && !is.null(data$data)) {
    # Input is list with data and labels
    X_trials <- data$data
    labels <- data$labels
  } else if (is.array(data) && length(dim(data)) == 4) {
    # Already in neurons x time x cond1 x cond2 format
    return(list(
      X = data,
      trial_data = NULL,
      condition_labels = NULL,
      marginalizations = list(
        time = list(dims = 2, type = "time"),
        condition1 = list(dims = 3, type = "condition"),
        condition2 = list(dims = 4, type = "condition"),
        interaction = list(dims = c(3, 4), type = "interaction")
      )
    ))
  } else if (is.array(data) && length(dim(data)) == 3 && !is.null(labels)) {
    X_trials <- data
  } else {
    stop("data must be 4D array or 3D array with labels")
  }

  # Convert trial data to condition-averaged data
  n_neurons <- dim(X_trials)[1]
  n_time <- dim(X_trials)[2]
  n_trials <- dim(X_trials)[3]

  # Get unique conditions
  if (is.data.frame(labels)) {
    condition_cols <- names(labels)
  } else {
    labels <- data.frame(condition = labels)
    condition_cols <- "condition"
  }

  # Create condition indices
  labels$trial_idx <- seq_len(n_trials)
  unique_conditions <- unique(labels[, condition_cols, drop = FALSE])
  n_conds <- nrow(unique_conditions)

  # Average over trials within conditions
  X_avg <- array(0, dim = c(n_neurons, n_time, n_conds))
  for (i in seq_len(n_conds)) {
    cond <- unique_conditions[i, , drop = FALSE]
    mask <- apply(labels[, condition_cols, drop = FALSE], 1, function(x) {
      all(x == cond)
    })
    trial_idx <- which(mask)
    if (length(trial_idx) > 0) {
      X_avg[, , i] <- apply(X_trials[, , trial_idx, drop = FALSE], c(1, 2), mean)
    }
  }

  # Build marginalization structure
  marginalizations <- list(
    time = list(dims = 2, type = "time")
  )

  for (col in condition_cols) {
    marginalizations[[col]] <- list(dims = 3, type = "condition", var = col)
  }

  # Add interaction term if multiple conditions
  if (length(condition_cols) > 1) {
    marginalizations$interaction <- list(dims = 3, type = "interaction")
  }

  list(
    X = X_avg,
    trial_data = X_trials,
    condition_labels = unique_conditions,
    marginalizations = marginalizations
  )
}

#' Compute marginalized covariance
#' @keywords internal
.dpca_marginalize_cov <- function(X, marg_info, marg_name) {
  n_neurons <- dim(X)[1]
  n_time <- dim(X)[2]

  if (marg_info$type == "time") {
    # Time marginalization: average over conditions, covariance over time
    X_time <- apply(X, c(1, 2), mean)
    C <- tcrossprod(X_time) / (n_time - 1)

  } else if (marg_info$type == "condition") {
    # Condition marginalization: average over time, covariance over conditions
    X_cond <- apply(X, c(1, 3), mean)
    n_cond <- ncol(X_cond)
    X_cond_centered <- X_cond - rowMeans(X_cond)
    C <- tcrossprod(X_cond_centered) / (n_cond - 1)

  } else if (marg_info$type == "interaction") {
    # Interaction: residual after removing time and condition effects
    X_time <- apply(X, c(1, 2), mean)
    X_cond <- apply(X, c(1, 3), mean)
    X_grand <- rowMeans(X)

    # Compute interaction residuals
    n_cond <- dim(X)[3]
    X_interaction <- array(0, dim(X))
    for (t in seq_len(n_time)) {
      for (c in seq_len(n_cond)) {
        X_interaction[, t, c] <- X[, t, c] - X_time[, t] - X_cond[, c] + X_grand
      }
    }

    X_flat <- matrix(X_interaction, n_neurons)
    C <- tcrossprod(X_flat) / (ncol(X_flat) - 1)

  } else {
    stop("Unknown marginalization type")
  }

  # Ensure symmetric
  C <- (C + t(C)) / 2
  C
}

#' Project New Data Through dPCA Model
#'
#' Project new trial data onto dPCA components.
#'
#' @param dpca_model Fitted dPCA model from fit_dpca()
#' @param new_data New data array (neurons x time x trials)
#'
#' @return List of projections per marginalization
#' @export
dpca_project <- function(dpca_model, new_data) {
  n_neurons <- dim(new_data)[1]
  n_time <- dim(new_data)[2]
  n_trials <- dim(new_data)[3]

  # Center with original mean
  new_centered <- sweep(new_data, 1, dpca_model$X_mean, "-")

  projections <- list()

  for (marg_name in dpca_model$marg_names) {
    decoder <- dpca_model$decoders[[marg_name]]
    n_comp <- ncol(decoder)

    proj <- array(0, dim = c(n_comp, n_time, n_trials))
    for (trial in seq_len(n_trials)) {
      proj[, , trial] <- t(decoder) %*% new_centered[, , trial]
    }
    projections[[marg_name]] <- proj
  }

  projections
}

#' Compute dPCA Decoding Accuracy
#'
#' Assess how well conditions can be decoded from each marginalization.
#'
#' @param dpca_model Fitted dPCA model
#' @param trial_data 3D array (neurons x time x trials)
#' @param labels Trial condition labels
#' @param n_shuffles Number of label shuffles for significance
#'
#' @return Data frame with decoding accuracy per marginalization and time
#' @export
dpca_decode <- function(dpca_model, trial_data, labels, n_shuffles = 100) {
  n_time <- dim(trial_data)[2]
  n_trials <- dim(trial_data)[3]

  # Project data
  projections <- dpca_project(dpca_model, trial_data)

  # Decode from each marginalization
  results <- list()

  for (marg_name in dpca_model$marg_names) {
    proj <- projections[[marg_name]]

    acc_by_time <- numeric(n_time)
    shuffle_acc <- matrix(0, n_time, n_shuffles)

    for (t in seq_len(n_time)) {
      X <- t(proj[, t, ])  # trials x components

      # Simple LDA-style decoder
      acc_by_time[t] <- .dpca_lda_cv(X, labels)

      # Shuffle test
      for (s in seq_len(n_shuffles)) {
        shuffle_acc[t, s] <- .dpca_lda_cv(X, sample(labels))
      }
    }

    results[[marg_name]] <- data.frame(
      marginalization = marg_name,
      time = seq_len(n_time),
      accuracy = acc_by_time,
      shuffle_mean = rowMeans(shuffle_acc),
      shuffle_sd = apply(shuffle_acc, 1, sd),
      p_value = rowMeans(shuffle_acc >= acc_by_time)
    )
  }

  do.call(rbind, results)
}

#' Simple LDA cross-validation
#' @keywords internal
.dpca_lda_cv <- function(X, labels, n_folds = 5) {
  n <- nrow(X)
  folds <- sample(rep(1:n_folds, length.out = n))

  correct <- 0
  for (fold in seq_len(n_folds)) {
    train_idx <- folds != fold
    test_idx <- folds == fold

    # Compute class means
    classes <- unique(labels)
    class_means <- sapply(classes, function(cl) {
      colMeans(X[train_idx & labels == cl, , drop = FALSE])
    })

    # Classify test points by nearest centroid
    for (i in which(test_idx)) {
      dists <- colSums((class_means - X[i, ])^2)
      pred <- classes[which.min(dists)]
      if (pred == labels[i]) correct <- correct + 1
    }
  }

  correct / n
}

#' Plot dPCA Results
#'
#' Visualize dPCA components and projections.
#'
#' @param dpca_model Fitted dPCA model
#' @param marginalization Which marginalization to plot
#' @param components Which components to plot
#' @param plot_type "trajectories", "variance", or "weights"
#'
#' @return ggplot object
#' @export
plot_dpca <- function(dpca_model, marginalization = NULL,
                      components = 1:3, plot_type = "trajectories") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }

  if (is.null(marginalization)) {
    marginalization <- dpca_model$marg_names[1]
  }

  if (plot_type == "trajectories") {
    proj <- dpca_model$projections[[marginalization]]
    n_comp <- min(length(components), dim(proj)[1])
    n_time <- dim(proj)[2]
    n_cond <- dim(proj)[3]

    # Build plot data
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

    ggplot2::ggplot(plot_data,
                    ggplot2::aes(x = time, y = value, color = factor(condition))) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~component, scales = "free_y",
                          labeller = ggplot2::labeller(
                            component = function(x) paste("PC", x))) +
      ggplot2::labs(x = "Time", y = "Projection",
                    color = "Condition",
                    title = paste(marginalization, "components")) +
      ggplot2::theme_minimal()

  } else if (plot_type == "variance") {
    var_data <- do.call(rbind, lapply(names(dpca_model$explained_var), function(m) {
      var_exp <- dpca_model$explained_var[[m]]
      data.frame(
        marginalization = m,
        component = seq_along(var_exp),
        variance = var_exp,
        cumulative = cumsum(var_exp)
      )
    }))

    ggplot2::ggplot(var_data,
                    ggplot2::aes(x = component, y = variance, fill = marginalization)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::labs(x = "Component", y = "Variance explained",
                    title = "Variance by marginalization") +
      ggplot2::theme_minimal()

  } else if (plot_type == "weights") {
    weights <- dpca_model$decoders[[marginalization]]
    n_neurons <- nrow(weights)
    n_comp <- min(length(components), ncol(weights))

    plot_data <- do.call(rbind, lapply(1:n_comp, function(c) {
      data.frame(
        component = c,
        neuron = seq_len(n_neurons),
        weight = weights[, components[c]]
      )
    }))

    ggplot2::ggplot(plot_data,
                    ggplot2::aes(x = neuron, y = weight)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::facet_wrap(~component) +
      ggplot2::labs(x = "Neuron", y = "Weight",
                    title = paste(marginalization, "decoder weights")) +
      ggplot2::theme_minimal()
  }
}

#' Summary Statistics for dPCA
#'
#' Compute summary statistics for dPCA decomposition.
#'
#' @param dpca_model Fitted dPCA model
#'
#' @return Data frame with variance explained by each marginalization
#' @export
dpca_summary <- function(dpca_model) {
  summary_data <- do.call(rbind, lapply(dpca_model$marg_names, function(m) {
    var_exp <- dpca_model$explained_var[[m]]
    data.frame(
      marginalization = m,
      total_var = sum(var_exp),
      pc1_var = var_exp[1],
      pc1_2_var = sum(var_exp[1:min(2, length(var_exp))]),
      pc1_3_var = sum(var_exp[1:min(3, length(var_exp))]),
      n_components = length(var_exp)
    )
  }))

  # Add percentage of total
  total <- sum(summary_data$total_var)
  summary_data$pct_of_total <- summary_data$total_var / total * 100

  summary_data
}
