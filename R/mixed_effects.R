#' Mixed-Effects Models for Hierarchical Neural Data
#'
#' Proper statistical inference accounting for nested structure:
#' cells within animals, trials within sessions, etc.
#'
#' @name mixed_effects
NULL

# ============================================================================
# Linear Mixed-Effects Models
# ============================================================================

#' Fit Linear Mixed-Effects Model for Neural Data
#'
#' Fit LMM to account for hierarchical structure in neural experiments.
#'
#' @param formula Model formula (e.g., response ~ condition + (1|animal/cell))
#' @param data Data frame with observations
#' @param REML Use restricted maximum likelihood
#' @param control Control parameters for optimization
#'
#' @return Mixed model fit object
#' @export
#'
#' @examples
#' \dontrun{
#' # Account for cell-within-animal structure
#' data <- data.frame(
#'   response = cell_responses,
#'   condition = conditions,
#'   cell_id = cell_ids,
#'   animal_id = animal_ids
#' )
#'
#' fit <- fit_lmm(response ~ condition + (1|animal_id/cell_id), data)
#' summary(fit)
#' }
fit_lmm <- function(formula, data, REML = TRUE, control = NULL) {

  # Check for lme4
  if (requireNamespace("lme4", quietly = TRUE)) {
    if (is.null(control)) {
      control <- lme4::lmerControl(optimizer = "bobyqa")
    }
    fit <- lme4::lmer(formula, data = data, REML = REML, control = control)
    class(fit) <- c("ca_lmm", class(fit))
    return(fit)
  }

  # Fallback to nlme
  if (requireNamespace("nlme", quietly = TRUE)) {
    # Convert lme4-style formula to nlme
    fit <- nlme::lme(formula, data = data, method = if (REML) "REML" else "ML")
    class(fit) <- c("ca_lmm", class(fit))
    return(fit)
  }

  # Manual implementation for simple random intercept
  fit_lmm_manual(formula, data, REML)
}

#' Manual LMM Implementation (Simple Random Intercept)
#' @keywords internal
fit_lmm_manual <- function(formula, data, REML = TRUE) {

  # Parse formula
  formula_parts <- as.character(formula)
  response_var <- formula_parts[2]
  predictors <- formula_parts[3]

  # Extract random effect
  random_match <- regmatches(predictors, regexpr("\\(1\\|([^)]+)\\)", predictors))
  if (length(random_match) == 0) {
    stop("Manual LMM requires (1|group) random intercept specification")
  }

  group_var <- gsub("\\(1\\|([^)]+)\\)", "\\1", random_match)

  # Extract fixed effects
  fixed_part <- gsub("\\s*\\+?\\s*\\(1\\|[^)]+\\)", "", predictors)
  fixed_part <- trimws(fixed_part)

  if (nchar(fixed_part) == 0) {
    fixed_formula <- as.formula(paste(response_var, "~ 1"))
  } else {
    fixed_formula <- as.formula(paste(response_var, "~", fixed_part))
  }

  # Get response and groups
  y <- data[[response_var]]
  groups <- as.factor(data[[group_var]])
  n_groups <- nlevels(groups)
  n <- length(y)

  # Design matrix for fixed effects
  X <- model.matrix(fixed_formula, data = data)
  p <- ncol(X)

  # Design matrix for random effects
  Z <- model.matrix(~ groups - 1)

  # EM algorithm for variance components
  # Initialize
  sigma2_e <- var(y) * 0.5
  sigma2_u <- var(y) * 0.5
  beta <- solve(t(X) %*% X) %*% t(X) %*% y

  for (iter in 1:100) {
    # E-step: Compute BLUPs
    V <- sigma2_e * diag(n) + sigma2_u * Z %*% t(Z)
    V_inv <- solve(V + 1e-6 * diag(n))

    # Fixed effects
    beta_new <- solve(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv %*% y

    # Random effects (BLUPs)
    resid <- y - X %*% beta_new
    u <- sigma2_u * t(Z) %*% V_inv %*% resid

    # M-step: Update variance components
    fitted <- X %*% beta_new + Z %*% u
    sigma2_e_new <- sum((y - fitted)^2) / (n - p)
    sigma2_u_new <- sum(u^2) / n_groups

    # Check convergence
    if (abs(sigma2_e_new - sigma2_e) < 1e-6 && abs(sigma2_u_new - sigma2_u) < 1e-6) {
      break
    }

    sigma2_e <- sigma2_e_new
    sigma2_u <- sigma2_u_new
    beta <- beta_new
  }

  # Compute standard errors
  V <- sigma2_e * diag(n) + sigma2_u * Z %*% t(Z)
  V_inv <- solve(V + 1e-6 * diag(n))
  vcov_beta <- solve(t(X) %*% V_inv %*% X)

  # Log-likelihood
  logLik <- -0.5 * (n * log(2 * pi) + determinant(V)$modulus + t(y - X %*% beta) %*% V_inv %*% (y - X %*% beta))

  # ICC
  icc <- sigma2_u / (sigma2_u + sigma2_e)

  structure(
    list(
      beta = as.vector(beta),
      beta_names = colnames(X),
      random_effects = as.vector(u),
      group_names = levels(groups),
      sigma2_e = sigma2_e,
      sigma2_u = sigma2_u,
      vcov = vcov_beta,
      se = sqrt(diag(vcov_beta)),
      logLik = as.numeric(logLik),
      icc = icc,
      n = n,
      n_groups = n_groups,
      REML = REML,
      formula = formula
    ),
    class = "ca_lmm_manual"
  )
}

#' @export
print.ca_lmm_manual <- function(x, ...) {
  cat("Linear Mixed-Effects Model\n")
  cat("==========================\n")
  cat(sprintf("Formula: %s\n", deparse(x$formula)))
  cat(sprintf("Observations: %d\n", x$n))
  cat(sprintf("Groups: %d\n", x$n_groups))
  cat("\nFixed Effects:\n")
  for (i in seq_along(x$beta)) {
    cat(sprintf("  %s: %.4f (SE: %.4f)\n", x$beta_names[i], x$beta[i], x$se[i]))
  }
  cat("\nVariance Components:\n")
  cat(sprintf("  Residual: %.4f\n", x$sigma2_e))
  cat(sprintf("  Random Intercept: %.4f\n", x$sigma2_u))
  cat(sprintf("  ICC: %.4f\n", x$icc))
  invisible(x)
}

#' @export
summary.ca_lmm_manual <- function(object, ...) {
  z_values <- object$beta / object$se
  p_values <- 2 * pnorm(-abs(z_values))

  cat("Linear Mixed-Effects Model Summary\n")
  cat("==================================\n\n")
  print(object)
  cat("\nSignificance:\n")
  for (i in seq_along(object$beta)) {
    sig <- if (p_values[i] < 0.001) "***" else if (p_values[i] < 0.01) "**" else if (p_values[i] < 0.05) "*" else ""
    cat(sprintf("  %s: z = %.2f, p = %.4f %s\n", object$beta_names[i], z_values[i], p_values[i], sig))
  }
  invisible(object)
}

# ============================================================================
# Generalized Linear Mixed Models
# ============================================================================

#' Fit Generalized Linear Mixed Model
#'
#' For count data (spikes) or binary outcomes.
#'
#' @param formula Model formula
#' @param data Data frame
#' @param family GLM family ("poisson", "binomial", "gaussian")
#' @param nAGQ Number of adaptive Gauss-Hermite quadrature points
#'
#' @return GLMM fit
#' @export
fit_glmm <- function(formula, data, family = "poisson", nAGQ = 1) {

  if (requireNamespace("lme4", quietly = TRUE)) {
    fam <- switch(family,
      "poisson" = poisson(),
      "binomial" = binomial(),
      "gaussian" = gaussian()
    )

    fit <- lme4::glmer(formula, data = data, family = fam, nAGQ = nAGQ)
    class(fit) <- c("ca_glmm", class(fit))
    return(fit)
  }

  stop("lme4 package required for GLMM. Install with: install.packages('lme4')")
}

# ============================================================================
# Convenience Functions for Neural Data
# ============================================================================

#' Prepare Data for Mixed Model Analysis
#'
#' @param traces Trace matrix (cells x time)
#' @param conditions Condition labels (length = time)
#' @param cell_info Data frame with cell metadata (cell_id, animal_id, etc.)
#' @param time_bins Optional time binning factor
#'
#' @return Long-format data frame ready for mixed modeling
#' @export
prepare_mixed_data <- function(traces, conditions, cell_info, time_bins = NULL) {

  n_cells <- nrow(traces)
  n_time <- ncol(traces)

  # Create long format
  data <- data.frame(
    cell_id = rep(rownames(traces) %||% seq_len(n_cells), each = n_time),
    time = rep(seq_len(n_time), n_cells),
    response = as.vector(t(traces))
  )

  # Add conditions
  data$condition <- rep(conditions, n_cells)

  # Add cell info
  if (!is.null(cell_info)) {
    data <- merge(data, cell_info, by = "cell_id", all.x = TRUE)
  }

  # Add time bins
  if (!is.null(time_bins)) {
    data$time_bin <- cut(data$time, breaks = time_bins, labels = FALSE)
  }

  data
}

#' Test Condition Effect with Mixed Model
#'
#' @param traces Trace matrix (cells x time)
#' @param conditions Condition vector
#' @param cell_groups Cell grouping variable (e.g., animal_id)
#' @param baseline_frames Frames to use as baseline
#' @param response_frames Frames to use as response
#'
#' @return Mixed model results
#' @export
test_condition_mixed <- function(traces, conditions, cell_groups,
                                  baseline_frames = NULL,
                                  response_frames = NULL) {

  n_cells <- nrow(traces)

  # Compute summary statistics per cell
  if (!is.null(baseline_frames) && !is.null(response_frames)) {
    baseline <- rowMeans(traces[, baseline_frames, drop = FALSE])
    response <- rowMeans(traces[, response_frames, drop = FALSE])
    cell_response <- response - baseline
  } else {
    cell_response <- rowMeans(traces)
  }

  # Create data frame
  data <- data.frame(
    response = cell_response,
    condition = conditions,
    group = cell_groups,
    cell_id = seq_len(n_cells)
  )

  # Fit mixed model
  fit <- fit_lmm(response ~ condition + (1|group), data = data)

  fit
}

#' Compare Multiple Conditions with Mixed Model
#'
#' @param traces_list List of trace matrices (one per condition)
#' @param condition_names Names of conditions
#' @param cell_groups Cell grouping variable
#'
#' @return Mixed model with post-hoc comparisons
#' @export
compare_conditions_mixed <- function(traces_list, condition_names, cell_groups) {

  # Combine data
  all_data <- lapply(seq_along(traces_list), function(i) {
    traces <- traces_list[[i]]
    data.frame(
      response = rowMeans(traces),
      condition = condition_names[i],
      group = cell_groups,
      cell_id = seq_len(nrow(traces))
    )
  })

  data <- do.call(rbind, all_data)
  data$condition <- factor(data$condition, levels = condition_names)

  # Fit model
  fit <- fit_lmm(response ~ condition + (1|group), data = data)

  # Extract contrasts
  if (requireNamespace("lme4", quietly = TRUE) && inherits(fit, "lmerMod")) {
    if (requireNamespace("emmeans", quietly = TRUE)) {
      emm <- emmeans::emmeans(fit, "condition")
      contrasts <- emmeans::pairs(emm)
      return(list(model = fit, emmeans = emm, contrasts = contrasts))
    }
  }

  list(model = fit)
}

#' Compute Intraclass Correlation Coefficient
#'
#' @param response Response variable
#' @param group Grouping variable
#'
#' @return ICC value and confidence interval
#' @export
compute_icc <- function(response, group) {
  data <- data.frame(response = response, group = as.factor(group))

  # Fit null model with only random intercept
  fit <- fit_lmm(response ~ 1 + (1|group), data = data)

  if (inherits(fit, "ca_lmm_manual")) {
    icc <- fit$icc

    # Bootstrap CI
    n_boot <- 100
    icc_boot <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      groups_sample <- sample(unique(data$group), replace = TRUE)
      data_boot <- do.call(rbind, lapply(groups_sample, function(g) {
        data[data$group == g, ]
      }))
      fit_boot <- tryCatch(
        fit_lmm(response ~ 1 + (1|group), data = data_boot),
        error = function(e) NULL
      )
      if (!is.null(fit_boot)) {
        icc_boot[b] <- fit_boot$icc
      }
    }

    ci <- quantile(icc_boot, c(0.025, 0.975), na.rm = TRUE)

  } else if (requireNamespace("lme4", quietly = TRUE)) {
    vc <- as.data.frame(lme4::VarCorr(fit))
    sigma2_u <- vc$vcov[vc$grp == "group"]
    sigma2_e <- vc$vcov[vc$grp == "Residual"]
    icc <- sigma2_u / (sigma2_u + sigma2_e)

    # Approximate CI
    se_icc <- sqrt(icc * (1 - icc) / length(unique(group)))
    ci <- c(icc - 1.96 * se_icc, icc + 1.96 * se_icc)
  }

  list(icc = icc, ci = ci)
}

#' Power Analysis for Mixed Model Design
#'
#' @param n_cells Number of cells per group
#' @param n_groups Number of groups (e.g., animals)
#' @param effect_size Expected effect size (Cohen's d)
#' @param icc Intraclass correlation
#' @param alpha Significance level
#' @param n_sim Number of simulations
#'
#' @return Estimated power
#' @export
power_mixed <- function(n_cells = 20, n_groups = 5, effect_size = 0.5,
                        icc = 0.3, alpha = 0.05, n_sim = 100) {

  # Design effect
  design_effect <- 1 + (n_cells - 1) * icc
  effective_n <- (n_cells * n_groups) / design_effect

  # Simulate power
  significant <- 0

  for (sim in seq_len(n_sim)) {
    # Generate data
    sigma2_total <- 1
    sigma2_u <- icc * sigma2_total
    sigma2_e <- (1 - icc) * sigma2_total

    group_effects <- rnorm(n_groups, 0, sqrt(sigma2_u))
    response <- numeric(n_cells * n_groups * 2)
    condition <- rep(c("control", "treatment"), each = n_cells * n_groups)
    group <- rep(rep(1:n_groups, each = n_cells), 2)

    for (g in 1:n_groups) {
      idx_ctrl <- ((g - 1) * n_cells + 1):(g * n_cells)
      idx_trt <- idx_ctrl + n_cells * n_groups

      response[idx_ctrl] <- group_effects[g] + rnorm(n_cells, 0, sqrt(sigma2_e))
      response[idx_trt] <- effect_size + group_effects[g] + rnorm(n_cells, 0, sqrt(sigma2_e))
    }

    data <- data.frame(response = response, condition = condition, group = factor(group))

    # Fit and test
    fit <- tryCatch({
      fit_lmm(response ~ condition + (1|group), data = data)
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      if (inherits(fit, "ca_lmm_manual")) {
        z <- fit$beta[2] / fit$se[2]
        p <- 2 * pnorm(-abs(z))
      } else if (requireNamespace("lme4", quietly = TRUE)) {
        coefs <- summary(fit)$coefficients
        p <- coefs[2, "Pr(>|t|)"]
      } else {
        p <- 1
      }

      if (p < alpha) significant <- significant + 1
    }
  }

  power <- significant / n_sim

  list(
    power = power,
    n_cells = n_cells,
    n_groups = n_groups,
    effect_size = effect_size,
    icc = icc,
    effective_n = effective_n,
    design_effect = design_effect
  )
}

# Helper for NULL default
`%||%` <- function(a, b) if (is.null(a)) b else a
