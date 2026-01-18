#' Package Utilities
#'
#' Common utility functions used across the package.
#'
#' @name utils
#' @keywords internal
NULL

# ============================================================================
# Null Coalescing Operator
# ============================================================================

#' Null Coalescing Operator
#'
#' Returns the left-hand side if not NULL, otherwise the right-hand side.
#' This is the canonical definition used throughout the package.
#'
#' @param x Left-hand side value
#' @param y Right-hand side value (default)
#'
#' @return x if not NULL, otherwise y
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' NULL %||% "default"  # Returns "default"
#' "value" %||% "default"  # Returns "value"
#' }
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# ============================================================================
# Safe Matrix Operations
# ============================================================================

#' Safe Matrix Inverse
#'
#' Compute matrix inverse with regularization for numerical stability.
#'
#' @param x Square matrix
#' @param eps Regularization parameter (added to diagonal)
#'
#' @return Inverse matrix
#' @keywords internal
safe_solve <- function(x, eps = 1e-8) {
  if (!is.matrix(x)) x <- as.matrix(x)
  n <- nrow(x)

  tryCatch({
    solve(x + eps * diag(n))
  }, error = function(e) {
    # Fall back to pseudoinverse via SVD
    s <- svd(x)
    d <- s$d
    d[d > eps] <- 1 / d[d > eps]
    d[d <= eps] <- 0
    s$v %*% diag(d) %*% t(s$u)
  })
}

#' Safe Determinant
#'
#' Compute matrix determinant with protection against numerical issues.
#'
#' @param x Square matrix
#' @param eps Minimum value to return
#' @param log Return log determinant
#'
#' @return Determinant or log-determinant
#' @keywords internal
safe_det <- function(x, eps = 1e-300, log = FALSE) {
  if (!is.matrix(x)) x <- as.matrix(x)

  d <- tryCatch({
    det(x)
  }, error = function(e) {
    # Use sum of log eigenvalues as fallback
    eig <- eigen(x, only.values = TRUE)$values
    prod(pmax(Re(eig), eps))
  })

  if (d <= 0) d <- eps

  if (log) log(d) else d
}

# ============================================================================
# Type Checking Utilities
# ============================================================================

#' Check if Object is from lme4 Package
#'
#' @param x Object to check
#' @return Logical
#' @keywords internal
is_lme4_fit <- function(x) {
  inherits(x, c("lmerMod", "glmerMod", "merMod"))
}

#' Check if Object is from nlme Package
#'
#' @param x Object to check
#' @return Logical
#' @keywords internal
is_nlme_fit <- function(x) {

  inherits(x, c("lme", "gls", "nlme"))
}

#' Safely Extract Predictions
#'
#' Extract predictions handling both vector and list returns.
#'
#' @param pred Prediction result from predict()
#' @return Numeric vector of predictions
#' @keywords internal
safe_extract_predictions <- function(pred) {
  if (is.list(pred)) {
    pred$predictions %||% pred$fit %||% pred[[1]]
  } else {
    as.numeric(pred)
  }
}
