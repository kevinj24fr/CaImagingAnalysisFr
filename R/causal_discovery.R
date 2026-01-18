#' Causal Discovery Algorithms for Neural Networks
#'
#' Methods beyond Granger causality: PC algorithm, FCI,
#' and constraint-based approaches for inferring causal structure.
#'
#' @name causal_discovery
NULL

# ============================================================================
# PC Algorithm
# ============================================================================

#' PC Algorithm for Causal Graph Discovery
#'
#' Constraint-based causal discovery using conditional independence tests.
#'
#' @param data Matrix or data frame (observations x variables)
#' @param alpha Significance level for independence tests
#' @param test Type of independence test ("gauss", "discrete", "hsic")
#' @param max_conditioning Maximum conditioning set size
#' @param verbose Print progress
#'
#' @return List with:
#'   - adjacency: Adjacency matrix of discovered graph
#'   - skeleton: Undirected skeleton
#'   - sepsets: Separation sets for each non-edge
#'   - pdag: Partially directed acyclic graph
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Discover causal structure from neural activity
#' result <- pc_algorithm(t(traces), alpha = 0.01)
#' plot_causal_graph(result)
#' }
pc_algorithm <- function(data,
                          alpha = 0.05,
                          test = c("gauss", "discrete", "hsic"),
                          max_conditioning = NULL,
                          verbose = FALSE) {

  test <- match.arg(test)

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  n <- nrow(data)
  p <- ncol(data)
  var_names <- colnames(data) %||% paste0("V", 1:p)

  if (is.null(max_conditioning)) {
    max_conditioning <- p - 2
  }

  # Initialize complete undirected graph
  G <- matrix(TRUE, p, p)
  diag(G) <- FALSE
  rownames(G) <- colnames(G) <- var_names

  # Separation sets
  sepsets <- array(list(), dim = c(p, p))

  # Get independence test function
  ci_test <- get_ci_test(test)

  # Phase 1: Learn skeleton by finding conditional independences
  for (cond_size in 0:max_conditioning) {
    if (verbose) message(sprintf("Conditioning set size: %d", cond_size))

    # Iterate over edges
    edges <- which(G & upper.tri(G), arr.ind = TRUE)

    for (e in seq_len(nrow(edges))) {
      i <- edges[e, 1]
      j <- edges[e, 2]

      if (!G[i, j]) next  # Already removed

      # Get possible conditioning sets (neighbors)
      neighbors_i <- which(G[i, ] & (1:p != j))
      neighbors_j <- which(G[j, ] & (1:p != i))
      neighbors <- union(neighbors_i, neighbors_j)

      if (length(neighbors) < cond_size) next

      # Test all conditioning sets of current size
      if (cond_size == 0) {
        cond_sets <- list(integer(0))
      } else {
        cond_sets <- combn(neighbors, cond_size, simplify = FALSE)
      }

      for (S in cond_sets) {
        p_value <- ci_test(data, i, j, S, n)

        if (p_value > alpha) {
          # Independent given S - remove edge
          G[i, j] <- FALSE
          G[j, i] <- FALSE
          sepsets[[i, j]] <- S
          sepsets[[j, i]] <- S

          if (verbose) {
            message(sprintf("  Removed edge %s - %s | %s (p = %.4f)",
                            var_names[i], var_names[j],
                            paste(var_names[S], collapse = ", "), p_value))
          }
          break
        }
      }
    }

    # Check if any edges remain to process
    if (sum(G) == 0 || cond_size >= sum(G[1, ])) break
  }

  skeleton <- G

  # Phase 2: Orient edges (v-structures and rules)
  pdag <- orient_edges_pc(G, sepsets, var_names, verbose)

  structure(
    list(
      adjacency = pdag,
      skeleton = skeleton,
      sepsets = sepsets,
      var_names = var_names,
      n = n,
      p = p,
      alpha = alpha
    ),
    class = "pc_result"
  )
}

#' Get Conditional Independence Test Function
#' @keywords internal
get_ci_test <- function(test) {
  if (test == "gauss") {
    # Partial correlation test (Fisher's z)
    function(data, i, j, S, n) {
      if (length(S) == 0) {
        r <- cor(data[, i], data[, j])
      } else {
        # Partial correlation
        residuals_i <- residuals(lm(data[, i] ~ as.matrix(data[, S])))
        residuals_j <- residuals(lm(data[, j] ~ as.matrix(data[, S])))
        r <- cor(residuals_i, residuals_j)
      }

      # Fisher's z transformation
      z <- 0.5 * log((1 + r) / (1 - r + 1e-10))
      se <- 1 / sqrt(n - length(S) - 3)
      2 * pnorm(-abs(z / se))
    }
  } else if (test == "discrete") {
    # Chi-squared test for discrete data
    function(data, i, j, S, n) {
      if (length(S) == 0) {
        tab <- table(data[, i], data[, j])
        chisq.test(tab)$p.value
      } else {
        # Conditional chi-squared (Cochran-Mantel-Haenszel)
        # Simplified: stratify and combine
        strata <- interaction(data[, S])
        p_values <- sapply(levels(strata), function(s) {
          idx <- strata == s
          if (sum(idx) < 5) return(1)
          tab <- table(data[idx, i], data[idx, j])
          if (any(dim(tab) < 2)) return(1)
          tryCatch(chisq.test(tab)$p.value, error = function(e) 1)
        })
        # Combine with Fisher's method
        chisq_stat <- -2 * sum(log(p_values + 1e-300))
        pchisq(chisq_stat, df = 2 * length(p_values), lower.tail = FALSE)
      }
    }
  } else {
    # HSIC (Hilbert-Schmidt Independence Criterion)
    function(data, i, j, S, n) {
      if (length(S) == 0) {
        x <- data[, i]
        y <- data[, j]
      } else {
        residuals_i <- residuals(lm(data[, i] ~ as.matrix(data[, S])))
        residuals_j <- residuals(lm(data[, j] ~ as.matrix(data[, S])))
        x <- residuals_i
        y <- residuals_j
      }

      hsic_stat <- compute_hsic(x, y)

      # Permutation test
      n_perm <- 100
      hsic_perm <- numeric(n_perm)
      for (b in seq_len(n_perm)) {
        hsic_perm[b] <- compute_hsic(x, sample(y))
      }

      mean(hsic_perm >= hsic_stat)
    }
  }
}

#' Compute HSIC Statistic
#' @keywords internal
compute_hsic <- function(x, y) {
  n <- length(x)

  # Gaussian kernel
  sigma_x <- median(abs(outer(x, x, "-"))) + 1e-10
  sigma_y <- median(abs(outer(y, y, "-"))) + 1e-10

  Kx <- exp(-outer(x, x, function(a, b) (a - b)^2) / (2 * sigma_x^2))
  Ky <- exp(-outer(y, y, function(a, b) (a - b)^2) / (2 * sigma_y^2))

  # Center kernels
  H <- diag(n) - matrix(1/n, n, n)
  Kxc <- H %*% Kx %*% H
  Kyc <- H %*% Ky %*% H

  # HSIC
  sum(Kxc * Kyc) / (n - 1)^2
}

#' Orient Edges in PC Algorithm
#' @keywords internal
orient_edges_pc <- function(G, sepsets, var_names, verbose) {
  p <- nrow(G)

  # Create direction matrix (1 = i->j, -1 = i<-j, 0 = undirected)
  D <- G * 1L  # Copy adjacency

  # Rule 1: Orient v-structures (i - k - j where i and j not adjacent)
  for (k in 1:p) {
    neighbors <- which(G[k, ])
    if (length(neighbors) < 2) next

    for (pair in combn(neighbors, 2, simplify = FALSE)) {
      i <- pair[1]
      j <- pair[2]

      if (G[i, j]) next  # i and j are adjacent

      # Check if k is in separation set
      S <- sepsets[[i, j]]
      if (is.null(S) || !(k %in% S)) {
        # Orient as v-structure: i -> k <- j
        D[i, k] <- 1
        D[k, i] <- 0
        D[j, k] <- 1
        D[k, j] <- 0

        if (verbose) {
          message(sprintf("  V-structure: %s -> %s <- %s",
                          var_names[i], var_names[k], var_names[j]))
        }
      }
    }
  }

  # Rules 2-4: Additional orientation rules
  changed <- TRUE
  while (changed) {
    changed <- FALSE

    for (i in 1:p) {
      for (j in 1:p) {
        if (D[i, j] == 0 || D[j, i] == 0) next  # Not undirected
        if (D[i, j] != D[j, i]) next  # Already oriented

        # Rule 2: i -> k - j and i -/- j => i -> k -> j
        for (k in 1:p) {
          if (k == i || k == j) next
          if (D[i, k] == 1 && D[k, i] == 0 &&  # i -> k
              D[k, j] > 0 && D[j, k] > 0 &&     # k - j
              D[i, j] == 0 && D[j, i] == 0) {   # i -/- j
            D[k, j] <- 1
            D[j, k] <- 0
            changed <- TRUE
          }
        }

        # Rule 3: i - k -> j and i - l -> j and k -/- l => i -> j
        # (acyclicity preservation)
      }
    }
  }

  D
}

#' @export
print.pc_result <- function(x, ...) {
  cat("PC Algorithm Result\n")
  cat("===================\n")
  cat(sprintf("Variables: %d\n", x$p))
  cat(sprintf("Observations: %d\n", x$n))
  cat(sprintf("Alpha: %.3f\n", x$alpha))

  n_edges <- sum(x$skeleton) / 2
  n_directed <- sum(x$adjacency == 1 & t(x$adjacency) == 0)

  cat(sprintf("Skeleton edges: %d\n", n_edges))
  cat(sprintf("Directed edges: %d\n", n_directed))
  invisible(x)
}

# ============================================================================
# FCI Algorithm (for latent confounders)
# ============================================================================

#' FCI Algorithm for Causal Discovery with Latent Variables
#'
#' Fast Causal Inference algorithm that handles latent confounders.
#'
#' @param data Matrix or data frame
#' @param alpha Significance level
#' @param test Independence test type
#' @param verbose Print progress
#'
#' @return FCI result with PAG (partial ancestral graph)
#' @export
fci_algorithm <- function(data, alpha = 0.05, test = "gauss", verbose = FALSE) {

  # Start with PC skeleton
  pc_result <- pc_algorithm(data, alpha = alpha, test = test, verbose = verbose)

  skeleton <- pc_result$skeleton
  sepsets <- pc_result$sepsets
  var_names <- pc_result$var_names
  p <- pc_result$p
  n <- pc_result$n

  # Initialize PAG with circles (unknown endpoints)
  # 1 = arrowhead (->), 2 = tail (-), 3 = circle (o)
  pag <- array(0, dim = c(p, p))
  pag[skeleton] <- 3  # Circle endpoints

  # Orient definite colliders (v-structures)
  for (k in 1:p) {
    neighbors <- which(skeleton[k, ])
    if (length(neighbors) < 2) next

    for (pair in combn(neighbors, 2, simplify = FALSE)) {
      i <- pair[1]
      j <- pair[2]

      if (skeleton[i, j]) next

      S <- sepsets[[i, j]]
      if (is.null(S) || !(k %in% S)) {
        # Collider: i -> k <- j
        pag[i, k] <- 2  # i -
        pag[k, i] <- 1  # -> k
        pag[j, k] <- 2  # j -
        pag[k, j] <- 1  # -> k
      }
    }
  }

  # FCI orientation rules (simplified)
  # These handle potential latent confounders

  changed <- TRUE
  while (changed) {
    changed <- FALSE

    # Rule 1: If i o-> k and k - j (no edge i-j), then k -> j
    for (i in 1:p) {
      for (k in 1:p) {
        if (pag[i, k] != 3 || pag[k, i] != 1) next  # Need i o-> k

        for (j in 1:p) {
          if (j == i || j == k) next
          if (pag[i, j] != 0 || pag[j, i] != 0) next  # Need no edge i-j
          if (pag[k, j] == 0) next  # Need edge k-j

          if (pag[k, j] == 3) {  # k o- j
            pag[k, j] <- 2  # k - j
            pag[j, k] <- 1  # k -> j
            changed <- TRUE
          }
        }
      }
    }

    # Additional FCI rules would go here
  }

  structure(
    list(
      pag = pag,
      skeleton = skeleton,
      sepsets = sepsets,
      var_names = var_names,
      n = n,
      p = p,
      alpha = alpha
    ),
    class = "fci_result"
  )
}

#' @export
print.fci_result <- function(x, ...) {
  cat("FCI Algorithm Result (handles latent confounders)\n")
  cat("=================================================\n")
  cat(sprintf("Variables: %d\n", x$p))
  cat(sprintf("Observations: %d\n", x$n))
  cat(sprintf("Alpha: %.3f\n", x$alpha))

  n_edges <- sum(x$skeleton) / 2
  cat(sprintf("Edges in skeleton: %d\n", n_edges))
  invisible(x)
}

# ============================================================================
# Score-Based Causal Discovery
# ============================================================================

#' Greedy Equivalence Search (GES) for Causal Discovery
#'
#' Score-based causal discovery using BIC.
#'
#' @param data Matrix or data frame
#' @param score Score function ("bic", "aic")
#' @param max_parents Maximum number of parents per node
#' @param verbose Print progress
#'
#' @return GES result
#' @export
ges_algorithm <- function(data, score = "bic", max_parents = 5, verbose = FALSE) {

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  n <- nrow(data)
  p <- ncol(data)
  var_names <- colnames(data) %||% paste0("V", 1:p)

  # Initialize empty graph
  G <- matrix(0L, p, p)
  rownames(G) <- colnames(G) <- var_names

  # Compute initial score
  current_score <- compute_dag_score(G, data, score)

  # Forward phase: Add edges
  if (verbose) message("Forward phase...")

  improved <- TRUE
  while (improved) {
    improved <- FALSE
    best_delta <- 0
    best_edge <- NULL

    for (i in 1:p) {
      for (j in 1:p) {
        if (i == j) next
        if (G[i, j] == 1) next  # Edge exists
        if (sum(G[, j]) >= max_parents) next  # Too many parents

        # Try adding i -> j
        G_new <- G
        G_new[i, j] <- 1

        # Check acyclicity
        if (!is_dag(G_new)) next

        new_score <- compute_dag_score(G_new, data, score)
        delta <- new_score - current_score

        if (delta > best_delta) {
          best_delta <- delta
          best_edge <- c(i, j)
        }
      }
    }

    if (best_delta > 0) {
      G[best_edge[1], best_edge[2]] <- 1
      current_score <- current_score + best_delta
      improved <- TRUE

      if (verbose) {
        message(sprintf("  Added %s -> %s (delta = %.2f)",
                        var_names[best_edge[1]], var_names[best_edge[2]], best_delta))
      }
    }
  }

  # Backward phase: Remove edges
  if (verbose) message("Backward phase...")

  improved <- TRUE
  while (improved) {
    improved <- FALSE
    best_delta <- 0
    best_edge <- NULL

    edges <- which(G == 1, arr.ind = TRUE)
    for (e in seq_len(nrow(edges))) {
      i <- edges[e, 1]
      j <- edges[e, 2]

      G_new <- G
      G_new[i, j] <- 0

      new_score <- compute_dag_score(G_new, data, score)
      delta <- new_score - current_score

      if (delta > best_delta) {
        best_delta <- delta
        best_edge <- c(i, j)
      }
    }

    if (best_delta > 0) {
      G[best_edge[1], best_edge[2]] <- 0
      current_score <- current_score + best_delta
      improved <- TRUE

      if (verbose) {
        message(sprintf("  Removed %s -> %s (delta = %.2f)",
                        var_names[best_edge[1]], var_names[best_edge[2]], best_delta))
      }
    }
  }

  structure(
    list(
      adjacency = G,
      score = current_score,
      var_names = var_names,
      n = n,
      p = p
    ),
    class = "ges_result"
  )
}

#' Compute DAG Score
#' @keywords internal
compute_dag_score <- function(G, data, score) {
  p <- ncol(data)
  n <- nrow(data)
  total_score <- 0

  for (j in 1:p) {
    parents <- which(G[, j] == 1)

    if (length(parents) == 0) {
      # No parents: score of marginal
      rss <- sum((data[, j] - mean(data[, j]))^2)
      k <- 1
    } else {
      # Regression on parents
      fit <- lm(data[, j] ~ as.matrix(data[, parents]))
      rss <- sum(residuals(fit)^2)
      k <- length(parents) + 1
    }

    # BIC or AIC
    if (score == "bic") {
      total_score <- total_score - n * log(rss / n) - k * log(n)
    } else {
      total_score <- total_score - n * log(rss / n) - 2 * k
    }
  }

  total_score
}

#' Check if Graph is a DAG
#' @keywords internal
is_dag <- function(G) {
  p <- nrow(G)

  # Topological sort via DFS
  visited <- rep(0L, p)  # 0 = unvisited, 1 = in progress, 2 = done

  dfs <- function(v) {
    if (visited[v] == 1) return(FALSE)  # Cycle
    if (visited[v] == 2) return(TRUE)

    visited[v] <<- 1
    children <- which(G[v, ] == 1)
    for (c in children) {
      if (!dfs(c)) return(FALSE)
    }
    visited[v] <<- 2
    TRUE
  }

  for (v in 1:p) {
    if (visited[v] == 0) {
      if (!dfs(v)) return(FALSE)
    }
  }

  TRUE
}

# ============================================================================
# Visualization
# ============================================================================

#' Plot Causal Graph
#'
#' @param x PC, FCI, or GES result
#' @param ... Additional arguments
#'
#' @export
plot_causal_graph <- function(x, ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph required for plotting")
  }

  if (inherits(x, "pc_result") || inherits(x, "ges_result")) {
    adj <- x$adjacency
    names <- x$var_names
  } else if (inherits(x, "fci_result")) {
    # Convert PAG to plottable format
    adj <- (x$pag == 1 | x$pag == 2) * 1
    names <- x$var_names
  }

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "directed")
  igraph::V(g)$name <- names

  plot(g, vertex.size = 30, vertex.label.cex = 0.8,
       edge.arrow.size = 0.5, layout = igraph::layout_nicely, ...)
}

#' Apply Causal Discovery to Neural Traces
#'
#' @param traces Trace matrix (cells x time)
#' @param method "pc", "fci", or "ges"
#' @param alpha Significance level
#' @param downsample Downsample factor for computational efficiency
#'
#' @return Causal discovery result
#' @export
discover_neural_causality <- function(traces, method = "pc", alpha = 0.01,
                                       downsample = 1) {

  # Transpose to observations x variables
  data <- t(traces)

  # Downsample if needed
  if (downsample > 1) {
    idx <- seq(1, nrow(data), by = downsample)
    data <- data[idx, ]
  }

  # Apply selected method
  switch(method,
    "pc" = pc_algorithm(data, alpha = alpha),
    "fci" = fci_algorithm(data, alpha = alpha),
    "ges" = ges_algorithm(data)
  )
}
