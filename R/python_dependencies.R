#' Check and Install Python Dependencies
#' 
#' Manages Python package dependencies for the CaImagingAnalysisFr package.
#' 
#' @param packages Vector of package names to check/install
#' @param install_missing Whether to install missing packages (default: TRUE)
#' @param verbose Whether to show progress messages (default: TRUE)
#' @return List of package status
#' @export
manage_python_dependencies <- function(packages = NULL, 
                                      install_missing = TRUE,
                                      verbose = TRUE) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for Python integration. Install with: install.packages('reticulate')")
  }
  
  config <- get_config()
  
  if (is.null(packages)) {
    packages <- names(config$python_packages)
  }
  
  if (verbose) {
    message("Checking Python dependencies...")
  }
  
  # Initialize Python if not already done
  if (!reticulate::py_available()) {
    if (verbose) message("  Initializing Python...")
    reticulate::use_python(reticulate::py_config()$python)
  }
  
  # Package name mapping for installation vs import
  package_mapping <- list(
    scikit_learn = "scikit-learn",  # Install as scikit-learn, import as sklearn
    umap_learn = "umap-learn",      # Install as umap-learn, import as umap
    openTSNE = "openTSNE"           # Install as openTSNE, import as openTSNE
  )
  
  results <- list()
  
  for (pkg in packages) {
    if (verbose) message("  Checking ", pkg, "...")
    
    # Determine import name and install name
    import_name <- pkg
    install_name <- pkg
    
    if (pkg %in% names(package_mapping)) {
      install_name <- package_mapping[[pkg]]
    }
    
    # Special handling for scikit-learn (imported as sklearn in Python)
    if (pkg == "scikit_learn") {
      import_name <- "sklearn"
    }
    
    # Check if package is available
    available <- reticulate::py_module_available(import_name)
    
    if (available) {
      # Check version if specified
      version_ok <- TRUE
      if (pkg %in% names(config$python_packages)) {
        required_version <- config$python_packages[[pkg]]
        tryCatch({
          pkg_module <- reticulate::import(import_name)
          if ("__version__" %in% names(pkg_module)) {
            current_version <- pkg_module$`__version__`
            version_ok <- compare_versions(current_version, required_version)
          }
        }, error = function(e) {
          version_ok <<- FALSE
        })
      }
      
      results[[pkg]] <- list(
        available = TRUE,
        version_ok = version_ok,
        message = ifelse(version_ok, "OK", "Version mismatch")
      )
    } else {
      results[[pkg]] <- list(
        available = FALSE,
        version_ok = FALSE,
        message = "Not installed"
      )
      
      if (install_missing) {
        if (verbose) message("    Installing ", install_name, "...")
        tryCatch({
          # Set environment variable to allow deprecated scikit-learn package if needed
          if (pkg == "scikit_learn") {
            Sys.setenv(SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL = "TRUE")
          }
          
          reticulate::py_install(install_name, pip = TRUE)
          
          # Verify installation
          if (reticulate::py_module_available(import_name)) {
            results[[pkg]]$available <- TRUE
            results[[pkg]]$version_ok <- TRUE
            results[[pkg]]$message <- "Installed successfully"
          } else {
            results[[pkg]]$message <- "Installation completed but module not available"
          }
        }, error = function(e) {
          results[[pkg]]$message <- paste("Installation failed:", e$message)
        })
      }
    }
  }
  
  if (verbose) {
    message("Python dependency check completed")
  }
  
  return(results)
}

#' Compare version strings
#' 
#' @param current Current version string
#' @param required Required version string with operator
#' @return TRUE if version requirement is met
#' @keywords internal
compare_versions <- function(current, required) {
  if (is.null(current) || is.null(required)) return(TRUE)
  
  # Parse version strings
  current_parts <- as.numeric(strsplit(current, "\\.")[[1]])
  required_parts <- as.numeric(strsplit(gsub("[^0-9.]", "", required), "\\.")[[1]])
  
  # Pad with zeros for comparison
  max_len <- max(length(current_parts), length(required_parts))
  current_parts <- c(current_parts, rep(0, max_len - length(current_parts)))
  required_parts <- c(required_parts, rep(0, max_len - length(required_parts)))
  
  # Compare versions
  for (i in seq_along(current_parts)) {
    if (current_parts[i] > required_parts[i]) return(TRUE)
    if (current_parts[i] < required_parts[i]) return(FALSE)
  }
  
  return(TRUE)
}

#' Get Python environment info
#' 
#' @return List of Python environment information
#' @export
get_python_info <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(list(error = "reticulate package not available"))
  }
  
  tryCatch({
    config <- reticulate::py_config()
    list(
      python = config$python,
      version = config$version,
      numpy = config$numpy,
      available = reticulate::py_available(),
      modules = reticulate::py_list_packages()
    )
  }, error = function(e) {
    list(error = e$message)
  })
}

#' Pre-install all required Python packages
#' 
#' @param verbose Whether to show progress messages (default: TRUE)
#' @return TRUE if successful
#' @export
install_python_dependencies <- function(verbose = TRUE) {
  if (verbose) {
    message("Installing all required Python packages...")
  }
  
  results <- manage_python_dependencies(install_missing = TRUE, verbose = verbose)
  
  # Check if all packages are available
  all_available <- all(sapply(results, function(x) x$available))
  all_versions_ok <- all(sapply(results, function(x) x$version_ok))
  
  if (all_available && all_versions_ok) {
    if (verbose) {
      message("All Python dependencies installed successfully!")
    }
    return(TRUE)
  } else {
    failed_packages <- names(results)[!sapply(results, function(x) x$available && x$version_ok)]
    warning("Failed to install/verify packages: ", paste(failed_packages, collapse = ", "))
    
    # Provide helpful error messages for common issues
    if ("scikit_learn" %in% failed_packages) {
      message("\nFor scikit-learn installation issues:")
      message("1. Try: pip install scikit-learn")
      message("2. If you get scikit-learn deprecation errors, use: pip install scikit-learn")
      message("3. Make sure you have a compatible Python version (>=3.8)")
    }
    
    if ("pymc" %in% failed_packages) {
      message("\nFor PyMC installation issues:")
      message("1. Try: pip install pymc")
      message("2. PyMC requires specific versions of numpy and scipy")
      message("3. Consider using conda: conda install -c conda-forge pymc")
    }
    
    if ("tensorflow" %in% failed_packages) {
      message("\nFor TensorFlow installation issues:")
      message("1. Try: pip install tensorflow")
      message("2. For Apple Silicon Macs: pip install tensorflow-macos")
      message("3. Consider using conda: conda install tensorflow")
    }
    
    return(FALSE)
  }
}

#' Install specific Python packages with error handling
#' 
#' @param packages Vector of package names to install
#' @param verbose Whether to show progress messages (default: TRUE)
#' @return TRUE if successful
#' @export
install_specific_python_packages <- function(packages, verbose = TRUE) {
  if (verbose) {
    message("Installing specific Python packages: ", paste(packages, collapse = ", "))
  }
  
  results <- manage_python_dependencies(packages = packages, install_missing = TRUE, verbose = verbose)
  
  # Check if all packages are available
  all_available <- all(sapply(results, function(x) x$available))
  
  if (all_available) {
    if (verbose) {
      message("All specified Python packages installed successfully!")
    }
    return(TRUE)
  } else {
    failed_packages <- names(results)[!sapply(results, function(x) x$available)]
    warning("Failed to install packages: ", paste(failed_packages, collapse = ", "))
    return(FALSE)
  }
} 