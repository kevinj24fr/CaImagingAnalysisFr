#' Example `targets` pipeline
#'
#' This defines a simple workflow for generating synthetic data,
#' performing correction and spike inference.
#'
#' @export
calcium_pipeline <- function() {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required for the pipeline")
  }
  list(
    targets::tar_target(raw, generate_synthetic_data()),
    targets::tar_target(corrected, calcium_correction(raw)),
    targets::tar_target(spikes, infer_spikes_advanced(corrected$Cell_1))
  )
}
