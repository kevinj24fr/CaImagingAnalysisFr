#' Interactive and Reproducible Reporting
#'
#' Generate interactive HTML reports with embedded plots, QC metrics, and code,
#' and export full analysis pipelines as reproducible scripts or notebooks.
#'
#' @name reporting
#' @docType package
NULL

#' Generate Comprehensive Analysis Report
#'
#' Create an interactive HTML report with all analysis results, plots, and code.
#'
#' @param analysis_results List containing all analysis results
#' @param output_file Path for output HTML file (default: "calcium_analysis_report.html")
#' @param report_title Title for the report (default: "Calcium Imaging Analysis Report")
#' @param include_code Whether to include R code in the report (default: TRUE)
#' @param include_plots Whether to include interactive plots (default: TRUE)
#' @param include_qc Whether to include quality control metrics (default: TRUE)
#' @param ... Additional arguments
#' @return Path to generated report
#' @export
generate_analysis_report <- function(analysis_results, output_file = "calcium_analysis_report.html", report_title = "Calcium Imaging Analysis Report", include_code = TRUE, include_plots = TRUE, include_qc = TRUE, ...) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required for report generation. Install with: install.packages('rmarkdown')")
  }
  
  # Create temporary Rmd file
  rmd_content <- create_report_rmd(analysis_results, report_title, include_code, include_plots, include_qc)
  temp_rmd <- tempfile(fileext = ".Rmd")
  writeLines(rmd_content, temp_rmd)
  
  # Generate report
  rmarkdown::render(
    input = temp_rmd,
    output_file = output_file,
    output_format = "html_document",
    params = list(
      analysis_results = analysis_results,
      include_code = include_code,
      include_plots = include_plots,
      include_qc = include_qc
    ),
    quiet = FALSE
  )
  
  # Clean up
  unlink(temp_rmd)
  
  message("Report generated: ", output_file)
  return(output_file)
}

#' Create RMarkdown Content for Report
#'
#' Generate RMarkdown content for the analysis report.
#'
#' @param analysis_results List containing analysis results
#' @param title Report title
#' @param include_code Whether to include code
#' @param include_plots Whether to include plots
#' @param include_qc Whether to include QC
#' @return RMarkdown content as character vector
#' @export
create_report_rmd <- function(analysis_results, title, include_code, include_plots, include_qc) {
  c(
    "---",
    paste0("title: \"", title, "\""),
    "date: \"`r Sys.Date()`\"",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    code_folding: show",
    "    theme: cosmo",
    "    highlight: tango",
    "params:",
    "  analysis_results: NULL",
    "  include_code: TRUE",
    "  include_plots: TRUE",
    "  include_qc: TRUE",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = params$include_code, warning = FALSE, message = FALSE)",
    "library(CaImagingAnalysisFr)",
    "library(ggplot2)",
    "library(knitr)",
    "```",
    "",
    "# Executive Summary",
    "",
    "This report presents the results of calcium imaging analysis performed using the CaImagingAnalysisFr package.",
    "",
    "## Analysis Overview",
    "",
    "```{r overview}",
    "cat(\"Analysis completed on:\", format(Sys.time(), \"%Y-%m-%d %H:%M:%S\"), \"\\n\")",
    "cat(\"Package version:\", packageVersion(\"CaImagingAnalysisFr\"), \"\\n\")",
    "```",
    "",
    "# Data Summary",
    "",
    "```{r data_summary}",
    "if (!is.null(params$analysis_results$data_summary)) {",
    "  kable(params$analysis_results$data_summary, caption = \"Data Summary\")",
    "}",
    "```",
    "",
    "# Quality Control",
    "",
    if (include_qc) {
      c(
        "## Signal Quality Metrics",
        "",
        "```{r qc_metrics}",
        "if (!is.null(params$analysis_results$qc_metrics)) {",
        "  kable(params$analysis_results$qc_metrics, caption = \"Quality Control Metrics\")",
        "}",
        "```",
        ""
      )
    } else NULL,
    "",
    "# Analysis Results",
    "",
    "## Spike Detection",
    "",
    "```{r spike_results}",
    "if (!is.null(params$analysis_results$spike_results)) {",
    "  cat(\"Number of detected spikes:\", sum(params$analysis_results$spike_results$spikes), \"\\n\")",
    "  cat(\"Mean spike rate:\", mean(params$analysis_results$spike_results$spike_rates), \"Hz\\n\")",
    "}",
    "```",
    "",
    if (include_plots) {
      c(
        "## Visualization",
        "",
        "```{r plots, fig.width=10, fig.height=6}",
        "if (!is.null(params$analysis_results$plots)) {",
        "  print(params$analysis_results$plots$main_plot)",
        "}",
        "```",
        ""
      )
    } else NULL,
    "",
    "# Statistical Analysis",
    "",
    "```{r stats}",
    "if (!is.null(params$analysis_results$statistical_analysis)) {",
    "  kable(params$analysis_results$statistical_analysis$summary, caption = \"Statistical Summary\")",
    "}",
    "```",
    "",
    "# Conclusions",
    "",
    "The analysis has been completed successfully. Please review the results and consult with your research team for interpretation.",
    "",
    "---",
    "",
    "*Report generated by CaImagingAnalysisFr package*"
  )
}

#' Export Analysis Pipeline as R Script
#'
#' Export the complete analysis pipeline as a reproducible R script.
#'
#' @param analysis_results List containing analysis results
#' @param output_file Path for output R script (default: "calcium_analysis_pipeline.R")
#' @param include_comments Whether to include detailed comments (default: TRUE)
#' @param include_data_loading Whether to include data loading code (default: TRUE)
#' @param ... Additional arguments
#' @return Path to generated script
#' @export
export_analysis_pipeline <- function(analysis_results, output_file = "calcium_analysis_pipeline.R", include_comments = TRUE, include_data_loading = TRUE, ...) {
  # Create script content
  script_content <- create_pipeline_script(analysis_results, include_comments, include_data_loading)
  
  # Write to file
  writeLines(script_content, output_file)
  
  message("Analysis pipeline exported: ", output_file)
  return(output_file)
}

#' Create Pipeline Script Content
#'
#' Generate R script content for the analysis pipeline.
#'
#' @param analysis_results List containing analysis results
#' @param include_comments Whether to include comments
#' @param include_data_loading Whether to include data loading
#' @return Script content as character vector
#' @export
create_pipeline_script <- function(analysis_results, include_comments, include_data_loading) {
  c(
    if (include_comments) {
      c(
        "# Calcium Imaging Analysis Pipeline",
        "# Generated by CaImagingAnalysisFr package",
        paste0("# Date: ", Sys.Date()),
        "#",
        "# This script reproduces the complete analysis pipeline.",
        "# Make sure to install required packages before running.",
        "#",
        ""
      )
    } else NULL,
    "# Load required packages",
    "library(CaImagingAnalysisFr)",
    "library(ggplot2)",
    "",
    if (include_data_loading) {
      c(
        "# Load your data",
        "# Replace 'your_data.csv' with your actual data file",
        "# data <- read.csv('your_data.csv')",
        "# raw_data <- data[, grep('^Cell_', colnames(data))]",
        "",
        "# For demonstration, we'll use synthetic data",
        "raw_data <- generate_synthetic_data(n_cells = 5, n_time = 1000)",
        ""
      )
    } else NULL,
    "# Step 1: Calcium Correction",
    "corrected_data <- calcium_correction(raw_data$calcium_traces)",
    "",
    "# Step 2: Spike Inference",
    "spike_results <- infer_spikes(corrected_data$Cell_1, method = 'oasis')",
    "",
    "# Step 3: Statistical Analysis",
    "statistical_analysis <- comprehensive_statistical_analysis(",
    "  raw_data$calcium_traces,",
    "  corrected_data,",
    "  spike_results",
    ")",
    "",
    "# Step 4: Generate Plots",
    "main_plot <- plot_cell_trace(corrected_data, 'Cell_1')",
    "",
    "# Step 5: Quality Control",
    "qc_metrics <- calculate_correction_quality(raw_data$calcium_traces, corrected_data)",
    "",
    "# Print Results",
    "cat('Analysis completed successfully!\\n')",
    "cat('Number of spikes detected:', sum(spike_results$spikes), '\\n')",
    "cat('Quality score:', qc_metrics$overall_quality, '\\n')",
    ""
  )
}

#' Export Analysis as Jupyter Notebook
#'
#' Export the analysis pipeline as a Jupyter notebook.
#'
#' @param analysis_results List containing analysis results
#' @param output_file Path for output notebook (default: "calcium_analysis.ipynb")
#' @param kernel R kernel name (default: "ir")
#' @param ... Additional arguments
#' @return Path to generated notebook
#' @export
export_jupyter_notebook <- function(analysis_results, output_file = "calcium_analysis.ipynb", kernel = "ir", ...) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required for notebook export")
  }
  
  # Create temporary Rmd file
  rmd_content <- create_notebook_rmd(analysis_results)
  temp_rmd <- tempfile(fileext = ".Rmd")
  writeLines(rmd_content, temp_rmd)
  
  # Convert to notebook
  knitr::pandoc(temp_rmd, output = output_file, format = "ipynb")
  
  # Clean up
  unlink(temp_rmd)
  
  message("Jupyter notebook exported: ", output_file)
  return(output_file)
}

#' Create Notebook RMarkdown Content
#'
#' Generate RMarkdown content for Jupyter notebook.
#'
#' @param analysis_results List containing analysis results
#' @return RMarkdown content as character vector
#' @export
create_notebook_rmd <- function(analysis_results) {
  c(
    "---",
    "title: Calcium Imaging Analysis",
    "jupyter: ir",
    "---",
    "",
    "```{r}",
    "library(CaImagingAnalysisFr)",
    "library(ggplot2)",
    "```",
    "",
    "## Data Loading",
    "",
    "```{r}",
    "# Load your data here",
    "raw_data <- generate_synthetic_data(n_cells = 5, n_time = 1000)",
    "head(raw_data$calcium_traces[, 1:5])",
    "```",
    "",
    "## Calcium Correction",
    "",
    "```{r}",
    "corrected_data <- calcium_correction(raw_data$calcium_traces)",
    "```",
    "",
    "## Spike Detection",
    "",
    "```{r}",
    "spike_results <- infer_spikes(corrected_data$Cell_1, method = 'oasis')",
    "cat('Number of spikes:', sum(spike_results$spikes))",
    "```",
    "",
    "## Visualization",
    "",
    "```{r}",
    "plot_cell_trace(corrected_data, 'Cell_1')",
    "```",
    ""
  )
}

#' Generate Interactive Dashboard
#'
#' Create an interactive Shiny dashboard for exploring analysis results.
#'
#' @param analysis_results List containing analysis results
#' @param output_dir Directory to save dashboard files (default: "dashboard")
#' @param ... Additional arguments
#' @return Path to dashboard directory
#' @export
generate_interactive_dashboard <- function(analysis_results, output_dir = "dashboard", ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required for dashboard generation. Install with: install.packages('shiny')")
  }
  
  # Create dashboard directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create app.R file
  app_content <- create_dashboard_app(analysis_results)
  writeLines(app_content, file.path(output_dir, "app.R"))
  
  # Create global.R file
  global_content <- create_dashboard_global(analysis_results)
  writeLines(global_content, file.path(output_dir, "global.R"))
  
  message("Interactive dashboard created in: ", output_dir)
  message("Run 'shiny::runApp('", output_dir, "')' to launch the dashboard")
  
  return(output_dir)
}

#' Create Dashboard App Content
#'
#' Generate Shiny app content for the dashboard.
#'
#' @param analysis_results List containing analysis results
#' @return App content as character vector
#' @export
create_dashboard_app <- function(analysis_results) {
  c(
    "library(shiny)",
    "library(ggplot2)",
    "library(CaImagingAnalysisFr)",
    "",
    "ui <- fluidPage(",
    "  titlePanel('Calcium Imaging Analysis Dashboard'),",
    "  ",
    "  sidebarLayout(",
    "    sidebarPanel(",
    "      selectInput('plot_type', 'Plot Type:',",
    "                  choices = c('Cell Traces', 'Spike Detection', 'Quality Metrics')),",
    "      selectInput('cell_id', 'Cell ID:',",
    "                  choices = if (!is.null(analysis_results$cell_ids)) analysis_results$cell_ids else 'Cell_1'),",
    "      sliderInput('time_range', 'Time Range:',",
    "                  min = 1, max = 1000, value = c(1, 1000))",
    "    ),",
    "    ",
    "    mainPanel(",
    "      plotOutput('main_plot'),",
    "      verbatimTextOutput('summary_stats')",
    "    )",
    "  )",
    ")",
    "",
    "server <- function(input, output) {",
    "  ",
    "  output$main_plot <- renderPlot({",
    "    # Generate plot based on user selection",
    "    if (input$plot_type == 'Cell Traces') {",
    "      plot_cell_trace(analysis_results$corrected_data, input$cell_id)",
    "    } else if (input$plot_type == 'Spike Detection') {",
    "      # Add spike detection plot",
    "      plot(1:100, type = 'l', main = 'Spike Detection')",
    "    } else {",
    "      # Add quality metrics plot",
    "      plot(1:100, type = 'l', main = 'Quality Metrics')",
    "    }",
    "  })",
    "  ",
    "  output$summary_stats <- renderPrint({",
    "    cat('Summary Statistics\\n')",
    "    cat('Selected cell:', input$cell_id, '\\n')",
    "    cat('Time range:', input$time_range[1], '-', input$time_range[2], '\\n')",
    "  })",
    "}",
    "",
    "shinyApp(ui = ui, server = server)"
  )
}

#' Create Dashboard Global Content
#'
#' Generate global.R content for the dashboard.
#'
#' @param analysis_results List containing analysis results
#' @return Global content as character vector
#' @export
create_dashboard_global <- function(analysis_results) {
  c(
    "# Global variables for Shiny dashboard",
    "library(CaImagingAnalysisFr)",
    "library(ggplot2)",
    "",
    "# Load analysis results",
    "analysis_results <- readRDS('analysis_results.rds')",
    "",
    "# Helper functions",
    "format_number <- function(x, digits = 3) {",
    "  round(x, digits)",
    "}",
    ""
  )
}

#' Generate Analysis Summary
#'
#' Create a concise summary of the analysis results.
#'
#' @param analysis_results List containing analysis results
#' @param output_file Path for output summary file (default: "analysis_summary.txt")
#' @param ... Additional arguments
#' @return Path to generated summary
#' @export
generate_analysis_summary <- function(analysis_results, output_file = "analysis_summary.txt", ...) {
  # Create summary content
  summary_content <- c(
    "=== Calcium Imaging Analysis Summary ===",
    "",
    paste("Analysis Date:", Sys.Date()),
    paste("Analysis Time:", Sys.time()),
    "",
    "DATA SUMMARY:",
    paste("  Number of cells:", if (!is.null(analysis_results$n_cells)) analysis_results$n_cells else "N/A"),
    paste("  Number of time points:", if (!is.null(analysis_results$n_time)) analysis_results$n_time else "N/A"),
    "",
    "SPIKE DETECTION:",
    paste("  Total spikes detected:", if (!is.null(analysis_results$total_spikes)) analysis_results$total_spikes else "N/A"),
    paste("  Mean spike rate:", if (!is.null(analysis_results$mean_spike_rate)) paste(round(analysis_results$mean_spike_rate, 3), "Hz") else "N/A"),
    "",
    "QUALITY CONTROL:",
    paste("  Overall quality score:", if (!is.null(analysis_results$quality_score)) round(analysis_results$quality_score, 3) else "N/A"),
    paste("  Passed QC checks:", if (!is.null(analysis_results$qc_passed)) analysis_results$qc_passed else "N/A"),
    "",
    "RECOMMENDATIONS:",
    if (!is.null(analysis_results$recommendations)) analysis_results$recommendations else "  No specific recommendations",
    "",
    "=== End of Summary ==="
  )
  
  # Write to file
  writeLines(summary_content, output_file)
  
  message("Analysis summary generated: ", output_file)
  return(output_file)
} 