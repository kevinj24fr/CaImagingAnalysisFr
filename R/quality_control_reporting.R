#' Batch Quality Control for Calcium Imaging Data
#'
#' Run quality control checks on a batch of data frames.
#'
#' @param data_list List of data frames (e.g., multiple recordings)
#' @param verbose Show progress messages
#' @return Data frame summarizing QC results
#' @export
batch_quality_control <- function(data_list, verbose = TRUE) {
  results <- lapply(seq_along(data_list), function(i) {
    if (verbose) message("QC for dataset ", i)
    df <- data_list[[i]]
    qc <- tryCatch({
      validate_data_frame(df)
      TRUE
    }, error = function(e) {
      FALSE
    })
    snr <- tryCatch({
      cell_cols <- grep(get_config()$cell_pattern, names(df), value = TRUE)
      mean(sapply(cell_cols, function(col) mean(df[[col]], na.rm = TRUE) / sd(df[[col]], na.rm = TRUE)), na.rm = TRUE)
    }, error = function(e) NA)
    data.frame(dataset = i, qc_pass = qc, mean_snr = snr)
  })
  do.call(rbind, results)
}

#' Flag Outliers in Calcium Imaging Data
#'
#' Flag outlier cells or time points using robust statistics.
#'
#' @param data Data frame or matrix (cells x time)
#' @param method Outlier detection method ("iqr", "zscore", "mad")
#' @param threshold Outlier threshold
#' @return Logical matrix indicating outliers
#' @export
flag_outliers <- function(data, method = "iqr", threshold = 1.5) {
  outlier_mat <- matrix(FALSE, nrow = nrow(data), ncol = ncol(data))
  for (i in 1:nrow(data)) {
    outlier_mat[i, ] <- detect_outliers(as.numeric(data[i, ]), method = method, threshold = threshold)
  }
  outlier_mat
}

#' Generate QC Report (RMarkdown)
#'
#' Generate an RMarkdown quality control report for a dataset.
#'
#' @param data Data frame
#' @param output_file Output HTML file
#' @param title Report title
#' @param open Whether to open the report after generation
#' @return Path to the generated report
#' @export
generate_qc_report <- function(data, output_file = "qc_report.html", title = "Quality Control Report", open = TRUE) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required for report generation. Install with install.packages('rmarkdown')")
  }
  temp_rmd <- tempfile(fileext = ".Rmd")
  rmd_content <- paste0(
    "---\ntitle: '", title, "'\noutput: html_document\n---\n",
    "\n## Data Summary\n",
    "```{r}\nstr(data)\nsummary(data)\n```\n",
    "\n## Quality Control\n",
    "```{r}\nvalidate_data_frame(data)\n```\n",
    "\n## Signal-to-Noise Ratio\n",
    "```{r}\ncell_cols <- grep(get_config()$cell_pattern, names(data), value = TRUE)\nsnr <- sapply(cell_cols, function(col) mean(data[[col]], na.rm = TRUE) / sd(data[[col]], na.rm = TRUE))\nbarplot(snr, main = 'Signal-to-Noise Ratio by Cell', ylab = 'SNR')\n```\n",
    "\n## Outlier Detection\n",
    "```{r}\noutlier_mat <- flag_outliers(data)\nimage(t(outlier_mat), main = 'Outlier Matrix', xlab = 'Time', ylab = 'Cell')\n```\n"
  )
  writeLines(rmd_content, temp_rmd)
  rmarkdown::render(temp_rmd, output_file = output_file, params = list(data = data), envir = new.env())
  if (open) browseURL(output_file)
  output_file
}

#' Launch Interactive QC Report (Shiny)
#'
#' Launch a Shiny app for interactive quality control exploration.
#'
#' @param data Data frame
#' @export
launch_qc_shiny <- function(data) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required for interactive QC. Install with install.packages('shiny')")
  }
  shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::titlePanel("Interactive Quality Control"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput("cell", "Cell", choices = grep(get_config()$cell_pattern, names(data), value = TRUE)),
          shiny::sliderInput("threshold", "Outlier Threshold", min = 0.5, max = 5, value = 1.5, step = 0.1)
        ),
        shiny::mainPanel(
          shiny::plotOutput("tracePlot"),
          shiny::plotOutput("outlierPlot")
        )
      )
    ),
    server = function(input, output) {
      output$tracePlot <- shiny::renderPlot({
        plot(data[[input$cell]], type = 'l', main = paste("Trace for", input$cell), ylab = "Fluorescence")
      })
      output$outlierPlot <- shiny::renderPlot({
        outliers <- detect_outliers(data[[input$cell]], threshold = input$threshold)
        plot(data[[input$cell]], type = 'l', main = paste("Outliers for", input$cell), ylab = "Fluorescence")
        points(which(outliers), data[[input$cell]][outliers], col = 'red', pch = 19)
      })
    }
  )
} 