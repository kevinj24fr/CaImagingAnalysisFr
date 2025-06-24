#' Launch Interactive Viewer for Calcium Traces
#' 
#' Provides a Shiny interface with adjustable parameters for exploring
#' calcium imaging data interactively.
#' 
#' @param raw_df Data frame containing raw traces
#' @param port Port number for the Shiny app (default: NULL for auto-assignment)
#' @param host Host address (default: "127.0.0.1")
#' @param launch_browser Whether to launch browser automatically (default: TRUE)
#' 
#' @return Shiny app object (invisibly)
#' 
#' @examples
#' # Basic usage
#' raw <- generate_synthetic_data(3, 500)
#' launch_interactive_viewer(raw)
#' 
#' # Custom port
#' launch_interactive_viewer(raw, port = 8080)
#' 
#' @export
launch_interactive_viewer <- function(raw_df, 
                                     port = NULL,
                                     host = "127.0.0.1",
                                     launch_browser = TRUE) {
  
  # Check required packages
  required_packages <- c("shiny", "plotly")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("Required packages not available: ", paste(missing_packages, collapse = ", "),
         "\nPlease install with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))")
  }
  
  # Validate input data
  validate_data_frame(raw_df)
  
  # Get cell columns
  config <- get_config()
  cells <- names(raw_df)[grepl(config$cell_pattern, names(raw_df))]
  
  if (length(cells) == 0) {
    stop("No cell columns found in data")
  }
  
  # UI definition
  ui <- shiny::fluidPage(
    shiny::titlePanel("Calcium Imaging Analysis Viewer"),
    
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        
        # Data info
        shiny::h4("Data Information"),
        shiny::verbatimTextOutput("data_info"),
        
        shiny::hr(),
        
        # Correction parameters
        shiny::h4("Correction Parameters"),
        shiny::selectInput("correction_method", "Correction Method:",
                          choices = c("modern", "legacy"),
                          selected = "modern"),
        shiny::sliderInput("span", "Loess Span:",
                          min = 0.01, max = 1.0, value = 0.45, step = 0.01),
        shiny::checkboxInput("normalize", "Normalize Traces", value = TRUE),
        
        shiny::hr(),
        
        # Cell selection
        shiny::h4("Cell Selection"),
        shiny::selectInput("cell", "Select Cell:",
                          choices = cells,
                          selected = cells[1]),
        
        shiny::hr(),
        
        # Spike inference parameters
        shiny::h4("Spike Inference"),
        shiny::selectInput("spike_method", "Method:",
                          choices = c("oasis", "caiman", "suite2p"),
                          selected = "oasis"),
        shiny::checkboxInput("show_spikes", "Show Spikes", value = TRUE),
        shiny::checkboxInput("show_deconvolved", "Show Deconvolved", value = TRUE),
        
        shiny::hr(),
        
        # Export options
        shiny::h4("Export"),
        shiny::downloadButton("download_plot", "Download Plot"),
        shiny::downloadButton("download_data", "Download Corrected Data")
      ),
      
      shiny::mainPanel(
        width = 9,
        
        # Main plot
        shiny::h3("Calcium Trace Visualization"),
        plotly::plotlyOutput("trace_plot", height = "500px"),
        
        # Statistics panel
        shiny::h3("Trace Statistics"),
        shiny::fluidRow(
          shiny::column(6, shiny::tableOutput("trace_stats")),
          shiny::column(6, shiny::tableOutput("spike_stats"))
        )
      )
    )
  )
  
  # Server logic
  server <- function(input, output, session) {
    
    # Data info
    output$data_info <- shiny::renderPrint({
      cat("Cells:", length(cells), "\n")
      cat("Time points:", nrow(raw_df), "\n")
      cat("Background regions:", length(names(raw_df)[grepl(config$background_pattern, names(raw_df))]), "\n")
    })
    
    # Reactive corrected data
    corrected_data <- shiny::reactive({
      shiny::req(input$correction_method, input$span, input$normalize)
      
      tryCatch({
        calcium_correction(
          raw_df, 
          method = input$correction_method,
          span = input$span,
          normalize = input$normalize,
          verbose = FALSE
        )
      }, error = function(e) {
        shiny::showNotification(
          paste("Correction failed:", e$message), 
          type = "error"
        )
        return(NULL)
      })
    })
    
    # Main plot
    output$trace_plot <- plotly::renderPlotly({
      shiny::req(corrected_data(), input$cell)
      
      tryCatch({
        p <- plot_cell_trace(
          corrected_data(),
          input$cell,
          method = input$spike_method,
          show_spikes = input$show_spikes,
          show_deconvolved = input$show_deconvolved,
          verbose = FALSE
        )
        
        plotly::ggplotly(p, tooltip = c("x", "y", "color"))
      }, error = function(e) {
        shiny::showNotification(
          paste("Plot generation failed:", e$message), 
          type = "error"
        )
        return(NULL)
      })
    })
    
    # Trace statistics
    output$trace_stats <- shiny::renderTable({
      shiny::req(corrected_data(), input$cell)
      
      trace <- corrected_data()[[input$cell]]
      stats <- data.frame(
        Statistic = c("Mean", "SD", "Min", "Max", "Range"),
        Value = c(
          round(mean(trace, na.rm = TRUE), 4),
          round(sd(trace, na.rm = TRUE), 4),
          round(min(trace, na.rm = TRUE), 4),
          round(max(trace, na.rm = TRUE), 4),
          round(max(trace, na.rm = TRUE) - min(trace, na.rm = TRUE), 4)
        )
      )
      stats
    }, striped = TRUE, hover = TRUE)
    
    # Spike statistics
    output$spike_stats <- shiny::renderTable({
      shiny::req(corrected_data(), input$cell, input$show_spikes)
      
      if (!input$show_spikes) {
        return(data.frame(Note = "Spike detection disabled"))
      }
      
      tryCatch({
        trace <- corrected_data()[[input$cell]]
        spikes <- infer_spikes(trace, method = input$spike_method, verbose = FALSE)
        
        n_spikes <- sum(spikes$spike > 0, na.rm = TRUE)
        spike_rate <- n_spikes / length(trace) * 100  # per 100 frames
        
        stats <- data.frame(
          Statistic = c("Total Spikes", "Spike Rate (%)", "Mean Amplitude"),
          Value = c(
            n_spikes,
            round(spike_rate, 2),
            ifelse(n_spikes > 0, 
                   round(mean(trace[spikes$spike > 0], na.rm = TRUE), 4),
                   "N/A")
          )
        )
        stats
      }, error = function(e) {
        data.frame(Error = paste("Spike analysis failed:", e$message))
      })
    }, striped = TRUE, hover = TRUE)
    
    # Download handlers
    output$download_plot <- shiny::downloadHandler(
      filename = function() {
        paste0("calcium_trace_", input$cell, "_", Sys.Date(), ".png")
      },
      content = function(file) {
        shiny::req(corrected_data(), input$cell)
        
        p <- plot_cell_trace(
          corrected_data(),
          input$cell,
          method = input$spike_method,
          show_spikes = input$show_spikes,
          show_deconvolved = input$show_deconvolved,
          verbose = FALSE
        )
        
        ggplot2::ggsave(file, p, width = 10, height = 6, dpi = 300)
      }
    )
    
    output$download_data <- shiny::downloadHandler(
      filename = function() {
        paste0("corrected_data_", Sys.Date(), ".csv")
      },
      content = function(file) {
        shiny::req(corrected_data())
        utils::write.csv(corrected_data(), file, row.names = FALSE)
      }
    )
  }
  
  # Launch app
  app <- shiny::shinyApp(ui = ui, server = server)
  
  if (launch_browser) {
    shiny::runApp(app, port = port, host = host, launch.browser = TRUE)
  } else {
    shiny::runApp(app, port = port, host = host, launch.browser = FALSE)
  }
  
  invisible(app)
}
