#' Launch interactive viewer for calcium traces
#'
#' Provides a Shiny interface with adjustable parameters.
#'
#' @param raw_df Data frame containing raw traces
#' @export
launch_interactive_viewer <- function(raw_df) {
  if (!requireNamespace("shiny", quietly = TRUE) ||
      !requireNamespace("plotly", quietly = TRUE)) {
    stop("Packages 'shiny' and 'plotly' are required for the interactive viewer")
  }
  cells <- names(raw_df)[grepl("^Cell_", names(raw_df))]
  ui <- shiny::fluidPage(
    shiny::sliderInput("span", "Loess span", min = 0.1, max = 1, value = 0.45, step = 0.05),
    shiny::selectInput("cell", "Cell", choices = cells),
    plotly::plotlyOutput("trace")
  )
  server <- function(input, output, session) {
    corrected <- shiny::reactive(calcium_correction(raw_df, span = input$span))
    output$trace <- plotly::renderPlotly({
      df <- corrected()
      plotly::ggplotly(plot_cell_trace(df, input$cell))
    })
  }
  shiny::shinyApp(ui, server)
}
