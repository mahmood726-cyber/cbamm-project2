
library(shiny)
library(cbamm)

ui <- fluidPage(
  titlePanel("CBAMM Simple Test"),
  sidebarLayout(
    sidebarPanel(
      actionButton("generate", "Generate Data"),
      br(), br(),
      actionButton("analyze", "Run Analysis")
    ),
    mainPanel(
      verbatimTextOutput("data_output"),
      verbatimTextOutput("analysis_output")
    )
  )
)

server <- function(input, output, session) {
  values <- reactiveValues(data = NULL, results = NULL)
  
  observeEvent(input$generate, {
    values$data <- simulate_cbamm_data()
    showNotification("Data generated!", type = "success")
  })
  
  observeEvent(input$analyze, {
    req(values$data)
    values$results <- cbamm(values$data)
    showNotification("Analysis complete!", type = "success")
  })
  
  output$data_output <- renderPrint({
    if (!is.null(values$data)) {
      cat("Data: ", nrow(values$data), "studies\n")
      head(values$data, 3)
    }
  })
  
  output$analysis_output <- renderPrint({
    if (!is.null(values$results)) {
      values$results$summary
    }
  })
}

shinyApp(ui = ui, server = server)

