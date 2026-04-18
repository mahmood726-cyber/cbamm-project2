
# CBAMM Shiny Dashboard (Fixed)
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(cbamm)

# UI
ui <- dashboardPage(
  dashboardHeader(title = "CBAMM Meta-Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data", tabName = "data", icon = icon("database")),
      menuItem("Analysis", tabName = "analysis", icon = icon("calculator")),
      menuItem("Results", tabName = "results", icon = icon("chart-line")),
      menuItem("Plots", tabName = "plots", icon = icon("chart-bar"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Data Tab
      tabItem(
        tabName = "data",
        h2("Data Input"),
        fluidRow(
          box(
            title = "Load Data",
            width = 12,
            radioButtons("data_source", "Choose data source:",
                        choices = c("Generate Example" = "example",
                                  "Upload CSV" = "upload")),
            
            conditionalPanel(
              condition = "input.data_source == 'upload'",
              fileInput("file", "Choose CSV File", accept = ".csv")
            ),
            
            conditionalPanel(
              condition = "input.data_source == 'example'",
              selectInput("example_type", "Example type:",
                         choices = c("Standard" = "standard",
                                   "IPD" = "ipd",
                                   "DTA" = "dta"))
            ),
            
            actionButton("load_data", "Load Data", class = "btn-primary"),
            br(), br(),
            DT::dataTableOutput("data_table")
          )
        )
      ),
      
      # Analysis Tab
      tabItem(
        tabName = "analysis",
        h2("Run Analysis"),
        fluidRow(
          box(
            title = "Analysis Settings",
            width = 12,
            
            selectInput("analysis_type", "Analysis Type:",
                       choices = c("Standard" = "standard",
                                 "IPD" = "ipd")),
            
            actionButton("run_analysis", "Run Analysis", 
                        class = "btn-success", icon = icon("play")),
            
            br(), br(),
            verbatimTextOutput("analysis_output")
          )
        )
      ),
      
      # Results Tab
      tabItem(
        tabName = "results",
        h2("Results"),
        fluidRow(
          valueBoxOutput("effect_size"),
          valueBoxOutput("heterogeneity"),
          valueBoxOutput("n_studies")
        ),
        fluidRow(
          box(
            title = "Summary",
            width = 12,
            verbatimTextOutput("results_summary")
          )
        )
      ),
      
      # Plots Tab
      tabItem(
        tabName = "plots",
        h2("Visualizations"),
        fluidRow(
          box(
            title = "Forest Plot",
            width = 12,
            plotlyOutput("forest_plot", height = "500px")
          )
        ),
        fluidRow(
          box(
            title = "Funnel Plot",
            width = 6,
            plotlyOutput("funnel_plot", height = "400px")
          ),
          box(
            title = "Diagnostic Plot",
            width = 6,
            plotOutput("diagnostic_plot", height = "400px")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values
  values <- reactiveValues(
    data = NULL,
    results = NULL
  )
  
  # Load data
  observeEvent(input$load_data, {
    if (input$data_source == "example") {
      if (input$example_type == "standard") {
        values$data <- simulate_cbamm_data()
      } else if (input$example_type == "ipd") {
        # Fix: simulate_ipd_data returns a list, convert to single data frame
        ipd_list <- simulate_ipd_data()
        if (is.list(ipd_list) && !is.data.frame(ipd_list)) {
          # Combine all data frames in the list
          values$data <- do.call(rbind, ipd_list)
        } else {
          values$data <- ipd_list
        }
      } else if (input$example_type == "dta") {
        values$data <- simulate_dta_data()
      }
      # Fix: Use correct notification syntax
      showNotification("Example data loaded!", duration = 3)
    } else if (input$data_source == "upload" && !is.null(input$file)) {
      values$data <- read.csv(input$file$datapath)
      showNotification("File uploaded successfully!", duration = 3)
    }
  })
  
  # Display data
  output$data_table <- DT::renderDataTable({
    req(values$data)
    DT::datatable(values$data, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Run analysis
  observeEvent(input$run_analysis, {
    req(values$data)
    
    withProgress(message = "Running analysis...", {
      if (input$analysis_type == "standard") {
        values$results <- tryCatch({
          cbamm(values$data)
        }, error = function(e) {
          showNotification(paste("Error:", e$message), duration = 5)
          NULL
        })
      } else if (input$analysis_type == "ipd") {
        values$results <- tryCatch({
          ipd_meta_analysis(ipd_data = values$data)
        }, error = function(e) {
          showNotification(paste("Error:", e$message), duration = 5)
          NULL
        })
      }
    })
    
    if (!is.null(values$results)) {
      showNotification("Analysis complete!", duration = 3)
    }
  })
  
  # Analysis output
  output$analysis_output <- renderPrint({
    req(values$results)
    print(values$results)
  })
  
  # Value boxes
  output$effect_size <- renderValueBox({
    val <- if(is.null(values$results)) "---" else {
      if(!is.null(values$results$summary$estimate)) {
        round(values$results$summary$estimate, 3)
      } else if(!is.null(values$results$core_results$estimate)) {
        round(values$results$core_results$estimate, 3)
      } else if(!is.null(values$results$pooled_effect)) {
        round(values$results$pooled_effect, 3)
      } else "---"
    }
    valueBox(val, "Effect Size", icon = icon("chart-line"), color = "green")
  })
  
  output$heterogeneity <- renderValueBox({
    val <- if(is.null(values$results)) "---" else {
      if(!is.null(values$results$summary$I2)) {
        paste0(round(values$results$summary$I2, 1), "%")
      } else if(!is.null(values$results$core_results$I2)) {
        paste0(round(values$results$core_results$I2, 1), "%")
      } else "---"
    }
    valueBox(val, "I² Heterogeneity", icon = icon("random"), color = "yellow")
  })
  
  output$n_studies <- renderValueBox({
    val <- if(is.null(values$data)) "---" else {
      if("study_id" %in% names(values$data)) {
        length(unique(values$data$study_id))
      } else {
        nrow(values$data)
      }
    }
    valueBox(val, "Studies", icon = icon("book"), color = "blue")
  })
  
  # Results summary
  output$results_summary <- renderPrint({
    req(values$results)
    if (!is.null(values$results$summary)) {
      values$results$summary
    } else {
      str(values$results)
    }
  })
  
  # Forest plot
  output$forest_plot <- renderPlotly({
    req(values$data)
    
    # Find effect and SE columns
    effect_col <- names(values$data)[grep("effect|yi|outcome", names(values$data), ignore.case = TRUE)[1]]
    se_col <- names(values$data)[grep("se|error|std", names(values$data), ignore.case = TRUE)[1]]
    
    if (!is.na(effect_col) && !is.na(se_col)) {
      # Get unique studies if study_id exists
      if ("study_id" %in% names(values$data)) {
        plot_data <- aggregate(values$data[[effect_col]], 
                              by = list(study = values$data$study_id), 
                              FUN = mean, na.rm = TRUE)
        plot_data$se <- aggregate(values$data[[se_col]], 
                                 by = list(study = values$data$study_id), 
                                 FUN = mean, na.rm = TRUE)$x
        names(plot_data)[2] <- "effect"
      } else {
        plot_data <- data.frame(
          study = 1:min(nrow(values$data), 20),
          effect = values$data[[effect_col]][1:min(nrow(values$data), 20)],
          se = values$data[[se_col]][1:min(nrow(values$data), 20)]
        )
      }
      
      plot_ly(plot_data) %>%
        add_trace(
          type = "scatter",
          mode = "markers",
          x = ~effect,
          y = ~study,
          name = "Studies",
          marker = list(size = 10, color = "blue"),
          error_x = list(
            array = ~se * 1.96,
            color = "blue"
          )
        ) %>%
        layout(
          title = "Forest Plot",
          xaxis = list(title = "Effect Size"),
          yaxis = list(title = "Study", autorange = "reversed"),
          hovermode = "closest"
        )
    } else {
      plot_ly() %>%
        layout(title = "No effect/SE columns found")
    }
  })
  
  # Funnel plot
  output$funnel_plot <- renderPlotly({
    req(values$data)
    
    effect_col <- names(values$data)[grep("effect|yi", names(values$data), ignore.case = TRUE)[1]]
    se_col <- names(values$data)[grep("se|error", names(values$data), ignore.case = TRUE)[1]]
    
    if (!is.na(effect_col) && !is.na(se_col)) {
      plot_ly(
        x = values$data[[effect_col]],
        y = values$data[[se_col]],
        type = "scatter",
        mode = "markers",
        marker = list(size = 10),
        name = "Studies"
      ) %>%
        layout(
          title = "Funnel Plot",
          xaxis = list(title = "Effect Size"),
          yaxis = list(title = "Standard Error", autorange = "reversed")
        )
    } else {
      plot_ly() %>%
        layout(title = "No effect/SE columns found")
    }
  })
  
  # Diagnostic plot
  output$diagnostic_plot <- renderPlot({
    req(values$data)
    
    par(mfrow = c(2, 2))
    
    # Find numeric columns for plotting
    numeric_cols <- sapply(values$data, is.numeric)
    if (sum(numeric_cols) > 0) {
      first_numeric <- values$data[, numeric_cols][, 1]
      hist(first_numeric, main = "Distribution", xlab = "Value")
      qqnorm(first_numeric)
      qqline(first_numeric)
      boxplot(first_numeric, main = "Boxplot")
      plot(density(first_numeric, na.rm = TRUE), main = "Density")
    }
    
    par(mfrow = c(1, 1))
  })
}

# Run app
shinyApp(ui = ui, server = server)

