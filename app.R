#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(nlme)
library(ggplot2)
library(yaml)
library(openxlsx)
library(dplyr)
library(knitr)    # For creating HTML tables
library(kableExtra) # For additional table styling (optional)

# Source functions
source("functions.R")


# Force rsconnect to detect openxlsx by using a dummy function
forceDependencyOpenxlsx <- function() {
  openxlsx::getNamedRegions  # Reference to a function in openxlsx
}

dummy_openxlsx <- openxlsx::getNamedRegions

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' is required but not installed.")
}



# Define UI
ui <- fluidPage(
  titlePanel("Photoelectric Effect Regression Analysis"),
  sidebarLayout(
    sidebarPanel(
      downloadButton('downloadTemplate', "Download Excel Template"),
      br(),
      fileInput('file1', 'Upload Excel Data File',
                accept = c(".xlsx")),

      radioButtons("energy_units", "Energy Units:",
                   choices = c("eV", "J"),
                   selected = "eV",
                   inline = TRUE),

      actionButton("run_analysis", "Run Analysis")
    ),
    mainPanel(
      plotOutput("scatterPlot"),
      uiOutput("modelSummary")  # Changed from tableOutput to uiOutput
    )
  )
)

# Define Server
server <- function(input, output) {
  PLANCK_EV <- 1.0
  PLANCK_J  <- 1.602176634e-19


  multiplier <- reactive({
    if (input$energy_units == "eV") PLANCK_EV else PLANCK_J
  })
  data_input <- eventReactive(input$run_analysis, {
    req(input$file1)

    # Define expected dimensions
    expected_rows <- 5
    expected_cols <- 7
    expected_indices <- 3:9

    # Initialize full data frame with NA
    df_raw <- as.data.frame(matrix(NA, nrow = expected_rows, ncol = expected_cols))

    # Read only what exists in the file
    df_partial <- tryCatch(
      openxlsx::read.xlsx(
        input$file1$datapath,
        sheet = 1,
        rows = 3:7,
        cols = expected_indices,
        colNames = FALSE
      ),
      error = function(e) {
        stop("Failed to read Excel file. Please make sure it has the expected structure.")
      }
    )

    # Copy existing data into the preallocated frame
    ncols_available <- min(ncol(df_partial), expected_cols)
    df_raw[, 1:ncols_available] <- df_partial[, 1:ncols_available]


    # Frequency = c / λ where λ is expressed in nm in the worksheet
    x <- 299792458/as.numeric(df_raw[,1]*1.0E-9)

    # σ_x = σ_λ/λ * x
    sigma_x <- as.numeric(df_raw[,2])/as.numeric(df_raw[,1])*x

    y_data <- suppressWarnings({
      data.frame(
        lapply(df_raw[, 3:7], function(x) as.numeric(as.character(x)))
      )
    })
    y_raw <- process_voltages(y_data)



    full_df <- data.frame(
      x = x,
      y = as.numeric(y_raw$y),
      sigma_x = sigma_x,
      sigma_y = as.numeric(y_raw$sigma_y)
    )

    na.omit(full_df)
  })

  # Update model_result()
  model_result <- reactive({
    df <- data_input()
    colnames(df) <- c("x", "y","sigma_x","sigma_y")
      df <- df %>%
        mutate(
          x = as.numeric(x),
          y = as.numeric(y),
          sigma_x = as.numeric(sigma_x),
          sigma_y = as.numeric(sigma_y)
        )

    # Perform Iterated GLS Regression
    result <- iterated_gls(df)
    return(result)
  })

  output$modelSummary <- renderUI({
    result <- model_result()
    model <- result$model

    tt <- summary(model)$tTable
    dof <- model$dims$N - model$dims$p

    est <- tt[, "Value"]
    se  <- tt[, "Std.Error"]
    t_crit <- qt(0.975, dof)
    conf.low  <- est - t_crit * se
    conf.high <- est + t_crit * se

    coef_table <- data.frame(
      `term` = rownames(tt),
      `Estimated Value` = est,
      `Standard Error` = se,
      `DF` = dof,
      `Statistic` = tt[, "t-value"],
      `P-value` = tt[, "p-value"],
      `Lower 95% CI` = conf.low,
      `Upper 95% CI` = conf.high
    )

    # Select relevant columns
    coef_table <- coef_table %>%
      select(all_of(c(
        "term", "Estimated.Value", "Standard.Error",
        "Lower.95..CI", "Upper.95..CI", "P.value"
      )))


    coef_table$term <- ifelse(
      coef_table$term == "(Intercept)",
      paste0("Intercept (", input$energy_units, ")"),
      paste0("Slope (", input$energy_units, "·s)")
    )

    coef_table[, 2:5] <- coef_table[, 2:5] * multiplier()

    coef_table[, 2:6] <- lapply(coef_table[, 2:6], smart_format)
    coef_table[, 6] <- formatC(as.numeric(coef_table[, 6]), format = "e", digits = 1)

    colnames(coef_table) <- c("Term", "Estimated Value", "Standard Error",
                              "Lower 95% CI", "Upper 95% CI", "P-value")


    # Generate HTML table using knitr::kable
    table_html <- kable(coef_table[,1:6], format = "html", table.attr = "class='table table-striped'", row.names=FALSE) %>%
      kable_styling(full_width = FALSE)
    HTML(paste0(table_html, "<br><p><b>Reminder:</b> In regression analysis, each p-value tests whether the corresponding coefficient is zero. That is, it assesses whether the predictor has a statistically significant effect on the response variable.</p>"))
    # Return the HTML code
    #HTML(table_html)
  })

  # Update output$scatterPlot
  output$scatterPlot <- renderPlot({
    df <- data_input()

    df <- df %>%
      mutate(
        x = as.numeric(x),
        y = as.numeric(multiplier()*y),
        sigma_x = as.numeric(sigma_x),
        sigma_y = as.numeric(multiplier()*sigma_y)
      )

    result <- model_result()

    # Create Plot
    p <- create_plot(df, result$model, units=input$energy_units, multiplier=multiplier())
    print(p)
  })

  output$downloadTemplate <- downloadHandler(
    filename = function() {
      "PhotoElectricEffect.xlsx"
    },
    content = function(file) {
      # Copy the template file from the 'www' directory to the temp file
      file.copy("www/PhotoElectricEffect.xlsx", file)
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )


}

# Run the app
shinyApp(ui = ui, server = server)
