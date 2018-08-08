# ui.R
# 
# Author: Gregory Brownson
#
# Description: frontend portion of the Shiny 
#              interface for loading data and
#              selecting location and scale
#              parameters

library(DT)
library(shiny)

# Some custom CSS code for formatting panels
CSS.format1 <- 
"
  /* Smaller font for preformatted text */
  pre, table.table {
    font-size: smaller;
  }

  body {
    min-height: 2000px;
  }

  .option-group {
    border: 1px solid #ccc;
    border-radius: 6px;
    padding: 0px 5px;
    margin: 5px -10px;
    background-color: #f5f5f5;
  }

  .option-header {
    color: #79d;
    text-transform: uppercase;
    margin-bottom: 5px;
  }
"

# JavaScript function to display slider input values in scientific notation
JS.log10 <-
"
// Function to compute logarithmic scale with base 10 for slider input
function log10 (sliderId) {
  $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return ('1e'+num); }
    })
}"

# JavaScript code needed to call above function
JS.onCall <-
"
$(document).ready(function() {
  // Wait for other scripts to execute
  setTimeout(function() {
    // Include call for each slider input
    log10('tolerance')
  }, 5)})
"

# Define UI for Shiny Application
shinyUI(navbarPage("Robust Statistics: Theory and Methods",
  
  # Tab to choose a data set
  tabPanel("Data",
    sidebarLayout(
      sidebarPanel(
        tags$head(tags$style(HTML(CSS.format1))),
        
        # Radio buttons to select data source
        radioButtons("source", "Data Source",
                     choices = c(Upload      = "upload",
                                 "R Package" = "library")),
        # Create panel for uploading data
        conditionalPanel(
          condition = "input.source == 'upload'",
          # Select a file
          fileInput("file", "Choose CSV File",
                    multiple = TRUE,
                    accept   = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
  
          tags$hr(),
    
          # Checkbox for header
          checkboxInput("header", "Header", TRUE),
    
          # Radio buttons for delimiting character
          radioButtons("sep", "Separator",
                       choices  = c(Comma     = ",",
                                    Semicolon = ";",
                                    Tab       = "\t"),
                       selected = ","),
    
          # Radio buttons for type of quotes used
          radioButtons("quote", "Quote",
                       choices  = c(None           = "",
                                    "Double Quote" = '"',
                                    "Single Quote" = "'"),
                       selected = '"')
        ),
        
        conditionalPanel(
          condition = "input.source == 'library'",
          # Selection of R packages
          selectInput("library", "Library Name",
                       choices = c(robustbase = "robustbase",
                                   RobStatTM  = "RobStatTM")),
          
          # Render UI given package selected
          uiOutput("select.data")
        ),
        
        # Button to display data table when pressed
        actionButton("display.table", "Load Data")
      ),
      
      # Panel for displaying the dataframe
      mainPanel(
        tags$head(tags$style(HTML(CSS.format1))),
        
        # Create Data table
        DT::dataTableOutput("contents.table")
      )
    )
  ),
  
  # Tab to choose parameters
  
  navbarMenu("Parameter Selection",
    tabPanel("Location",
      sidebarLayout(
        sidebarPanel(
          tags$head(tags$style(HTML(CSS.format1))),
          tags$head(tags$script(HTML(JS.log10))),
          tags$head(tags$script(HTML(JS.onCall))),
          
          # Renders selection of univariate vectors chosen from
          # dataset
          uiOutput("select.variable"),
          
          # Selection of score function family
          selectInput("psi", "Score Function (Psi)",
                      choices = c("bi-square" = "Bis",
                                  "Huber"     = "Hub")),
          
          # Radio buttons for desired asymptotic efficiency
          radioButtons("efficiency", "Asymptotic Efficiency",
                       choices  = c("0.85", "0.9", "0.95"),
                       selected = "0.9"),
          
          # Slider for maximum number of iterations
          sliderInput("max.iter", "Maximum Iterations",
                      min = 10, max = 100, value = 50),
          
          # Slider for tolerance level of convergence
          sliderInput("tolerance", "Error Tolerance",
                      min = -16, max = -3, value = -4),
           # Button to display data table when pressed
          actionButton("display.estimates", "View Output")
        ),
        
        mainPanel(
          tags$head(tags$style(HTML(CSS.format1))),
          # Display values for location and scale estimators
          verbatimTextOutput("estimates")
        )
      )
    ),
    
    tabPanel("Scale")
  ),
  
  navbarMenu("Linear Regression",
    tabPanel("Model"),
    tabPanel("Plots"),
    tabPanel("Predict")
  )
))