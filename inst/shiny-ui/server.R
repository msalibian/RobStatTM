# server.R
# 
# Author: Gregory Brownson
#
# Description: Backend portion of the Shiny 
#              interface for loading data and
#              selecting location and scale
#              parameters

library(DT)
library(RobStatTM)
library(shiny)



# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Need this to store reactive objects
  values <- reactiveValues()
  
  # Display correct options to obtain data given
  # method selected by user
  output$select.data <- renderUI({
    # If no package is chosen, then do nothing
    if(is.null(input$library)) {
      return()
    }
    
    # List of datasets
    list <- data(package=input$library)$results[,"Item"]
    
    # Get info on datasets
    values$data.info <- vcdExtra::datasets(c("robustbase"))
    
    # Create selection for datasets
    selectInput("dataset",
                label   = "Select Dataset",
                choices = list)
  })
  
  # On-click, load the data and return the data frame
  contents_table <- eventReactive(input$display.table, {
    if (input$source == "upload") {
      req(input$file)
      
      # Read data from file
      values$dat <- read.csv(input$file$datapath,
                      header = input$header,
                      sep    = input$sep,
                      quote  = input$quote)
      
      # If there are no headers, give data headers
      if (input$header == FALSE)
      {
        colnames(dat)[-1] <- paste0('X', 1:(ncol(dat) - 1))
      }
      
      
    } else {
      # If no dataset exists, return nothing
      if (is.null(input$dataset)) {
        return()
      }
      
      # Obtain specified dataset from specified package
      data(list = input$dataset, package = input$library)
      
      values$dat <- as.data.frame(get(input$dataset))
    }
    
    # Get variable names
    values$dat.variables <- colnames(values$dat)
      
    return (values$dat)
  })
  
  # Create data table
  output$contents.table <- DT::renderDataTable({
    contents_table()
  })
  
  # Render variable input list
  output$select.variable <- renderUI({
    # If there is no data, do nothing
    if (is.null(dim(values$dat))) {
      return()
    }
    
    # Render select input for variables
    selectInput("variable", "Variable",
                choices = values$dat.variables)
  })
  
  # On-click, find the estimators and create string object of results
  contents_estimators <- eventReactive(input$display.estimates, {
    # Get values for location and scale using 'MLocDis' function from
    # 'RobStatTM' package
    est <- MLocDis(x     = as.numeric(values$dat[,input$variable]),
                   psi   = input$psi,
                   eff   = input$efficiency,
                   maxit = input$max.iter,
                   tol   = input$tolerance)
    
    # Store results in string objects
    line1 <- paste0("Location estimate: ", round(est$mu, 4))
    line2 <- paste0("Scale estimate:    ", round(est$disper, 4))
    line3 <- paste0("Standard deviation of location: ", round(est$std.mu, 4))

    return(paste(line1, line2, "", line3, sep = "\n"))
  })
  
  # Display output
  output$estimates <- renderText({
    contents_estimators()
  })
})