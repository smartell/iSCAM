library(shiny)
source("helpers.R")


# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  # Subset Dataframe based on User Interface Selection.
  data <- reactive({
    a  <- subset(mse.DF,
                 Year      %in% input$years[1]:input$years[2] &
                 Scenario  %in% input$scenario                &
                 Procedure %in% input$procedure
                 )
  })


  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot

  output$distPlot <- renderPlot({
    x    <- faithful[, 2]  # Old Faithful Geyser data
    bins <- seq(min(x), max(x), length.out = input$bins + 1)

    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
  })

  output$funnelPlot <- renderPlot({
    funnel.plot(data(),input)
  })
})