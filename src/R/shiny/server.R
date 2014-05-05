library(shiny)
source("helpers.R")


# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  # Subset Dataframe based on User Interface Selection.
  data <- reactive({
      
    if(input$plotType=="Spawning biomass" || input$plotType=="Depletion")
    {
      DF <- BIO.DF
    }
    if(input$plotType=="Catch" || input$plotType=="Fishing mortality")
    {
      DF <- CAT.DF
    }
    if(   input$plotType=="Sub-legal Catch" || input$plotType=="Wastage")
    {
      DF <- SUB.DF
    }
    if(input$plotType=="AAV in Catch" )
    {
      DF <- AAV.DF
    }
    if(input$plotType=="Efficiency")
    {
      DF <- MSE.DF
    }

    a  <- subset(DF,
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
    tulip.plot(data(),input)
  })

  output$googleVisPlot <- renderGvis({
    motionChart(MOT.DF,input)
  })


  # MEDIAN DEPLETION TABLE
  output$viewDepletionTable <- renderTable({

    cat("Depletion table \n")
    df  <- subset(mse.data$biomass.df,
             Year      %in% input$years[1]:input$years[2] &
             Scenario  %in% input$scenario                &
             Procedure %in% input$procedure
             )
    mdf <- melt(df,id=c("Scenario","Procedure","Year"))
    tmp <- dcast(mdf,Procedure~Scenario,mean,na.rm=TRUE,margins="Scenario",
                 subset=.(variable=="t.Dt0.5"))
    return(tmp)
 
  })

  # MEDIAN CATCH TABLE
  output$viewCatchTable <- renderTable({

    cat("Catch table \n")
    df  <- subset(mse.data$catch.df,
             Year      %in% input$years[1]:input$years[2] &
             Scenario  %in% input$scenario                &
             Procedure %in% input$procedure
             )
    mdf <- melt(df,id=c("Scenario","Procedure","Year"))
    tmp <- dcast(mdf,Procedure~Scenario,mean,na.rm=TRUE,margins="Scenario",
                 subset=.(variable=="ct50"))
    return(tmp)
 
  })

  # PROBABILITY OF GOING BELOW THE LIMIT REFERENCE POINT: P(SB<0.20)
  output$viewSSBlimit <- renderTable({

    cat("SSB Limit table \n")
    df  <- subset(mse.data$biomass.df,
             Year      %in% input$years[1]:input$years[2] &
             Scenario  %in% input$scenario                &
             Procedure %in% input$procedure
             )
    mdf <- melt(df,id=c("Scenario","Procedure","Year"))
    tmp <- dcast(mdf,Procedure~Scenario,mean,na.rm=TRUE,margins="Scenario",
                 subset=.(variable=="P.SSB.0.20."))
    return(tmp)

  })

  # PROBABILITY OF GOING BELOW THE THRESHOLD REFERENCE POINT: P(SB<0.30)
  output$viewSSBthreshold <- renderTable({

    cat("SSB Threshold table \n")
    df  <- subset(mse.data$biomass.df,
             Year      %in% input$years[1]:input$years[2] &
             Scenario  %in% input$scenario                &
             Procedure %in% input$procedure
             )
    mdf <- melt(df,id=c("Scenario","Procedure","Year"))
    tmp <- dcast(mdf,Procedure~Scenario,mean,na.rm=TRUE,margins="Scenario",
                 subset=.(variable=="P.SSB.0.30."))
    return(tmp)

  })

})

