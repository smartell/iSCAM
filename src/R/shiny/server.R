library(shiny)
source("helpers.R")

paramNames <- c("size_limit",
                "discard_mortality_rate",
                "selex_fishery",
                "selex_bycatch",
                "num_bycatch")
# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

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

  # TULIP PLOTS FOR MSE PROCEDURES AND SCENARIOS
  output$funnelPlot <- renderPlot({
    tulip.plot(data(),input)
  })

  # GVIS GRAPHICS FOR COMPARING SCENARIOS DYNAMICALLY
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

  # OPERATING MODEL INTERFACE FUNCTIONS
  output$omiPlot <- renderPlot({
    cat(input$omiplotType)
    switch(input$omiplotType,
           "Spawning biomass"   = .plotSpawnBiomass(M),
           "Depletion"          = .plotDepletion(M),
           "Recruitment"        = .plotRecruitment(M),
           "Stock Recruitment"  = .plotStockRecruit(M),
           "Relative abundance" = .plotSurveyFit(M),
           "Mortality"          = .plotMortality(M)
           )
  })


  # EQUILIBRIUM MODEL INTERFACE 
  getParams <- function(prefix) {
    print("Hello")
    input[[paste0(prefix, "_recalc")]]

    params <- lapply(paramNames, function(p) {
      input[[paste0(prefix, "_", p)]]
    })
    names(params) <- paramNames
    params <- c(params,prefix=prefix)
    print(params)
    params
  }

  # output$a_selex <- renderPlot({
  #   x = 0:80
  #   par(mar=c(4,3,0,0))
  #   plot(x,.plogis95(x,input$a_selex_fishery[1],input$a_selex_fishery[2]),
  #        type="l",las=1,ylab=NA,xlab="Length (in.)",bty="l")
  # })
  # output$b_selex <- renderPlot({
  #   x = 0:80
  #   par(mar=c(4,3,0,0))
  #   plot(x,.plogis95(x,input$b_selex_fishery[1],input$b_selex_fishery[2]),
  #        type="l",las=1,ylab=NA,xlab="Length (in.)",bty="l")
  # })

  scnA <- reactive(do.call(equilibrium_model, getParams("a")))
  scnB <- reactive(do.call(equilibrium_model, getParams("b")))
  

  output$a_equilPlot <- renderPlot({
    # A <- scnA()
    switch(input$a_chartType,
           "Equilibrium Yield" = .plotEquilYield(scnA()),
           "Performance Metrics" = .plotPerformanceMetrics(scnA())
           )
    
  })
  output$b_equilPlot <- renderPlot({
    B <- rbind(scnA(),scnB())
    switch(input$b_chartType,
           "Equilibrium Yield"   = .plotEquilYield(B),
           "Performance Metrics" = .plotPerformanceMetrics(scnB())
           )
  })


})



.plotEquilYield <- function(Scenario)
{
  # qplot(fe,Ye,data=Scenario,geom_line)
  p <- ggplot(Scenario,aes(x=fe,y=Ye,color=prefix)) + geom_line()
  p <- p + labs(x="Fishing intensity",y="Directed yield (Mlb)",col="Scenario")
  print(p + mytheme())

}


.plotPerformanceMetrics <- function(Scenario)
{
  ir <- which.max(Scenario$Ye)
  print(ir)
  barplot(as.matrix(Scenario[ir,]))
  barplot(as.matrix(Scenario[ir,]))
}







