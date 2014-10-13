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
    AB <<- rbind(scnA(),scnB())
    switch(input$a_chartType,
           "Equilibrium Yield"          = .plotEquilYield(AB),
           "Performance Metrics at MSY" = .plotPerformanceMSY(scnA())
           )
  })

  output$b_table <- renderTable({
    cat("Equilibrium MSY summary table")
    AB <- rbind(scnA(),scnB())
    
    test <- ddply(AB,.(prefix),plyr::summarize,
                  Fmsy=fe[which.max(AB$Ye)],
                  MSY =Ye[which.max(AB$Ye)]
                  )
    return(test)
    # return(DD)
  })
  # output$b_equilPlot <- renderPlot({
  #   # B <- rbind(scnA(),scnB())
  #   switch(input$b_chartType,
  #          "Equilibrium Yield"          = .plotEquilYield(scnB()),
  #          "Performance Metrics at MSY" = .plotPerformanceMSY(scnB())
  #          )
  # })


})



.plotEquilYield <- function(Scenario)
{
  

  mdf<-melt(Scenario,id.vars=1:5)
  sdf<-subset(mdf,variable %in% c("Ye","De","We","wbar_m"))
  levels(sdf$variable)[levels(sdf$variable)=="Ye"] <- "Directed Fishery Yield (Mlb)"
  levels(sdf$variable)[levels(sdf$variable)=="De"] <- "Directed Fishery Discard (Mlb)"
  levels(sdf$variable)[levels(sdf$variable)=="We"] <- "Directed Fishery Wastage (Mlb)"
  
  p <- ggplot(sdf,(aes(fe,value,col=prefix))) +geom_line()
  p <- p + facet_wrap(~variable,scales="free")
  p <- p + labs(x="Fishing Intensity",col="Scenario",y=NA)
  print(p + theme_bw(14))

}


.plotPerformanceMSY<- function(Scenario)
{
  # Row in which fisheries yield is maximized.
  ir <- which.max(Scenario$"Fishery Yield")

  out <- as.vector(Scenario[ir,])
  mdf <- melt(out,measure.vars=c("Ye","De","bycatch"))

  YIELD    <- out$"Ye"
  DISCARD  <- out$"De"
  BYCATCH  <- out$"bycatch"
  TOTAL    <- YIELD + BYCATCH + DISCARD



  p <- ggplot(melt(out),aes(variable,value)) + geom_bar(stat="identity")
  # print(ir)\
  # barplot(as.matrix(Scenario[ir,]))
  # barplot(as.matrix(Scenario[ir,]))
}







