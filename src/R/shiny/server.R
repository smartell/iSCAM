library(shiny)
source("helpers.R")

paramNames <- c("size_limit",
                "discard_mortality_rate",
                "spr_target",
                "selex_fishery",
                "selex_bycatch",
                "num_bycatch",
                "five",
                "ten",
                "twenty",
                "forty")

## ------------------------------------------------------------------------------------ ##
## Define server logic required to run scripts
## ------------------------------------------------------------------------------------ ##
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
    
    # input[[paste0(prefix, "_recalc")]]

    params <- lapply(paramNames, function(p) {
      input[[paste0(prefix, "_", p)]]
    })
    names(params) <- paramNames
    params <- c(params,prefix=prefix)
    # print(params)
    params
  }

  ## Run equilibrium models
  scnA <- reactive(do.call(equilibrium_model, getParams("a")))
  scnB <- reactive(do.call(equilibrium_model, getParams("b")))


  output$a_equilPlot <- renderPlot({
    AB <<- rbind(scnA(),scnB())
    switch(input$selChartType,
           "Equilibrium Yield"          = .plotEquilYield(AB),
           "Performance Metrics at MSY" = .plotPerformanceMSY(AB),
           "Equilibrium Value"          = .plotEquilValue(AB),
           "Value at MSY"               = .plotPerformanceValue(AB)
           )
  })

  output$msytable <- renderTable({
    AB <- rbind(scnA(),scnB())
    .equilibriumTables(AB)
  })
  
  output$sprtable <- renderTable({
    AB <- rbind(scnA(),scnB())
    .sprTables(AB)
  })

  output$u26ratio <- renderTable({
    AB <- rbind(scnA(),scnB())
    .u26Table(AB)
  })

})  # End of ShinyServer
## ------------------------------------------------------------------------------------ ##

.u26Table <- function(Scenario)
{
  test <- ddply(Scenario,.(prefix),plyr::summarize,
                  "Fishery @ MSY"  =f26[which.max(Ye)],
                  "Bycatch @ MSY"  =b26[which.max(Ye)],
                  "Fishery @ SPR"  =f26[which.min(SPR>spr_target)],
                  "Bycatch @ SPR"  =b26[which.min(SPR>spr_target)]
                  )
  colnames(test)[1]="Scenario"
  print(test) 
}


.sprTables <- function(Scenario)
{
  test <- ddply(Scenario,.(prefix),plyr::summarize,
                  FSPR  =fe[which.min(SPR>spr_target)],
                  Catch =Ye[which.min(SPR>spr_target)],
                  BMSY  =Be[which.min(SPR>spr_target)],
                  Depl  =depletion[which.min(SPR>spr_target)],
                  DMSY  =De[which.min(SPR>spr_target)],
                  WMSY  =We[which.min(SPR>spr_target)],
                  EFF   =OE[which.min(SPR>spr_target)]
                  )
  colnames(test)[1]="Scenario"
  print(test)
}

.equilibriumTables <- function(Scenario)
{
  test <- ddply(Scenario,.(prefix),plyr::summarize,
                  Fmsy=fe[which.max(Ye)],
                  MSY =Ye[which.max(Ye)],
                  BMSY=Be[which.max(Ye)],
                  Depl=Be[which.max(Ye)]/max(Be),
                  DMSY=De[which.max(Ye)],
                  WMSY=We[which.max(Ye)],
                  EFF =OE[which.max(Ye)]
                  )
    colnames(test)=c("Scenario","F","MSY","Biomass (♀)",
                     "Depletion","Discards","Wastage","Efficiency")
    rownames(test)=NULL
    print(test)
}


.plotEquilYield <- function(Scenario)
{
  
  helpText("These are the equilibrium yield versus fishing effort.")

  mdf<-melt(Scenario,id.vars=1:5)
  sdf<-subset(mdf,variable %in% c("Ye","De","We","OE"))
  levels(sdf$variable)[levels(sdf$variable)=="Ye"] <- "Directed Fishery Yield (Mlb)"
  levels(sdf$variable)[levels(sdf$variable)=="De"] <- "Directed Fishery Discard (Mlb)"
  levels(sdf$variable)[levels(sdf$variable)=="We"] <- "Directed Fishery Wastage (Mlb)"
  levels(sdf$variable)[levels(sdf$variable)=="OE"] <- "Operational Efficiency (%)"
  
  
  p <- ggplot(sdf,(aes(fe,value,col=prefix))) +geom_line()
  # p <- p + geom_vline(xintercept=0.3, subset = .(variable == "Directed Fishery Yield (Mlb)"))
  # p <- p + geom_vline(data=sdg,aes(xintercept=fe[which.min("SPR">spr_target)],col="Scenario"),size=2,alpha=0.5)
  p <- p + facet_wrap(~variable,scales="free")
  p <- p + labs(x="Fishing Intensity",col="Scenario",y="")
  print(p + theme_bw(14))

}

.plotEquilValue <- function(Scenario)
{
  

  mdf<-melt(Scenario,id.vars=1:5)
  sdf<-subset(mdf,variable %in% c("YEv","DEv","BYv","WEv"))
  levels(sdf$variable)[levels(sdf$variable)=="YEv"] <- "Landed value"
  levels(sdf$variable)[levels(sdf$variable)=="DEv"] <- "Value of discards"
  levels(sdf$variable)[levels(sdf$variable)=="BYv"] <- "Value of bycatch mortality"
  levels(sdf$variable)[levels(sdf$variable)=="WEv"] <- "Value of wastage"
  
  p <- ggplot(sdf,(aes(fe,value,col=prefix))) +geom_line()
  p <- p + facet_wrap(~variable,scales="free")
  p <- p + labs(x="Fishing Intensity",col="Scenario",y="Millions of dollars")
  print(p + theme_bw(14))

}


.plotPerformanceMSY<- function(Scenario)
{
  x<-ddply(AB,.(prefix),plyr::summarize,
           "Fishing Intensity @ MSY"=fe[which.max(Ye)],
           "Maximum Fishery Yield"  =Ye[which.max(Ye)],
           "Discards"               =De[which.max(Ye)],
           "Wastage"                =We[which.max(Ye)]
           )

  p <- ggplot(melt(x,id.vars=1),aes(variable,value,fill=prefix))
  p <- p + geom_bar(stat="identity",position="dodge")
  p <- p + labs(x="Variable",y="Value (million lbs or fishing intensity)",col="Scenario")
  p <- p +facet_wrap(~variable,scales="free")
  print(p + theme_bw(14))
}

.plotPerformanceValue<- function(Scenario)
{
  x<-ddply(AB,.(prefix),plyr::summarize,
         "Landed Value @ MSY"           =YEv[which.max(Ye)],
         "Value of Wastage"             =WEv[which.max(Ye)],
         "Total Value of all mortality" =YEv[which.max(Ye)]+WEv[which.max(Ye)]+BYv[which.max(Ye)],
         "Value of losses"              =(YEv[which.max(Ye)]+WEv[which.max(Ye)]+BYv[which.max(Ye)])-YEv[which.max(Ye)]
         )

  p <- ggplot(melt(x,id.vars=1),aes(variable,value,fill=prefix))
  p <- p + geom_bar(stat="identity",position="dodge")
  p <- p + labs(x="Variable",y="Value (million $)",col="Scenario")
  p <- p +facet_wrap(~variable,scales="free")
  print(p + theme_bw(14))
}





