source("helpers.R")




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
    print(params)
    params
  }

  # getProcedure <- function(prefix) {

  # }


  ## Run equilibrium models
  scnA <- reactive(do.call(equilibrium_model_cpp, getParams("A")))
  scnB <- reactive(do.call(equilibrium_model_cpp, getParams("B")))

  output$plot_equil <- renderPlot({
    AB <- rbind(scnA(),scnB())
    xx <- input$selEquilPlot
    print(xx)
    .plotObject(AB,xx)
  })

  # output$a_equilPlot <- renderPlot({
  #   AB <- rbind(scnA(),scnB())
  #   switch(input$selChartType,
  #          "Equilibrium Yield"          = .plotEquilYield(AB),
  #          "Performance Metrics at MSY" = .plotPerformanceMSY(AB),
  #          "Equilibrium Value"          = .plotEquilValue(AB),
  #          "Value at MSY"               = .plotPerformanceValue(AB)
  #          )
  # })

  output$table_biological <- renderTable({
    AB <- rbind(scnA(),scnB())
    .biologicalTable(AB)
  })

  output$table_fishery <- renderTable({
    AB <- rbind(scnA(),scnB())
    .fisheryTable(AB)
  })

  output$table_economics <- renderTable({
    AB <- rbind(scnA(),scnB())
    .economicTable(AB)
  })

  output$msytable <- renderTable({
    AB <- rbind(scnA(),scnB())
    .equilibriumTables(AB)
  })
  
  # output$sprtable <- renderTable({
  #   AB <- rbind(scnA(),scnB())
  #   # .sprTables(AB)
  # })

  # output$u26ratio <- renderTable({
  #   AB <- rbind(scnA(),scnB())
  #   # .u26Table(AB)
  # })

})  # End of ShinyServer
## ------------------------------------------------------------------------------------ ##

.plotObject <- function(Scenario,objs)
{
  mdf<-melt(Scenario,id.vars=1:7)
  # print(tail(Scenario))
  SS <- Scenario  %>% group_by(prefix) %>% filter(Yield == max(Yield))
  # print(SS)
  mSS <- subset(melt(SS,id.vars=1:7),variable %in% objs)
  sdf<-subset(mdf,variable %in% objs)

  p <- ggplot(sdf,(aes(fe,value,col=prefix))) + geom_line()
  
  p <- p + geom_segment(aes(x=fe,xend=fe,y=value,yend=0,col=prefix),data=mSS, arrow = arrow(length = unit(0.5,"cm")))
  p <- p + geom_segment(aes(x=fe,xend=0,y=value,yend=value,col=prefix),data=mSS, arrow = arrow(length = unit(0.5,"cm")))
  p <- p + facet_wrap(~variable,scales="free")
  p <- p + labs(x="Fishing Intensity",col="Procedure",y="")
  p <- p + theme_bw(18) + theme(legend.position="top") 
  print(p)
  

}

.biologicalTable <- function(Scenario)
{
    # x <- runif(1e5)
    # cat("Mean of x = ",meanC(x))

   test <- ddply(Scenario,.(prefix),plyr::summarize,
                  BMSY    = Be[which.max(Ye)],
                  DMSY    = depletion[which.max(Ye)],
                  BTarget = Be[which.min(depletion>ssb_threshold)],
                  BLimit  = Be[which.min(depletion>ssb_limit)],
                  SPR     = SPR[which.min(depletion>ssb_threshold)],
                  Fmsy    = fe[which.max(Ye)],
                  Ftarget = fe[which.min(depletion>ssb_threshold)],
                  Flimit  = fe[which.min(depletion>ssb_limit)]
                  )
    colnames(test)=c("Procedure",
                     "♀ SSB_MSY",
                     "♀ Depletion",
                     "♀ SSB Threshold",
                     "♀ SSB Limit",
                     "SPR Threshold",
                     "F MSY",
                     "F Target",
                     "F Limit")
    rownames(test)=NULL
    print(test)
}

.fisheryTable <- function(Scenario)
{
  
  test <- ddply(Scenario,.(prefix),plyr::summarize,
                  "TCEY"    = TCEY[which.min(depletion>ssb_threshold)],#+O26[which.min(depletion>ssb_threshold)]+We[which.min(depletion>ssb_threshold)]+Ye[which.min(depletion>ssb_threshold)],
                  "O26 Bycatch" = O26[which.min(depletion>ssb_threshold)],
                  "Wastage" = We[which.min(depletion>ssb_threshold)],
                  "FCEY"    = Ye[which.min(depletion>ssb_threshold)],
                  "FCEY/TCEY (%)" = OE[which.min(depletion>ssb_threshold)]*100,
                  "U26 Bycatch" = U26[which.min(depletion>ssb_threshold)]

                  )
    colnames(test)[1]="Procedure"
    print(test)
}

.economicTable <- function(Scenario)
{
  test <- ddply(Scenario,.(prefix),plyr::summarize,
                  "Landed value"= YEv[which.min(depletion>ssb_threshold)],
                  "Wastage value" = WEv[which.min(depletion>ssb_threshold)],
                  "Bycatch value" = BYv[which.min(depletion>ssb_threshold)],
                  "Total value of all mortality" = YEv[which.min(depletion>ssb_threshold)]+WEv[which.min(depletion>ssb_threshold)]+BYv[which.min(depletion>ssb_threshold)]
                  )
    colnames(test)[1]="Procedure"
    print(test)
}

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
  colnames(test)[1]="Procedure"
  print(test)
}

.equilibriumTables <- function(Scenario)
{
  test <- ddply(Scenario,.(prefix),plyr::summarize,
                  Fmsy=fe[which.max(Ye)],
                  MSY =Ye[which.max(Ye)],
                  BMSY=Be[which.max(Ye)],
                  Depl=depletion[which.max(Ye)],
                  DMSY=De[which.max(Ye)],
                  WMSY=We[which.max(Ye)]
                  )
    colnames(test)=c("Procedure","F","MSY","Biomass (♀)",
                     "Depletion","Discards","Wastage")
    rownames(test)=NULL
    print(test)
}


.plotEquilYield <- function(Scenario)
{
  helpText("These are the equilibrium yield versus fishing effort.")

  mdf<-melt(Scenario,id.vars=1:7)
  sdf<-subset(mdf,variable %in% c("Ye","De","We","SPR"))
  levels(sdf$variable)[levels(sdf$variable)=="Ye"] <- "Directed Fishery Yield (Mlb)"
  levels(sdf$variable)[levels(sdf$variable)=="De"] <- "Directed Fishery Discard (Mlb)"
  levels(sdf$variable)[levels(sdf$variable)=="We"] <- "Directed Fishery Wastage (Mlb)"
  levels(sdf$variable)[levels(sdf$variable)=="SPR"] <- "SPR"
  
  print(head(sdf))
  
  p <- ggplot(sdf,(aes(fe,value,col=prefix))) + geom_line()
  p <- p + facet_wrap(~variable,scales="free")
  p <- p + labs(x="Fishing Intensity",col="Procedure",y="")
  p <- p + theme_bw(14) + theme(legend.position="top") 
  print(p)
  

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
  p <- p + labs(x="Fishing Intensity",col="Procedure",y="Millions of dollars")
  p <- p + theme_bw(14) + theme(legend.position="top") 
  print(p)

}


.plotPerformanceMSY<- function(Scenario)
{
  x<-ddply(Scenario,.(prefix),plyr::summarize,
           "Fishing Intensity @ MSY"=fe[which.max(Ye)],
           "Maximum Fishery Yield"  =Ye[which.max(Ye)],
           "Discards"               =De[which.max(Ye)],
           "Wastage"                =We[which.max(Ye)]
           )

  p <- ggplot(melt(x,id.vars=1),aes(variable,value,fill=prefix))
  p <- p + geom_bar(stat="identity",position="dodge")
  p <- p + labs(x="Variable",y="Value (million lbs or fishing intensity)",fill="Procedure")
  p <- p +facet_wrap(~variable,scales="free")
  print(p + theme_bw(14))
}

.plotPerformanceValue<- function(Scenario)
{
  x<-ddply(Scenario,.(prefix),plyr::summarize,
         "Landed Value @ MSY"           =YEv[which.max(Ye)],
         "Value of Wastage"             =WEv[which.max(Ye)],
         "Total Value of all mortality" =YEv[which.max(Ye)]+WEv[which.max(Ye)]+BYv[which.max(Ye)],
         "Value of losses"              =(YEv[which.max(Ye)]+WEv[which.max(Ye)]+BYv[which.max(Ye)])-YEv[which.max(Ye)]
         )

  p <- ggplot(melt(x,id.vars=1),aes(variable,value,fill=prefix))
  p <- p + geom_bar(stat="identity",position="dodge")
  p <- p + labs(x="Variable",y="Value (million $)",fill="Procedure")
  p <- p +facet_wrap(~variable,scales="free")
  print(p + theme_bw(14))
}





