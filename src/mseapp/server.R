source('./www/globals.R')
source('helpers.R')


shinyServer(function(input,output,session){

    ## ------------------------------------------------------------ ##
    ## MANAGEMENT STRATEGY EVALUATION MODEL (May 15, 2015)
    ## ------------------------------------------------------------ ##

    #
    # plotMSE
    # 
    output$plotMSE <- renderPlot({
      .funnelPlot( input )
    })

    # 
    # Median depletion table
    # 
    output$viewDepletionTable <- renderTable({
      .tablePeformanceMetric(input,"t.Dt0.5")
    })

    # 
    # Probability of falling below SB 20%
    # 
    output$viewSSBLimitTable <- renderTable({
      .tablePeformanceMetric(input,"P.SSB.0.20.")
    })

    # 
    # Probability of falling below SB 30%
    # 
    output$viewSSBThresholdruTable <- renderTable({
      .tablePeformanceMetric(input,"P.SSB.0.30.")
    })

    # 
    # Median catch table
    # 
    output$viewCatchTable <- renderTable({
      .tablePeformanceMetric(input,"ct50")
    })

    # 
    # Median annual variation in catch
    # 
    output$viewAAVTable <- renderTable({
      .tablePeformanceMetric(input,"AAV50")
    })

	  ## ------------------------------------------------------------ ##
    ## EQUILIBRIUM MODEL (May 7, 2015)
    ## ------------------------------------------------------------ ##
    

    ## ------------------------------------------------------------ ##
    ## Run equilibrium models
    ## ------------------------------------------------------------ ##
    scnA <- reactive(do.call(equilibrium_model_cpp, getParams("A",input)))
    scnB <- reactive(do.call(equilibrium_model_cpp, getParams("B",input)))


    ## ------------------------------------------------------------ ##
    ## Plot Equilibrium values versus fishing mortality
    ## ------------------------------------------------------------ ##
    output$plot_equil <- renderPlot({
      AB <- rbind(scnA(),scnB())
      xx <- input$chkEquilPlot
      if(length(xx) != 0)
      {
        .plotEquilFe(AB,xx)
      }
    })		

    ## ------------------------------------------------------------ ##
    ## Run Selex plots
    ## ------------------------------------------------------------ ##
    output$plotSelex <-renderPlot({
      pars <- list(getParams("A",input),getParams("B",input))
      .plotBycatchSelex(pars)
    })

    output$plotFishSelex <-renderPlot({
      pars <- list(getParams("A",input),getParams("B",input))
      .plotFishSelex(pars)
    })

    ## ------------------------------------------------------------ ##
    ## Print Equilibrium MSY Table
    ## ------------------------------------------------------------ ##
    output$msyTable <- renderTable({
      AB <- rbind(scnA(),scnB())
      xx <- input$chkMSYTable

      .msyTable(AB,xx)

    })

    ## ------------------------------------------------------------ ##
    ## Print Equilibrium SPR Table
    ## ------------------------------------------------------------ ##
    output$sprTable <- renderTable({
      AB <- rbind(scnA(),scnB())
      xx <- input$chkMSYTable

      .sprTable(AB,xx)

    })

    ## ------------------------------------------------------------ ##
    ## Print Equilibrium MEY Table
    ## ------------------------------------------------------------ ##
    output$meyTable <- renderTable({
      AB <- rbind(scnA(),scnB())
      xx <- input$chkMSYTable

      .meyTable(AB,xx)

    })
})