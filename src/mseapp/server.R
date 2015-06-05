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
    output$viewSSBThresholdTable <- renderTable({
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



    # BYCATCH
    observe({
      A_dmr = input$A_bycatch_dmr
      B_dmr = input$A_bycatch_dmr
      print(A_dmr)
      A_bycatch = input$A_num_bycatch_total
      B_bycatch = input$B_num_bycatch_total

      A_bcm = A_dmr * A_bycatch
      B_bcm = B_dmr * B_bycatch
      updateNumericInput(session, paste0("A","_","num_bycatch"), value = A_bcm)
      updateNumericInput(session, paste0("B","_","num_bycatch"), value = B_bcm)
    })

    # output$num_bycatch <- renderText({
    #     bycatch = input$bycatch_dmr * input$num_bycatch_total
    # })
  #   observe({
  #   pars <- list(getParams("A",input),getParams("B",input))
  #   print(pars)
  #   dmr <- input$bycatch_dmr
  #   bct <- input$num_bycatch_total
  #   bcm <- 8
    
  #   updateNumericInput(session, paste0("A","_","num_bycatch"), value = bcm)
  #   updateNumericInput(session, paste0("B","_","num_bycatch"), value = bcm)

  #   # updateNumericInput(session, "inNumber2",
  #   #   label = paste("Number label ", x),
  #   #   value = x, min = x-10, max = x+10, step = 5)
  # })









})
# End of shinyServer