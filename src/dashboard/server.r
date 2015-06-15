if(!require("shinydashboard"))      install.packages("shinydashboard")
source('./www/globals.R')
source('helpers.R')
server <- 
shinyServer(function(input, output, session){

    # 
    # DASHBOARD OUTPUT
    # 
    output$byr_ratio = renderText(1.32)









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
})