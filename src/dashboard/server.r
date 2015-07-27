if(!require("shinydashboard"))      install.packages("shinydashboard")
# source('./www/globals.R')
source('helpers.R')
server <- 
shinyServer ( function(input, output, session){
    cat(input$selexSldr)
    # 
    # DASHBOARD OUTPUT
    # 
    #output$mitigationBox <- renderInfoBox({
    #  infoBox(
    #    "Mitigation", paste0(25 , input$something, "%"), icon = icon("list"),
    #    color = "purple"
    #  )
    #  print("mitigationBox")
    #})

    # observe({
    #   print("Fleets")
    #   if (input$nfleets == 0)
    #     return()
    #   isolate({
    #     output$allocationTable <-renderTable({
    #       # num.inputs.col1 <- paste0("<input id='c1n", 1:input$test, "' class='shiny-bound-input' type='number' value='2'>")
    #       # num.inputs.col2 <- paste0("<input id='c2n", 1:input$test, "' class='shiny-bound-input' type='number' value='2'>")
    #       # data.frame(num.inputs.col1, num.inputs.col2)
    #       num.inputs.col1 <- paste0("<input id='c1n", 1:input$nfleets, "' class='shiny-bound-input' type='number' value='2'>")
    #       data.frame(num.inputs.col1)
    #     }, sanitize.text.function = function(x) x)
    #   })
    # })
    
    ## ------------------------------------------------------------ ##
    ## TOTAL MORTALITY ALLOCATION (July 24, 2015)
    ## ------------------------------------------------------------ ##
    getTMAparams <- function(input)
    {
  
      params <- lapply(TMA_PARNAMES, function(p) {
        input[[p]]
      })
      names(params) <- TMA_PARNAMES
      # params <- c(params,prefix=prefix)
      cat(params)
      params   
    }

    #getF <- reactive(do.call(getFs,getTMAparams(input)))

    # spr  <- reactiveValues()
    output$TMAtable <- renderTable({
      data.frame(spr = input$sprTarget )
    })

    # output$TMAtable <- renderTable({
    #   print("TMA Table")
    #   sprTarget = input$sprTarget

    #   # df <- getF()
    #   dff <-data.frame(a=c(1,2,3,sprTarget))
    #   return(dff)
    # })

















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

