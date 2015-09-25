# source('./www/globals.R')
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
    ## Total mortality allocation - 2 (Sep 3, 2015)
    ## ------------------------------------------------------------ ##


    #
    # Allocation Input-output function
    #
    validate <- function(tbl,prefix){           

            if(input[[paste0(prefix,"_","Dist_type")]]=="fixed PSC"){

                vld<-c(1,3,4)
                sumAlloc<-as.numeric(tbl[vld,"proportion"])
                    
                updateTableStyle(session, paste0(prefix,"_","tbl"), "invalid", 2, 1)
                updateTableStyle(session, paste0(prefix,"_","tbl"), "valid", 2, 2)
                updateTableStyle(session, paste0(prefix,"_","tbl"), "invalid", c(1,3,4), 2)

            }else{

                vld<-1:4
                sumAlloc<-as.numeric(tbl[,"proportion"])
                updateTableStyle(session, paste0(prefix,"_","tbl"), "invalid", 1:(length(glbl)), 2)
            }   


            if(sum(sumAlloc)!=1.0){
                updateTableStyle(session, paste0(prefix,"_","tbl"), "warning", vld, 1) 
            }else{
                updateTableStyle(session, paste0(prefix,"_","tbl"), "valid", vld, 1)
            }              
                                      
    }



    output$A_tbl <- renderHtable({
        do.call(buildOuttbl,list(prefix="A"))
    })

    output$B_tbl <- renderHtable({
        do.call(buildOuttbl,list(prefix="B"))
    })

    buildOuttbl<-function(prefix){

        if(is.null(input[[paste0(prefix,"_","tbl")]])){

            tbl<-data.frame(list(proportion=c(0.80,0.00,0.17,0.03),cap=c(0.00,7.75,0.00,0.00)),row.names = c("IFQ","PSC","SPT","PER"))

            validate(tbl,prefix)
            
            return(tbl)
     
        }else{

        
            if(input[[paste0(prefix,"_","Dist_type")]]=="yield per recruit"){
                 
                 tbl<-data.frame(input[[paste0(prefix,"_","tbl")]],row.names = c("IFQ","PSC","SPT","PER"))
            
            }else if(input[[paste0(prefix,"_","Dist_type")]]=="mortality per recruit"){
                
                 tbl<-data.frame(input[[paste0(prefix,"_","tbl")]],row.names = c("IFQ","PSC","SPT","PER"))
            
            }else{ 
                
                tbl<-data.frame(input[[paste0(prefix,"_","tbl")]],row.names = c("IFQ","PSC","SPT","PER") )

            }  

            
            validate(tbl,prefix)
      
          return(tbl)
              
        }

     }   
       
  
    getAllocA<-reactive(do.call(getResultAllocation2,getArgsTMA2(input,"A")))
    getAllocB<-reactive(do.call(getResultAllocation2,getArgsTMA2(input,"B")))

    

    output$res_alloc <-renderTable({
        
        Atab<-getAllocA()
        Btab<-getAllocB()

        allocTable(Atab,Btab)
    })


    output$res_plot <-renderPlot({
        
        Atab<-getAllocA()
        Btab<-getAllocB()
        
        plotResultAllocation2(Atab,Btab)

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


  output$ui_Animation <- renderUI({
    sex <- .SEXS[[input$si_sex]]
    age <- input$si_age
    # get html file name
    fn <- paste0("./www/",sex,"_",age,".html")
    print(fn)
    return(includeHTML(fn))
  })






})
# End of shinyServer
