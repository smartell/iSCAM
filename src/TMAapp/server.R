# source('./www/globals.R')

shinyServer(function(input,output,session){

    validate <- function(tbl){

        if(sum(as.numeric(tbl$allocation))>1){
            #tbl$allocation<-as.numeric(tbl$allocation)/sum(as.numeric(tbl$allocation))
            updateTableStyle(session, "tbl", "invalid", 1:ng, 1)        
        }else if(sum(as.numeric(tbl$allocation))<1){
            updateTableStyle(session, "tbl", "warning", 1:ng, 1)
        }else{
            updateTableStyle(session, "tbl", "valid", 1:ng, 1)
        }

    }
    
    output$tbl <- renderHtable({
    	if (is.null(input$tbl)){
    		rows <- ng
    		nomes=NULL
			for (i in 1:ng){
    			nomes[i]=paste0("fisheries ",i) 
    		}     		
    		if(input$Allocation_type=="yield per recruit"){

    			tbl<-data.frame(list(allocation=rep(0,ng)))	
			}
    		if(input$Allocation_type=="mortality per recruit"){
		
    			tbl<-data.frame(list(allocation=rep(0,ng)))
			}

        rownames(tbl) <- nomes  
        validate(tbl)     
      	return(tbl)
    	}else{
      
    		tbl <- input$tbl
            validate(tbl)
      		return(tbl)
    	}
    	  
    	print(input$tbl)
	})  
    
    runTMA <- reactive(do.call(getFdataframe, getArgs(input)))  

    output$to_fishSpecs<-renderTable({
        runTMA()
    })

    getAlloc<-reactive(do.call(getResultAllocation,getArgs(input)))

    output$res_alloc <-renderTable({
        getAlloc()
    })

    output$plotTMA <- renderPlot({
      .plotfs( getArgs(input) )
    })

})
# End of shinyServer

