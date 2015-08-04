# source('./www/globals.R')

shinyServer(function(input,output,session){

    validate <- function(tbl){

        if(sum(as.numeric(tbl$allocation))>1){
            tbl$allocation<-as.numeric(tbl$allocation)/sum(as.numeric(tbl$allocation))
            updateTableStyle(session, "tbl", "invalid", 1:ng, 1)
        }else if(sum(as.numeric(tbl$allocation))<1){
            updateTableStyle(session, "tbl", "warning", 
                            1:ng, 1)
        }else{
            updateTableStyle(session, "tbl", "valid", 
                            1:ng, 1)
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

    getAlloc<-reactive(do.call(make_no,getArgs(input)))

    output$res_alloc <-renderTable({
        getAlloc()
    })

    output$plotTMA <- renderPlot({
      .plotfs( getArgs(input) )
    })

})
# End of shinyServer


getFdataframe<- function(Alloc_type,sprtarget=0.4,akdf=data.frame(list(allocation=rep(0,ng)))){

    print("calculating Fdataframe")

    target_spr<<- sprtarget
    ak<<-as.numeric(akdf$allocation)
    
    fitB <- optim(fe,fnB,method="BFGS",hessian=TRUE)
    fitC <- optim(fe,fnC,method="BFGS",hessian=TRUE)

    print(cbind(fitC$par,fitB$par))

}

getArgs <-function(input){

    print("in getargs")

    args <- list(Alloc_type=input$Allocation_type,sprtarget=input$ni_sprTarget,akdf=input$tbl)
    return(args)
}



make_no <- function(Alloc_type,sprtarget=0.4,akdf=data.frame(list(allocation=rep(0,ng)))){

	print("in make_no")
    target_spr<<- sprtarget
    ak<<-as.numeric(akdf$allocation)
    
    if(Alloc_type=="mortality per recruit"){

        fitC  	<- optim(fe,fnC,method="BFGS",hessian=TRUE)
        result 	<- equilibriumModel(fitC$par)
        tmp1 	<- grep("yield", names(result))
        YPR_allocation 	<- result[tmp1]/sum(result[tmp1])
        out 	<- data.frame(YPR_allocation)   

    }
    if(Alloc_type=="yield per recruit"){
        
        fitB 	<- optim(fe,fnB,method="BFGS",hessian=TRUE)
        result 	<- equilibriumModel(fitB$par)
        tmp1 	<- grep("pmort", names(result))
        MPR_allocation 	<- result[tmp1 ]/sum(result[tmp1])
        out 	<- data.frame(MPR_allocation)

    }

    print(out)
}
