# source('./www/globals.R')

shinyServer(function(input,output,session){

    
    

    output$tbl <- renderHtable({
    	if (is.null(input$tbl)){
    		rows <- ng
    		nomes=NULL
			for (i in 1:ng){
    			nomes[i]=paste0("fisheries ",i) 
    		}  
   		
    		if(input$Allocation_type=="yield per recruit"){

    			tbl<-data.frame(list(ak=rep(1/ng,ng)))
      			rownames(tbl) <- nomes		
			}
    		if(input$Allocation_type=="mortality per recruit"){
		
    			tbl<-data.frame(list(ak=rep(.99/ng,ng)))
      			rownames(tbl) <- nomes	
			}      
      	return(tbl)
    	}else{
      
    		tbl <- input$tbl
      		return(tbl)
    	}
    	print("cheguei aqui")
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

})
# End of shinyServer


getFdataframe<- function(Alloc_type,sprtarget=0.4,akdf=data.frame(list(ak=c(0.4,0.6)))){

    print("calculating Fdataframe")

    target_spr<<- sprtarget
    ak<<-as.numeric(akdf$ak)
    
    fitB <- optim(fe,fnB,method="BFGS",hessian=TRUE)
    fitC <- optim(fe,fnC,method="BFGS",hessian=TRUE)

    print(cbind(fitC$par,fitB$par))

}

getArgs <-function(input){

    print("in getargs")

    args <- list(Alloc_type=input$Allocation_type,sprtarget=input$ni_sprTarget,akdf=input$tbl)
    return(args)
}



make_no <- function(Alloc_type,sprtarget=0.4,akdf=data.frame(list(ak=c(0.5,0.5)))){

	print("in make_no")
	print(ak)
    target_spr<<- sprtarget
    ak<<-as.numeric(akdf$ak)
    
    if(Alloc_type=="mortality per recruit"){

        fitC  	<- optim(fe,fnC,method="BFGS",hessian=TRUE)
        result 	<- equilibriumModel(fitC$par)
        tmp1 	<- grep("yield", names(result))
        tmp2 	<- result[tmp1]/sum(result[tmp1])
        out 	<- data.frame(tmp2)   

    }
    if(Alloc_type=="yield per recruit"){
        
        fitB 	<- optim(fe,fnB,method="BFGS",hessian=TRUE)
        result 	<- equilibriumModel(fitB$par)
        tmp1 	<- grep("pmort", names(result))
        tmp2 	<- result[tmp1 ]/sum(result[tmp1])
        out 	<- data.frame(tmp2)

    }

    print(out)
}
