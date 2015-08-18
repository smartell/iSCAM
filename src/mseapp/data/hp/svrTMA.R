

getArgsTMA <- function(input){

    print("in getargs")
    argsTMA <- list(Alloc_type=input$Allocation_type,sprtarget=input$ni_sprTarget,akdf=input$tbl)
    print(argsTMA)
    
    return(argsTMA)
}



getResultAllocation <- function(Alloc_type,sprtarget,akdf){

    print("in getResultAllocation")

    hpSTQ$target_spr <<- sprtarget
    hpSTQ$allocation  <<- as.numeric(akdf[,1])
        
    if(Alloc_type=="mortality per recruit"){

        hpSTQ$type<<-"MPR"
        tmp<-runModel(theta,hpSTQ)
        YPR_allocation  <- tmp$ypr/sum(tmp$ypr)
        out <- data.frame(YPR_allocation)   
   
    }else if(Alloc_type=="yield per recruit"){


        hpSTQ$type <<- "YPR"
        tmp <- runModel(theta,hpSTQ) 
        print("here?")      
        MPR_allocation  <- tmp$mpr/sum(tmp$mpr)
        out <- data.frame(MPR_allocation)

    }
    print(out)
    return(out)
}



