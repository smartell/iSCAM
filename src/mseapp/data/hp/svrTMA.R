

getArgsTMA <- function(input){

    print("in getargs")
    argsTMA <- list(Dist_type=input$Dist_type,sprtarget=input$ni_sprTarget,akdf=input$tbl)
    print(argsTMA)
    
    return(argsTMA)
}



getResultAllocation <- function(Dist_type,sprtarget,akdf){

    print("in getResultAllocation")

    MP0$sprTarget <<- sprtarget

    #hpSTQ$allocation  <<- as.numeric(akdf[,1])
        
    if(Dist_type=="mortality per recruit"){

        MP<-MP0

        MP$pMPR<-as.numeric(akdf[,1])
        MP$type<-"MPR"
 
        MP$fstar    <- exp(getFspr(MP)$par)

        rtmp     <- run(MP) 
        
        out <- data.frame(YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,f=rtmp$fe) 

   
    } else if(Dist_type=="yield per recruit"){

        MP<-MP0

        MP$pYPR<-as.numeric(akdf[,1])
        MP$type<-"YPR"

        MP0$fstar  <- exp(getFspr(MP)$par)

        rtmp     <- run(MP)
        
        out <- data.frame(YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,f=rtmp$fe)   
   
    } else if(Dist_type=="fixed PSC"){

        EM<-run(MP0)
        
        out <- data.frame(YPR=EM$ypr, MPR=EM$mpr, yield=EM$ye,f=EM$fe) 
    }
    print(out)
    return(out)
}



