

getArgsTMA <- function(input){

    print("in getargs")
    argsTMA <- list(Dist_type=input$Dist_type,sprtarget=input$ni_sprTarget,akdf=input$tbl, limPsc=input$pscLim)
    print(argsTMA)
    
    return(argsTMA)
}



getResultAllocation <- function(Dist_type,sprtarget,akdf,limPsc){

    print("in getResultAllocation")

    MP0$sprTarget <<- sprtarget

    #hpSTQ$allocation  <<- as.numeric(akdf[,1])
        
    if(Dist_type=="mortality per recruit"){

        MP<-MP0

        MP$pMPR<-as.numeric(akdf[,1])
        MP$type<-"MPR"

        #MP$pscLimit  <-c(NA,as.numeric(limPsc),NA,NA)
 
        MP$fstar    <- exp(getFspr(MP)$par)

        rtmp     <- run(MP) 
        
        out <- data.frame(YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,f=rtmp$fe) 

   
    } else if(Dist_type=="yield per recruit"){

        MP<-MP0

        MP$pYPR<-as.numeric(akdf[,1])
        MP$type<-"YPR"

        #MP$pscLimit  <-c(NA,as.numeric(limPsc),NA,NA)

        MP$fstar  <- exp(getFspr(MP)$par)

        rtmp     <- run(MP)
        
        out <- data.frame(YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,f=rtmp$fe)   
   
    } else if(Dist_type=="fixed PSC"){

        MP<-MP0

        MP$type<-"YPR"

        MP$pscLimit  <-c(NA,as.numeric(limPsc),NA,NA) 


        ak    <- MP$pYPR
        bGear <- !is.na(MP$pscLimit)
        iGear <- which(!is.na(MP$pscLimit))
        pk    <- ak[!bGear]/sum(ak[!bGear])

        fs<-getFsprPSC(MP)
        
        tmp        <- ak
        tmp[bGear] <- fs$par[-1]
        tmp[!bGear]<- (1-sum(tmp[bGear]))*pk

        
        
        

        MP$fstar  <- exp(fs$par[1])

        MP$pYPR<- tmp

        rtmp  <- run(MP)

        
        out <- data.frame(YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,f=rtmp$fe) 
    }
    print(out)
    return(out)
}



