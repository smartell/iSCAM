

getArgsTMA <- function(input){

    print("in getargs")
    argsTMA <- list(Dist_type=input$Dist_type,sprtarget=input$ni_sprTarget,intbl=input$tbl)
    print(argsTMA)
    
    return(argsTMA)
}



getResultAllocation <- function(Dist_type,sprtarget,intbl,limPsc){

    print("in getResultAllocation")

    MP0$sprTarget <<- sprtarget

    #if(is.null(akdf)){
    #    akdf<-data.frame(proportion=list(rep(0,length(glbl))),row.names = c("IFQ","PSC","SPT","PER"))
    #}
        
    if(Dist_type=="mortality per recruit"){

        MP<-MP0

        MP$pMPR<-as.numeric(intbl$proportion)
        MP$type<-"MPR"
 
        tmpfs<-getFspr(MP)$par        
        MP$fstar  <- exp(tmpfs)

        rtmp     <- run(MP) 
        
        out <- data.frame(sector=c("IFQ","PSC","SPT","PER"),YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,f=rtmp$fe) 

   
    }else if(Dist_type=="yield per recruit"){

        MP<-MP0
       
        MP$pYPR<-as.numeric(intbl$proportion)
        MP$type<-"YPR"
       
        tmpfs<-getFspr(MP)$par
        MP$fstar  <- exp(tmpfs)
        
        rtmp     <- run(MP)
        
        out <- data.frame(sector=c("IFQ","PSC","SPT","PER"),YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,f=rtmp$fe)   
        
    }else if(Dist_type=="fixed PSC"){

        MP<-MP0

        MP$type<-"YPR"
        MP$pYPR<-as.numeric(intbl$proportion)
        print("batata")

        print(MP$pscLimit)
        print(intbl$cap[which(intbl$cap>0.0)])
        
        MP$pscLimit[which(intbl$cap>0.0)] <- as.numeric(intbl$cap[which(intbl$cap>0.0)])

        MP$pscLimit<-as.numeric(MP$pscLimit)
        print(as.numeric(MP$pscLimit))


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

        

        out <- data.frame(sector=c("IFQ","PSC","SPT","PER"),YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,f=rtmp$fe) 
    }
    print(out)
    return(out)
}

plotResultAllocation<-function(df){

    mdf<-melt(df)
    mdf$sector <- factor(mdf$sector , levels = c("IFQ","PSC","SPT","PER"))

    p<-ggplot(data=mdf, aes(x=sector, y=value, fill=variable)) 
    p<- p+geom_bar(stat="identity") +guides(fill=FALSE)
    p<- p + facet_wrap(~ variable, ncol=2, scales="free")
    p <- p +.THEME 
    print(p)


}

