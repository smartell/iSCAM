



getArgsTMA2 <- function(input, prefix){

    print("in getargs")

    argsTMA <- list(Dist_type=input[[paste0(prefix,"_","Dist_type")]],sprtarget=input[[paste0(prefix,"_","ni_sprTarget")]],intbl=input[[paste0(prefix,"_","tbl")]],ck_mortRate=input[[paste0(prefix,"_","ck_mortRate")]])
    print(argsTMA)

    
    return(argsTMA)
}



getResultAllocation2 <- function(Dist_type,sprtarget,intbl,ck_mortRate){

    #Dist_type,sprtarget,intbl,limPsc
    print("in getResultAllocation2")

    if(ck_mortRate==TRUE){
        cm<<-c(0.1,0.1)
    }else{
        cm<<-c(0.0,0.0)
    }


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
        
        out <- data.frame(sector=c("IFQ","PSC","SPT","PER"),YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,effort=rtmp$fe) 

   
    }else if(Dist_type=="yield per recruit"){

        MP<-MP0
       
        MP$pYPR<-as.numeric(intbl$proportion)
        MP$type<-"YPR"
       
        tmpfs<-getFspr(MP)$par
        MP$fstar  <- exp(tmpfs)
        
        rtmp     <- run(MP)
        
        out <- data.frame(sector=c("IFQ","PSC","SPT","PER"),YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,effort=rtmp$fe)   
        
    }else if(Dist_type=="fixed PSC"){

        MP<-MP0

        MP$type<-"YPR"
        MP$pYPR<-as.numeric(intbl$proportion)

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

        out <- data.frame(sector=c("IFQ","PSC","SPT","PER"),YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,effort=rtmp$fe) 
    
    }
    print(out)
    return(out)
}

allocTable<-function(A,B){

        nomes<-c("YPR", "MPR", "yield", "effort")
        mps<-c(" A"," B")

        nomes1<-NULL
        for(i in 1:length(nomes)){
            tmp<-c(paste0(nomes[i],mps))
            nomes1<-c(nomes1,tmp)
        }

        A<-round(A[,-1],2)
        B<-round(B[,-1],2)
        

        sector=c("IFQ","PSC","SPT","PER")

        AB<-data.frame(sector,A[,1],B[,1],A[,2],B[,2],A[,3],B[,3],A[,4],B[,4])
        names(AB)<-c("sector",nomes1)            
        
        return(AB)
}

plotResultAllocation2<-function(df,df2){

    procedure<-rep(c("A","B"),each=nrow(df))
    df0<-cbind(rbind(df,df2),procedure)

    mdf<-melt(df0)
    print(mdf)

    mdf$sector <- factor(mdf$sector , levels = c("IFQ","PSC","SPT","PER"))

    p<-ggplot(data=mdf, aes(x=sector, y=value, fill=procedure)) 
    p<- p+geom_bar(stat="identity",position="dodge") 
    p<- p + facet_wrap(~ variable, ncol=2, scales="free")
    p <- p +.THEME 
    print(p)


}


