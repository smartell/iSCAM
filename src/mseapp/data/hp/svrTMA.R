



getArgsTMA2 <- function(input, prefix){

    print("in getargs")

    argsTMA <- list(Dist_type=input[[paste0(prefix,"_","Dist_type")]],intbl=input[[paste0(prefix,"_","tbl")]],
        sl_sizLim=input[[paste0(prefix,"_","sl_sizeLim")]],sl_50sel=input[[paste0(prefix,"_","sl_50sel")]],
        sl_mortRate=input[[paste0(prefix,"_","sl_mortRate")]],rb_excluder=input[[paste0(prefix,"_","Excluder")]])
    
    print(argsTMA)
    
    if(is.null(argsTMA$intbl)){
        
        argsTMA$intbl<-data.frame(list(proportion=c(0.80,0.0,0.17,0.03),cap=c(0.0,7.75,0.0,0.0)),row.names = c("IFQ","PSC","SPT","PER"))
        print(prefix)
        print(argsTMA$intbl)
    }
    
    return(argsTMA)
}



getResultAllocation2 <- function(Dist_type,intbl,sl_sizLim,sl_50sel,sl_mortRate,rb_excluder){

    
    print("in getResultAllocation2")

    MP0$slim<<-c(sl_sizLim[1]/0.393701,00,sl_sizLim[1]/0.393701,00)
    MP0$ulim<<-c(sl_sizLim[2]/0.393701,sl_sizLim[2]/0.393701,sl_sizLim[2]/0.393701,sl_sizLim[2]/0.393701)

    MP0$slx$slx1[1]<<- sl_50sel
    

    cm<<-rep(sl_mortRate,2)

    if(rb_excluder=="no excluder"){
        MP0$slx$slx3[2] <<- 0.072
        
    }else if(rb_excluder=="moderate excluder"){
        
        MP0$slx$slx3[2] <<- 0.1
    
    }else{
        
        MP0$slx$slx3[2] <<- 0.2
    }


    #if(is.null(akdf)){
    #    akdf<-data.frame(proportion=list(rep(0,length(glbl))),row.names = c("IFQ","PSC","SPT","PER"))
    #}
        
    if(Dist_type=="mortality per recruit"){

        MP<-MP0

        MP$pMPR<-as.numeric(intbl$proportion)
        MP$type<-"MPR"
 
        MP<-getFspr(MP) 
                   
        rtmp     <- run(MP) 
        
        out <- data.frame(sector=c("IFQ","PSC","SPT","PER"),YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,effort=rtmp$fe*100) 
        return(out)
   
    }else if(Dist_type=="yield per recruit"){

        MP<-MP0
       
        MP$pYPR<-as.numeric(intbl$proportion)
        MP$type<-"YPR"
       
        MP<-getFspr(MP)               
        rtmp     <- run(MP) 
        
        out <- data.frame(sector=c("IFQ","PSC","SPT","PER"),YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,effort=rtmp$fe*100)   
        return(out)
        
    }else if(Dist_type=="fixed PSC"){

        MP<-MP0

        MP$type<-"YPR"
        MP$pYPR<-as.numeric(intbl$proportion)
       
        MP$pscLimit[which(intbl$cap>0.0)] <- as.numeric(intbl$cap[which(intbl$cap>0.0)])

        MP$pscLimit<-as.numeric(MP$pscLimit)
        
        print("fantastic")
        print(MP)
        fs<-getFsprPSC(MP)
        print("allons'y")        
        rtmp  <- run(fs)

        out <- data.frame(sector=c("IFQ","PSC","SPT","PER"),YPR=rtmp$ypr, MPR=rtmp$mpr, yield=rtmp$ye,effort=rtmp$fe*100) 
        
        return(out)

    }
    print(out)
    #return(out)
}

allocTable<-function(A,B){

        nomes<-c("YPR", "MPR", "yield", "effort")
        mps<-c(" A"," B")

        nomes1<-NULL
        for(i in 1:length(nomes)){
            tmp<-c(paste0(nomes[i],mps))
            nomes1<-c(nomes1,tmp)
        }

        #A<-round(A[,-1],2)
        #B<-round(B[,-1],2)
        
        A<-A[,-1]
        B<-B[,-1]

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


