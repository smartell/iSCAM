# source('./www/globals.R')
source('helpers.R')



shinyServer(function(input,output,session){

    runTMA <- reactive(do.call(getFdataframe, getArgs(input)))  

    output$to_fishSpecs<-renderTable({
        runTMA()
    })

    output$ui <- renderUI({
        do.call(mainPanel,make_ni(input$Allocation_type))
    })   
    
        

})
# End of shinyServer


getFdataframe<- function(sprtarget=0.4){

    print("calculating Fdataframe")

    target_spr<<- sprtarget
    #tma_allocation<<-ak
    
    fitB <- optim(fe,fnB,method="BFGS",hessian=TRUE)
    fitC <- optim(fe,fnC,method="BFGS",hessian=TRUE)

    print(cbind(fitC$par,fitB$par))

}

getArgs <-function(input){

    print("in getargs")

    args <- list(sprtarget=input$ni_sprTarget)
    return(args)
}




make_ni <- function(input){

    ni_list=list()
    aks=NULL
    if(input=="yield per recruit"){
        for (i in 1:ng){
             id=paste0("ni_ypr_F",i)
             #aks[i]=id
             lb=paste0("YPR fisheries ",i)
             ni_list[[i]]=numericInput(id, lb,value = 0.3, min=0, max=1, step=0.01)
            if(i==ng){
                 ni_list[[i]]=numericInput(id, lb,value = (1-0.3*(ng-1)) , min=0, max=1, step=0.01)
            }
        }
    }
    if(input=="mortality per recruit"){
        for (i in 1:ng){
             id=paste0("ni_mpr_F",i)
             aks[i]=id
             lb=paste0("MPR fisheries ",i)
             ni_list[[i]]=numericInput(id, lb,value = 0.2, min=0, max=1)
             if(i==ng){
                 ni_list[[i]]=numericInput(id, lb,value = (1-0.2*(ng-1)) , min=0, max=1, step=0.01)
            }
        }
    }
    return(ni_list)
}


make_no <- function(Alloc_type){

    target_spr<<- sprtarget
    tma_allocation<<-ak
    
    if(Alloc_type=="mortality per recruit"){

        fitC <- optim(fe,fnC,method="BFGS",hessian=TRUE)
        result<- equilibriumModel(fitC$par)
        y_pos<-grep("yield", names(result))
        y_allocation<-result[y_pos]/sum(result[y_pos])
        out <- matrix(y_allocation,ncol=1)   

    }
    if(Alloc_type=="yield per recruit"){
        
        fitB <- optim(fe,fnB,method="BFGS",hessian=TRUE)
        result<- equilibriumModel(fitB$par)
        m_pos<-grep("pmort", names(result))
        m_allocation<-result[m_pos]/sum(result[m_pos])
        out <- matrix(m_allocation,ncol=1)

    }
    return(out)
}
