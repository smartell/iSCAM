
## RENDER TMA USER INTERFACE
renderTMA <- function(){
	tabPanel("Total Mortality Allocation",
	   buildTMAGui()
	)
}



buildTMAGui  <- function(){

   fluidPage( 
        sidebarLayout(
            sidebarPanel("Fisheries footprint",

                fluidRow(
                    
                    wellPanel( "SPR Target",
                        numericInput("ni_sprTarget", "Enter SPR target",0.3,0,1,0.1)
                    )
                )     
            
            ),
            
            mainPanel(


                fluidRow(
                    radioButtons("Dist_type", "Enter Distribution Method:",
                        c("yield per recruit", "mortality per recruit","fixed PSC")),
                        
                    wellPanel(                  
                        #allocation table
                        "Fisheries footprint",
                        htable("tbl", colHeaders=c("provided"),rowNames = c("provided")),

                        "PSC cap:",
                        htable("pscLim", colHeaders="provided")

                        ),
                                
                    
                        wellPanel( 
                            tableOutput("res_alloc")
                        ) 

                )
            )           
        )
   ) 
                
}