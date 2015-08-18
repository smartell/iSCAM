
## RENDER TMA USER INTERFACE
renderTMA <- function(){
	tabPanel("Total Mortality Allocation",
	   buildTMAGui()
	)
}



buildTMAGui  <- function(){

   fluidPage( 
        sidebarLayout(
            sidebarPanel("Total Mortality allocation",
                fluidRow(
                    selectInput("Allocation_type", "Enter Allocation Method:",c("yield per recruit", "mortality per recruit")),
                    
                    wellPanel( "Enter Allocation:",
                                        
                        htable("tbl", colHeaders="provided")
                    ),
                    
                    wellPanel( 

                        tableOutput("res_alloc")
                    ) 

                )
            ),
            mainPanel(

                fluidRow(
                    wellPanel( "Harvest specs",
                        numericInput("ni_sprTarget", "Enter SPR target",0.4,0,1,0.1)
                    )
                )            
            )
        )
   ) 
                
}