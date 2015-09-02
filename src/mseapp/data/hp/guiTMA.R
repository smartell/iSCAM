
## RENDER TMA USER INTERFACE
renderTMA <- function(){
	tabPanel("Footprint",
	   buildTMAGui()
	)
}



buildTMAGui  <- function(){

   fluidPage( 
        sidebarLayout(
            sidebarPanel("Fisheries footprint",

                fluidRow(
                    
                    wellPanel( "SPR Target",
                        numericInput("ni_sprTarget", "Enter SPR target",0.3,0,1,0.1),

                        radioButtons("Dist_type", "Enter Distribution Method:",
                        c("yield per recruit", "mortality per recruit","fixed PSC"),selected = "yield per recruit")),
                        
                    wellPanel(                  
                        #allocation table
                        "Cells highlighted in red are ignored by the program",
                        htable("tbl", colHeaders="provided",rowNames = "provided")


                        
                    )
                )     
            
            ),
            
            mainPanel(


                fluidRow(
                    
                        #actionButton("actionButtonID","apply table edits"),
                                
                    
                        wellPanel( 
                            tableOutput("res_alloc"),
                            plotOutput("res_plot")
                        ) 

                )
            )           
        )
   ) 
                
}