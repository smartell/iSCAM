
## RENDER TMA USER INTERFACE
renderTMA <- function(){
	tabPanel("Footprint",
	   buildTMAGui()
	)
}



buildTMAGui  <- function(){

   fluidPage( 
       
        column(5,

            #buildTMAInterface()   
             wellPanel( 
            buildTMAinpProc()    
            )
            ),
            
            column(7,


                fluidRow(
                    
                        #actionButton("actionButtonID","apply table edits"),
                                
                    
                        wellPanel( 
                            tableOutput("res_alloc"),
                            plotOutput("res_plot")
                        ) 

                )
            )           
        
   ) 
                
}



buildTMAInputs<-function(prefix){
    fluidRow(
                    
                wellPanel( 

                    radioButtons(paste0(prefix,"_","Dist_type"), "Enter Distribution Method:",
                    c("yield per recruit", "mortality per recruit","fixed PSC"),selected = "fixed PSC")
                ),
                        
                wellPanel("Distribution table:",                
                    
                    htable(paste0(prefix,"_","tbl"), colHeaders="provided",rowNames = "provided"),
                    tags$li("Proportion column should add up to 1."),
                    tags$li("Red cells are ignored by the program.")
                       
                ),

                wellPanel("Directed fisheries controls",                  

                    sliderInput(paste0(prefix,"_","sl_sizeLim"), "Size Limits (in)", 2,100,c(32,100),2),

                    sliderInput(paste0(prefix,"_","sl_50sel"), "50% selectivity size (in)", 10,80,26.5,0.5)
                    #radioButtons(paste0(prefix,"_","Excluder"), "Excluder Option:",
                    #c("no excluder", "moderate excluder","intensive excluder"),selected = "no excluder")                        
                        
                ),

                wellPanel("Bycatch fisheries controls",                  

                    sliderInput(paste0(prefix,"_","sl_mortRate"), "Juvenile Mortality", -1.0,1.0,0.0,0.1),

                    radioButtons(paste0(prefix,"_","Excluder"), "Excluder Option:",
                    c("no excluder", "moderate excluder","intensive excluder"),selected = "no excluder")                        
                        
                )
            )     
}

buildTMAinpProc <- function(){
    
    fluidRow(
    
        column(6,
            tags$h4("Procedure A"),
            buildTMAInputs("A")
            ), 

        column(6,
            tags$h4("Procedure B"),
            buildTMAInputs("B")
            )       
    )

}
