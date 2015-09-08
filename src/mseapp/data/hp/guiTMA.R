
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
                    
                    wellPanel( "SPR Target",
                        numericInput(paste0(prefix,"_","ni_sprTarget"), "Enter SPR target",0.3,0,1,0.1),

                        radioButtons(paste0(prefix,"_","Dist_type"), "Enter Distribution Method:",
                        c("yield per recruit", "mortality per recruit","fixed PSC"),selected = "yield per recruit")),
                        
                    wellPanel(                  
                        #allocation table

                        htable(paste0(prefix,"_","tbl"), colHeaders="provided",rowNames = "provided"),
                        tags$span(style="color:red", "Proportion column should add up to 1."),
                        tags$br(),
                        tags$span(style="color:red", "Cells highlighted in red are ignored by the program.")
                        
                    ),

                    wellPanel(                  

                        #checkboxInput(paste0(prefix,"_","ck_mortRate"), "High Juvenile Mortality", FALSE)

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
