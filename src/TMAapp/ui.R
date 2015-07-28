# ----------------------------------------#
# MAIN USER INTERFACE FOR THE APPLICATION #
# ----------------------------------------#

source("helpers.R")

shinyUI(
    fluidPage(
        sidebarLayout(
            sidebarPanel("Total Mortality allocation",
                fluidRow( 
                    selectInput("Allocation_type", "Allocation method",c("yield per recruit", "mortality per recruit")),
                    
                    uiOutput("ui")

                   
                    )
                 ),
            mainPanel(

                fluidRow(
                    wellPanel( "Harvest specs",
                        numericInput("ni_sprTarget", "Enter SPR target",0.4,0,1,0.1)
                    ),

                    wellPanel("fisheries Specs",
                        tableOutput("to_fishSpecs")),

                    wellPanel("resulting allocation",
                        tableOutput("res_alloc"))    
                )            
            )
        )
    )
)                


