# ----------------------------------------#
# MAIN USER INTERFACE FOR THE APPLICATION #
# ----------------------------------------#

source("helpers.R")

library(shiny)
library(shinyTable)

shinyUI(
    fluidPage(
        sidebarLayout(
            sidebarPanel("Total Mortality allocation",
                fluidRow( 
                    selectInput("Allocation_type", "Allocation method",c("yield per recruit", "mortality per recruit")),
                    
                    wellPanel(
                                        
                        htable("tbl", colHeaders="provided")),
                    
                    wellPanel(
                        tableOutput("res_alloc")) 

                   
                    )
                 ),
            mainPanel(

                fluidRow(
                    wellPanel( "Harvest specs",
                        numericInput("ni_sprTarget", "Enter SPR target",0.4,0,1,0.1)
                    ),

                    wellPanel("fisheries Specs",
                        tableOutput("to_fishSpecs"))

                       
                )            
            )
        )
    )
)                


