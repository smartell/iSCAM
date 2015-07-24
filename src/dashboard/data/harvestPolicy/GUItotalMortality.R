#GUItotalMortality.R
buildUItotalMortality <- function()
{

	print("renderTMA")
	fluidRow(
	    box(title = "Harvest Specification", width = 3, status = "primary", solidHeader = TRUE, 
	        # nfleets
	        numericInput("nfleets","No. fleets",1,1,10,1),

	        # sprTarget
			numericInput("sprTarget","Target SPR",0.40,0,1,0.01),

			# allocation
			tableOutput("allocationTable")

	    ),

	    box(title="Selectivity",width=5,status="primary",solidHeader=TRUE,
	        sliderInput("selexSldr","50% and 95% percentiles",min=0,max=1,value=c(0.3,0.9),step=0.01)
	    ),

	    box(title="Results",width=5,status="info",solidHeader=TRUE,icon=icon("cog"),
	        tableOutput("TMAtable")
	    )
	)

}



# shiny::runApp(list(
#   ui = pageWithSidebar(
#       headerPanel("test"),
#       sidebarPanel(
#         actionButton("test","add a row")),
#       mainPanel(
#         tableOutput("value"))
#   ),   
#   server = function(input,output){
#     observe({
#       if (input$test == 0) 
#         return()
#       isolate({
#         output$value <-renderTable({
#           num.inputs.col1 <- paste0("<input id='c1n", 1:input$test, "' class='shiny-bound-input' type='number' value='2'>")
#           num.inputs.col2 <- paste0("<input id='c2n", 1:input$test, "' class='shiny-bound-input' type='number' value='2'>")
#           data.frame(num.inputs.col1, num.inputs.col2)
#         }, sanitize.text.function = function(x) x)
#       })
#     })
#   }
# ))