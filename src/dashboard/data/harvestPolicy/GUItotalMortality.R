#GUItotalMortality.R
renderTMA <- function()
{

	print("renderTMA")

    box(title = "Harvest Specification", width = 3, status = "primary", solidHeader = TRUE, 
        # nfleets
        numericInput("nfleets","No. fleets",1,1,10,1),

        # sprTarget
		numericInput("sprTarget","Target SPR",0.40,0,1,0.01),

		# allocation
		tableOutput("allocationTable")


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