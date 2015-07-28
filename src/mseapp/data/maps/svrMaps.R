

.animationMap <- function(input,output)
{
	output$ui_Animation <- renderUI({
	  	sex <- .SEXS[[input$si_sex]]
  		age <- input$si_age
  		# get html file name
  		fn <- paste0("./www/",sex,"_",age,".html")
  		print()
  		return(includeHTML(fn))
  })
}