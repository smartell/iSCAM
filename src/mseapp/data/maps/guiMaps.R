renderMaps <- function()
{
	tabPanel("Maps",
		buildMapGui()
	)
}

.AGES = c(10,15,20)
.SEXS = list(male="M",female="F")
buildMapGui <- function()
{
	fluidPage(
	    sidebarLayout(
	      sidebarPanel(
	        selectInput("si_sex","Choose Sex",
	                    c("male","female"),selected=1,
	                    multiple=FALSE),
	        selectInput("si_age","Choose Age",
	                    .AGES,multiple=FALSE,selected=15)
	    ),
	      mainPanel(
	        # tags$h1("That was easy"),
	        # includeHTML("www/F15animate.html")
	        htmlOutput("ui_Animation")
	      )
	    )
  )

}