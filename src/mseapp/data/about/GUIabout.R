# ---------------------------------------- #
# INFORMATION INTERFACE (NEEDS TOC)
# ---------------------------------------- #

renderAbout <- function()
{
	tabPanel("About",
	  fluidRow(
	    	navlistPanel(widths=c(2,10),
	    		"Navigation",
	    		tabPanel("About",
					includeMarkdown("www/About.md")
	    		),
	    		tabPanel("Equilibrium",
	    			includeMarkdown("www/Equilibrium.md")
	    		),
	    		tabPanel("MSE Dashboard",
	    		 	includeMarkdown("www/MSEDashboard.md")       
	    		)
	    		# tabPanel("OMI",
	    		# 	includeMarkdown("www/OMI.md")
	    		# ),
	    		# "----",
	    		# tabPanel("MAP")
	    	)
		)

	)
}