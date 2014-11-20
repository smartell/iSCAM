# NOTATION
# For R-scripts that deal with user interface, prefix file with "gui_*.R"
library(leaflet)
source("helpers.R")





renderOMI <- function()
{
	fluidRow(
	  column(4,
			wellPanel("Operating Model Interface",
				selectInput('omiplotType',"Choose Graphic",
			  	          c("Spawning biomass",
			    	          "Depletion",
			    	          "Recruitment",
			    	          "Stock Recruitment",
			    	          "Relative abundance",
			    	          "Mortality"),
			      	      selected="Spawning biomass")
			)
		),
		column(8,
		  plotOutput("omiPlot")
		)
	)

}



# ----------------------------------------#
# MAIN USER INTERFACE FOR THE APPLICATION #
# ----------------------------------------#
shinyUI(fluidPage(
        navbarPage("IPHC MSE TOOL",
                   id="nav",
                   footer=img(src="iphclogo.png",  height = 60, width = 60),
	
	
    # ---------------------------------------- #
	# EQUILIBRIUM INTERFACE
    # ---------------------------------------- #
	tabPanel("Equilibrium",
	    buildEquilibriumGui()
	),


	# ---------------------------------------- #
	# MANAGEMENT STRATEGY EVALUATION INTERFACE
	# ---------------------------------------- #
	tabPanel("MSE",
		buildMSEGui()
	),


	# ---------------------------------------- #
	# OPERATING MODEL INTERFACE
	# ---------------------------------------- #
	tabPanel("OMI",
		renderOMI()
	),

	

	# ---------------------------------------- #
	# INFORMATION INTERFACE (NEEDS TOC)
	# ---------------------------------------- #
	tabPanel("About",
	  fluidRow(
	    	navlistPanel(widths=c(3,9),
	    		"Navigation",
	    		tabPanel("About",
					includeMarkdown("www/About.md")
	    		),
	    		tabPanel("Equilibrium",
	    			includeMarkdown("www/Equilibrium.md")
	    		),
	    		tabPanel("MSE"),
	    		tabPanel("OMI",
	    			includeMarkdown("www/OMI.md")
	    		),
	    		"----",
	    		tabPanel("MAP")
	    	)
		)

	),

	# ---------------------------------------- #
	# MAPS
	# ---------------------------------------- #
	tabPanel("MAP",
	 	leafletMap("map",width=1600,height=800,
	 	    # initialTileLayer = "//{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
	 	   	initialTileLayer = "http://oatile3.mqcdn.com/tiles/1.0.0/sat/{z}/{x}/{y}.jpg",
		    options = list(
			center = c(52, -155),
			zoom = 5
		))
	)


)))
# ---------------------------------------- #
# END OF SHINY UI
# ---------------------------------------- # 


