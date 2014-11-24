# NOTATION
# For R-scripts that deal with user interface, prefix file with "gui_*.R"

source("helpers.R")


# RENDER FILTER FOR CHOOSING SCENARIOS AND PROCEDURES
renderFilterInputs <- function()
{
	wellPanel(
		h4("Filter"),
		sliderInput("years", "Years:",
              			min = min(BIO.DF$Year), 
              			max = max(BIO.DF$Year), 
              			value = range(BIO.DF$Year),
              			format= "####"),

		selectInput('scenario',"Secnario",
					levels(BIO.DF$Scenario),
					selected = levels(BIO.DF$Scenario)[1],
					multiple = TRUE),

		selectInput('procedure',"Procedure",
              		levels(BIO.DF$Procedure),
              		selected =  levels(BIO.DF$Procedure)[1],
              		multiple = TRUE)
	)
}

# RENDER LAYOUT CONTROLS FOR TUPLIP PLOTS
renderLayoutInputs <- function()
{
	wellPanel(
	  fluidRow(
		  column(3,
			selectInput("icolor","Facet color",
								c(None = ".",
								  "Scenario",
								  "Procedure",
								  "gear",
								  "area",
								  "sex",
								  "group"),
								selected=".")
			),

			column(3,
			selectInput("facet_row","Facet row",
		            c(None = ".",
		              "Scenario",
		              "Procedure",
		              "gear",
		              "area",
		              "sex",
		              "group"),
		            selected="Procedure")
			),

			column(3,
			selectInput("facet_col","Facet column",
								c(None = ".",
								  "Scenario",
								  "Procedure",
								  "gear",
								  "area",
								  "sex",
								  "group"),
								selected="Scenario")
			),

			column(3,
			selectInput('plotType',"Facet variable",
								c( "Spawning biomass",
								   "Depletion",
								   "Catch",
								   "Sub-legal Catch",
								   "AAV in Catch",
								   "Wastage",
								   "Efficiency",
								   "Fishing mortality"),
								selected="Spawning biomass")	       
			)
		)
	)
}

renderMSEtabs <- function()
{
	tabsetPanel(type="tabs",id="tab",

    # Tulip plots
	  tabPanel("Tulip plots", 
	    renderLayoutInputs(),
		  plotOutput("funnelPlot")
		),

		# Google Vis plots
		tabPanel("Motion chart",
		  htmlOutput("googleVisPlot")
		),

		# Summary tables
		tabPanel("Performance Metrics",
			h4("Median depletion"),
			tableOutput("viewDepletionTable"),
			h4("Median catch"),
			tableOutput("viewCatchTable"),
			h4("Probability of falling below limit reference point P(SB<0.20)"),
			tableOutput("viewSSBlimit"),
			h4("Probability of falling below threshold reference point P(SB<0.30)"),
			tableOutput("viewSSBthreshold")
		)

	)
}


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

renderBanner <- function()
{
	wellPanel(
		# Logo image
	  column(9),
	  column(3,
	  	img(src="iphclogo.png",  height = 80, width = 80),
			img(src="iscamLogo.png", height = 80, width = 80)
	  )
	)
}

# ----------------------------------------#
# MAIN USER INTERFACE FOR THE APPLICATION #
# ----------------------------------------#
shinyUI(fluidPage(navbarPage("IPHC MSE TOOL",
	



	
    # ---------------------------------------- #
	# EQUILIBRIUM INTERFACE
    # ---------------------------------------- #
	tabPanel("Equilibrium",
	    includeCSS("styles.css"),
		fluidRow(
			tags$h4("Equilibrium Model")
		),

		renderEquilriumInterface(),
		renderShutter(),	  
		fluidRow(
			renderBanner()
		)

	),

	# ---------------------------------------- #
	# MANAGEMENT STRATEGY EVALUATION INTERFACE
	# ---------------------------------------- #
	tabPanel("MSE",

	  fluidRow(
			column(3,
				renderFilterInputs()
			),

			column(9,
		    renderMSEtabs()
			)
		),

		fluidRow(
			renderBanner()
		)
	),

	# ---------------------------------------- #
	# OPERATING MODEL INTERFACE
	# ---------------------------------------- #
	tabPanel("OMI",
		renderOMI(),
	fluidRow(
			renderBanner()
		)

	),

	# ---------------------------------------- #
	# MAPS
	# ---------------------------------------- #
	tabPanel("MAP",
	    htmlOutput("map")
		# leafletMap("map", width="100%", height="100%",
	 #        initialTileLayer = "//{s}.tiles.mapbox.com/v3/jcheng.map-5ebohr46/{z}/{x}/{y}.png",
	 #        initialTileLayerAttribution = HTML('Maps by <a href="http://www.mapbox.com/">Mapbox</a>'),
	 #        options=list(
	 #          center = c(37.45, -93.85),
	 #          zoom = 4,
	 #          maxBounds = list(list(15.961329,-129.92981), list(52.908902,-56.80481)) # Show US only
	 #        )
  #     	)
	),


	# ---------------------------------------- #
	# INFORMATION INTERFACE (NEEDS TOC)
	# ---------------------------------------- #
	tabPanel("About",
	  fluidRow(
			includeMarkdown("About.md")
		),

		fluidRow(
			renderBanner()
		)
	)


)))
# ---------------------------------------- #
# END OF SHINY UI
# ---------------------------------------- # 


