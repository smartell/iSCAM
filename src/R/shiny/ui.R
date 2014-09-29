library(shiny)
library(markdown)

# RENDER EQUILIBRIUM USER INTERFACE
renderEquilInputs <- function()
{
	fluidRow(

	  column(3,
			wellPanel("Equilibrium model inputs",
			  br(),
				sliderInput("sldr_fe","Fishing rate (fe)",
				            min = 0,
				            max = 1.0,
				            value = 0.10,
				            format="#.##"),

				 checkboxGroupInput("chck_plot","View plot",
				                    c("Yield","Biomass","Recruitment"),selected ="Yield")
			) 
		),
		column(9)
	)

}


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

renderBanner <- function()
{
	wellPanel(
		# Logo image
	  column(10),
	  column(2,
	  	img(src="iphclogo.png",  height = 80, width = 80),
			img(src="iscamLogo.png", height = 80, width = 80)
	  )
	)
}

# ----------------------------------------#
# MAIN USER INTERFACE FOR THE APPLICATION #
# ----------------------------------------#
shinyUI(fluidPage(navbarPage("IPHC MSE TOOL",
	
	tabPanel("About",
	  fluidRow(
			includeMarkdown("About.md")
		),

		fluidRow(
			renderBanner()
		)
	),


  tabPanel("MSE",
	  # titlePanel("IPHC Management Strategy Evaluation App"),

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

	tabPanel("Equilibrium",

	  fluidRow(
	  	renderEquilInputs()
	  ),

		fluidRow(
			renderBanner()
		)
	)



	
)))


