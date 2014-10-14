library(shiny)
library(markdown)

# RENDER EQUILIBRIUM USER INTERFACE
renderEquilInputs <- function(prefix)
{

	wellPanel(
		fluidRow(
			# column(6,
				sliderInput(paste0(prefix,"_","selex_fishery"),"Fishery: 50% & 95% selectivity (inches)",min=15,max=60,value=c(34,40),step=1),
				sliderInput(paste0(prefix,"_","size_limit"),"Min & Max size limit (inches)",min=0,max=100,value=c(32,100),step=1),
				sliderInput(paste0(prefix,"_","discard_mortality_rate"),"Discard mortality rate",min=0,max=1,value=0.16,step=0.01),
				# sliderInput(paste0(prefix,"_","max_size_limit"),"Maximum size limit (inches)",min=0,max=100,value=100,step=1)
			# ),
			# column(6,
				sliderInput(paste0(prefix,"_","selex_bycatch"),"Bycatch: 50% & 95% selectivity (inches)",min=15,max=60,value=c(24,40),step=1),
				numericInput(paste0(prefix,"_","num_bycatch"), label = "Bycatch cap (Mlb)", value = 8),
				

				# Economic inputs
				#fluidRow(
				tags$p("Price per pound ($)"),
				column(2,
					numericInput(paste0(prefix,"_","five"), label = "5-10", value = 0, step=0.10)
				,offset=0),
				column(2,
					numericInput(paste0(prefix,"_","ten"), label = "10-20", value = 5, step=0.10)
				,offset=1),
				column(2,
					numericInput(paste0(prefix,"_","twenty"), label = "20-40", value = 5, step=0.10)
				,offset=1),
				column(2,
					numericInput(paste0(prefix,"_","forty"), label = "40+", value = 5, step=0.10)
				,offset=1)
				#)
				# selectInput(paste0(prefix,"_",'chartType'),"Model Output",
				#             c("Equilibrium Yield",
				#               "Performance Metrics at MSY",
				#               "Reference Points",
				#               "Selectivity curves"))

				# absolutePanel(id=paste0(prefix,"_","selex_panel"),class="modal",
				#               draggable=TRUE,fixed=FALSE,cursor="auto",
				#               top="auto",left="auto",right="auto",bottom="auto",
				#               width=300,
				#   h4("Selectivity curves"),
				#   plotOutput(paste0(prefix,"_","selex"), height = "200px"),
				#   style = "opacity: 0.60"
				# )
			# ),

			# p(actionButton(paste0(prefix, "_", "recalc"),
   #    "Re-run scenario", icon("random")
   #  	))
		)
	)

}

  # absolutePanel(id = "controls", class = "modal", fixed = TRUE, draggable = TRUE,
  #       top = 60, left = "auto", right = 20, bottom = "auto",
  #       width = 330, height = "auto",
        

renderInputs <- function(prefix) {
  wellPanel(
    fluidRow(
      column(6,
        sliderInput(paste0(prefix, "_", "n_obs"), "Number of observations (in Years):", min = 0, max = 40, value = 20),
        sliderInput(paste0(prefix, "_", "start_capital"), "Initial capital invested :", min = 100000, max = 10000000, value = 2000000, step = 100000, format="$#,##0", locale="us"),
        sliderInput(paste0(prefix, "_", "annual_mean_return"), "Annual investment return (in %):", min = 0.0, max = 30.0, value = 5.0, step = 0.5),
        sliderInput(paste0(prefix, "_", "annual_ret_std_dev"), "Annual investment volatility (in %):", min = 0.0, max = 25.0, value = 7.0, step = 0.1)
      ),
      column(6,
        sliderInput(paste0(prefix, "_", "annual_inflation"), "Annual inflation (in %):", min = 0, max = 20, value = 2.5, step = 0.1),
        sliderInput(paste0(prefix, "_", "annual_inf_std_dev"), "Annual inflation volatility. (in %):", min = 0.0, max = 5.0, value = 1.5, step = 0.05),
        sliderInput(paste0(prefix, "_", "monthly_withdrawals"), "Monthly capital withdrawals:", min = 1000, max = 100000, value = 10000, step = 1000, format="$#,##0", locale="us",),
        sliderInput(paste0(prefix, "_", "n_sim"), "Number of simulations:", min = 0, max = 2000, value = 200)
      )
    ),
    p(actionButton(paste0(prefix, "_", "recalc"),
      "Re-run simulation", icon("random")
    ))
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
	
	# INFORMATION INTERFACE (NEEDS TOC)
	tabPanel("About",
	  fluidRow(
			includeMarkdown("About.md")
		),

		fluidRow(
			renderBanner()
		)
	),

	# EQUILIBRIUM INTERFACE
	tabPanel("Equilibrium",
	  fluidRow(
	    column(6,
		  	tags$h3("Equilibrium Model: reference points")
	   	),
	   	column(6,
      	selectInput('selChartType',"Model Output",
            c("Equilibrium Yield",
              "Performance Metrics at MSY",
              "Equilibrium Value",
              "Value at MSY"))
	   	)
	  ),
	  fluidRow(
	  	column(3,
	  	  tags$h4("Scenario A"),
	  	  renderEquilInputs("a")
	  	),
	  	column(3,
	  	  tags$h4("Scenario B"),
	  	  renderEquilInputs("b")
	  	),
	  	column(6,
	      plotOutput("a_equilPlot", height = "550px")
	    )
	    
	  ),
	  # fluidRow(
   #  	column(6, tags$h4("Scenario A")),
   #  	column(6, tags$h4("Scenario B"))
  	# ),
	  # fluidRow(
	  # 	column(6,renderEquilInputs("a")),
	  # 	column(6,renderEquilInputs("b"))
	  # ),

	  # fluidRow(
	  #   column(6,
	  #     plotOutput("a_equilPlot", height = "500px")
	  #   ),
	  #   column(6,
	  #     # plotOutput("b_equilPlot", height = "500px")
	  #   	tableOutput("b_table")
	  #   )
	  # ),
		fluidRow(
			renderBanner()
		)
	),

	# MANAGEMENT STRATEGY EVALUATION INTERFACE
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

  # OPERATING MODEL INTERFACE
  tabPanel("OMI",
   	renderOMI(),
    fluidRow(
  		renderBanner()
  	)

  )



	
)))


