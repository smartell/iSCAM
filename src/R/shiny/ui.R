source("helpers.R")


# RENDER EQUILIBRIUM USER INTERFACE
renderEquilInputs <- function(prefix)
{

	wellPanel(
		fluidRow(
			# column(6,
				sliderInput(paste0(prefix,"_","selex_fishery"),"Fishery: 50% & 95% selectivity (inches)",min=15,max=60,value=c(34,40),step=1),
				sliderInput(paste0(prefix,"_","size_limit"),"Min & Max size limit (inches)",min=0,max=100,value=c(32,100),step=1),
				sliderInput(paste0(prefix,"_","discard_mortality_rate"),"Discard mortality rate",min=0,max=1,value=0.16,step=0.01),
				sliderInput(paste0(prefix,"_","spr_target"),"SSB limit-threshold reference",min=0,max=1,value=c(0.2,0.3),step=0.05),
			# ),
			# column(6,
				sliderInput(paste0(prefix,"_","selex_bycatch"),"Bycatch: 50% & 95% selectivity (inches)",min=15,max=60,value=c(24,40),step=1),
				numericInput(paste0(prefix,"_","selex_asymptote"),"Asymptote",value=0.65,step=0.05),
				numericInput(paste0(prefix,"_","num_bycatch"), label = "Bycatch mortality cap (Mlb)", value = 8),
				

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
		)
	)

}

renderStockInputs <- function(prefix)
{
	wellPanel(
		fluidRow(
		    sliderInput(paste0(prefix,"_","linf_dev"),"Deviation in max length (%)",min=-40,max=40,value=0,step=5),
		    sliderInput(paste0(prefix,"_","vonk_dev"),"Deviation in metabolic rate (%)",min=-20,max=20,value=0,step=5),
			sliderInput(paste0(prefix,"_","maternal_effect"),"Maternal effect",min=0,max=3,value=1,step=0.05)
		)
	)
}


renderEquilriumInterface <- function()
{
	
	

	fluidRow(
	  
	  column(6,
		  tabsetPanel(type="pills",id="eqpar",
		    tabPanel("Procedures",
					column(6,
					  tags$h4("Procedure A"),
					  renderEquilInputs("A")
					), 

					column(6,
					  tags$h4("Procedure B"),
					  renderEquilInputs("B")
					)
				),
				tabPanel("Scenarios",
					column(6,
						tags$h4("Scenario A"),
						renderStockInputs("A")
					),
					column(6,
						tags$h4("Scenario B"),
						renderStockInputs("B")
					)
				)
			)
		),


		column(6,
		  tabsetPanel(type="pills",id="eqtab",
		    tabPanel("Plots",
		      # selectInput('selChartType',"Model Output",
							 #      c("Equilibrium Yield",
							 #        "Performance Metrics at MSY",
							 #        "Equilibrium Value",
							 #        "Value at MSY")),
		      selectInput('selEquilPlot',"Equilibrium Vs Fishing Intensity",
		                  c("Yield","Discard","Waste","SPR","Spawning.Biomass","Recruitment","Avg.Weight"),
		                  multiple=TRUE,selected="Yield"),
	      	# plotOutput("a_equilPlot", height = "500px"),
	      	plotOutput("plot_equil",height = "500px")

	      ),
	      tabPanel("Tables",
	        tags$p("Biological sustainability"),
	        tableOutput('table_biological'),
	        tags$p("Fisheries sustainability at SSB-threshold. TCEY=(O26 bycatch)+(Wastage)+(FCEY)"),
	        tableOutput('table_fishery'),
	        tags$p("Economic performance at SSB-threshold (million $)"),
	        tableOutput('table_economics'),
	        tags$p("MSY-based reference points"),
	     		tableOutput('msytable')
	   		)
	    )
	  )
	  
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
			tags$h3("Equilibrium Model: reference points")
	  ),

	  fluidRow(
			renderEquilriumInterface()
	  ),

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


