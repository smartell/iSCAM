#gui_Equilibrium.R

.PLOTVARS = c(
			"Yield"                    = "Yield",
			"Total Mortality"          = "Total.Mortality",
			"Discard"                  = "Discard",
			"Waste"                    = "Waste",
			"Fbycatch"                 = "Fbycatch",
			"Spawning Biomass"         = "Spawning.Biomass",
			"SPR"                      = "SPR",
			"YPR"                      = "YPR",
			"Recruitment"              ="Recruitment",
			"Average weight"           ="Avg.Weight",
			"Landed value (million $)" ="Landed.Value"
			)

.TABLEVARS = c(
			"Fishing Rate"             = "Fishing.Mortality",
			"Total Mortality"          = "Total.Mortality",
			"FCEY"                     = "FCEY",
			"Discard"                  = "Discard",
			"Waste"                    = "Waste",
			"SSB"                      = "Spawning.Biomass",
			"SPR"                      = "SPR",
			"YPR"                      = "YPR",
			"Recruitment"              ="Recruitment",
			"Avg. weight"              ="Avg.Weight",
			"Landed value (million $)" ="Landed.Value",
			"Legal numbers caught"     = "Legal.numbers",
			"Sublegal numbers caught"  = "Sublegal.numbers",
			"Sublegals per 100"        = "Sublegal.100"
	           )


# RENDER EQUILIBRIUM USER INTERFACE
renderEquil <- function()
{
	tabPanel("Equilibrium",
	    buildEquilibriumGui()
	)
}


buildEquilibriumGui <- function()
{
	print("Start")
	#includeCSS("styles.css")
		# fluidRow(
		# 	tags$h4("Equilibrium Model")
		# )
		
		renderEquilriumInterface()
		
}


renderEquilInputs <- function(prefix)
{

	wellPanel(
		# fluidRow(

		    # Directed fishery controls.
		    wellPanel( fluidRow(
		    tags$h5("Directed fishery"),
		    
			sliderInput(paste0(prefix,"_","selex_fishery"),tags$h6("Fishery: 50% & 95% selectivity (in.)"),min=15,max=60,value=c(34,40),step=1),
			sliderInput(paste0(prefix,"_","size_limit"),tags$h6("Min & Max size limit (in.)"),min=0,max=100,value=c(32,100),step=1),
			sliderInput(paste0(prefix,"_","discard_mortality_rate"),tags$h6("Discard mortality rate"),min=0,max=1,value=0.16,step=0.01)
			)),
			

		    # Bycatch panel
			wellPanel( fluidRow(
		      tags$h5("Bycatch controls"),
			  column(6,
			    sliderInput(paste0(prefix,"_","selex_bycatch"),tags$h6("Acesending 50% & 95% (in.)"),min=15,max=40,value=c(24,30),step=1)
			  ),
			  column(6,
			    sliderInput(paste0(prefix,"_","selex_bycatch_desc"),tags$h6("Descending 95% & 50% (in.)"),min=40,max=80,value=c(60,80),step=1)
			  ),

			  # bycatch DMR
			  sliderInput(paste0(prefix,"_","bycatch_dmr"),tags$h6("Bycatch discard mortality rate"),min=0,max=1,value=0.09,step=0.01),

			  # bycatch amounts
			  column(6, offset=0,
			  	numericInput(paste0(prefix,"_","num_bycatch_total"), label = tags$h6("Bycatch (Mlb)"), value = 39)
			  ),
			  column(6, offset=0,
			    numericInput(paste0(prefix,"_","num_bycatch"), label = tags$h6("Mortality (Mlb)"), value = 8)
			  )
			  ) 
			  	
			),



			# Economic inputs
			wellPanel( fluidRow(
			tags$h5("Price per pound ($)"),
			column(2,
				numericInput(paste0(prefix,"_","five"), label = h6("5-10"), value = 0, step=0.10)
			,offset=0),
			column(2,
				numericInput(paste0(prefix,"_","ten"), label = h6("10-20"), value = 5, step=0.10)
			,offset=1),
			column(2,
				numericInput(paste0(prefix,"_","twenty"), label = h6("20-40"), value = 5, step=0.10)
			,offset=1),
			column(2,
				numericInput(paste0(prefix,"_","forty"), label = h6("40+"), value = 5, step=0.10)
			,offset=1)
			))


			
		# )
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
	  	# checkboxInput("chk_viewSelex","View selectivity curves"),
			tabsetPanel(type="pills",id="eqtab",
				tabPanel("Plots",
				 
					# selectInput('selEquilPlot',"Equilibrium Vs Fishing Intensity",
				 #              c("Yield","Discard","Waste","SPR","Spawning.Biomass",
				 #                "Recruitment","Avg.Weight","Landed.Value"),
				 #              multiple=TRUE,selected="Yield"),
					
					checkboxGroupInput("chkEquilPlot", tags$b("Select Variable to display in plot"),
							.PLOTVARS,
	                   	selected = "Yield", inline=FALSE),

					tags$hr(),

					plotOutput("plot_equil",height = "750px")

				),

				tabPanel("Tables",
				    # selectInput('selMSYTable',"Select input columns",
				    #             c("Spawning.Biomass","SPR","Recruitment","Yield",
				    #               "Fe","Fbycatch","Numbers","Landed.Value"),
				    #             multiple=TRUE,
				    #             selected=c("Fe","Spawning.Biomass","Yield","SPR")),
				    # column(4,
				    checkboxGroupInput("chkMSYTable", tags$b("Select Variable to display in tables"),
	                   .TABLEVARS,
	                   	selected = c("Yield","Spawning.Biomass"), inline=FALSE),

				    tags$hr(),
				    # ),
				    # column(8,
				    tags$h5("Equilibrium values at MSY"),
					tableOutput('msyTable'),

					tags$h5("Eqiulibrium values at SPR=30%"),
					tableOutput('sprTable'),

					tags$h5("Equilibrium values at MEY"),
					tableOutput('meyTable')		
					# ),
				
				),
				
				tabPanel("Selectivity",
					h5("Fishery selectiity"),
			        plotOutput("plotFishSelex",height=300),

			        h5("Bycatch selectivity"),
			        plotOutput("plotSelex",height=300)
				),

				tabPanel("Size-at-age",
				    h5("Size-at-age for a given area"),
				    plotOutput("plotSAA",height=500)
				)
			)
		)
	)
	  
	


}



renderShutter <-function()
{
	
	conditionalPanel(" input.chk_viewSelex == true ",
	absolutePanel(id = "controls", class = "modal", fixed = TRUE, draggable = TRUE,
        top = 40, left = "auto", right = 20, bottom = "auto",
        width = 330, height = "auto",
        
        h4("Selectivity Explorer"),
        
        h5("Fishery selectiity"),
        plotOutput("plotFishSelex",height=200),

        h5("Bycatch selectivity"),
        plotOutput("plotSelex",height=200)
        
        # selectInput("color", "Color", vars),
        # selectInput("size", "Size", vars, selected = "adultpop"),
        # conditionalPanel("input.color == 'superzip' || input.size == 'superzip'",
        #   # Only prompt for threshold when coloring or sizing by superzip
        #   numericInput("threshold", "SuperZIP threshold (top n percentile)", 5)
        # ),
        
        # plotOutput("histCentile", height = 200),
        # plotOutput("scatterCollegeIncome", height = 250)
	))
	
}


renderStockInputs <- function(prefix)
{
	wellPanel(
		fluidRow(
		    selectInput(paste0(prefix,"_","regArea"),"Regulatory Area 2014 Growth Parameters",
		                .REGAREA),

		    sliderInput(paste0(prefix,"_","linf_dev"),"Deviation in max length (%)",min=-40,max=40,value=0,step=5),
		    sliderInput(paste0(prefix,"_","vonk_dev"),"Deviation in metabolic rate (%)",min=-20,max=20,value=0,step=5),
			sliderInput(paste0(prefix,"_","maternal_effect"),"Maternal effect",min=0,max=3,value=1,step=0.05)

		)
	)
}




