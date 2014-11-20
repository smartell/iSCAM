#gui_Equilibrium.R

# RENDER EQUILIBRIUM USER INTERFACE
buildEquilibriumGui <- function()
{
	print("Start")
	includeCSS("styles.css")
		fluidRow(
			tags$h4("Equilibrium Model")
		)
		
		renderEquilriumInterface()
		
}


renderEquilInputs <- function(prefix)
{

	wellPanel(
		# fluidRow(

		    # Directed fishery controls.
		    wellPanel( fluidRow(
		    tags$h5("Directed fishery"),
		    
			sliderInput(paste0(prefix,"_","selex_fishery"),"Fishery: 50% & 95% selectivity (inches)",min=15,max=60,value=c(34,40),step=1),
			sliderInput(paste0(prefix,"_","size_limit"),"Min & Max size limit (inches)",min=0,max=100,value=c(32,100),step=1),
			sliderInput(paste0(prefix,"_","discard_mortality_rate"),"Discard mortality rate",min=0,max=1,value=0.16,step=0.01)
			)),
			

			# Economic inputs
			wellPanel( fluidRow(
			tags$h5("Price per pound ($)"),
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
			)),


			# Bycatch panel
			wellPanel( fluidRow(
		      tags$h5("Bycatch controls"),
			  column(5,
			    sliderInput(paste0(prefix,"_","selex_bycatch"),"Acesending 50% & 95% (in.)",min=15,max=40,value=c(24,30),step=1)
			  ),
			  column(5,
			    sliderInput(paste0(prefix,"_","selex_bycatch_desc"),"Descending 95% & 50% (in.)",min=40,max=80,value=c(60,80),step=1)
			  ,offset=1),

			  # bycatch DMR
			  sliderInput(paste0(prefix,"_","bycatch_dmr"),"Bycatch discard mortality rate",min=0,max=1,value=0.06,step=0.01),

			  # bycatch amounts
			  column(5, offset=0,
			  	numericInput(paste0(prefix,"_","num_bycatch_total"), label = "Bycatch (Mlb)", value = 39)
			  ),
			  column(5,offset=1,
			    numericInput(paste0(prefix,"_","num_bycatch"), label = "Mortality (Mlb)", value = 8)
			  )
			  ) 
			  	
			)
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
					
					checkboxGroupInput("chkEquilPlot", "Plot Variable:",
	                   c("Yield" = "Yield",
	                     "Discard" = "Discard",
	                     "Waste" = "Waste",
	                     "Spawning Biomass" = "Spawning.Biomass",
	                     "SPR" = "SPR",
	                     "Recruitment"="Recruitment",
	                     "Average weight"="Avg.Weight",
	                     "Landed value (million $)"="Landed.Value"),
	                   	selected = "Yield", inline=TRUE),

					plotOutput("plot_equil",height = "550px")

				),

				tabPanel("Tables",
				    # selectInput('selMSYTable',"Select input columns",
				    #             c("Spawning.Biomass","SPR","Recruitment","Yield",
				    #               "Fe","Fbycatch","Numbers","Landed.Value"),
				    #             multiple=TRUE,
				    #             selected=c("Fe","Spawning.Biomass","Yield","SPR")),

				    checkboxGroupInput("chkMSYTable", "Plot Variable:",
	                   c("Fishing Rate" = "Fe",
	                     "Yield" = "Yield",
	                     "Discard" = "Discard",
	                     "Waste" = "Waste",
	                     "SSB" = "Spawning.Biomass",
	                     "SPR" = "SPR",
	                     "Recruitment"="Recruitment",
	                     "Avg. weight"="Avg.Weight",
	                     "Landed value (million $)"="Landed.Value"),
	                   	selected = c("Yield","Spawning.Biomass"), inline=TRUE),


				    tags$h5("Equilibrium values at MSY"),
					tableOutput('msyTable'),

					tags$h5("Eqiulibrium values at SPR=30%"),
					tableOutput('sprTable'),

					tags$h5("Equilibrium values at MEY"),
					tableOutput('meyTable')		

				
				),
				
				tabPanel("Selectivity",
					h5("Fishery selectiity"),
			        plotOutput("plotFishSelex",height=200),

			        h5("Bycatch selectivity"),
			        plotOutput("plotSelex",height=200)
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
		    sliderInput(paste0(prefix,"_","linf_dev"),"Deviation in max length (%)",min=-40,max=40,value=0,step=5),
		    sliderInput(paste0(prefix,"_","vonk_dev"),"Deviation in metabolic rate (%)",min=-20,max=20,value=0,step=5),
			sliderInput(paste0(prefix,"_","maternal_effect"),"Maternal effect",min=0,max=3,value=1,step=0.05)
		)
	)
}




