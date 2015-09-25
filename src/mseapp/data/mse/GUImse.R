 renderMSE <- function()
{
	tabPanel("MSE Dashboard",
	    buildMseGUI()
	)
}


buildMseGUI <- function()
{
	lblSldr_minSSB = tags$h5("Spawning biomass limit (10% of Bo)")
	lblSldr_limSSB = tags$h5("Spawning biomass threshold (20% of Bo)")
	lblSldr_thrSSB = tags$h5("Spawning biomass target (30% of Bo)")

	lblSldr_fshOPEN = tags$h5("Maintain directed fishery (x million pounds)")
	lblSldr_fshTARG = tags$h5("Maintain threshold average catch (x million pounds)")
	lblSldr_fshTARG = tags$h5("Maintain target average catch (x million pounds)")

	fluidRow(
		column(6,
		navlistPanel(widths=c(3,9),
	    # sidebarLayout(
			"Navigation",
			tabPanel("Objectives",
			    column(12,     
					tags$h3("Biological Objectives"),
					sliderInput(inputId="sldr_minSSB",lblSldr_minSSB,min=0,max=1,value=0.99,step=0.01),
					sliderInput(inputId="sldr_limSSB",lblSldr_limSSB,min=0,max=1,value=0.95,step=0.01),
					sliderInput(inputId="sldr_thrSSB",lblSldr_thrSSB,min=0,max=1,value=0.75,step=0.01)
				)
			 #    column(6,     
				# 	tags$h3("Fishery Objectives"),
				# 	sliderInput(inputId="sldr_fshOPEN",lblSldr_fshOPEN,min=0,max=1,value=0.99,step=0.01),
				# 	sliderInput(inputId="sldr_limSSB",lblSldr_limSSB,min=0,max=1,value=0.95,step=0.01),
				# 	sliderInput(inputId="sldr_thrSSB",lblSldr_thrSSB,min=0,max=1,value=0.75,step=0.01)
				# )
				# column(4,     
				# 	tags$h3("Biological Objectives"),
				# 	sliderInput(inputId="sldr_minSSB",lblSldr_minSSB,min=0,max=1,value=0.99,step=0.01),
				# 	sliderInput(inputId="sldr_limSSB",lblSldr_limSSB,min=0,max=1,value=0.95,step=0.01),
				# 	sliderInput(inputId="sldr_thrSSB",lblSldr_thrSSB,min=0,max=1,value=0.75,step=0.01)
				# )
			),
			tabPanel("Procedures",
				h3("Management Procedures"),

				
				selectInput("siProcedure", h5("Procedure:"),
				    choices  = prc,
				    selected = prc[1],
            		multiple=TRUE),

				# harvest control rules
				sliderInput("sldr_Urate",h5("Exploitation rate"),
							min   =0.1,
							max   =0.3,
							value =0.2,
							step  =0.01),

				sliderInput("sldr_Brefp",h5("Depletion (limit, threshold)"),
							min   =0,
							max   =1,
							value =c(0.2,0.3),
							step  =0.05),

				# size limits
				sliderInput('sldr_Slims',h5("Size limits (min,max) - inches"),
				            min   = 0,
				            max   = 100,
				            value = c(32,100),
				            step  = 2)
				# wellPanel(fluidRow(
				#     column(6,
				# 		checkboxGroupInput("dataSelect", "DATA:",
	   #                 	 c("Setline CPUE" = "cpue",
	   #                     "Commercial WPUE" = "wpue",
	   #                     "Age Composition" = "comps"),inline=FALSE),
						
				# 		radioButtons("methodSelect", "ASSESSMENT METHOD:",
			 #             c("Simple" = "simple",
			 #               "SCA" = "sca",
			 #               "VPA" = "vpa",
			 #               "Kalman Filter" = "kf"),inline=FALSE)
				# 	),
				# 	column(6,
				# 		radioButtons("hcrSelect","HARVEST CONTROL RULE:",
			 #             c("Fixed F" = "fixedF",
			 #               "Sloping F"="slopingF"),inline=FALSE)
				# 	)
				# ))
			),

			tabPanel("Scenarios",
			    h3("Operating Model Scenarios"),

			    
				selectInput("siScenario", h5("Scenarios:"),
				    choices  = scn,
				    selected = scn[1],
            		multiple=TRUE),

				
				selectInput("siRecruitment",h5("Recruitment:"),
				    choices  = rcn,
				    selected = rcn,
				    multiple = TRUE),

				
				selectInput("siGrowth",h5("Size-at-age:"),
				    choices  = gcn,
				    selected = gcn[1],
				    multiple = TRUE)
			),

			tabPanel("Performance",
			    h3("Performance Metrics")
			),
			"Dashboard"
		)),

		# Right side devoted to output.
		column(6,
			tabsetPanel(type="pills",id="mseTab",
				tabPanel("Figures",

				    selectInput("siPlotType","Select variable",
				                c("Spawning biomass" = "ssb",
				                  "SSB Depletion"    = "sbd",
				                  "Fishery landings" = "fct"),
				                selected="ssb"),

					plotOutput("plotMSE"),

					hr(),

					
					sliderInput("sldr_year_range",h5("Trim years"),
					            min=syr,max=nyr,value=range(yrs),
					            sep="",width=600)
				),

				tabPanel("Performance Metrics",
					h5("Median depletion"),
					tableOutput("viewDepletionTable"),

					h5("Probability of falling below 20% unfished"),
					tableOutput("viewSSBLimitTable"),

					h5("Probability of falling below 30% unfished"),
					tableOutput("viewSSBThresholdTable"),					

					h5("Median catch"),
					tableOutput("viewCatchTable"),

					h5("5-year Average Annual Variation in Catch"),
					tableOutput("viewAAVTable")
				)
			)
		)
		
	)

}