library(shiny)

# Define UI for application that draws a histogram
#shinyUI(navbarPage(
shinyUI(fluidPage(
	  titlePanel("IPHC MSE-tool"),

	fluidRow(
		column(3,
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
		              		multiple = TRUE),

				selectInput('plotType',"Select variable to plot",
			        		c("Spawning biomass",
			        		  "Depletion",
			        	  	   "Catch",
			        	  	   "Sub-legal Catch",
			        	  	   "AAV in Catch",
			        	  	   "Wastage",
			        	  	   "Efficiency",
			        	  	   "Fishing mortality"),
			        		selected="Spawning biomass")
			),
			wellPanel(
				selectInput("icolor","Color Variable",
				            c(None = ".",
				              "Scenario",
				              "Procedure",
				              "gear",
				              "area",
				              "sex",
				              "group"),
				            selected="."),

				selectInput("facet_row","Facet row",
				            c(None = ".",
				              "Scenario",
				              "Procedure",
				              "gear",
				              "area",
				              "sex",
				              "group"),
				            selected="Procedure"),

				selectInput("facet_col","Facet column",
				            c(None = ".",
				              "Scenario",
				              "Procedure",
				              "gear",
				              "area",
				              "sex",
				              "group"),
				            selected="Scenario")
			)
		),
		column(9,
		    tabsetPanel(type="tabs",
			    # Funnel plots
			    tabPanel("Tulip plots", 
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
			),
			wellPanel(
				# Logo image
			    img(src="iphclogo.png", height = 100, width = 100),
				img(src="iscamLogo.png", height = 100, width = 100)
			)
		)
	)
	# tabPanel HELP
	# tabPanel
	# (
	#  	"Help",
	#  	h3("Introduction"),
	#  	p("Welcome to the IPHC Managment Strategy Evaluation Toolbox. The purpose of this
	#  	  tool is to explore alternative managment procedures for Pacific halibut and 
	#  	  evaluating the performance of each procedure with respect to the management
	#  	  objectives."),

	#  	h4("Toolbox menu"),
	#  	p("Under the toolbox menu item, there are two submenus: Proceudres and Scenarios.
	#  	  Default management procedures and scenarios are already selected, but you may 
	#  	  wish to select or unselect alternative procedures and scenarios to explore. By 
	#  	  definition, procedures are things that you can manage (e.g., size-limits), and 
	#  	  scenarios are things that you cannot manage (e.g., recruitment variability)."),

	#  	h4("Plots menu"),
	#  	p("Select the plots menu to view the chosen combinations of management procedures 
	#  	  and scenarios.  The sidebar allows you to explore a variety of different 
	#  	  plotting options."),

	#  	h4("Tables menu"),
	#  	p("Select the table menu to view summary performance metrics of the chosen 
	#  	  combinations of management procedures and scenarios. Again the sidebar will 
	#  	  allows you to customize the summary statistics.")
	# ),

	
))


