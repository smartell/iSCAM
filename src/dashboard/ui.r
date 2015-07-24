# --------------------------------------------------------------------------- #
# Steve Martell
# May 31, 2015
#
# --------------------------------------------------------------------------- #

if(!require("shinydashboard"))      install.packages("shinydashboard")
if(!require("dygraphs"))            install.packages("dygraphs")
source('helpers.R')





header <- dashboardHeader(
  title = "MSAB Dashboard"
)

sidebar <- dashboardSidebar(
  buildSideBarMenu,
  hr()

)

body <- dashboardBody(
  tabItems(
    tabItem(tabName="dashboard",
    	fluidRow(
    		valueBoxOutput("mitigationBox")
    	)
    ),

    # tabItem(tabName="harvest_policy"),



    tabItem(tabName="equilibrium",
      fluidRow(
        box( title="Graphics" ,width=12, status="info", solidHeader=TRUE )
      ),
      fluidRow(
        box( width = 12, status = "primary", solidHeader = TRUE, title="Equilibrium Model",
    		## -------------- ##
		    # Equil Interface  #
		    ## -------------- ##
		    renderEquil()
    	))
    ),

    tabItem(tabName="total_mortality",
      buildUItotalMortality()
    ),



    tabItem(tabName="ui",
    	box( width = NULL, status = "primary", solidHeader = TRUE, title="ui.R",
                 downloadButton('downloadData2', 'Download'),
                 br(),br(),
                 pre(includeText("ui.R"))
            )
    ),
    tabItem(tabName="server",
    	box( width = NULL, status = "primary", solidHeader = TRUE, title="ui.R",
                 downloadButton('downloadData2', 'Download'),
                 br(),br(),
                 pre(includeText("server.R"))
            )   
    ),
    tabItem(tabName="code")
  )
)

ui <- dashboardPage(skin = "yellow",
  header,
  sidebar,
  body
)

