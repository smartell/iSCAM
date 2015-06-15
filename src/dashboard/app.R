## app.R ##
library(shinydashboard)
library(dygraphs)

header <- 
dashboardHeader(
  title = "MSAB Dashboard",
  dropdownMenu(type = "messages",
  messageItem(
    from = "Sales Dept",
    message = "Sales are steady this month."
  ),
  messageItem(
    from = "New User",
    message = "How do I register?",
    icon = icon("question"),
    time = "13:45"
  ),
  messageItem(
    from = "Support",
    message = "The new server is ready.",
    icon = icon("life-ring"),
    time = "2014-12-01"
  )
))




sidebar <- ## Sidebar content
dashboardSidebar(
  hr(),
  sidebarMenu(
    menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
    menuItem("Widgets", tabName = "widgets", icon = icon("th")),
    menuItem("Harvest Policy", tabName = "harvestpolicy", icon = icon("bar-chart-o")),
    menuItem("About", tabName = "about", icon = icon("question"),
      menuSubItem("ui.R", tabName = "ui", icon = icon("angle-right")),
      menuSubItem("server.R", tabName = "server", icon = icon("angle-right"))
    )
  )
)


body <-   ## Body content
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
      # valueBoxes
        fluidRow(
          valueBox(
            uiOutput("valBoxDepletion"), "Median depletion", icon = icon("cog", lib = "glyphicon")
          ),
          valueBox(
            uiOutput("progress"), "Progress", icon = uiOutput("progressIcon"),
            color = "purple"
          ),
          # An entire box can be in a uiOutput
          uiOutput("approvalBox")
        ),
        fluidRow(
          box(title="Histogram",
              # box(
            plotOutput("plot1", height = 250)),
              box(
            plotOutput("plot2", height = 250))
          ),

          box(
            title = "Controls",
            sliderInput("slider", "Number of observations:", 0, 500, 150)
          ),
          box(
          	title = "Funnel Plots",
          	status = "primary",
          	solidHeader = TRUE,
          	collapsible = TRUE,
          	dygraphOutput("dygraph", height = 250)
          )
        )
      ),

      # Second tab content
      tabItem(tabName = "widgets",
        h2("Widgets tab content")
      ),

      # UI code
      tabItem(tabName = "ui",
            box( width = NULL, status = "primary", solidHeader = TRUE, title="ui.R",
                 downloadButton('downloadData2', 'Download'),
                 br(),br(),
                 pre(includeText("app.R"))
            )
      )
    )







ui <- dashboardPage(skin = "yellow",
  header,
  sidebar,
  body
)



server <- function(input, output) {
  set.seed(122)
  histdata <- rnorm(500)
  output$plot1 <- renderPlot({
    data <- histdata[seq_len(input$slider)]
    hist(data)
    output$valBoxDepletion <- renderText({median(data)})
  })

  # Predicted values for dygraph
  predicted <- reactive({
    hw <- HoltWinters(ldeaths)
    predict(hw, n.ahead = input$slider,
            prediction.interval = TRUE,
            level = as.numeric(0.95))
  })

  output$dygraph <- renderDygraph({
  	lungDeaths <- cbind(mdeaths, fdeaths)
    dygraph(predicted(), main = "Predicted Deaths/Month") %>%
      dySeries(c("lwr", "fit", "upr"), label = "Deaths") %>%
      dyOptions(drawGrid = TRUE)
  })
}

shinyApp(ui, server)


