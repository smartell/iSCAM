#sideBarMenu.R


# Sidebar Menu
buildSideBarMenu <- sidebarMenu(
	
	# Main Dashboard
	menuItem("Dashboard",tabName="dashboard",icon=icon("dashboard")),


	# Harvest Policy interfaces
	menuItem("Harvest Policy",tabName="harvest_policy",icon=icon("cog"),
		menuSubItem("Equilibrium",tabName="equilibrium",icon=icon("line-chart")),
		menuSubItem("Total Mortality",tabName="total_mortality",icon=icon("bar-chart"))
	),

	# Jane Sullivan's Size-at-age 
	menuItem("Size-at-age",tabName="saa"),

	# Source code for geeks
	menuItem("Code",tabName="code",icon=icon("file-text-o"),
		menuSubItem("ui.R", tabName = "ui", icon = icon("angle-right")),
		menuSubItem("server.R", tabName = "server", icon = icon("angle-right"))
	)
)