# ----------------------------------------#
# MAIN USER INTERFACE FOR THE APPLICATION #
# ----------------------------------------#
shinyUI(fluidPage(
        navbarPage("IPHC Shiny App (MSEApp)",
                   id="nav", collapsible=TRUE,
                   footer=img(src="iphclogo.png",  height = 30, width = 70),

    ## -------------- ##
    # About Interface  #
    ## -------------- ##
    renderAbout(),

    ## -------------- ##
    # Equil Interface  #
    ## -------------- ##
    renderEquil(),

    ## -------------- ##
    # MSE Interface  #
    ## -------------- ##
    renderMSE()


)))