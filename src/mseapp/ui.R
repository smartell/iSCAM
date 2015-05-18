# ----------------------------------------#
# MAIN USER INTERFACE FOR THE APPLICATION #
# ----------------------------------------#
shinyUI(fluidPage(theme = shinytheme("cosmo"),
        navbarPage("IPHC Shiny App (ISA)",
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