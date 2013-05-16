#|-----------------------------------------------------------------------------------|#
#| R4iSCAM.R                                                                         |#
#| AUTHOR: Steven Martell                                                            |#
#|-----------------------------------------------------------------------------------|#
#| 


# |----------------------------------------------------------------------------------|
# | DEFINITIONS
# |----------------------------------------------------------------------------------|
# | .PWD       <- Global Parent Working Directory for R-scripts
# | .FIGUREDIR <- Directory for saving figures.
# | .RFILES    <- List of R functions to source from the lib directory.
.PWD        <- "~/Documents/iSCAM-project/src/R"
.FIGUREDIR  <- "../FIGS/"
.RFILES     <- list.files("./lib/",pattern="\\.[Rr]$")
.VIEWTRCK   <- "iSCAMViewTracker.txt"
.BOOLREADFN <- TRUE
.THEME      <- theme_bw(12)

# |----------------------------------------------------------------------------------|
# | guiView: Main function for starting the iSCAM Gui
# |----------------------------------------------------------------------------------|
# |
guiView  <- function()
 {
	setwd(.PWD)
	for(nm in .RFILES) source(file.path("./lib", nm), echo=FALSE)
	.gui2("iSCAMView") 	
 }


# |----------------------------------------------------------------------------------|
# | .gui2: creates gui interface and reads .VIEWTRCK file
# |----------------------------------------------------------------------------------|
# | win    <- name of X11 window to close if one already exists.
# | [ ]    <- TODO, if .VIEWTRCK does not exist write header for file.
.gui2	<- function(win)
{
	require(PBSmodelling)
	
	trckExists <- file.exists( .VIEWTRCK )
	
	if (trckExists)
	{
		cat( "MSG (.hCamViewSetup): Viewer tracking file ",.VIEWTRCK, " found.\n" )
		tmp    <- read.table( file = .VIEWTRCK,  as.is=TRUE,  header=TRUE,  sep="," )
		ifiles <- tmp
	}
	else
	{
		cat( .VIEWTRCK," does not exist, please check file name. ")
	}
	if(exists(".PBSmod")) closeWin(win)
	createWin("iscamWin2.txt")
	
	setWinVal(list(txtFigDir=.FIGUREDIR))
}



cat("Type: \n guiView()\n to start the iSCAM gui")