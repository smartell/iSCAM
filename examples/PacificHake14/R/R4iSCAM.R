#|-----------------------------------------------------------------------------------|#
#| R4iSCAM.R                                                                         |#
#| AUTHOR: Steven Martell                                                            |#
#|-----------------------------------------------------------------------------------|#
#| NOTES, if using sublime text 3.0
#|    use ⌘ + \ to set working directory to the directory of the current file.
#|    use ⌘ + enter to source a single line or highligted lines.
#|    use ⌘ + B to source entire file.

# |----------------------------------------------------------------------------------|
# | DEFINITIONS
# |----------------------------------------------------------------------------------|
# | .PWD       <- Global Parent Working Directory for R-scripts
# | .FIGUREDIR <- Directory for saving figures.
# | .RFILES    <- List of R functions to source from the lib directory.
.PWD        <- "~/Documents/iSCAM/examples/PacificHake14/R"
setwd(.PWD)
# .FIGUREDIR  <- "../FIGS/"
.FIGUREDIR  <- "../logo/"
.RFILES     <- list.files("./lib/",pattern="\\.[Rr]$")
.VIEWTRCK   <- "iscamViewTracker.txt"
.BOOLREADFN <- TRUE
require(ggplot2)
.THEME      <- theme_bw(10)
.UNITS      <- "(mlb)"

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
	ifiles     <- NULL
	if (trckExists)
	{
		cat( "MSG (.hCamViewSetup): Viewer tracking file ",.VIEWTRCK, " found.\n" )
		ifiles <- read.table( file = .VIEWTRCK,  as.is=TRUE,  header=TRUE,  sep="," )
		# ifiles <- tmp
		print(ifiles)
		# browser()
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