#|-----------------------------------------------------------------------------------|#
#| R4iSCAM.R                                                                         |#
#| AUTHOR: Steven Martell                                                            |#
#|-----------------------------------------------------------------------------------|#
#| NOTES, if using sublime text 3.0
#|    use ? + \ to set working directory to the directory of the current file.
#|    use ? + enter to source a single line or highligted lines.
#|    use ? + B to source entire file.

# |----------------------------------------------------------------------------------|
# | DEFINITIONS
# |----------------------------------------------------------------------------------|
# | .PWD       <- Global Parent Working Directory for R-scripts
# | .FIGUREDIR <- Directory for saving figures.
# | .RFILES    <- List of R functions to source from the lib directory.

# .FIGUREDIR  <- "../FIGS/"
require(ggplot2)
.PWD        <- getwd()
.FIGUREDIR  <- "../FIGS/"
.RFILES     <- list.files("./lib/",pattern="\\.[Rr]$")
.VIEWTRCK   <- "iSCAMViewTracker.txt"
.BOOLREADFN <- TRUE
.THEME      <- theme_bw(12)
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
	if (trckExists)
	{
		cat( "MSG (.hCamViewSetup): Viewer tracking file ",.VIEWTRCK, " found.\n" )
	}
		# if(exists(".PBSmod")) closeWin(win)
	tmp    <- read.table( file = .VIEWTRCK,  as.is=TRUE,  header=TRUE,  sep="," )
	ifiles <<- as.data.frame(tmp)

	createWin("iscamWin2.txt")


	setWinVal(list(txtFigDir=.FIGUREDIR))
}



cat("Type: \n guiView()\n to start the iSCAM gui \n")
