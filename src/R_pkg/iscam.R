##                                      ##
##  R-script for the iSCAM project                                           ##
##  Author: Steven Martell,  University of British Columbia                  ##
##                                                                           ##
##  List of source files:                                                    ##
##                       - read.admb.R                                       ##
##                       - summary.R                                         ##
##                       - plot.ft.R                                         ##
##                       - plot.resid.R                                      ##
##                                                                           ##
##                                      ##

## Global Variables:
.VIEWTRCK	<- "iSCAMViewTracker.txt"	# filename for viewer file.

## Dependencies:
  require(PBSmodelling)
  require(Hmisc)
  require(ggplot2)


## Source files:
  source("read.admb.R")
  source("summary.R")
  source("plot.ft.R")
  source("plot.resid.R")
  source("iscamView.R")
  source("plotCatch.R")


## Instructions:
.inst <-
function()
{
	cat("=====================\n")
	cat("= type: guiView()   =\n")
	cat("= to launch viewer  =\n")
	cat("=====================\n")
}
.inst()