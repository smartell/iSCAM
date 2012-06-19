##                                      ##
##  R-script for the iSCAM project                                           ##
##  Author: Steven Martell,  University of British Columbia                  ##
##  LAST MODIFIED: June 18, 2012                                             ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
##  SOURCE FILES           STATUS              TODO                          ##                
##  - read.admb.R          - working                                         ##
##  - summary.R            - working                                         ##
##  - plot.ft.R                                                              ##
##  - plot.resid.R                                                           ##
##  - iscamView.R          - working           - fix x-axis labels           ##
##  - plotCatch.R                                                            ##
##  - plotAgeComs.R                                                          ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
##                                                                           ##
##                                      ##

## Global Variables:
.VIEWTRCK	<- "iSCAMViewTracker.txt"	# filename for viewer file.

## Dependencies:
  require(PBSmodelling)
  require(Hmisc)
  require(ggplot2)
  require(reshape)
  require(grid)
  


## Source files:
  source("read.admb.R")
  source("theme_iscam.R")
  source("summary.R")
  source("plot_bt.R")
  source("plot.ft.R")
  source("plot.resid.R")
  source("iscamView.R")
  source("plotCatch.R")
  source("plotAgeComps.R")
  source("plotMeanWt.R")
  source("plotIt.R")


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