# -------------------------------------------------------------------------- ##
# HaliSim.R                                                                  ##
# Author: Steven Martell                                                     ##
# Date: March 23, 2012                                                       ##
#                                                                            ##
# A GUI INTERFACE FOR THE HALIBUT SIMULATION MODEL FOR BYCATCH RESEARCH      ##
#                                                                            ##
# -------------------------------------------------------------------------- ##
.WDIR               <- '/Users/stevenmartell/Documents/iSCAM-project/fba/Halibut/R/HalibutModel'
.MODEL_DIRECTORY    <- "/Users/stevenmartell/Documents/iSCAM-project/fba/Halibut/DATA/"
.SIMULATION_FILE    <- "Halibut_2sex_develop.sim"
.HARVESTPOLICY_FILE <- "iphcHP.txt"
.TONNES2LBS      <- 2204.62262



guiView	<- function()
{
	.guiSetUp()
}

.guiSetUp	<- function()
{
	setwd(.WDIR)
	
	#Required libraries
	require(PBSmodelling)
	
	#Read Harvest Policy controls
	hpFile <- read.table(.HARVESTPOLICY_FILE, header=TRUE)
	
	#Close any open graphics devices
	graphics.off()
	closeWin()
	
	#Creat window based on HaliSimWin.txt
	createWin("HaliSimWin.txt")
}


.runSimulation	<- function()
{
	#Read simulation controls from GUI and write Simulation file
	# Get the guiPerf parameters so that plot controls available.
	guiInfo <- getWinVal(scope="L")
	
	.writeSimulationFile(guiInfo)
	print(guiInfo)
	
	#Required libraries
	wdir<- "/Users/stevenmartell/Documents/iSCAM-project/fba/Halibut/DATA/"
	arg <- "make"
	setwd(wdir)
	system(arg)
}

.writeSimulationFile	<- function(gI)
{
	fn <-paste( .MODEL_DIRECTORY ,  .SIMULATION_FILE, sep="")
	with(gI, {
		
		write("# Controls for Halibut simulation model from R-GUI ", fn)
		write(spnNyrs, fn,  append=TRUE);
		
		write("# Area based harvest policy from R-GUI", fn, append=TRUE)
		write.table(t(hpFile[, -1]), fn, append=TRUE, row.names=FALSE, col.names=FALSE)
	})
}













