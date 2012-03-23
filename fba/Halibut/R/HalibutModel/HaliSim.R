# -------------------------------------------------------------------------- ##
# HaliSim.R                                                                  ##
# Author: Steven Martell                                                     ##
# Date: March 23, 2012                                                       ##
#                                                                            ##
# A GUI INTERFACE FOR THE HALIBUT SIMULATION MODEL FOR BYCATCH RESEARCH      ##
#                                                                            ##
# -------------------------------------------------------------------------- ##
.MODEL_DIRECTORY <- "../../DATA/"
.SIMULATION_FILE <- "Halibut_2sex_develop.sim"
.TONNES2LBS      <- 2204.62262



guiView	<- function()
{
	.guiSetUp()
}

.guiSetUp	<- function()
{
	#Required libraries
	require(PBSmodelling)
	
	
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
	})
}













