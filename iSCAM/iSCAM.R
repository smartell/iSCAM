#-------------------------------------------------------------------------------#
#   iSCAM Viewer: A gui based viewer for iscam inputs and outputs               #
#                                                                               #
#                                                                               #
#   Authors: Steven Martell (with lots of borrowed code from A.R.Kronlund)      #
#            A.R. Kronlund (Pacific Biological Station, Nanaimo, B.C.)          #
#   Date: Nov. 22,  2010                                                        #
#                                                                               #
#                                                                               #
#                                                                               #
#                                                                               #
#                                                                               #
# NOTES:                                                                        #
# 1. requires PBSmodelling                                                      #
#                                                                               #
#                                                                               #
#-------------------------------------------------------------------------------#




.iscamViewSetup <- function(win)
{
	#Required libraries
	require(PBSmodelling)
	
	
	#Close any open graphics devices
	graphics.off()
	closeWin()
	
	#Create new window based on iscamWin.txt
	createWin("iscamWin.txt")
}

.iscamViewSetup("test")