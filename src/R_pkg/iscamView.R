## The following is the R-script for generating the iscamViewer.

# guiView
# Purpose:    wrapper for .isamView    
# Parameters: non 
# Returns:    NULL 
# Source:     
guiView <-
function()
{
	.iscamView("iSCAM")
}




# .iscamView
# Purpose:    Open an x11 gui    
# Parameters: winName is a string containing the name of the window to set up. 
# Returns:    NULL
# Source:     Steve Martell
.iscamView <-
function( winName )
{
	# Requires libraries
	require(PBSmodelling)
	
	# Close any open graphics devices or X11 devices
	graphics.off()
	#closeWin()
	
	# Check for global .VIEWTRCK and read header files
	trckExists <- file.exists( .VIEWTRCK )
	if (trckExists)
	{
		tmp <- read.table( file = .VIEWTRCK,  as.is=TRUE,  header=TRUE,  sep="," )
		cat( "MSG (.iscamView): Viewer tracking file ",.VIEWTRCK, " found.\n" )
		ifiles <- tmp
		gomenu <- TRUE
	}
	else
	{
		cat( "Error MSG (.iscamView): missing", .VIEWTRCK, "file.\n")
		gomenu <- FALSE
	}
	
	if(gomenu)
	createWin("iscamWin.txt")
	
}

# .viewPlot
# Purpose:    Determine which plot device was selected from the gui and call
#    		  the appropriate function to plot.
# Parameters:  
# Returns:     
# Source:     
.viewPlot <-
function()
{
	print(".mpdView")
	
	# Get the guiPerf parameters so that plot controls available.
	guiInfo <- getWinVal(scope="L")
	
	# List of gui changes
	guiChanges <- list()
	
	# Determine which files have been selected
	# and read files and store into model object (M)
	hdr	<- ifiles[ ifiles$Select, ]
	fn	<- hdr$Control.File
	M	<- lapply(fn, read.admb)
	names(M) <- hdr$Model
	
	# use plotType to determine which function to call.
	switch(plotType, 
		catch={
			print("catch")
			.plotCatch(M)
		}, 
		survey={
			print("survey")
			.plotIt(M)
		}, 
		agecomp={
			print("agecomp")
			.plotAgeComps(M)
		},
		agehist={
			print("agehist")
			.plotAgeHist(M)
		},  
		meanwt={
			print("meanwt")
			.plotMeanWt(M)
		}, 
		vbio={
			print("vbio")
			.plot_bt(M, "bt")
		}, 
		sbio={
			print("sbio")
			.plot_bt(M, "sbt")
		}, 
		urate={
			print("urate")
			.plot_bt(M, "ut")
		}
		)
	
	# Return a global model object that are in play
	iMod <<- M
}





 