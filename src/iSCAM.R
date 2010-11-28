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

require(Riscam)	#custom library built specifically for iscam.


# Graphics defaults.
.VIEWCEX    <- 1            # Generic default cex for axis labels etc.
.VIEWPANCEX <- 1            # Default cex for panLab.

.VIEWMAR  <- c(2, 2, 1, 1)  # Multi-panel plots: plot margin sizes c(b,l,t,r).
.VIEWOMA  <- c(2, 2, 1, 1)  # Multi-panel plots: outer margin sizes c(b,l,t,r).
.VIEWLAS  <- 2









repObj	<- read.rep("iscam.rep")



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

.mpdView	<- function()
{
	print(".mpdView")
	
	# Get the guiPerf parameters so that plot controls available.
	guiInfo <- getWinVal(scope="L")
	
	
	#Conditional statements for radio button selections of plots
	if ( plotType=="catch" )
	{
		.plotCatch( repObj, annotate=FALSE )
	}
	
	if ( plotType=="biomass" )
    {
    	.plotBiomass( repObj, annotate=TRUE )
    }
	
	if ( plotType=="depletion" )
	{
		.plotDepletion( repObj, annotate=TRUE )
	}
}


.plotCatch	<- function( repObj, annotate=FALSE )
{
	#barplot of the observed catch
	with(repObj, {
		barplot(obs_ct, names.arg=yr,axes=FALSE, 
			xlab="Year", ylab="Catch (t)")
		axis( side=2, las=.VIEWLAS )
		
	})
}

.plotDepletion	<- function( repObj, annotate=FALSE )
{
	#plot the spawning biomass depletion level & reference points
	with(repObj, {
		xx=yrs
		yy=sbt/bo
		yrange=c(0,1.1*max(yy, na.rm=TRUE))
		
		plot(xx, yy, type="n", axes=FALSE,
			xlab="Year", ylab="Depletion", 
			ylim=yrange)
		lines(xx, yy)
		rlvl=c(1.0, 0.8, 0.4)
		abline(h=rlvl*bmsy/bo,lty=2,lwd=rlvl)
		
		
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		grid()
		
		if ( annotate )
		{
			mfg <- par( "mfg" )
			if ( mfg[1]==1 && mfg[2]==1 )
			legend( "top",legend=c( "Spawning biomass","MSY depletion level",
				"Upper stock reference","Limit reference point"),
				bty='n',lty=c(1,2,2,2),lwd=c(1,rlvl),pch=-1,ncol=2 )
		}
	})
}

.plotBiomass	<- function( repObj, annotate=FALSE )
{
	#plot total biomass & spawning biomass 
	with(repObj, {
		xx=yrs
		yy=cbind(bt, sbt)
		
		yrange=c(0, 1.2*max(yy, na.rm=TRUE))
		
		matplot(xx, yy, type="n",axes=FALSE,
				xlab="Year", ylab="Biomass (t)", 
				ylim=yrange)
		
		matlines(xx,yy,
			type="l", col="black",
			ylim=c(0,max(yy,na.rm=T)))
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		
		if ( annotate )
		{
			mfg <- par( "mfg" )
			if ( mfg[1]==1 && mfg[2]==1 )
			legend( "top",legend=c( "Pre-fishery biomass","Spawning biomass"),
				bty='n',lty=c(1,2),lwd=c(1,1),pch=c(-1,-1),ncol=1 )
		}
	})
	
}

