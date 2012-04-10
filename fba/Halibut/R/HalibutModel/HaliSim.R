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
.FILE_PREFIX        <- "SIMREP2"
.HARVESTPOLICY_FILE <- "iphcHP.txt"
.TONNES2LBS         <- 2204.62262
.REP_FILE           <- paste(.MODEL_DIRECTORY, list.files(.MODEL_DIRECTORY,"^SIMREP"), sep="")
.NEW_RUN            <- TRUE

# DEPENDENCIES
source("../read.admb.R")
#WSQ <- read.rep("/Users/stevenmartell/Documents/iSCAM-project/fba/Halibut/MISC/wobblesq.rep")
#RF   <- lapply(.REP_FILE, read.rep)


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
	
	.NEW_RUN <<- TRUE
	cat("New Run ", .NEW_RUN)
	
	guiInfo$plotType="ebio"
	setWinVal(guiInfo)
	.mpdView()
}

.writeSimulationFile	<- function(gI)
{
	fn <-paste( .MODEL_DIRECTORY ,  .SIMULATION_FILE, sep="")
	with(gI, {
		
		write("# Controls for Halibut simulation model from R-GUI ", fn)
		write(spnNyrs, fn, append=TRUE);
		write(spnDdgr, fn, append=TRUE);
		
		write("# Area based harvest policy from R-GUI", fn, append=TRUE)
		write.table(t(hpFile[, -1]), fn, append=TRUE, row.names=FALSE, col.names=FALSE)
	})
}

.readReportFile <- function()
{
	# Read the report files and make a global object
	#if(!exists("RF"))
	RF   <- lapply(.REP_FILE, read.rep)
	if(!exists("WSQ"))
	WSQ  <<- read.rep("/Users/stevenmartell/Documents/iSCAM-project/fba/Halibut/MISC/wobblesq.rep")
	return(RF)
}

.mpdView	<- function()
{
	print(".mpdView")
	
	# Get the guiPerf parameters so that plot controls available.
	guiInfo <- getWinVal(scope="L")
	
	# List of gui changes
	guiChanges <- list()
	
	# Read input & output files
	if(.NEW_RUN)
	{
		RF       <<- .readReportFile()
		.NEW_RUN <<- FALSE
	
		repfn    <- paste(.MODEL_DIRECTORY, .FILE_PREFIX,".rep",  sep="")
		repObj   <<- read.rep(repfn)
	}
	
	if ( plotType=="ebio" )
	{
		.plotEBio( RF,xarg="yrsim", yarg="EBio",  yscale=1e6 )
	}
	if ( plotType=="enum" )
	{
		.plotENum( repObj, annotate=TRUE )
	}
	if ( plotType=="n8plus" )
	{
		.plotN8plus( RF, annotate=TRUE )
	}
	if ( plotType=="agesel" )
	{
		.plotAgeSel( repObj )
	}
	if ( plotType=="catch" )
	{
		.plotCatch( repObj )
	}
	if ( plotType=="waste" )
	{
		.plotWaste( repObj )
	}
	if ( plotType=="yieldloss" )
	{
		.plotYieldLoss( RF )
	}
	if ( plotType=="ylr" )
	{
		.plotYieldLossRatio( RF )
	}
	if ( plotType=="catchage" )
	{
		.plotCatchAge( repObj )
	}
}
.plotCatchAge	<- function( repObj, ...)
{
	# Bubble plots of Catch-at-age (coast-wide)
	require(PBSmodelling)
	with(repObj, {
		sex  <- unique(Chat[, 1])
		gear <- unique(Chat[, 2])
		for(h in sex)
			for(k in gear)
			{
				z <- subset(Chat, Chat[, 1]==h)
				z <- subset(z   , z[, 2]==k)
				yrs <- z[, 3]
				plotBubbles(t(z[, -1:-3]), yrs, age , xlab="Year", ylab="Age", 
					cpro=TRUE, prettyaxis=TRUE, hide0=TRUE, frange=0.01, fill=TRUE)
					
				
				y2=max(yrs)-2007+1
				polygon(x=c(2006.5, max(yrs)+.5, max(yrs)+.5, 2006.5), 
					y=c(1, 1, y2, 1), col=colr(4, 0.2), border=NA)
			}
		
	})
}

.plotWaste	<- function( repObj, ...)
{
	# Matrix plot of Catch 
	require(hacks)
	with(repObj, {
		matplot(yrsim, t(ht), type="l", xlab="Year", ylab="Wastage (million lbs)", lwd=2)
		abline(h=33.135, lwd=20, col=colr("grey", 0.2))
		text(1996, 33.135, "2012 Coastwide CEY", adj=0, col=colr(1, 0.6))
		
		gears <- c("Setline", "U32", "O32", "Recreational", "Personal")
		legend("top", gears, lty=1:5, col=1:5, lwd=2, ncol=2, bty="n")
	})
}

.plotYieldLoss <- function( repObj, ...)
{
	# plot lost yield
	n <- length(repObj)
	for(ii in 1:n)
	{
		with(repObj[[ii]], {
			xx = yrsim
			yy = yieldLoss
		
			yy[yy<=0] = NA
		
			if(ii==1)
			matplot(xx, t(yy), type="l", xlab="Year", ylab="Yield Loss (million lb)", xlim=c(2010, max(xx)))
			if(ii>1)
			matlines(xx, t(yy))
		})
		
	}
}

.plotYieldLossRatio <- function( repObj, ...)
{
	# plot lost yield
	n <- length(repObj)
	for(ii in 1:n)
	{
		with(repObj[[ii]], {
			xx = yrsim
			yy = as.vector(yieldLoss[1, ]/colSums(ct[2:3, ]))
			yy[yy<=0] = NA
			
			if(ii == 1)
			plot(xx, yy, type="l", xlab="Year", ylab="Yield Loss Ratio", ylim=c(0, 2), xlim=c(2010, max(xx)))
			if(ii > 1)
			lines(xx,yy, lty=ii, col=ii)
		})
		
	}
}



.plotCatch	<- function( repObj, ...)
{
	# Matrix plot of Catch 
	require(hacks)
	with(repObj, {
		matplot(yrsim, t(ct), type="l", xlab="Year", ylab="Catch (million lbs)", lwd=2)
		abline(h=33.135, lwd=20, col=colr("grey", 0.2))
		text(1996, 33.135, "2012 Coastwide CEY", adj=0, col=colr(1, 0.6))
		
		gears <- c("Setline", "U32", "O32", "Recreational", "Personal")
		legend("top", gears, lty=1:5, col=1:5, lwd=2, ncol=2, bty="n")
	})
}

.plotAgeSel	<- function( repObj, ...)
{
	# plot first row of age-selectivity for each sex/gear
	opar<- par(no.readonly=TRUE)
	par(mfcol=c(5, 2), mar=c(2, 2.5, 1, 1), oma=c(2, 2, 1, 1), las=1)
	with(repObj, {
		hsex <- unique(log_sel[, 1])
		str_sex=c("Female", "Male")
		for(h in hsex)
		{
			D     <- subset(log_sel, log_sel[, 1]==h)
			kgear <- unique(D[, 2])
			
			for(k in kgear)
			{
				E <- subset(D, D[, 2]==k)
				plot(age, exp(E[1, -1:-2]), type="l",xlab="", ylab="", bty="l", ylim=c(0, 1.2))
				if(k==1) title(main=str_sex[h])
				gletter(k)
			}
		}
		mtext(c("Age", "Selectivity"), c(1, 2), outer=TRUE, las=0)
		par(opar)
	})
	dev.copy2pdf(file="../../FIGURES/fig:AgeSel.pdf")
	
}
# S3 Classes for iscam plots
plot.iscam	<- 
function( object, 
	xarg="yrsim", yarg="EBio", 
	xscale=1, yscale=1,
	xlab=NULL, ylab=NULL,
	xlim=NULL, ylim=NULL,   
	type="l", ... )
{
	n  <- length(object)
	xx <- NULL
	yy <- NULL
	for(ii in 1:n)
	{
		ix <- which(names(object[[ii]]) == xarg)
		xx <- cbind(xx, object[[ii]][[ix]]/xscale)
		iy <- which(names(object[[ii]]) == yarg)
		yy <- cbind(yy, object[[ii]][[iy]]/yscale)
	}
	if(is.null(xlab)) xlab=xarg
	if(is.null(ylab)) ylab=yarg
	xlim <- if(is.null(xlim)) range(xx[is.finite(xx)])
	if(is.null(ylim)) ylim <- range(yy[is.finite(yy)])
	
	matplot(xx, yy, type="l", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
}

.plotEBio	<- function( repObj, ..., xarg, yarg, xscale=1, yscale=1 )
{
	#Determine length of repObj &:  
	require(hacks)
	nL  <- length(repObj)

	if ( nL >= 5)  # A single report object
	{
		ix  <- which( names(repObj) == xarg )
		iy  <- which( names(repObj) == yarg )
		xx  <- repObj[[ix]]/xscale;
		yy  <- repObj[[iy]]/yscale;       
		
		names(xx) <- "Year"
		names(yy) <- "Expoitable biomass (million pounds net)"
		
		plot (xx, yy, type="l", xlab=names(xx), ylab=names(yy))
	}
	
	if ( nL < 5) # Multiple report objects use matrix plot
	{
		xx  <- yy <- NULL
		for (ii in 1:nL)
		{
			ix  <- which( names(repObj[[ii]]) == xarg )
			iy  <- which( names(repObj[[ii]]) == yarg )
			print(ix)
			xx  <- cbind(xx, repObj[[ii]][[ix]]/xscale);
			yy  <- cbind(yy, repObj[[ii]][[iy]]/yscale);
		}
		xlbl <- "Year"
		ylbl <- "Expoitable biomass (million pounds net)"
		matplot (xx, yy, type="l", xlab=xlbl, ylab=ylbl)
		
		
	}
	# Wobbleseq references
	lines(repObj[[1]]$yrs,WSQ$EBio,lwd=15, col=colr("grey", 0.5))
	
	

}


.plotENum	<- function( repObj, annotate=FALSE )
{
	#plot total biomass & spawning biomass 
	require(hacks)
	with(repObj, {
		xx=yrsim; names(xx)="Year"
		yy=ENum/1e6;  names(yy)="Expoitable Numbers (million)"
		
		plot(xx, yy, type="line", xlab=names(xx), ylab=names(yy), ylim=c(0, 1.2*max(yy)) )
		#abline(h=260, lwd=20, col=colr("grey", 0.2))
		#text(1996, 260, "IPHC Assessment of EBio (260 mlb)", adj=0, col=colr(1, 0.6))
		
	})	
}

.plotN8plus	<- function( repObj, annotate=FALSE )
{
	#plot total biomass & spawning biomass 
	require(hacks)
	n=length(repObj)
	for(ii in 1:n)
	{
		with(repObj[[ii]], {
			xx=yrsim; names(xx)="Year"
			yy=N8plus/1e6;  names(yy)="Number of 8+ (million)"
		
			if(ii==1)
			{
				plot(xx, yy, type="line", xlab=names(xx), ylab=names(yy), ylim=c(0, 1.4*max(yy)) )
		
				xx=yrs
				yy=(rowSums(WSQ$N.F[, 8:30]) + rowSums(WSQ$N.M[, 8:30]))/1e6
				lines(xx, yy, lwd=20, col=colr("grey", 0.5))
				#abline(h=260, lwd=20, col=colr("grey", 0.2))
				text(2006, max(yy), "IPHC N 8+ Total", adj=0, col=colr(1, 0.6))
			}
			if(ii>1)
			lines(xx, yy)
		
		})	
		
	}
}






