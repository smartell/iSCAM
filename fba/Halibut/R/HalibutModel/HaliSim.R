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
.HARVESTPOLICY_FILE <- c("iphcHP.txt")
.TONNES2LBS         <- 2204.62262
.REP_FILE           <- paste(.MODEL_DIRECTORY, list.files(.MODEL_DIRECTORY,"^SIMREP"), sep="")
.NEW_RUN            <- TRUE
.iCLR               <- c("red","black","blue")
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
		
		hptable <- hpFile[, -1]
		hptable[6:8, 4:5] <- spnBSAI*hptable[6:8, 4:5]
		hptable[4:5, 4:5] <- spnGULF*hptable[4:5, 4:5]
		
		write("# Controls for Halibut simulation model from R-GUI ", fn)
		write(spnNyrs, fn, append=TRUE);
		write(spnDdgr, fn, append=TRUE);
		
		write("# Area based harvest policy from R-GUI", fn, append=TRUE)
		write.table(t(hptable), fn, append=TRUE, row.names=FALSE, col.names=FALSE)
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
	
	# Graphical parameters
	par(bg="white", cex.lab=1.5, cex.axis=1.2, mar=c(5.1,  4.6,  4.1,  2.1))
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
		.plotEBio( RF,xarg="yrsim", yarg="EBio",  yscale=1e6, ylim=c(0, 1200), lwd=5, lty=1 )
		# Wobbleseq references
		lines(RF[[1]]$yrs,WSQ$EBio,lwd=15, col=colr("grey", 0.5))
		grid(col="grey", lwd=2 )
	}
	if ( plotType=="sbio" )	
	{
		.plotEBio( RF,xarg="yrsim", yarg="SBio",  yscale=1e6, ylim=c(0, 800), lwd=5, lty=1, 
		ylab="Spawning biomass (million lb)" )
		grid(col="grey", lwd=2 )
		abline(h=0.3*RF[[1]]$bo/1e6, lwd=14, col=colr("grey", 0.5))
	}
	if ( plotType=="enum" )
	{
		.plotENum( repObj, annotate=TRUE, ylim=c(max()) )
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
		#.plotCatch( RF, ylim=c(0, 200) )
		.plotLandedValue(RF, yarg="ct", ylab="Catch (million lb)", lwd=5, lty=1, ylim=c(0, 200) )
		grid(col="grey", lwd=2)
	}
	if ( plotType=="waste" )
	{
		.plotWaste( RF, ylim=c(0, 5), lwd=5, lty=1 )
		grid(col="grey", lwd=2 )
	}
	if ( plotType=="valcatch" )
	{
		.plotLandedValue( RF, lwd=5, lty=1, ylim=c(0, 1200) )
		grid(col="grey", lwd=2 )
	}
	if ( plotType=="valwaste" )
	{
		.plotLandedValue( RF,yarg="val_ht", ylab="Value of wastage (million $)", lwd=5, lty=1, ylim=c(0, 30) )
		grid(col="grey", lwd=2 )
	}
	if ( plotType=="valdisc" )
	{
		.plotLandedValue( RF,yarg="val_dt", ylab="Value of commercial discards (million $)", lwd=5, lty=1, ylim=c(0, 150) )
		grid(col="grey", lwd=2 )
	}
	if ( plotType=="discfrac" )
	{
		.plotLandedValue( RF,yarg="dt", ylab="Fraction of catch discarded (numbers)", lwd=5, lty=1, ylim=c(0, 1.2) )
		grid(col="grey", lwd=2 )
	}
	if ( plotType=="yieldloss" )
	{
		.plotYieldLoss( RF , ylim=c(0, 15))
	}
	if ( plotType=="ylr" )
	{
		.plotYieldLossRatio( RF, xlim=c(2012, 2026), lty=1, lwd=5, ylim=c(0.7, 1.0))
		grid(col="grey", lwd=2 )
	}
	if ( plotType=="catchage" )
	{
		.plotCatchAge( repObj )
	}
	if ( plotType=="fishrate" )
	{
		par(mfcol=c(1, 2))
		.plotFishRate( repObj, sex="Female" )
		.plotFishRate( repObj, sex="Male" )
		dev.copy2pdf(file="../../FIGURES/fig:FishRate.pdf", width=8, height=4)
	}
}
.plotFishRate	<- function( repObj, annotate=TRUE, irow=1:5, sex="Female",  ...)
{
	with(repObj, {
		xx=yrsim
		if(sex=="Female") ir=1:5
		if(sex=="Male")   ir=6:10
		if(is.matrix(ft))
		{
			yy=t(as.matrix(ft[ir, ]))
		}
		else
		{
			yy=as.matrix(ft)
		}
		#n=nrow(t(as.matrix(yy)))
		icol=apply(yy,2,function(x){sum(cumsum(x))!=0.0})
		ng=length(icol[icol==T])
		
		#yy = cbind( rowMeans(M_tot), yy[,icol] )
		yy = yy[, icol]
		
		csyy = t(apply(yy,1, cumsum))	#cumulative sum
			
		yrange=c(0, max(csyy, na.rm=TRUE))
		lw = c(1, rep(1,ng))
		lt = 1
		iclr = colr(1:(ng+1),0.5)
		
		matplot(xx, csyy, type="n", axes=FALSE, log="y",  
			xlab="Year", ylab="Fishing mortality rate", main=paste(sex))
			
		
		#lines(xx, yy[,1], lwd=2, lty=1, col=1)
		matlines(xx, csyy, log="y", col=iclr, lwd=lw, lty=lt)
		axis( side=1 )
		axis( side=2, las=1 )
		box()
		grid()
		
		ry=cbind(1.e-30,csyy)
		for(i in 1:dim(csyy)[2])
		{
			x2=c(xx, rev(xx))
			y2=c(ry[,i+1], rev(ry[,i]))
			#browser()
			polygon(x2, y2, border=NA,col=colr(i,0.2), log="y")
		}
		
		
		
		if ( annotate )
		{
			txt=c("Commercial", "U32 Bycatch", "O32 Bycatch", "Recreational", "Personal")
			#txt = c("Natural mortality", paste("Gear",1:ng))
			#txt = c(paste("Gear",1:ng))
			#mfg <- par( "mfg" )
			#if ( mfg[1]==1 && mfg[2]==1 )
			legend( "topleft",legend=txt,col=iclr, 
				bty='n',lty=lt,lwd=5,pch=-1,ncol=2, cex=0.7)
		}
	})
	
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


.plotYieldLoss <- function( repObj, ...)
{
	# plot lost yield
	n <- length(repObj)
	for(ii in 1:n)
	{
		with(repObj[[ii]], {
			xx = yrsim
			yy = yieldLoss[1, ]
		
			yy[yy<=0] = NA
		
			if(ii==1)
			plot(xx, (yy), type="l", xlab="Year", ylab="Yield Loss (million lb)", 
			xlim=c(2010, max(xx)), lwd=5, ...)
			if(ii>1)
			lines(xx, (yy), col=ii, lwd=5)
		})
		
	}
}





.plotCatch	<- function( repObj, ...)
{
	# Matrix plot of Catch 
	require(hacks)
	n <- length(repObj)

	for(ii in 1:n)
	{
		with(repObj[[ii]], {
			if(ii==1)
			{
				plot(yrsim, ct[1, ], type="l", xlab="Year", ylab="Catch (million lb)", col=.iCLR[ii], lwd=3, ...)
				abline(h=60, lwd=20, col=colr("grey", 0.2))
				text(1996, 60, "60 Million Pounds", adj=0, col=colr(1, 0.6))
			}
			else{
				lines(yrsim, ct[1, ], col=.iCLR[ii], lwd=3)
			}
		
			#gears <- c("Setline", "U32", "O32", "Recreational", "Personal")
			#legend("top", gears, lty=1:5, col=1:5, lwd=2, ncol=2, bty="n")
			if(ii==2)
			text(1996, 100, paste(round(mean(ct[1, 25:30]), 1)), cex=3.5, pos=4)
		})
				
	}
	
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

.plotEBio	<- function( repObj, xarg, yarg, xscale=1, yscale=1, xlab=NULL, ylab=NULL, ... )
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
		
		if(is.null(xlab)) xlab="Year"
		if(is.null(ylab)) ylab="Exploitable biomass (million pounds net)"
		
		plot (xx, yy, type="l", xlab=xlab, ylab=ylab, ...)
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
		if(is.null(xlab)) xlab="Year"
		if(is.null(ylab)) ylab="Exploitable biomass (million pounds net)"
		matplot (xx, yy, type="l", xlab=xlab, ylab=ylab, col=.iCLR, ...)
		
		mu <- round(mean(yy[25:30, ]), 0)
		text(1996, mu, paste(mu) , cex=3.5, pos=4)
	}
	
}

.plotWaste	<- function( repObj,yarg=NULL , ...)
{
	# Matrix plot of Catch 
	require(hacks)
	n <- length(repObj)
	xx <- yy <- NULL
	yarg="ht"
	for(ii in 1:n)
	{
		ix <- which( names(repObj[[ii]])=="yrsim" )
		iy <- which( names(repObj[[ii]])==yarg )
		xx <- cbind(xx, repObj[[ii]][[ix]])
		yy <- cbind(yy, repObj[[ii]][[iy]][1, ])
	}
	xlab="Year"
	ylab="Wastage (million lbs)"
	matplot (xx, yy, type="l", xlab=xlab, ylab=ylab, col=.iCLR, ...)
	mu <- round(mean(yy[25:30, ]), 3)
	text(1996, 2.5, paste(mu), cex=3.5, pos=4)

}

.plotYieldLossRatio <- function( repObj, ...)
{
	# plot lost yield
	n <- length(repObj)
	xx <- yy <- dd <- NULL
	for(ii in 1:n)
	{
		ix <- which( names(repObj[[ii]])=="yrsim")
		iy <- which( names(repObj[[ii]])=="yieldLoss")
		id <- which( names(repObj[[ii]])=="ct")
		xx <- cbind(xx, repObj[[ii]][[ix]])
		yy <- cbind(yy, repObj[[ii]][[iy]][1, ])
		dd <- cbind(dd, colSums(repObj[[ii]][[id]][2:3, ]))
	}
	iyr <- which(xx[,1]==2012)
	ii  <- iyr:dim(xx)[1]
	matplot(xx[ii, ], yy[ii, ]/dd[ii, ], type="l", xlab="Year", ylab="Yield Loss Ratio",col=.iCLR,  ...)
	mu <- round(mean(yy[25:30, ]/dd[25:30, ]), 2)
	text(2018, mu, paste(mu), cex=3.5, pos=4)
		##with(repObj[[ii]], {
		##	xx = yrsim
		##	yy = as.vector(yieldLoss[1, ]/colSums(ct[2:3, ]))
		##	yy[yy<=0] = NA
		##	
		##	if(ii == 1)
		##	plot(xx, yy, type="l", xlab="Year", ylab="Yield Loss Ratio", xlim=c(2010, max(xx)))
		##	if(ii > 1)
		##	lines(xx,yy, lty=ii, col=ii)
		##})
		##
	
}


.plotLandedValue <- function(repObj, xarg=NULL,yarg=NULL, xlab=NULL, ylab=NULL, ...)
{
	# Matrix plot of Catch 
	require(hacks)
	n <- length(repObj)
	xx <- yy <- NULL
	if(is.null(xarg)) xarg="yrsim"
	if(is.null(yarg)) yarg="val_ct"
	for(ii in 1:n)
	{
		ix <- which( names(repObj[[ii]])==xarg )
		iy <- which( names(repObj[[ii]])==yarg )
		xx <- cbind(xx, repObj[[ii]][[ix]])
		yy <- cbind(yy, repObj[[ii]][[iy]][1, ])
	}
	if(is.null(xlab)) xlab="Year"
	if(is.null(ylab)) ylab="Landed value ($ million)"
	matplot (xx, yy, type="l", xlab=xlab, ylab=ylab, col=.iCLR, ...)
	mu <- round(mean(yy[25:30, ]), 1)
	text(1996, mu, paste(mu), cex=3.5, pos=4)
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

plotSelRet<-function(ii=1, A)
{
	with(A, {
		sj=exp(log_sel[ii, -1:-2])
		rj=exp(log_ret[ii, -1:-2])
		age=1:30
		plot(age, sj*rj, xlab="Age", ylab="P(retention) & P(discard)", type="l", col=1, lwd=5)
		lines(age, sj*(rj*(1-rj)), lwd=5, col=2)
		lines(age, sj, lty=2)
		legend("topleft",c("Retention", "Discard"), lwd=5, col=c(1, 2), bty="n")
	})
}

plotWa <- function()
{
	par(mfcol=c(1, 2))
	require(hacks)
	wa.f = cbind(RF[[1]]$wt_obs[31, ],RF[[2]]$wt_obs[31, ],RF[[3]]$wt_obs[31, ],RF[[2]]$wt_obs[16, ])
	wa.m = cbind(RF[[1]]$wt_obs[63, ],RF[[2]]$wt_obs[63, ],RF[[3]]$wt_obs[63, ],RF[[2]]$wt_obs[48, ])

	matplot(1:20, wa.f[1:20, ], xlab="Age", ylab="Weight (lb)", type="l",lty=1,lwd=c(2, 2, 2, 5), col=c(.iCLR, colr("grey", 0.5)), main="Female")
	
	matplot(1:20, wa.m[1:20, ], xlab="Age", ylab="Weight (lb)", type="l",lty=1,lwd=c(2, 2, 2, 5), col=c(.iCLR, colr("grey", 0.5)), main="Male")
	
}




