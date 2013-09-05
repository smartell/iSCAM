##                                       ##
#-------------------------------------------------------------------------------#
#   iSCAM Viewer: A gui based viewer for iscam inputs and outputs               #
#                                                                               #
#                                                                               #
#   Authors: Steven Martell (with lots of borrowed code from A.R.Kronlund)      #
#            A.R. Kronlund (Pacific Biological Station, Nanaimo, B.C.)          #
#   Date: Nov. 22,  2010                                                        #
#   Date: Aug. 28,  2012                                                        #
#                                                                               #
#                                                                               #
#                                                                               #
#   DIRECTORY TREE                                                              #
#    .                                                                          #
#    |____.RData                                                                #
#    |____.Rhistory                                                             #
#    |____iSCAM.R                                                               #
#    |____iSCAMequil_soln.R                                                     #
#    |____iSCAMViewTracker.txt                                                  #
#    |____iSCAMWin.txt                                                          #
#    |____logo                                                                  #
#    | |____iscamLogo.eps                                                       #
#    | |____iscamLogo.gif                                                       #
#    | |____iscamLogo.png                                                       #
#    | |____iscamLogoSmall.png                                                  #
#    | |____logo.r                                                              #
#    |____R                                                                     #
#    | |____.RData                                                              #
#    | |____.Rhistory                                                           #
#    | |____plotAgeCompResiduals.R                                              #
#    | |____plotAgeComps.R                                                      #
#    | |____plotBiomass.R                                                       #
#    | |____plotCatch.R                                                         #
#    | |____plotCatchResiduals.R                                                #
#    | |____plotDepletion.R                                                     #
#    | |____plotIndex.R                                                         #
#    | |____plotMarginalPosteriors.R                                            #
#    | |____plotMeanWt.R                                                        #
#    | |____plotMortality.R                                                     #
#    | |____plotNaturalMortality.R                                              #
#    | |____plotRecruitment.R                                                   #
#    | |____plotRecruitmentResiduals.R                                          #
#    | |____plotReferencePoints.R                                               #
#    | |____plotRiskTable.R                                                     #
#    | |____plotSelectivity.R                                                   #
#    | |____plotSSBretrospective.R                                              #
#    | |____plotStockRecruitment.R                                              #
#    | |____plotStockStatus.R                                                   #
#    | |____plotSurveyFit.R                                                     #
#    | |____plotSurveyResiduals.R                                               #
#    | |____read.admb.R                                                         #
#                                                                               #
#                                                                               #
#                                                                               #
#                                                                               #
#                                                                               #
#                                                                               #
#-------------------------------------------------------------------------------#
setwd("/Users/stevenmartell1/Documents/iSCAM-project/src/r-code/")
require(Hmisc)
.RFILES     <- list.files("./R/",pattern="\\.[Rr]$")
for(nm in .RFILES) source(file.path("./R", nm), echo=FALSE)


# Graphics defaults.
.VIEWCEX    <- 1            # Generic default cex for axis labels etc.
.VIEWPANCEX <- 1            # Default cex for panLab.

.VIEWMAR  <- c(3, 3, 1, 1)  # Multi-panel plots: plot margin sizes c(b,l,t,r).
.VIEWOMA  <- c(2, 2, 1, 1)  # Multi-panel plots: outer margin sizes c(b,l,t,r).
.VIEWLAS  <- 1

#.REPFILES   <- list.files(pattern="\\.rep")
.VIEWTRCK   <- "iSCAMViewTracker.txt"  # File containing list of report files.
.TABLEDIR   <- "../TABLES/"
.FIGUREDIR  <- "../FIGS/"









.iscamViewSetup <- function(win)
{
	#Required libraries
	require(PBSmodelling)
	
	
	#Close any open graphics devices
	graphics.off()
	closeWin(win)
	
	#Create a file list object for selection
	trckExists <- file.exists( .VIEWTRCK )
	if (trckExists)
	{
		tmp <- read.table( file = .VIEWTRCK,  as.is=TRUE,  header=TRUE,  sep="," )
		cat( "MSG (.hCamViewSetup): Viewer tracking file ",.VIEWTRCK, " found.\n" )
		ifiles <- tmp
	}
	else
	{
		ifiles=data.frame("Report Files"=.REPFILES,Select=TRUE)
	}
	
	#FIXME need to depricate this,  causes error if cannot find ifiles b/c directory is wrong	
	#dummy data frame for parameter controls for loading
	A=(matrix(1,nrow=7,ncol=7))
	dumfile = paste(ifiles[1, 4],".rep", sep="")
	#A=read.rep(dumfile)$ctrl
	
	print("OK")
	#Build data frame
	colhdr=c("ival", "lb", "ub", "phz", "prior", "mu\nshape","SD\nrate")
	rownme=c("log(Ro)","steepness","log(M)","log(Rbar)","log(Rinit)","rho","precision")
	ctrlDF<<-as.data.frame(A)
	rownames(ctrlDF)<<-rownme
	colnames(ctrlDF)<<-colhdr
	ctrlDF<<-cbind(ctrlDF,View=TRUE)
	
	#Create new window based on iscamWin.txt
	createWin("iscamWin.txt")
	
	
	#Default Graphic directory
	wdir = paste(getwd(), sep="")
	setWinVal(list(graphicDirectory=wdir))
	
}

.gui2	<- function(win)
{
	require(PBSmodelling)
	
	#Create a file list object for selection
	trckExists <- file.exists( .VIEWTRCK )
	cat(trckExists)
	if (trckExists)
	{
		tmp <- read.table( file = .VIEWTRCK,  as.is=TRUE,  header=TRUE,  sep="," )
		cat( "MSG (.hCamViewSetup): Viewer tracking file ",.VIEWTRCK, " found.\n" )
		ifiles <- tmp
	}
	else
	{
		ifiles=data.frame("Report Files"=.REPFILES,Select=TRUE)
	}
	
	
	browser()
	
	closeWin(win)
	createWin("iscamWin2.txt")
	
	# setWinVal(list(txtFigDir=.FIGUREDIR))
}

guiView	<- function()
{
	#.iscamViewSetup("iscam")
	.gui2("iSCAMView")
}


.iscamTable	<- function()
{
	print(".iscamTable")
	
	# Get the guiPerf parameters so that plot controls available.
	guiInfo <- getWinVal(scope="L")
	
	# Determine which files have been selected
	hdr	<- ifiles[ ifiles$Select, ]
	#print(hdr)
	nRuns <- nrow(hdr)
	
	if ( tableType == "catch" )
	{
		.tableCatch( hdr )
	}
	if ( tableType == "survey" )
	{
		.tableSurvey( hdr )
	}
	if ( tableType == "refpoints" )
	{
		.tableRefpoints( hdr )
	}
	if ( tableType == "forecast" )
	{
		.tableForecast( hdr )
	}
}

.tableForecast	<- function( hdr )
{
	#This function is currently set up to reproduce the
	#decision table for the Pacific herring assessments.
	#Colum headers are:
	#Stock 2010SSB 2011Age-4 PreFishBiomass Cuttoff AvailHarvest
	#Use poor average good recruitment for PrefishBio and AvailHarvest
	#All results are based on median values.
	nRuns <- nrow(hdr)
	
	cutoff = c(10700, 12100, 17600, 21200, 18800)
	xTable = NULL
	for ( i in 1:nRuns )
	{
		repObj	<- read.rep(paste(hdr$Control.File[i],".rep", sep = ""))
		mcfile=paste(hdr$Control.File[i],".mcmc", sep="")
		if(file.exists(mcfile))
			repObj$mcmc = read.table( mcfile, header=TRUE )
		else
			cat("NB. MCMC file missing.")
		
		Bo = quantile(repObj$mcmc$bo, prob=0.5)*1000
		#cutoff[i] = 0.25*Bo
		SSB = quantile(repObj$mcmc$SSB, prob=0.5)*1000
		Bt4 = quantile(repObj$mcmc$Age.4, prob=0.5)*1000
		Btpoor = quantile(repObj$mcmc$Poor, prob=0.5)*1000
		Btaverage = quantile(repObj$mcmc$Average, prob=0.5)*1000
		Btgood = quantile(repObj$mcmc$Good, prob=0.5)*1000
		tmp1 = c(hdr$Stock[i],SSB, Bt4, Btpoor, Btaverage, Btgood,cutoff[i])
		
		
		fhr <-function(bt, eps, hr=0.2)
		{
			if(bt < eps)
				return(0)
			if(bt-hr*bt>eps)
				return(hr)
			else if(bt>eps && bt-hr*bt<=eps)
				return((bt-eps)/bt)
		}
		
		Ctpoor = fhr(Btpoor, cutoff[i], 0.2)*Btpoor
		Ctaverage = fhr(Btaverage, cutoff[i], 0.2)*Btaverage
		Ctgood = fhr(Btgood, cutoff[i], 0.2)*Btgood
		
		tmp1 = c(tmp1, Ctpoor, Ctaverage, Ctgood)
		tmp1[-1] = prettyNum(round(as.numeric(tmp1[-1]), 0), big.mark=",")
		xTable=rbind(xTable, tmp1)
		
	}
	cgrp = c("", "Pre-fishery forecast biomass", "", "Available harvest")
	ncgrp = c(3, 3, 1, 3)
	colnames(xTable) = c("Stock", "SSB","4+ Biomass","Poor","Average","Good"
					, "Cutoff", "Poor","Average","Good")
	print(xTable)
	cap = "Estimated spawning stock biomass,  age-4+ biomass and pre-fishery
			biomass for poor average and good recruitment,  cutoffs,  and 
			available harvest."
	fn=paste(.TABLEDIR, "DecisionTable.tex", sep="")
	tmp <- latex(xTable, file=fn, rowname=NULL, longtable=FALSE
		, landscape=FALSE, cgroup=cgrp, n.cgroup=ncgrp
		, caption=cap, label="TableCatchAdvice", na.blank=TRUE, vbar=FALSE
		, size="small")

	cat("Latex table saved as", fn)
}

.tableRefpoints	<- function( hdr )
{
	#Creat a latex table with stock-specific reference points.
	nRuns <- nrow(hdr)
	rpTable = NULL
	for( i in 1:nRuns )
	{
		repObj	<- read.rep(paste(hdr$Control.File[i],".rep", sep = ""))
		repObj$fit	<- read.fit(paste(hdr$Control.File[i],"", sep = ""))
		tmp = c(round(repObj$fmsy, 2)
			,round(c(repObj$msy, repObj$bo, 0.25*repObj$bo
			, repObj$bmsy,  0.8*repObj$bmsy, 0.4*repObj$bmsy)*1000, 0)
			, round(repObj$sbt[length(repObj$sbt)-1]/repObj$bo, 2))
		tmp = c(hdr$Stock[i], round(repObj$fit$nopar, 0), prettyNum(tmp, big.mark=","))
		
		rpTable=rbind(rpTable, tmp)
	}
	colnames(rpTable) = c("Stock", "No.", "\\fmsy","MSY","$B_0$", "0.25$B_0$", 
						"\\bmsy","0.8\\bmsy", "0.4\\bmsy", "Spawn depletion")
	cap <- "Reference points"

	fn=paste(.TABLEDIR, "RefPointsTable.tex", sep="")
	tmp <- latex(rpTable, file=fn, rowname=NULL, longtable=FALSE
		, landscape=FALSE, cgroup=NULL, n.cgroup=NULL
		, caption=cap, label="TableRefPoints", na.blank=TRUE, vbar=FALSE
		, size="small")

	cat("Latex table saved as", fn)
}

.tableSurvey	<- function( hdr )
{
	#Create a latex table with all the relative abundance data in it.
	#hdr contains thie checked files for extracting the survey data
	nRuns <- nrow(hdr)
	ytTable=cgrp=ncgrp=NULL
	for( i in 1:nRuns )
	{
		repObj	<- read.rep(paste(hdr$Control.File[i],".rep", sep = ""))
		it 		<- t(as.matrix(repObj$it))
		yr		<- t(as.matrix(repObj$iyr))
		nCols	<- dim(it)[2]
		for(j in 1:nCols)
		{
			tmp = cbind(yr[, j], it[, j])
			colnames(tmp) <- c("Year", paste("Survey", j))
			ytTable = cbind(ytTable, tmp)
		}
		cgrp = c(cgrp, hdr$Stock[i])
		ncgrp = c(ncgrp, 2*nCols)
	}
	print(ytTable)
	cap <- "Abundance for each survey by year for each stock."
	
	fn=paste(.TABLEDIR, "SurveyTable.tex", sep="")
	tmp <- latex(ytTable, file=fn, rowname=NULL, longtable=TRUE
		, landscape=TRUE, lines.page=60, cgroup=cgrp, n.cgroup=ncgrp
		, caption=cap, label="TableSurvey", na.blank=TRUE, vbar=TRUE
		, size="tiny")
	
	cat("Latex table saved as", fn)	
}




getRepObj   <- function(fileName)
{
	print("getRepObj")
	cat(fileName)
	repObj	<- read.admb(fileName)
	#repObj$stock = hdr$Stock[i]
	repObj$Control.File = fileName
	mcfile=paste(fileName,".mcmc", sep="")
	if(file.exists(mcfile))
	{
		repObj$mcmc = read.table( mcfile, header=TRUE )
		mcfile=paste(fileName,".mcst", sep="")
		repObj$mcsbt = read.table( mcfile, header=FALSE )
	}
	else
		cat("NB. MCMC file missing.")
	
	return(repObj)
}

.mpdView	<- function()
{
	print(".mpdView")
	
	# Get the guiPerf parameters so that plot controls available.
	guiInfo <- getWinVal(scope="L")
	
	# List of gui changes
	guiChanges <- list()
	
	# Determine which files have been selected
	hdr	<- ifiles[ ifiles$Select, ]
	print(hdr)
	nRuns <- nrow(hdr)
	#hdr$Report.File contains the vector of report files to examine.

	## TODO auto layout stuff.
	# if ( autoLayout )
	# 	{
	# 		if (nRuns < 4 )
	# 		{
	# 			winCols <- 1
	# 			winRows <- nRuns
	# 		}
	# 		else if( nRuns >3 && nRuns<=6 )
	# 		{
	# 			winCols <- 2
	# 			winRows <- 3
	# 		}
	# 		
	# 		if (plotType=="selectivity")
	# 		{
	# 			winCols <- 1
	# 			winRows <- 1
	# 		}
	# 	}
	# 	
	# 	newPlot <- TRUE
	# 	if ( !autoLayout )
	# 	{
	# 		winCols <- ncols
	# 		winRows <- nrows
	# 		
	# 		## Determin if newPlot should be set to TRUE
	# 		# This sets newPlot to TRUE if current plot panel is the last row and
	# 	    # last column of the layout.
	# 		newPlot <- FALSE
	# 	    mfg <- par( "mfg" )
	# 	    if ( (mfg[1]==mfg[3]) && (mfg[2]==mfg[4]) )
	# 	      newPlot <- TRUE
	# 	    
	# 	}
	# 	
	# 	## Update the gui
	# 	guiChanges$ncols = winCols
	# 	guiChanges$nrows = winRows
	# 	setWinVal( guiChanges )
	# 	
	# 	
	# 	## set graphical parameters
	# 	if (newPlot)
	# 	{
	# 		if ( plotbyrow )
	# 			par( oma=.VIEWOMA, mar=.VIEWMAR, mfrow=c(winRows,winCols) )
	# 		else
	# 			par( oma=.VIEWOMA, mar=.VIEWMAR, mfcol=c(winRows,winCols) )
	# 	}
	
	## TODO Loop over runs and plot corresponding graph
	par(las=1)
	for(i in 1:nRuns)
	{
		# Read the report file
		#repObj	<- read.rep(hdr$Report.File[i])
		repObj	<- read.admb(hdr$Control.File[i])
		repObj$stock = hdr$Stock[i]
		repObj$Control.File = hdr$Control.File[i]
		mcfile=paste(hdr$Control.File[i],".mcmc", sep="")
		if(file.exists(mcfile))
			repObj$mcmc = read.table( mcfile, header=TRUE )
		else
			cat("NB. MCMC file missing.")
		
		if(plotType=="sbmcmc" || plotType=="depletionmcmc")
		{
			mcfile=paste(hdr$Control.File[i],".mcst", sep="")
			repObj$mcsbt = read.table( mcfile, header=FALSE )
		}
	
		# Conditional statements for radio button selections of plots
		# Add new calls to radio buttons here.
		if ( plotType=="catch" )
		{
			.plotCatch( repObj, legend.txt=NULL )
		}
	
		if ( plotType=="survey" )
		{
			.plotIndex( repObj, annotate=TRUE )
		}
	
		if ( plotType=="biomass" )
	    {
	    	.plotBiomass( repObj, annotate=TRUE )
	    }
	
		if ( plotType=="depletion" )
		{
			.plotDepletion( repObj, annotate=TRUE )
		}
	
		if ( plotType=="surveyfit" )
		{
			cat("plotSurveyFit\n")
			.plotSurveyfit( repObj, annotate=TRUE )
		}
	
		if ( plotType=="mortality" )
		{
			.plotMortality( repObj, annotate=TRUE )
		}
	
		if ( plotType=="agecomps" )
		{
			.plotAgecomps( repObj, meanAge=TRUE )
		}
	
		if ( plotType=="catchresid" )
		{
			.plotCatchResiduals( repObj, annotate=TRUE )
		}
	
		if ( plotType=="surveyresid" )
		{
			.plotSurveyResiduals( repObj, annotate=TRUE )
		}
	
		if ( plotType=="recresid" )
		{
			.plotRecruitmentResiduals( repObj )
		}
	
		if ( plotType=="agecompsresid" )
		{
			.plotAgeCompResiduals( repObj )
		}
	
		if ( plotType=="recruitment" )
		{
			.plotRecruitment( repObj )
		}
	
		if ( plotType=="meanwt" )
		{
			.plotMeanwt( repObj )
		}
	
		if ( plotType=="selectivity" )
		{
			.plotSelectivity( repObj )
		}
	
		if ( plotType=="ssbretro" )
		{
			.plotSSBretrospective( repObj )
		}
		
		if ( plotType=="stockrecruit" )
		{
			.plotStockRecruit( repObj )
		}
	
		if ( plotType=="productivity" )
		{
			.plotRecruitsPerSpawner( repObj )
		}
		
		if ( plotType=="parameters" )
		{
			.plotMarginalPosteriors( repObj )
		}
	
		if ( plotType=="refpoints" )
		{
			.plotReferencePoints( repObj )
		}
	
		if ( plotType=="kobeplot" )
		{
			#admbObj = read.admb( "iscam" )
			#admbObj$mcmc = read.table( "iscam.mcmc", header=TRUE )
			.plotStockStatus( repObj )
		}
	
		if ( plotType=="sbmcmc" )
		{
			#admbObj = read.admb( "iscam" )
			#admbObj$mcsbt = read.table( "sbt.mcmc" )
			#admbObj$mcmc = read.table( "iscam.mcmc", header=TRUE )
			
			.plotSbtPosterior( repObj )
		}
		
		if ( plotType=="risktable" )
		{
			.plotRiskTable( hdr )
		}
	
		if ( plotType=="depletionmcmc" )
		{
			#admbObj = read.admb( "iscam" )
			#admbObj$mcsbt = read.table( "sbt.mcmc" )
			#admbObj$mcmc = read.table( "iscam.mcmc", header=TRUE )
			.plotSbtPosterior( repObj, TRUE, annotate=TRUE )
		}
	
		if ( plotType=="traceplot" )
		{
			#admbObj = read.admb( "iscam" )
			#admbObj$mcmc = read.table( "iscam.mcmc", header=TRUE )
			.plotMCMCtrace( repObj )
		}
		
		if ( plotType=="mcmcpairs" )
		{
			.plotMCMCpairs( repObj )
		}
	
		if ( plotType=="simplot" )
		{
			admbObj = read.admb( "iscam" )
			admbObj$sim = read.rep( "iscam.sim" )
			
			.plotSimulationSummary( admbObj )
		}
	
		# if(gLetter)
		# {
		# 	mfg=par("mfg")
		# 	if(mfg[1]==1 && mfg[2]==1) ix=1
		# 	if(mfg[1]==2 && mfg[2]==1) ix=2
		# 	if(mfg[1]==1 && mfg[2]==2) ix=3
		# 	if(mfg[1]==2 && mfg[2]==2) ix=4
		# 	gletter(ix)
		# }
	}  #end of loop over nRuns
	
	## Saving PDF of figure in current device
	# if(savePDF)
	# 	{
	# 		#if(graphicFileName=="Graphics File Name or File Prefix")
	# 		#	graphicFileName = "iSCAM"
	# 		#filePrefix=graphicFileName
	# 		#dev.copy2pdf( file=paste( filePrefix,"fig",plotType,".pdf",sep="" ) )
	# 		setWinVal(c(figFileName=plotType))
	# 		.saveGraphic()
	# 	}
}




.plotSimulationSummary	<- function( admbObj )
{
	print("	.plotSimulationSummary")
	op=par(no.readonly=T)
	with( admbObj, {
		par(las=1,mar=c(5, 5, 1, 1), oma=c(1, 1, 0, 0), mfcol=c(2, 2))
		
		#Spawing biomass
		plot(yr, sim$sbt[1:length(yr)], type="l", xlab="Year", ylab="Spawning biomass (t)")
		lines(yr, sbt[1:length(yr)], lwd=5, col=colr(1, 0.25))
		gletter(1)
		
		#Fishing mortality rates
		if(is.matrix(ft)){
			yy=t(as.matrix(ft))
			y2=t(as.matrix(sim$ft))
		}
		else{
			yy=ft
			y2=sim$ft
		}
		matplot(yr,y2,type="l",lty=1,ylab="Fishing mortality",xlab="Year")
		matlines(yr,yy,lwd=5,col=colr(1:3,0.25),lty=1)
		gletter(2)
		
		#Survey abundance
		.plotSurveyfit( admbObj )
		gletter(3)
		
		#Survey residuals
		.plotSurveyResiduals( admbObj )
		gletter(4)
	})
	par(op)
}
# Deprecated
# .plotStockStatus	<- function( admbObj )
# {
# 	print("	.plotStockStatus")
# 	
# 	require(MASS)
# 	require( KernSmooth)
# 	fried.egg=function(xx,yy,...)
# 	{
# 		bw=25
# 		bwx=diff(extendrange(xx))/bw; bwy=diff(extendrange(yy))/bw
# 		#bwx=(max(xx)-min(xx))/bw
# 		#bwy=(max(yy)-min(yy))/bw
# 		est <- bkde2D(cbind(xx,yy),bandwidth=c(bwx,bwy),gridsize=c(81, 81))
# 		est$fhat=est$fhat/max(est$fhat)
# 		#plot(xx,yy,pch=".",col="dark grey",xlab=NA,ylab=NA,type="n")
# 		#text(max(xx),max(yy),labels="D",adj=c(1,1))
# 		lvs=c(0.05,0.25,0.75,0.95)
# 		maxct=max(lvs)
# 		nlvs=length(lvs)
# 		thelines=contourLines(est$x1,est$x2,est$fhat,levels=lvs)
# 		iclr=colr("khaki", 0.9)
# 		polygon(thelines[[nlvs-3]]$x,thelines[[nlvs-3]]$y,col=iclr,border=iclr,lwd=1)
# 		iclr=colr("snow", 0.9)
# 		polygon(thelines[[nlvs-2]]$x,thelines[[nlvs-2]]$y,col=iclr,border=iclr,lwd=2)
# 		iclr=colr("yellow", 0.9)
# 		polygon(thelines[[nlvs-1]]$x,thelines[[nlvs-1]]$y,col=iclr,border=iclr,lwd=3)
# 		polygon(thelines[[nlvs]]$x,thelines[[nlvs]]$y,col="lightyellow",border="yellow",lwd=1)
# 		#contour(est$x1,est$x2,est$fhat,drawlabels=T,add=T,levels=lvs,lty=1,lwd=1,labcex= 0.7)
# 		#Add salt and pepper
# 		#xi=sample(1:length(xx),300)
# 		#points(xx[xi],yy[xi],pch=".",col=grey(0:10/10))
# 	}
# 	
# 	#KOBE plots
# 	## This routine needs some work to accomodate multiple
# 	## fleets. Also need to address the Fmsy calculation in iscam.
# 	with(admbObj, {
# 		xx = sbt[1:length(yr)]/bmsy
# 		yy = ft/fmsy  #yy can be a matrix
# 		yy[yy==0]=NA; ii=!is.na(yy)
# 
# 		matplot(xx, (yy[ii]), type="l", xlim=c(0,max(2,xx)), 
# 		ylim=c(0,max(2,yy[ii])),xlab="Spawning biomass/SBmsy", ylab="Ft/Fmsy")
# 		rect(0, 0, 1, 1, col=colr("yellow", 0.5), border=NA)
# 		rect(1, 1, max(2, xx),max(2, yy),col=colr("yellow", 0.5),border=NA)
# 		rect(1, 0, max(2, xx), 1, col=colr("green", 0.5), border=NA)
# 		rect(0, 1, 1, max(2, yy), col=colr("red", 0.5), border=NA)
# 		## add bayesian fried egg 
# 		## need to get marginal samples for ft and sbt
# 		## to correctly plot the fried egg uncertainty
# 		xxx=sbt[length(yr)]/mcmc$bmsy
# 		yyy=ft[length(yr)]/mcmc$fmsy
# 		fried.egg(xxx, yyy)
# 		lines(xx, yy[ii], type="l", col=colr(1, 1))
# 		text(xx, yy[ii], yr, cex=0.75, col=colr(1, 1))
# 	})
# }

# Deprecated
# .plotReferencePoints	<- function( admbObj )
# {
# 	print("	.plotReferencePoints")
# 	op=par(no.readonly=T)
# 	with(admbObj, {
# 		par(las=1,mar=c(5, 5, 1, 1), oma=c(1, 1, 1, 0))
# 		par(mfcol=c(2, 2))
# 		for(i in 8:11)
# 		{
# 			ps=mcmc[, i]
# 			xl=range(ps)
# 			hist(ps,xlab=colnames(mcmc[i]),prob=T, 
# 				main="", ylab="",
# 				xlim=xl)#, ...)
# 		}
# 		mtext(c("MSY reference points", "Probability density",paste(stock)), 
# 		c(1, 2, 3),outer=T, line=-1, las=0)
# 	})
# 	
# 	par(op)
# 	
# }

# .plotSbtPosterior	<- function( admbObj, depletion=FALSE, annotate=FALSE )
# {
# 	## Median and 95% CI for spawning biomass or depletion
# 	print("	.plotSbtPosterior")
# 	with( admbObj, {
# 		xx <- yr
# 		yy <- t(apply( mcsbt, 2, quantile, probs=c(0.5, 0.025, 0.975) ))
# 		
# 		bo <- quantile( mcmc$bo, probs=c(0.5, 0.025, 0.975) )
# 		bmsy <- quantile( mcmc$bmsy, probs=c(0.5, 0.025, 0.975) )
# 		yl="Spawning biomass (1000 t)"
# 		
# 		if(depletion)
# 		{
# 			yy <- t(apply( mcsbt/mcmc$bo, 2, quantile, probs=c(0.5, 0.025, 0.975) ))
# 			yl = "Spawning biomass depletion"
# 		}
# 			
# 
# 		matplot(xx, yy, type="n", xlab="Year", 
# 			ylab=yl, axes=FALSE, 
# 			ylim=c(0, 1.0*max(yy)), main=paste(stock) )
# 		
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		
# 		polygon(c(xx, rev(xx)), 
# 			c(yy[,2],rev(yy[,3])), 
# 			col=colr("red",0.25), border=NA)
# 		print(yy)
# 		lines(xx, yy[,1], lwd=1.5)
# 		#matlines(xx, yy)
# 		
# 		if(!depletion)
# 		{
# 			points(min(xx)-0.4,bo[1])
# 			arrows(min(xx)-0.4,bo[2], min(xx)-0.4, bo[3], length=0)
# 		}
# 		
# 		if ( annotate && depletion )
# 		{
# 			#mfg <- par( "mfg" )
# 			#if ( mfg[1]==1 && mfg[2]==1 )
# 			#legend( "top",legend=c( "Spawning biomass","MSY depletion level",
# 			#	"Upper stock reference","Limit reference point"),
# 			#	bty='n',lty=c(1,2,2,2),lwd=c(1,rlvl),pch=-1,ncol=2 )
# 			
# 			#Delinate critical zone,  cautious zone, healthy zone.
# 			rect(min(xx)-5,-0.5,max(xx)+5,0.4*bmsy[1]/bo[1],col=colr("red",0.1), border=NA)
# 			rect(min(xx)-5,0.4*bmsy[1]/bo[1],max(xx)+5,0.8*bmsy[1]/bo[1],col=colr("yellow",0.1),border=NA)
# 			rect(min(xx)-5,0.8*bmsy[1]/bo[1],max(xx)+5,1.5,col=colr("green",0.1), border=NA)
# 		}
# 		
# 		
# 	})
# }

# .plotMCMCpairs	<- function( admbObj )
# {
# 	## pairs plot of mcmc samples for theta
# 	## code kindly stolen from ARK.
# 	panel.mcmc <- function( x,y,z=modes )
#   {
# 	xMean <- mean( x,na.rm=T )
#     yMean <- mean( y,na.rm=T )
# 	points( x,y,pch=16,cex=0.5,col=colr("black", 0.15) )
# 	abline( h=yMean,v=xMean,col="blue",lty=3 )
# 	points( xMean,yMean, bg="cyan", pch=21,cex=1.25 )
# 	if ( !is.null(modes) )
#    	{
#       # This is logic to figure out what "pair" is being plotted.
#       # The modal estimates are the first row of the mcmcObj.
#       # The par()$mfg calls finds the current row and column indices of
#       # the panel being plotted.
# 
# 		xMode <- z[ par()$mfg[2] ]
#     	yMode <- z[ par()$mfg[1] ]
# 		points( xMode,yMode, bg="red", pch=22, cex=1.25 )
#     }
#   }
# 
#   panel.hist <- function( x,... )
#   {
#     # Histograms for diagonal of pairs plot (from PBS Modelling CCA).
# 	  usr <- par("usr")
#       on.exit( par(usr) )
# 	  h <- hist( x, breaks="Sturges", plot=FALSE )
# 	  breaks <- h$breaks
#       nB <- length(breaks)
# 	  y <- h$counts
#       y <- y / sum(y)
# 	  par( usr = c(usr[1:2], 0, max(y)*1.5) )
# 	  rect( breaks[-nB], 0, breaks[-1], y, col="#FFD18F" )
#       box()
#   }
# 
#   # Find the active parameters.  If the chain is all equal, then the parameter
#   # was fixed in the model configuration.  This gets a Boolean vector that
#   # indicates which columns have fixed values.
#   mcmcObj <- admbObj$mcmc[, 1:11]
#   iPars <- apply( mcmcObj,2,function(x) { sum(diff(x))!=0.0 } )
#   nPars <- sum( iPars )     # Number of active parameters in mcmc output.
# 
#   tmp <- mcmcObj[ ,iPars ]
#   tmpNames <- names( tmp )
# 
#   modes <- mcmcObj[1,]
#   pairs( tmp, panel=panel.mcmc, diag.panel=panel.hist, gap=0 )
#   mtext(admbObj$stock, side=3, outer=T, line=-1.5)
#   
# }

# Deprecated
# .mcmcTrace	<- function( admbObj, label=NULL )
# {
# 	## this function examines the trace plots for the
# 	## estimated leading parameters
# 	print("	.mcmcTrace")
# 	op=par(no.readonly=T)
# 	guiInfo <- getWinVal(scope="L")
# 	if ( !autoLayout )
# 	{
# 		winCols <- ncols
# 		winRows <- nrows
# 	}
# 	else
# 	{
# 		winCols <- 3
# 		winRows <- 4
# 	}
# 	
# 	## set graphical parameters
# 	if ( plotbyrow )
# 		par( oma=.VIEWOMA, mar=.VIEWMAR, mfrow=c(winRows,winCols) )
# 	else
# 		par( oma=.VIEWOMA, mar=.VIEWMAR, mfcol=c(winRows,winCols) )
# 	
# 	#par(las=1,mar=c(5, 4, 1, 1), oma=c(1, 2, 1, 0), mfcol=c(3, 3))
# 	plotTrace <- function( obj )
# 	{
# 	  # Input "obj" is a VECTOR of MCMC samples.
# 	  # Produces one panel trace plot.
# 
# 	  nSample <- length( obj )
# 	  plot( c(1:nSample), obj, type="n", axes=FALSE, xlab="", ylab="" )
# 	  points( c(1:nSample),obj, cex=0.5, pch=20, col="darkgray" )
# 
# 	  lines( lowess( c(1:nSample),obj,f=1/4), lty=1, lwd=1 )
# 	  abline( h=mean(obj), lty=2 )
# 
# 	  # Plot MPD point (1st element).
# 	  points( 1,obj[1], cex=1.0, pch=16, col="green" )
# 	  points( 1,obj[1], cex=1.0, pch=1 )    
# 
# 	  axis( side=1 )
# 	  axis( side=2 )
# 	  box()
# 	}
#   
# 	with(admbObj, {
# 	  # Find the active parameters.  If the chain is all equal, then the parameter
# 	  # was fixed in the model configuration.  This gets a Boolean vector that
# 	  # indicates which columns have fixed values.
# 	  mcmcObj=mcmc[, 1:11]
# 	  iPars <- apply( mcmcObj,2,function(x) { sum(diff(x))!=0.0 } )
# 	  nPars <- sum( iPars )     # Number of active parameters in mcmc output.
# 
# 	  tmp <- mcmcObj[ ,iPars ]
# 	  tmpNames <- names( tmp )
# 	
# 	  for ( i in 1:ncol(tmp) )
# 	  {
# 	    plotTrace( tmp[,i] )
# 		title(ylab=tmpNames[i], line=2)
# 		print(tmpNames[i])
# 	  }
# 
# 	  mtext(c("Sample", "Parameter",paste(stock)), c(1, 2, 3), 
# 			outer=T, line=c(-1,0, -1), las=0)
# 	  
# 	})
# 	par(op)
# }


# Deprecated
# .plotMarginalPosteriors	<- function( admbObj )
# {
# 	print("	.plotMarginalPosteriors")
# 	#Marginal distributions & priors for theta
# 	op=par(no.readonly=T)
# 	par(las=1,mar=c(5, 4, 1, 1), oma=c(1, 1, 1, 0))
# 
# 	## Read control file to get bounds and priors for theta
# 	## ctrl=read.table(A$control.file, header=F, skip=13, nrow=6)
# 
# 	with(admbObj, {
# 		std=apply(mcmc[,1:7],2,sd)
# 		nr=length(std[std!=0])
# 		if(nr > 6)
# 		{
# 			nRow=3; nCol=3
# 		}
# 		else if(nr>4)
# 		{
# 			nRow=3; nCol=2
# 		}
# 		else
# 		{
# 			nRow=2; nCol=2
# 		}
# 		par(mfcol=c(nRow, nCol))
# 		for(i in 1:7){
# 			if(std[i]!=0){
# 				ps = mcmc[, i]  #posterior samples
# 				xl=range(ps)
# 
# 				hist(ps,xlab=colnames(mcmc[i]),prob=T, 
# 					main="", ylab="", col="lightgrey",breaks=30, 
# 					xlim=xl)#, ...)
# 
# 				## Add priors
# 				nfn=c("dunif","dnorm","dlnorm","dbeta","dgamma")
# 				pt = ctrl[i, 5]+1
# 				fn=match.fun(nfn[pt])
# 				p1=ctrl[i, 6]; p2=ctrl[i, 7]
# 				#browser()
# 				if(pt!=4)
# 					curve(unlist(lapply(x,fn,p1,p2)),
# 						xl[1],xl[2],add=T, col=colr(4, 0.7), lwd=2)
# 				else
# 					curve(unlist(lapply((x-ctrl[i,2])/
# 						 (ctrl[i,3]-ctrl[i,2])
# 						,fn,p1,p2)),xl[1],xl[2],add=T, col=colr(4, 0.7), lwd=2)
# 			}
# 		}
# 		mtext(c("Parameter", "Probability density",paste(stock)), c(1, 2, 3), 
# 			outer=T, line=-1, las=0)
# 	})
# 	par(op)
# }
# DEPRECATED
# .plotStockRecruit	<- function( repObj )
# {
# 	with(repObj, {
# 		xx = sbt[1:(length(yr)-min(age))]
# 		yy = rt
# 		
# 		plot(xx, yy, type="n",ylim=c(0, max(yy, ro)),xlim=c(0, max(xx,bo)), 
# 			xlab="Spawning biomass", 
# 			ylab=paste("Age-",min(age)," recruits", sep=""), 
# 			main=paste(stock))
# 			
# 		points(xx, yy)
# 		points(xx[1],yy[1], pch=20, col="green")
# 		points(xx[length(xx)], yy[length(xx)], pch=20, col=2)
# 		
# 		st=seq(0, max(sbt, bo), length=100)
# 		if(rectype==1)
# 		{
# 			#Beverton-Holt
# 			rrt=kappa*ro*st/(bo+(kappa-1)*st)*exp(-0.5*tau^2)  
# 		}
# 		if(rectype==2)
# 		{
# 			#Ricker
# 			rrt=kappa*ro*st*exp(-log(kappa)*st/bo)/bo *exp(-0.5*tau^2) 
# 		}
# 		lines(st, rrt)
# 		ro=ro*exp(-0.5*tau^2)
# 		points(bo, ro, pch="O", col=2)
# 		points(bo, ro, pch="+", col=2)
# 	})
# }

# Deprecated
# .plotSelectivity	<- function( repObj )
# {
# 	#plot the selectivity curves (3d plots)
# 	with(repObj, {
# 		#par(mgp=c(3, 3, 5))
# 		plot.sel<-function(x, y, z, ...)
# 		{
# 			#z=exp(A$log_sel)*3
# 			#x=A$yr
# 			#y=A$age
# 			z <- z/max(z)
# 			z0 <- 0#min(z) - 20
# 			z <- rbind(z0, cbind(z0, z, z0), z0)
# 			x <- c(min(x) - 1e-10, x, max(x) + 1e-10)
# 			y <- c(min(y) - 1e-10, y, max(y) + 1e-10)
# 			clr=colorRampPalette(c("honeydew","lawngreen"))
# 			nbcol=50
# 			iclr=clr(nbcol)
# 			nrz <- nrow(z)
# 			ncz <- ncol(z)
# 			zfacet <- z[-1, -1]+z[-1, -ncz]+z[-nrz, -1]+z[-nrz, -ncz]
# 			facetcol <- cut(zfacet, nbcol)
# 			fill <- matrix(iclr[facetcol],nr=nrow(z)-1,nc=ncol(z)-1)
# 			fill[ , i2 <- c(1,ncol(fill))] <- "white"
# 			fill[i1 <- c(1,nrow(fill)) , ] <- "white"
# 
# 			par(bg = "transparent")
# 			persp(x, y, z, theta = 35, phi = 25, col = fill, expand=5, 
# 				shade=0.75,ltheta=45 , scale = FALSE, axes = TRUE, d=1,  
# 				xlab="Year",ylab="Age",zlab="Selectivity", 
# 				ticktype="simple", ...)
# 			
# 			#require(lattice)
# 			#wireframe(z, drap=TRUE, col=fill)
# 		}
# 		ix=1:length(yr)
# 		for(k in 1:ngear){
# 			plot.sel(yr, age, exp(log_sel[log_sel[,1]==k,-1]), 
# 			main=paste(stock, "Gear", k))
# 			#file.name=paste(prefix, "Fig9",letters[k],".eps", sep="")
# 			#if(savefigs) dev.copy2eps(file=file.name, height=8, width=8)
# 		}
# 		
# 	})
# }

# DEPRECATED
# .plotMeanwt	<- function( repObj )
# {
# 	#plot mean weight-at-age by cohort
# 	with(repObj, {
# 		xx = yr		## xaxis labels
# 		yy = age	## yaxis labels
# 		nage=length(age)
# 		
# 		if(sum(par("mfcol"))==2)
# 		{
# 			xl = "Cohort year";xlm=""
# 			yl = "Weight-at-age (kg)";ylm=""
# 		}
# 		else
# 		{
# 			xlm = "Cohort year";xl=""
# 			ylm = "Weight-at-age (kg)";yl=""
# 		}
# 		
# 		plot(range(xx), range(wt_obs), type="n", axes=FALSE,
# 		xlab=xl, ylab=yl, main=paste(stock))
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		grid()
# 		
# 		for(i in 1:(dim(wt_obs)[1]-1))
# 		{
# 			yy = (diag(as.matrix(wt_obs[0:-i, ]))) 
# 			xx = 1:length(yy)+yr[i]-min(age)+1
# 			
# 			yy[yy==0]=NA;xx[yy==NA]=NA
# 			lines(xx,yy)
# 			
# 			points(xx[1],yy[1],pch=20,col="steelblue",cex=0.5)
# 			points(xx[nage],yy[nage],pch=20,col="salmon",cex=0.5)
# 		}
# 		for(i in 1:dim(wt_obs)[2]-1)
# 		{
# 			yy = diag(as.matrix(wt_obs[,-1:-i]))
# 			n = length(yy)
# 			xx = yr[1]:(yr[1]+n-1)
# 			lines(xx, yy)
# 			points(xx[n], yy[n], pch=20, col="salmon", cex=0.5)
# 		}
# 		
# 		mtext(xlm, side=1, outer=T, line=0)
# 		mtext(ylm, side=2, outer=T, line=0)
# 	})
# }

# DEPRECATED
# .plotRecruitment	<- function( repObj )
# {
# 	#plot age-a recruits.
# 	with(repObj, {
# 		xx = yr
# 		yy = exp(ln_rt)
# 		yy=yy
# 		yrange=c(0, max(yy, na.rm=T))
# 		
# 		plot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
# 			xlab="Year", main=paste(stock), 
# 			ylab=paste("Age-", min(age), " recruits", sep=""))
# 		
# 		lines(xx, yy, type="h")
# 		
# 		#add 0.33 and 0.66 quantile lines
# 		qtl = quantile(yy, prob=c(1/3, 2/3))
# 		abline(h=qtl,lty=1, col=colr("darkgrey",0.75))
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 	})
# }

# DEPRECATED
# .plotCatchResiduals		<- function( repObj, annotate=FALSE )
# {
# 	#Plot residuals between observed and predicted catches
# 	#residuals (epsilon=log(obs_ct)-log(ct))
# 	with(repObj, {
# 		epsilon=log(obs_ct)-log(ct)
# 		if(is.matrix(epsilon)){
# 			xx = yr
# 			yy = t(epsilon)
# 			t1 = colSums(yy,na.rm=T)
# 			ng = length(t1[t1!=0])
# 		}else{
# 			xx = yr
# 			yy = epsilon
# 			ng = 1
# 		}
# 		#browser()
# 		absmax = abs(max(yy, na.rm=TRUE))
# 		if(absmax<=1e-3)absmax=1
# 		yrange=c(-absmax, absmax)
# 		
# 		matplot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
# 			xlab="Year", ylab="Catch residual", main=paste(stock))
# 		
# 		matlines(xx, yy, type="h", col="black",lty=1)
# 		matpoints(xx, yy,pch=1:ng, cex=0.75, col=1)
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		if ( annotate )
# 		{
# 			#n=dim(yy)[2]
# 			txt=paste("Gear",1:ng)
# 			
# 			mfg <- par( "mfg" )
# 			if ( mfg[1]==1 && mfg[2]==1 )
# 			legend( "top",legend=txt,cex=0.75, 
# 				bty='n',pch=1:ng,lwd=1,lty=-1,ncol=ng )
# 		}
# 	})
# }

# DEPRECATED
# .plotSurveyResiduals	<- function( repObj, annotate=FALSE )
# {
# 	#Plot residuals between observed and predicted relative abundance
# 	#indicies (epsilon)
# 	with(repObj, {
# 		if(is.matrix(epsilon)){
# 			xx = t(iyr)
# 			yy = t(epsilon)
# 		}else{
# 			xx = iyr
# 			yy = epsilon
# 		}
# 		
# 		absmax = abs(max(yy, na.rm=TRUE))
# 		if(absmax< 1e-3) absmax=1
# 		yrange=c(-absmax, absmax)
# 		
# 		matplot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
# 			xlab="Year", ylab="Survey residual", main=paste(stock))
# 		
# 		matlines(xx, yy, type="h", col="black")
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		if ( annotate )
# 		{
# 			n=dim(xx)[2]
# 			txt=paste("Survey",1:n)
# 			
# 			mfg <- par( "mfg" )
# 			if ( mfg[1]==1 && mfg[2]==1 )
# 			legend( "top",legend=txt,
# 				bty='n',lty=1:n,lwd=1,pch=-1,ncol=1 )
# 		}
# 	})
# }

# DEPRECATED
#.plotRecruitmentResiduals	<- function( repObj )
#{
#	#Plot the log residuals between the estimated recruits and
#	#those obtained from the recruitment model (delta)
#	with(repObj, {
#		ii = 1:min(age)
#		xx = yr[-ii]
#		yy = delta
#		absmax = abs(max(yy, na.rm=TRUE))
#		if(absmax< 1e-3) absmax=1
#		yrange=c(-absmax, absmax)
#		
#		plot(xx, yy, type="n",  axes=FALSE, ylim=yrange, 
#			xlab="Year", ylab="Recruitment residuals", main=paste(stock))
#			
#		lines(xx, yy, type="h", col="black")
#		axis( side=1 )
#		axis( side=2,  las=.VIEWLAS )
#		box()
#	})
#}
#
# DEPRECATED
# .plotIndex	<- function( repObj, annotate=FALSE )
# {
# 	#line plot for relative abundance indices
# 	with(repObj, {
# 		if(is.matrix(it)){
# 			xx=t(iyr)
# 			yy=t(it)
# 		}else{
# 			xx=iyr
# 			yy=it
# 		}
# 		n=nrow(t(as.matrix(yy)))
# 		yrange=c(0, max(yy, na.rm=TRUE))
# 		
# 		matplot(xx, yy, type="n", axes=FALSE,
# 			xlab="Year", ylab="Relative abundance", 
# 			ylim=yrange , main=paste(stock))
# 		
# 		matlines(xx, yy, col="black",type="o", pch=1:n)
# 		
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		
# 		if( annotate )
# 		{
# 			txt=paste("Survey",1:n)
# 			legend("top", txt, lty=1:n, pch=1:n, bty="n")
# 		}
# 	})
# }

# DEPRECATED
# .plotCatch	<- function( repObj, legend.txt=NULL )
# {
# 	#barplot of the observed catch
# 	with(repObj, {
# 		tmp = obs_ct
# 		iRows <- apply( tmp,1,function(x) { sum(diff(x))!=0.0 } )
# 		barplot( tmp[iRows,], names.arg=yr,axes=FALSE,  
# 			xlab="Year", ylab="Catch (1000 t)",main=paste(stock),  
# 			legend.text = legend.txt )
# 		axis( side=2, las=.VIEWLAS )
# 	})
# }

# DEPRECATED
# .plotDepletion	<- function( repObj, annotate=FALSE )
# {
# 	#plot the spawning biomass depletion level & reference points
# 	with(repObj, {
# 		xx=yr
# 		yy=sbt[1:length(xx)]/bo
# 		yrange=c(0,1.1*max(yy, na.rm=TRUE))
# 		
# 		plot(xx, yy, type="n", axes=FALSE,
# 			xlab="Year", ylab="Spawning depletion",main=paste(stock), 
# 			ylim=yrange)
# 		lines(xx, yy)
# 		rlvl=c(1.0, 0.8, 0.4)
# 		abline(h=rlvl*bmsy/bo,lty=2,lwd=0.5)
# 		
# 		
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		grid()
# 		
# 		if ( annotate )
# 		{
# 			#mfg <- par( "mfg" )
# 			#if ( mfg[1]==1 && mfg[2]==1 )
# 			#legend( "top",legend=c( "Spawning biomass","MSY depletion level",
# 			#	"Upper stock reference","Limit reference point"),
# 			#	bty='n',lty=c(1,2,2,2),lwd=c(1,rlvl),pch=-1,ncol=2 )
# 			
# 			#Delinate critical zone,  cautious zone, healthy zone.
# 			rect(min(xx)-5,-0.5,max(xx)+5,0.4*bmsy/bo,col=colr("red",0.1), border=NA)
# 			rect(min(xx)-5,0.4*bmsy/bo,max(xx)+5,0.8*bmsy/bo,col=colr("yellow",0.1),border=NA)
# 			rect(min(xx)-5,0.8*bmsy/bo,max(xx)+5,1.5,col=colr("green",0.1), border=NA)
# 		}
# 	})
# }

# DEPRECATED
#.plotBiomass	<- function( repObj, annotate=FALSE )
#{
#	#plot total biomass & spawning biomass 
#	with(repObj, {
#		xx=yr
#		yy=cbind(bt[1:length(xx)], sbt[1:length(xx)])
#		
#		yrange=c(0, 1.2*max(yy, na.rm=TRUE))
#		
#		matplot(xx, yy, type="n",axes=FALSE,
#				xlab="Year", ylab="Biomass (1000 t)",main=paste(stock), 
#				ylim=yrange)
#		
#		matlines(xx,yy,
#			type="l", col="black",
#			ylim=c(0,max(yy,na.rm=T)))
#		axis( side=1 )
#		axis( side=2, las=.VIEWLAS )
#		box()
#		
#		if ( annotate )
#		{
#			mfg <- par( "mfg" )
#			if ( mfg[1]==1 && mfg[2]==1 )
#			legend( "top",legend=c( "Pre-fishery biomass","Spawning biomass"),
#				bty='n',lty=c(1,2),lwd=c(1,1),pch=c(-1,-1),ncol=1 )
#		}
#	})	
#}

# Deprecated
# .plotSurveyfit	<- function( repObj, annotate=FALSE)
# {
# 	with(repObj, {
# 		if(is.matrix(it)){
# 			xx = t(iyr)
# 			m = apply(it,1,max, na.rm=T)
# 			yy = t(pit/m)
# 			y2 = t(it/m)
# 		}else{
# 			xx = iyr
# 			yy = pit
# 			y2 = it
# 		}
# 		n=nrow(t(as.matrix(yy)))
# 		#n=dim(xx)[2]
# 		yrange=c(0, 1.15*max(yy, y2, na.rm=TRUE))
# 		
# 		matplot(xx, yy, type="n",axes=FALSE,ylim=yrange, 
# 			xlab="Year", ylab="Relative abundance", main=paste(stock))
# 		
# 		matlines(xx, yy, col=1:n, lty=1)
# 		matpoints(xx, y2, col=1:n, pch=1:n)
# 		
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		grid()
# 		
# 		if ( annotate )
# 		{
# 			
# 			txt=paste("Survey ",1:n,", q=",round(q, 3), sep="")
# 			
# 			mfg <- par( "mfg" )
# 			#if ( mfg[1]==1 && mfg[2]==1 )
# 			legend( "top",legend=txt,
# 				bty='n',lty=1,lwd=1,pch=1:n,ncol=n, col=1:n )
# 				
# 			#print(q)
# 		}
# 	})
# }

# Deprecated
# .plotMortality	<- function( repObj, annotate=FALSE )
# {
# 	# SJDM June 5, 2011 Changed to plot Average M and sum of average Fs by gear
# 	#plot average total mortality,  fishing mortality & natural mortality
# 	with(repObj, {
# 		xx=yr
# 		if(is.matrix(ft))
# 		{
# 			yy=t(as.matrix(ft))
# 		}
# 		else
# 		{
# 			yy=as.matrix(ft)
# 		}
# 		#n=nrow(t(as.matrix(yy)))
# 		icol=apply(yy,2,function(x){sum(cumsum(x))!=0.0})
# 		ng=length(icol[icol==T])
# 		
# 		yy = cbind( rowMeans(M_tot), yy[,icol] )
# 		
# 		csyy = t(apply(yy,1, cumsum))	#cumulative sum
# 			
# 		yrange=c(0, max(csyy, na.rm=TRUE))
# 		lw = c(1, rep(1,ng))
# 		lt = 1
# 		iclr = colr(1:(ng+1),0.5)
# 		
# 		matplot(xx, csyy, type="n", axes=FALSE, log="y",  
# 			xlab="Year", ylab="Mortality rate", main=paste(stock))
# 			
# 		
# 		#lines(xx, yy[,1], lwd=2, lty=1, col=1)
# 		matlines(xx, csyy, log="y", col=iclr, lwd=lw, lty=lt)
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		grid()
# 		
# 		ry=cbind(1.e-30,csyy)
# 		for(i in 1:dim(csyy)[2])
# 		{
# 			x2=c(xx, rev(xx))
# 			y2=c(ry[,i+1], rev(ry[,i]))
# 			#browser()
# 			polygon(x2, y2, border=NA,col=colr(i,0.2), log="y")
# 		}
# 		
# 		
# 		
# 		if ( annotate )
# 		{
# 			txt = c("Natural mortality", paste("Gear",1:ng))
# 			#mfg <- par( "mfg" )
# 			#if ( mfg[1]==1 && mfg[2]==1 )
# 			legend( "topright",legend=txt,col=iclr, 
# 				bty='n',lty=lt,lwd=5,pch=-1,ncol=2)
# 		}
# 	})
# }

# DEPRECATED
# .plotAgecomps	<- function(repObj, meanAge = FALSE )
# {
# 	#Bubble plot of age-composition data
# 	#A is the observed age-comps
# 	#Ahat is the predicted age-comps (proportions)
# 	with( repObj, {
# 		if(!is.null(repObj$A)){
# 			nagear = unique(A[, 2])
# 			xrange = range(A[, 1])
# 			#par(mfcol=c(length(nagear), 1))
# 			for(i in nagear)
# 			{
# 				ac = subset(A, A[, 2]==i)
# 				xx = ac[, 1]
# 				zz = t(ac[, -1:-2])
# 				
# 				
# 				# plot proportions-at-age (cpro=TRUE)
# 				plotBubbles(zz, xval = xx, yval = age, cpro=TRUE, hide0=TRUE,  
# 					las=.VIEWLAS, xlab="Year", ylab="Age", frange=0.0, size=0.1, 
# 					bg=colr("steelblue", 0.5),main=paste(stock, "Gear", i), 
# 					xlim=xrange)
# 				
# 				if( meanAge )
# 				{
# 					tz = t(zz)
# 					p = t(tz/rowSums(tz))
# 					abar = colSums(t(tz/rowSums(tz))*age)
# 					sbar = sqrt(colSums(p*(1-p)*age))
# 					sbar = 1.96*colSums(sqrt(p*(1-p))/sqrt(age))
# 				
# 					lines( xx, abar, col=colr("steelblue", 0.75), lwd=2 )
# 				
# 					yy = c(exp(log(abar)+log(sbar)), rev(exp(log(abar)-log(sbar))))
# 					polygon(c(xx, rev(xx)),yy,border=NA,col=colr("steelblue",0.25))
# 				}
# 			}
# 		}
# 		else{print("There is no age-composition data")}
# 	})
# }

# DEPRECATED
#.plotAgeCompResiduals	<- function(repObj)
#{
#	#Bubble plot of age-composition data
#	#A is the observed age-comps
#	#Ahat is the predicted age-comps (proportions)
#	with( repObj, {
#		if(!is.null(repObj$A)){
#			nagear = unique(A[, 2])
#			xrange = range(A[, 1])
#			#par(mfcol=c(length(nagear), 1))
#			j=0
#			for(i in nagear)
#			{
#				j=j+1
#				ac = subset(A_nu, A_nu[, 2]==i)
#				xx = ac[, 1]
#				zz = t(ac[, -1:-2])
#			
#				# plot residuals
#				plotBubbles(zz, xval = xx, yval = age, rres=FALSE, hide0=TRUE,  
#					las=.VIEWLAS, xlab="Year", ylab="Age", frange=0.0, size=0.5*age_tau2[j],
#					bg=colr("white", 0.5), xlim=xrange,main=paste(stock, "Gear", i))
#				title(main=paste("Variance=",round(age_tau2[j],3)),line=-1,cex.main=0.75)
#			}
#			
#		}
#		else{print("There is no age-composition data")}
#	})
#}
#
.plotPriors	<-function()
{
	## This function plots the prior distributions for the parameters
	## specified in the control file (actually read from the report file)
	## after the model has been run.  The user can also change the mean, 
	## SD, shape or rate parameters,  and overlay alternative priors.
	guiInfo <- getWinVal(scope="L")
	print(".plotPriors")
	op <- par(no.readonly=TRUE)
	par(las=1, mfcol=c(3, 2))
	df <- guiInfo$ctrlDF			#Control data frame
	np <- dim(df)[1]		#number of parameters
	
	## loop over each parameter, set x-limits
	for(i in 1:np)
	{
		x	<- seq(df[i,2], df[i, 3], length=200)
		p1	<- df[i, 6]
		p2	<- df[i, 7]
		view<- df[i, 8]
		pt	<- df[i, 5]+1	#prior type
		xl	<- range(x)
		if(view)
		{
			## Add priors
			nfn=c("dunif","dnorm","dlnorm","dbeta","dgamma")
			fn=match.fun(nfn[pt])
			y = unlist(lapply(x,fn,p1,p2))
			
			if(i==2 && pt==4)  #steepness parameter BH model
				y= unlist(lapply((x-0.2)/0.8,fn,p1,p2))
			
			plot(x,y,type="n",
				xlab=rownames(df[i,]),ylab="Density")
			
			polygon(c(x,rev(x)),c(y,rep(0,length=length(x))),
					col=colr(df[i,4], 0.25))
			
			title(main=nfn[pt])
			#if(pt!=4)
			#	curve(unlist(lapply(x,fn,p1,p2)),
			#		xl[1],xl[2],add=T, col=4, lty=2)
			#else
			#	curve(unlist(lapply((x-ctrl[i,2])/
			#		 (ctrl[i,3]-ctrl[i,2])
			#		,fn,p1,p2)),xl[1],xl[2],add=T, col=4, lty=2)
			
		}
	}
	par(op)
	#browser()
}

#-------------------------------------------------------------------------------#
#   Miscellaneous functions called by above routines                            #
#                                                                               #
#-------------------------------------------------------------------------------#

.plotBubbles <- function (z, xval = FALSE, yval = FALSE, dnam = FALSE, rpro = FALSE, 
    cpro = FALSE, rres = FALSE, cres = FALSE, powr = 1, size = 0.2, 
    lwd = 1, fg.clrs = c("black", "red", "blue"), bg.clrs = c("white","white","white"), 
    hide0 = FALSE, debug = FALSE, ...) 
{
    dz <- dim(z)
    ny <- dz[1]
    nx <- dz[2]
    xval1 <- 1:nx
    yval1 <- 1:ny
    nx1 <- nx
    ny1 <- ny
    if (mode(xval) == "logical") {
        if (xval[1]) {
            xval1 <- z[1, ]
            ny1 <- ny - 1
        }
    }
    if (mode(yval) == "logical") {
        if (yval[1]) {
            yval1 <- z[, 1]
            nx1 <- nx - 1
        }
    }
    xind <- (nx - nx1 + 1):nx
    x2 <- xval1[xind]
    yind <- (ny - ny1 + 1):ny
    y2 <- yval1[yind]
    if ((mode(xval) == "numeric") & (length(xval) == nx1)) 
        x2 <- xval
    if ((mode(yval) == "numeric") & (length(yval) == ny1)) 
        y2 <- yval
    zz <- array(z[yind, xind], dim = c(length(yind), length(xind)), 
        dimnames = dimnames(z))
    if (dnam & !is.null(dimnames(zz))) {
        if (!is.null(dimnames(zz)[[2]])) {
            xnam <- as.numeric(dimnames(zz)[[2]])
            if (!any(is.na(xnam)) && (length(xnam) == 1 || all(diff(xnam) > 
                0 | all(diff(xnam) < 0)))) 
                x2 <- xnam
        }
        if (!is.null(dimnames(zz)[[1]])) {
            ynam <- as.numeric(dimnames(zz)[[1]])
            if (!any(is.na(ynam)) && (length(ynam) == 1 || all(diff(ynam) > 
                0 | all(diff(ynam) < 0)))) 
                y2 <- ynam
        }
    }
    xx <- rep(x2, each = length(y2))
    yy <- rep(y2, length(x2))
    minz <- min(zz, na.rm = TRUE)
    maxz <- max(zz, na.rm = TRUE)
    if (rpro | cpro) {
        if (minz < 0) {
            zz <- zz - minz
            minz <- 0
            maxz <- max(zz, na.rm = TRUE)
        }
    }
    if (rpro) {
        zs <- apply(zz, 1, sum, na.rm = TRUE)
        zz <- sweep(zz, 1, zs, "/")
    }
    if (cpro) {
        zs <- apply(zz, 2, sum, na.rm = TRUE)
        zz <- sweep(zz, 2, zs, "/")
    }
    if (rres) {
        zm <- apply(zz, 1, mean, na.rm = TRUE)
        zz <- sweep(zz, 1, zm, "-")
    }
    if (cres) {
        zm <- apply(zz, 2, mean, na.rm = TRUE)
        zz <- sweep(zz, 2, zm, "-")
    }
    zNA <- is.na(zz) | is.nan(zz) | is.infinite(zz)
    zz[zNA] <- 0
    z0 <- sign(zz) * abs(zz)^abs(powr)
    z1 <- z3 <- z0
    z1[z0 <= 0] <- NA
    z3[z0 < 0 | z0 > 0] <- NA
    z2 <- -z0
    z2[z0 >= 0] <- NA
    za <- max(z0, na.rm = TRUE)
    zb <- min(z0, na.rm = TRUE)
    zM <- max(abs(z0))
    sz1 <- max(za * size/zM, 0.001)
    sz2 <- max(-zb * size/zM, 0.001)
    if (debug) 
        browser()
    symbols(xx, yy, circles = as.vector(abs(z0)), inches = size, 
        fg = 0, ...)
    if (debug) 
        browser()
    if (!hide0 && !all(is.na(z3))) {
        symbols(xx, yy, circles = as.vector(z3), inches = 0.001, 
            fg = fg.clrs[1], bg = bg.clrs[3], lwd = lwd, add = TRUE, ...)
    }
    if (!all(is.na(z2))) {
        symbols(xx, yy, circles = as.vector(z2), inches = sz2, 
            fg = fg.clrs[1], bg = bg.clrs[2], lwd = lwd, add = TRUE, ...)
    }
    if (!all(is.na(z1))) {
        symbols(xx, yy, circles = as.vector(z1), inches = sz1, 
            fg = fg.clrs[1], bg = bg.clrs[1], lwd = lwd, add = TRUE, ...)
    }
    invisible(z0)
}


.saveGraphic	<- function()
{
	#This function saves the current graphic
	#with the specified file name as a pdf file.
	print("MSG (.saveGraphic)")
	
	# Get the guiPerf parameters so that plot controls available.
	guiInfo <- getWinVal(scope="L")
	
	if(graphicDirectory=="Working Directory")
		graphicDirectory = getwd()
	if(figFileName=="Figure File Name")
	{
		figFileName = plotType
		setWinVal(c(figFileName=plotType))
	}
	figFile = figFileName
	
	fileName = paste(graphicDirectory, "/iscam_fig_", figFile,".pdf", sep="")
	dev.copy2pdf(file=fileName)
	cat("Graphic saved as:", fileName, "\n")
}

.selectDirectory <- function()
{
	#This function selects the working directory for saving 
	#graphical images.
	print("MSG (.selectDirctory)")
	selectDir(usewidget="graphicDirectory")
}

.runSimulationTrials <- function()
{
	#This function runs the simulation trials 
	#given nTrials and randomSeed from the gui.
	
	#Output contains log2 ratios for estimated theta values
	#Also compares bias ratios in estimated msy-based reference points.
	
	# Get the guiPerf parameters so that plot controls available.
	guiInfo <- getWinVal(scope="L")
	
	# First read iscamMC.rda if it exists.
	if(file.exists("iscamMC.rda"))
	{
		MC=dget("iscamMC.rda")
	}
	else
	{
		## Run the model once to determine # of estimated parameters
		arg = paste("./iscam -nox -sim", 111)
		system(arg)
		admbObj <- read.admb("iscam")
	
		itheta  <- grep("theta", admbObj$fit$names)
		ix		<- as.numeric(substr(admbObj$fit$names[itheta],7,7))
		theta <- admbObj$ctrl[ix, 1]
		theta_p <- NULL
		rPoints	<- NULL
		for(i in 1:nTrials)
		{
			seed = randomSeed + 2*(i-1)
			arg = paste("./iscam -nox -sim", seed)
			system(arg)
		
			admbObj <- read.admb("iscam")
			admbObj$sim = read.rep( "iscam.sim" )
			theta_p = rbind( theta_p, admbObj$fit$est[itheta] )
			rp = c(admbObj$fmsy/admbObj$sim$fmsy, 
				admbObj$msy/admbObj$sim$msy, 
				admbObj$bmsy/admbObj$sim$bmsy)
			rPoints = rbind( rPoints, log2(rp) )
			.plotSimulationSummary( admbObj )
			print(paste("Trial ", i))
		}
		pn	<- c("log(Ro)","h","log(m)","log(Rbar)","log(Rinit)","rho","vartheta")
		theta_dev<-log2(t(t(theta_p)/theta))
		#theta_dev<-cbind(theta_dev, rPoints)
		#browser()

		colnames(theta_dev)= (pn[ix])#(pn[itheta])
		colnames(rPoints) = c("Fmsy","MSY","Bmsy")
		MC=list(theta_dev=theta_dev, rPoints=rPoints)
	}
	
	
	op = par(no.readonly=T)
	par(mfcol=c(2, 1), las=1)
	boxplot(MC$theta_dev, ylab="log2(theta/theta')", 
			ylim=c(-1, 1) )
			#ylim=c(-max(c(1, abs(theta_dev))),max(c(1, abs(theta_dev)))) )
	abline(h=0, col="grey")
	
	boxplot(MC$rPoints,ylab="log2 bias", 
			ylim=c(-1, 1) )
			#ylim=c(-max(c(1, abs(rPoints))), max(c(1, abs(rPoints)))) )
	abline(h=0, col="grey")	
	par(op)

	dput(MC,file="iscamMC.rda")
	
	#TODO	write theta_dev and rPoints to files so you don't have to repeat.
}


## The following has been deprecated and is scheduled for deletion
##.subView	<- function()
##{
##	guiInfo <- getWinVal(scope="L")
##	
##	
##	##Graphics options
##	# Check graphics options
##	if( autoLayout && ! plotbyrow )
##	{
##		par(mar=.VIEWMAR, oma=.VIEWOMA, las=.VIEWLAS, mfrow=c(1, 1))
##	} 
##	
##	if( plotbyrow )
##	{
##		winCols <- ncols
##		winRows <- nrows
##		par(mar=.VIEWMAR, oma=.VIEWOMA, las=.VIEWLAS, mfrow=c(winRows, winCols))
##	}
##	
##	if( !plotbyrow )
##	{
##		winCols <- ncols
##		winRows <- nrows
##		par(mar=.VIEWMAR, oma=.VIEWOMA, las=.VIEWLAS, mfcol=c(winRows, winCols))
##	}
##	#call the plotting function again if user change plotting options.
##	.mpdView()
##}
##

cat("Type: \n guiView()\n to start the iSCAM gui")
#######################
#Type: guiView()
#to start gui
#######################