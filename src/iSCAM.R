##                                       ##
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
# TODO: Finish  test priors widget to examine priors for parameters             #
# TODO: Save input List as an .rda file (dput & dget) for saving scenarios      #
#                                                                               #
#                                                                               #
#                                                                               #
#                                                                               #
#-------------------------------------------------------------------------------#

require(hacks)	#transparent colors using the function colr("color",tranparency)
require(Riscam)	#custom library built specifically for iscam.


# Graphics defaults.
.VIEWCEX    <- 1            # Generic default cex for axis labels etc.
.VIEWPANCEX <- 1            # Default cex for panLab.

.VIEWMAR  <- c(3, 3, 1, 1)  # Multi-panel plots: plot margin sizes c(b,l,t,r).
.VIEWOMA  <- c(2, 2, 1, 1)  # Multi-panel plots: outer margin sizes c(b,l,t,r).
.VIEWLAS  <- 1

.REPFILES <- list.files(pattern="\\.rep")











.iscamViewSetup <- function(win)
{
	#Required libraries
	require(PBSmodelling)
	
	
	#Close any open graphics devices
	graphics.off()
	closeWin()
	
	#Create a file list object for selection
	ifiles=data.frame("Report Files"=.REPFILES,Select=TRUE)
		
	#Read iscam report file and get controls for priors tab
	A=read.rep("iscam.rep")
	#Build data frame
	colhdr=c("ival", "lb", "ub", "phz", "prior", "mu\nshape","SD\nrate")
	rownme=c("log(Ro)","steepness","log(M)","log(Rbar)","rho","precision")
	ctrlDF<<-as.data.frame(A$ctrl)
	rownames(ctrlDF)<<-rownme
	colnames(ctrlDF)<<-colhdr
	ctrlDF<<-cbind(ctrlDF,View=TRUE)
	
	#Create new window based on iscamWin.txt
	createWin("iscamWin.txt")
	
	#Default Graphic directory
	wdir = paste(getwd(), sep="")
	setWinVal(list(graphicDirectory=wdir))
	
}

guiView	<- function()
{
	.iscamViewSetup("iscam")
}

.subView	<- function()
{
	guiInfo <- getWinVal(scope="L")
	
	##Graphics options
	# Check graphics options
	if( autolayout && ! plotbyrow )
	{
		par(mar=.VIEWMAR, oma=.VIEWOMA, las=.VIEWLAS, mfrow=c(1, 1))
	} 
	
	if( plotbyrow )
	{
		winCols <- ncols
		winRows <- nrows
		par(mar=.VIEWMAR, oma=.VIEWOMA, las=.VIEWLAS, mfrow=c(winRows, winCols))
	}
	
	#call the plotting function again if user change plotting options.
	.mpdView()
}


.mpdView	<- function()
{
	print(".mpdView")
	
	# Get the guiPerf parameters so that plot controls available.
	guiInfo <- getWinVal(scope="L")
	
	# Determine which files have been selected
	hdr	<- ifiles[ ifiles$Select, ]
	#hdr$Report.Files contains the vector of report files to examine.
	
	# Read the report file
	repObj	<- read.rep(hdr$Report.Files)
	#repObj	<- read.rep("iscam.rep")
	
	
	
	
	
	
	
	# Conditional statements for radio button selections of plots
	if ( plotType=="catch" )
	{
		.plotCatch( repObj, legend.txt=NULL )
	}
	
	if ( plotType=="survey" )
	{
		.plotIndex( repObj, annotate=annotate )
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
		.plotAgecompresiduals( repObj )
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
	
	if ( plotType=="stockrecruit" )
	{
		.plotStockRecruit( repObj )
	}
	
	if ( plotType=="parameters" )
	{
		admbObj = read.admb( "iscam" )
		admbObj$mcmc = read.table( "iscam.mcmc", header=TRUE )
		.plotMarginalPosteriors( admbObj )
	}
	
	if ( plotType=="refpoints" )
	{
		admbObj = read.admb( "iscam" )
		admbObj$mcmc = read.table( "iscam.mcmc", header=TRUE )
		.plotReferencePoints( admbObj )
	}
	
	if ( plotType=="kobeplot" )
	{
		admbObj = read.admb( "iscam" )
		admbObj$mcmc = read.table( "iscam.mcmc", header=TRUE )
		.plotStockStatus( admbObj )
	}
	
	if ( plotType=="sbmcmc" )
	{
		admbObj = read.admb( "iscam" )
		admbObj$mcsbt = read.table( "sbt.mcmc" )
		admbObj$mcmc = read.table( "iscam.mcmc", header=TRUE )
		.plotSbtPosterior( admbObj )
	}
	
	if ( plotType=="depletionmcmc" )
	{
		admbObj = read.admb( "iscam" )
		admbObj$mcsbt = read.table( "sbt.mcmc" )
		admbObj$mcmc = read.table( "iscam.mcmc", header=TRUE )
		.plotSbtPosterior( admbObj, TRUE )
	}
	
	if ( plotType=="traceplot" )
	{
		admbObj = read.admb( "iscam" )
		admbObj$mcmc = read.table( "iscam.mcmc", header=TRUE )
		.mcmcTrace( admbObj )
	}
	
	if ( plotType=="simplot" )
	{
		admbObj = read.admb( "iscam" )
		admbObj$sim = read.rep( "iscam.sim" )
		.plotSimulationSummary( admbObj )
	}
	
	## Saving PDF of figure in current device
	if(savePDF)
	{
		if(graphicFileName=="Graphics File Name or File Prefix")
			graphicFileName = "iSCAM"
		filePrefix=graphicFileName
		dev.copy2pdf( file=paste( filePrefix,"fig",plotType,".pdf",sep="" ) )
	}
}

.plotSimulationSummary	<- function( admbObj )
{
	print("	.plotSimulationSummary")
	op=par(no.readonly=T)
	with( admbObj, {
		par(las=1,mar=c(5, 5, 1, 1), oma=c(1, 1, 0, 0), mfcol=c(2, 2))
		
		#Spawing biomass
		plot(yrs, sim$sbt, type="l", xlab="Year", ylab="Spawning biomass (t)")
		lines(yrs, sbt, lwd=5, col=colr(1, 0.25))
		gletter(1)
		
		#Fishing mortality rates
		matplot(yr,t(sim$ft[1:3,]),type="l",lty=1,ylab="Fishing mortality",xlab="Year")
		matlines(yr,t(ft[1:3,]),lwd=5,col=colr(1:3,0.25),lty=1)
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

.plotStockStatus	<- function( admbObj )
{
	print("	.plotStockStatus")
	
	require(MASS)
	require( KernSmooth)
	fried.egg=function(xx,yy,...)
	{
		bw=25
		bwx=diff(extendrange(xx))/bw; bwy=diff(extendrange(yy))/bw
		#bwx=(max(xx)-min(xx))/bw
		#bwy=(max(yy)-min(yy))/bw
		est <- bkde2D(cbind(xx,yy),bandwidth=c(bwx,bwy),gridsize=c(81, 81))
		est$fhat=est$fhat/max(est$fhat)
		#plot(xx,yy,pch=".",col="dark grey",xlab=NA,ylab=NA,type="n")
		#text(max(xx),max(yy),labels="D",adj=c(1,1))
		lvs=c(0.05,0.25,0.75,0.95)
		maxct=max(lvs)
		nlvs=length(lvs)
		thelines=contourLines(est$x1,est$x2,est$fhat,levels=lvs)
		iclr=colr("khaki", 0.9)
		polygon(thelines[[nlvs-3]]$x,thelines[[nlvs-3]]$y,col=iclr,border=iclr,lwd=1)
		iclr=colr("snow", 0.9)
		polygon(thelines[[nlvs-2]]$x,thelines[[nlvs-2]]$y,col=iclr,border=iclr,lwd=2)
		iclr=colr("yellow", 0.9)
		polygon(thelines[[nlvs-1]]$x,thelines[[nlvs-1]]$y,col=iclr,border=iclr,lwd=3)
		polygon(thelines[[nlvs]]$x,thelines[[nlvs]]$y,col="lightyellow",border="yellow",lwd=1)
		#contour(est$x1,est$x2,est$fhat,drawlabels=T,add=T,levels=lvs,lty=1,lwd=1,labcex= 0.7)
		#Add salt and pepper
		#xi=sample(1:length(xx),300)
		#points(xx[xi],yy[xi],pch=".",col=grey(0:10/10))
	}
	
	#KOBE plots
	## This routine needs some work to accomodate multiple
	## fleets. Also need to address the Fmsy calculation in iscam.
	with(admbObj, {
		xx = sbt[1:length(yr)]/bmsy
		yy = ft/fmsy  #yy can be a matrix
		yy[yy==0]=NA; ii=!is.na(yy)

		matplot(xx, (yy[ii]), type="l", xlim=c(0,max(2,xx)), 
		ylim=c(0,max(2,yy[ii])),xlab="Spawning biomass/SBmsy", ylab="Ft/Fmsy")
		rect(0, 0, 1, 1, col=colr("yellow", 0.5), border=NA)
		rect(1, 1, max(2, xx),max(2, yy),col=colr("yellow", 0.5),border=NA)
		rect(1, 0, max(2, xx), 1, col=colr("green", 0.5), border=NA)
		rect(0, 1, 1, max(2, yy), col=colr("red", 0.5), border=NA)
		## add bayesian fried egg 
		## need to get marginal samples for ft and sbt
		## to correctly plot the fried egg uncertainty
		xxx=sbt[length(yr)]/mcmc$bmsy
		yyy=ft[length(yr)]/mcmc$fmsy
		fried.egg(xxx, yyy)
		lines(xx, yy[ii], type="l", col=colr(1, 1))
		text(xx, yy[ii], yr, cex=0.75, col=colr(1, 1))
	})
}


.plotReferencePoints	<- function( admbObj )
{
	print("	.plotReferencePoints")
	op=par(no.readonly=T)
	with(admbObj, {
		par(las=1,mar=c(5, 5, 1, 1), oma=c(1, 1, 0, 0))
		par(mfcol=c(2, 2))
		for(i in 7:10)
		{
			ps=mcmc[, i]
			xl=range(ps)
			hist(ps,xlab=colnames(mcmc[i]),prob=T, 
				main="", ylab="",
				xlim=xl)#, ...)
		}
		mtext(c("MSY reference points", "Probability density"), 
		c(1, 2),outer=T, line=-1, las=0)
	})
	
	par(op)
	
}

.plotSbtPosterior	<- function( admbObj, depletion=FALSE )
{
	#To do.  Add Sbo and posterior at syr-1
	#for reference.
	print("	.plotSbtPosterior")
	with( admbObj, {
		xx <- yrs
		yy <- t(apply( mcsbt, 2, quantile, probs=c(0.5, 0.025, 0.975) ))
		
		bo <- quantile( mcmc$bo, probs=c(0.5, 0.025, 0.975) )
		
		if(depletion) yy <- t(t(yy)/bo)

		matplot(yrs, yy, type="n", xlab="Year", 
			ylab="Spawning biomass (t)", axes=FALSE, 
			ylim=c(0, max(yy)) )
		
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		
		polygon(c(xx, rev(xx)), 
			c(yy[,2],rev(yy[,3])), 
			col=colr("red",0.5), border=NA)
		
		lines(xx, yy[,1], lwd=1.5)
		#matlines(xx, yy)
		
		if(!depletion)
		{
			points(min(xx)-0.4,bo[1])
			arrows(min(xx)-0.4,bo[2], min(xx)-0.4, bo[3], length=0)
		}
		
		
		
	})
}

.mcmcTrace	<- function( admbObj, label=NULL )
{
	## this function examines the trace plots for the
	## estimated leading parameters
	print("	.mcmcTrace")
	plotTrace <- function( obj )
	{
	  # Input "obj" is a VECTOR of MCMC samples.
	  # Produces one panel trace plot.

	  nSample <- length( obj )
	  plot( c(1:nSample), obj, type="n", axes=FALSE, xlab="", ylab="" )
	  points( c(1:nSample),obj, cex=0.5, pch=20, col="darkgray" )

	  lines( lowess( c(1:nSample),obj,f=1/4), lty=1, lwd=1 )
	  abline( h=mean(obj), lty=2 )

	  # Plot MPD point (1st element).
	  points( 1,obj[1], cex=1.0, pch=16, col="green" )
	  points( 1,obj[1], cex=1.0, pch=1 )    

	  axis( side=1 )
	  axis( side=2 )
	  box()
	}
  
	with(admbObj, {
	  # Find the active parameters.  If the chain is all equal, then the parameter
	  # was fixed in the model configuration.  This gets a Boolean vector that
	  # indicates which columns have fixed values.
	  mcmcObj=mcmc
	  iPars <- apply( mcmcObj,2,function(x) { sum(diff(x))!=0.0 } )
	  nPars <- sum( iPars )     # Number of active parameters in mcmc output.

	  tmp <- mcmcObj[ ,iPars ]
	  tmpNames <- names( tmp )
	
	  for ( i in 1:ncol(tmp) )
	  {
	    plotTrace( tmp[,i] )
	    #panLab( 0.5, 0.9, cex=1.0, tmpNames[i] )  
	  }

	  if ( !is.null(label) )
	    mtext( side=3, line=-0.5, cex=1.0, outer=T, label )
	  mtext( side=1, line=0.5, cex=1.0, outer=T, "Sample" )
	  
	})
}



.plotMarginalPosteriors	<- function( admbObj )
{
	print("	.plotMarginalPosteriors")
	#Marginal distributions & priors for theta
	op=par(no.readonly=T)
	par(las=1,mar=c(5, 4, 1, 1), oma=c(1, 1, 0, 0))

	## Read control file to get bounds and priors for theta
	## ctrl=read.table(A$control.file, header=F, skip=13, nrow=6)

	with(admbObj, {
		std=apply(mcmc[,1:6],2,sd)
		nr=length(std[std!=0])/2
		par(mfcol=c(nr, 2))
		for(i in 1:6){
			if(std[i]!=0){
				ps = mcmc[, i]  #posterior samples
				xl=range(ps)

				hist(ps,xlab=colnames(mcmc[i]),prob=T, 
					main="", ylab="", col="lightgrey",breaks=30, 
					xlim=xl)#, ...)

				## Add priors
				nfn=c("dunif","dnorm","dlnorm","dbeta","dgamma")
				pt = ctrl[i, 5]+1
				fn=match.fun(nfn[pt])
				p1=ctrl[i, 6]; p2=ctrl[i, 7]
				if(pt!=4)
					curve(unlist(lapply(x,fn,p1,p2)),
						xl[1],xl[2],add=T, col=4, lty=2)
				else
					curve(unlist(lapply((x-ctrl[i,2])/
						 (ctrl[i,3]-ctrl[i,2])
						,fn,p1,p2)),xl[1],xl[2],add=T, col=4, lty=2)
			}
		}
		mtext(c("Parameter", "Probability density"), c(1, 2), 
			outer=T, line=-1, las=0)
	})
	par(op)
}

.plotStockRecruit	<- function( repObj )
{
	with(repObj, {
		xx = sbt[1:(length(yr)-min(age))]
		yy = rt
		
		plot(xx, yy, type="n",ylim=c(0, max(yy)),xlim=c(0, max(xx,bo)), 
			xlab="Spawning biomass", ylab=paste("Age-",min(age)," recruits", sep=""))
			
		points(xx, yy)
		
		st=seq(0, max(sbt, bo), length=100)
		if(rectype==1)
		{
			#Beverton-Holt
			rrt=kappa*ro*st/(bo+(kappa-1)*st)*exp(-0.5*tau^2)  
		}
		if(rectype==2)
		{
			#Ricker
			rrt=kappa*ro*st*exp(-log(kappa)*st/bo)/bo *exp(-0.5*tau^2) 
		}
		lines(st, rrt)
		ro=ro*exp(-0.5*tau^2)
		points(bo, ro, pch="O", col=2)
		points(bo, ro, pch="+", col=2)
	})
}

.plotSelectivity	<- function( repObj )
{
	#plot the selectivity curves (3d plots)
	with(repObj, {
		#par(mgp=c(3, 3, 5))
		plot.sel<-function(x, y, z)
		{
			#z=exp(A$log_sel)*3
			#x=A$yr
			#y=A$age
			z0 <- 0#min(z) - 20
			z <- rbind(z0, cbind(z0, z, z0), z0)
			x <- c(min(x) - 1e-10, x, max(x) + 1e-10)
			y <- c(min(y) - 1e-10, y, max(y) + 1e-10)
			clr=colorRampPalette(c("honeydew", "lawngreen"))
			nbcol=50
			iclr=clr(nbcol)
			nrz <- nrow(z)
			ncz <- ncol(z)
			zfacet <- z[-1, -1]+z[-1, -ncz]+z[-nrz, -1]+z[-nrz, -ncz]
			facetcol <- cut(zfacet, nbcol)
			fill <- matrix(iclr[facetcol],nr=nrow(z)-1,nc=ncol(z)-1)
			fill[ , i2 <- c(1,ncol(fill))] <- "white"
			fill[i1 <- c(1,nrow(fill)) , ] <- "white"

			par(bg = "transparent")
			persp(x, y, z, theta = 35, phi = 25, col = fill, expand=3, 
				shade=0.75,ltheta=45 , scale = FALSE, axes = TRUE, d=1,  
				xlab="Year",ylab="Age",zlab="Selectivity",
				ticktype="detailed")
			
			#require(lattice)
			#wireframe(z, drap=TRUE, col=fill)
		}
		ix=1:length(yr)
		for(k in 1:ngear){
			plot.sel(yr, age, exp(log_sel[log_sel[,1]==k,-1]))
			#file.name=paste(prefix, "Fig9",letters[k],".eps", sep="")
			#if(savefigs) dev.copy2eps(file=file.name, height=8, width=8)
		}
		
	})
}

.plotMeanwt	<- function( repObj )
{
	#plot mean weight-at-age by cohort
	with(repObj, {
		xx = yr		## xaxis labels
		yy = age	## yaxis labels
		nage=length(age)
		
		plot(range(xx), range(wt_obs), type="n", axes=FALSE,
		xlab="Cohort year", ylab="Weight-at-age (kg)")
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		
		for(i in 1:dim(wt_obs)[1])
		{
			#ir = (age-min(age))+i
			#xx = yr[i]+(age-min(age))
			#yy = (diag(as.matrix(wt_obs[ir, ])))
			yy = (diag(as.matrix(wt_obs[0:-i, ]))) 
			xx = 1:length(yy)+yr[i]-min(age)+1
			
			yy[yy==0]=NA;xx[yy==NA]=NA
			lines(xx,yy)

			points(xx[1],yy[1],pch=20,col="steelblue",cex=0.5)
			points(xx[nage],yy[nage],pch=20,col="salmon",cex=0.5)
			
		}
	})
}

.plotRecruitment	<- function( repObj )
{
	#plot age-a recruits.
	with(repObj, {
		xx = yr
		yy = exp(ln_rt)
		yy=yy
		yrange=c(0, max(yy, na.rm=T))
		
		plot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
			xlab="Year", 
			ylab=paste("Age-", min(age), " recruits", sep=""))
		
		lines(xx, yy, type="h")
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
	})
}

.plotCatchResiduals		<- function( repObj, annotate=FALSE )
{
	#Plot residuals between observed and predicted catches
	#residuals (epsilon=log(obs_ct)-log(ct))
	with(repObj, {
		epsilon=log(obs_ct)-log(ct)
		if(is.matrix(epsilon)){
			xx = yr
			yy = t(epsilon)
		}else{
			xx = yr
			yy = epsilon
		}
		
		absmax = abs(max(yy, na.rm=TRUE))
		if(absmax<=1e-3)absmax=1
		yrange=c(-absmax, absmax)
		
		matplot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
			xlab="Year", ylab="Catch residual")
		
		matlines(xx, yy, type="h", col="black")
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		if ( annotate )
		{
			n=dim(xx)[2]
			txt=paste("Gear",1:n)
			
			mfg <- par( "mfg" )
			if ( mfg[1]==1 && mfg[2]==1 )
			legend( "top",legend=txt,
				bty='n',lty=1:n,lwd=1,pch=-1,ncol=1 )
		}
	})
}

.plotSurveyResiduals	<- function( repObj, annotate=FALSE )
{
	#Plot residuals between observed and predicted relative abundance
	#indicies (epsilon)
	with(repObj, {
		if(is.matrix(epsilon)){
			xx = t(iyr)
			yy = t(epsilon)
		}else{
			xx = iyr
			yy = epsilon
		}
		absmax = abs(max(yy, na.rm=TRUE))
		if(absmax< 1e-3) absmax=1
		yrange=c(-absmax, absmax)
		
		matplot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
			xlab="Year", ylab="Survey residual")
		
		matlines(xx, yy, type="h", col="black")
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		if ( annotate )
		{
			n=dim(xx)[2]
			txt=paste("Survey",1:n)
			
			mfg <- par( "mfg" )
			if ( mfg[1]==1 && mfg[2]==1 )
			legend( "top",legend=txt,
				bty='n',lty=1:n,lwd=1,pch=-1,ncol=1 )
		}
	})
}

.plotRecruitmentResiduals	<- function( repObj )
{
	#Plot the log residuals between the estimated recruits and
	#those obtained from the recruitment model (delta)
	with(repObj, {
		ii = 1:min(age)
		xx = yr[-ii]
		yy = delta
		absmax = abs(max(yy, na.rm=TRUE))
		if(absmax< 1e-3) absmax=1
		yrange=c(-absmax, absmax)
		
		plot(xx, yy, type="n",  axes=FALSE, ylim=yrange, 
			xlab="Year", ylab="Recruitment residuals")
			
		lines(xx, yy, type="h", col="black")
		axis( side=1 )
		axis( side=2,  las=.VIEWLAS )
		box()
	})
}


.plotIndex	<- function( repObj, annotate=FALSE )
{
	#line plot for relative abundance indices
	with(repObj, {
		if(is.matrix(it)){
			xx=t(iyr)
			yy=t(it)
		}else{
			xx=iyr
			yy=it
		}
				yrange=c(0, max(yy, na.rm=TRUE))
		
		matplot(xx, yy, type="n", axes=FALSE,
			xlab="Year", ylab="Relative abundance", 
			ylim=yrange )
		
		matlines(xx, yy, col="black",type="o", pch=1:ncol(yy))
		
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		
		if( annotate )
		{
			n=nrow(t(as.matrix(yy)))
			txt=paste("Survey",1:n)
			legend("top", txt, lty=1:n, pch=1:n, bty="n")
		}
	})
}

.plotCatch	<- function( repObj, legend.txt=NULL )
{
	#barplot of the observed catch
	with(repObj, {
		barplot(obs_ct, names.arg=yr,axes=FALSE, 
			xlab="Year", ylab="Catch (1000 t)", 
			legend.text = legend.txt)
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

.plotSurveyfit	<- function( repObj, annotate=FALSE)
{
	with(repObj, {
		if(is.matrix(it)){
			xx = t(iyr)
			yy = t(pit)
			y2 = t(it)
		}else{
			xx = iyr
			yy = pit
			y2 = it
		}
		yrange=c(0, max(yy, y2, na.rm=TRUE))
		
		matplot(xx, yy, type="n",axes=FALSE,ylim=yrange, 
			xlab="Year", ylab="Relative abundance")
		
		matlines(xx, yy, col="black")
		matpoints(xx, y2, col="black")
		
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		
		if ( annotate )
		{
			n=dim(xx)[2]
			txt=rep(c( "Predicted","Observed"),1)
			
			mfg <- par( "mfg" )
			if ( mfg[1]==1 && mfg[2]==1 )
			legend( "top",legend=txt,
				bty='n',lty=c(1,-1),lwd=1,pch=c(-1,"1"),ncol=1 )
		}
	})
}


.plotMortality	<- function( repObj, annotate=FALSE )
{
	#plot average total mortality,  fishing mortality & natural mortality
	with(repObj, {
		xx=yr
		if(is.matrix(ft))
			yy=t(as.matrix(ft))
		else
			yy=ft
		
		yy = cbind( yy, rowMeans(M_tot) )	
		yrange=c(0, max(yy, na.rm=TRUE))
		lw = c(rep(1,ngear),2)
		lt = c(rep(1,ngear),1)
		
		matplot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
			xlab="Year", ylab="Mortality rate")
			
		matlines(xx, yy, col="black", lwd=lw, lty=lt)
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		grid()
		
		if ( annotate )
		{
			txt = c(paste("Gear",1:ngear),"Natural mortality")
			mfg <- par( "mfg" )
			if ( mfg[1]==1 && mfg[2]==1 )
			legend( "right",legend=txt,
				bty='n',lty=c(1:ngear,1),lwd=lw,pch=-1,ncol=1)
		}
	})
}


.plotAgecomps	<- function(repObj, meanAge = FALSE )
{
	#Bubble plot of age-composition data
	#A is the observed age-comps
	#Ahat is the predicted age-comps (proportions)
	with( repObj, {
		if(!is.null(repObj$A)){
			nagear = unique(A[, 2])
			#par(mfcol=c(length(nagear), 1))
			for(i in nagear)
			{
				ac = subset(A, A[, 2]==i)
				xx = ac[, 1]
				zz = t(ac[, -1:-2])
				
				
				# plot proportions-at-age (cpro=TRUE)
				plotBubbles(zz, xval = xx, yval = age, cpro=TRUE, hide0=TRUE,  
					las=.VIEWLAS, xlab="Year", ylab="Age", frange=0.0, size=0.1, 
					bg=colr("steelblue", 0.5))
				
				if( meanAge )
				{
					tz = t(zz)
					p = t(tz/rowSums(tz))
					abar = colSums(t(tz/rowSums(tz))*age)
					sbar = sqrt(colSums(p*(1-p)*age))
					sbar = 1.96*colSums(sqrt(p*(1-p))/sqrt(age))
				
					lines( xx, abar, col=colr("steelblue", 0.75), lwd=2 )
				
					yy = c(exp(log(abar)+log(sbar)), rev(exp(log(abar)-log(sbar))))
					polygon(c(xx, rev(xx)),yy,border=NA,col=colr("steelblue",0.25))
				}
			}
		}
		else{print("There is no age-composition data")}
	})
}

.plotAgecompresiduals	<- function(repObj)
{
	#Bubble plot of age-composition data
	#A is the observed age-comps
	#Ahat is the predicted age-comps (proportions)
	with( repObj, {
		if(!is.null(repObj$A)){
			nagear = unique(A[, 2])
			#par(mfcol=c(length(nagear), 1))
			for(i in nagear)
			{
				ac = subset(A_nu, A_nu[, 2]==i)
				xx = ac[, 1]
				zz = t(ac[, -1:-2])
			
				# plot residuals
				plotBubbles(zz, xval = xx, yval = age, rres=FALSE, hide0=TRUE,  
					las=.VIEWLAS, xlab="Year", ylab="Age", frange=0.0, size=0.1,
					bg=colr("white", 0.5))
			}
		}
		else{print("There is no age-composition data")}
	})
}

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
	
	fileName = paste(graphicDirectory, "/", plotType,".pdf", sep="")
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
#######################
#Type: guiView()
#to start gui
#######################