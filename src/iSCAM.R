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
#-------------------------------------------------------------------------------#

require(Riscam)	#custom library built specifically for iscam.


# Graphics defaults.
.VIEWCEX    <- 1            # Generic default cex for axis labels etc.
.VIEWPANCEX <- 1            # Default cex for panLab.

.VIEWMAR  <- c(2, 2, 1, 1)  # Multi-panel plots: plot margin sizes c(b,l,t,r).
.VIEWOMA  <- c(2, 2, 1, 1)  # Multi-panel plots: outer margin sizes c(b,l,t,r).
.VIEWLAS  <- 2

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
	
	#Create new window based on iscamWin.txt
	createWin("iscamWin.txt")
	
}

.iscamViewSetup("test")

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
		.plotIndex( repObj, annotate=FALSE )
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
		.plotAgecomps( repObj )
	}
	
	if ( plotType=="surveyresid" )
	{
		.plotSurveyResiduals( repObj, annotate=TRUE )
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
}

.plotStockRecruit	<- function( repObj )
{
	with(repObj, {
		xx = sbt[1:(length(yr)-min(age))]
		yy = rt
		
		plot(xx, yy, type="n",ylim=c(0, max(yy)),xlim=c(0, max(xx)), 
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
		yrange=c(-absmax, absmax)
		
		matplot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
			xlab="Year", ylab="Residual")
		
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
		
		matlines(xx, yy, col="black",type="o")
		
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
	})
}

.plotCatch	<- function( repObj, legend.txt=NULL )
{
	#barplot of the observed catch
	with(repObj, {
		barplot(obs_ct, names.arg=yr,axes=FALSE, 
			xlab="Year", ylab="Catch (t)", 
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
		
		yy = cbind( rowMeans(M_tot), yy )	
		yrange=c(0, max(yy, na.rm=TRUE))
		
		matplot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
			xlab="Year", ylab="Mortality rate")
			
		matlines(xx, yy, col="black")
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		grid()
		
		if ( annotate )
		{
			txt = paste("Gear",1:ngear)
			mfg <- par( "mfg" )
			if ( mfg[1]==1 && mfg[2]==1 )
			legend( "top",legend=txt,
				bty='n',lty=1:ngear,lwd=1,pch=-1,ncol=1)
		}
	})
}

.plotAgecomps	<- function(repObj)
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
					las=.VIEWLAS, xlab="Year", ylab="Age", frange=0.0, size=0.2, bg="honeydew")
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
					las=.VIEWLAS, xlab="Year", ylab="Age", frange=0.0, size=0.2, bg="white")
			}
		}
		else{print("There is no age-composition data")}
	})
}

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
