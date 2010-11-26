#Rcode for SCAM.R
## ---------------------------------------------------------------- ##
## ---------------------------------------------------------------- ##
## flag for saving figures
   savefigs = FALSE
   prefix = "iscam User Guide/iscamFigs/phake"     #figure prefix
   tprefix = "iscam User Guide/iscamTables/phake"  #table prefix
## list of figures:
## fig1()		#Age-schedule information Length, Growth/maturity
## fig2()		#Total biomass and spawning biomass & Bo Bmsy
## fig3()		#Observed & predicted relative abundance data
## fig4()		#Estimated fishing mortality rate
## fig5()		#Marginal distributions & priors for theta
## fig6()		#Marginal distributions for MSY reference points.
## fig7()		#KOBE plot (need marginals for Sbt and Ft)
## fig8()		#DFO FMF plot
## fig9()		#3D Plots for time-varying selectivity curves
## fig10()		#Stock recruitment (fig10a) & St vs ln(Rt/St) (fig10b)
## fig11()		#Residuals in abundance indices, catch, & surveys
## fig12()		#Spawning biomass depletion & reference points
## fig13()		#Bubble plots for catch-at-age & residuals
## fig14()		#Estimates of sage recruits
## fig15()		#Observed landings barplot
## fig16()		#Spawning biomass & depletion with credible intervals
## fig17()		#Retrospective plots for spawning biomass
##
##
## list of tables:
## table1()		#Recent trends in spawning stock biomass & CI
## table2()		#MLE and Median estimates of theta[1:6]
## ________________________________________________________________ ##
require(Hmisc)
require(lattice)
require(Riscam)	#custom library built specifically for iscam.
source("reptoRlist.R")

A=read.admb('iscam')



#A=reptoRlist("iSCAM.rep")
#need to dump the control file to the Report file. DONE
#A$control.file=scan("iSCAM.rep", skip=1, nline=1, what="character")
#A$theta=read.fit("iSCAM")
B=reptoRlist("iSCAM.sim")
#if(!exists("A$mc")) A$mc=read.psv("iscam.psv")
if(!exists("A$mcmc")) A$mcmc=read.table("iscam.mcmc", header=T)
if(!exists("A$mc.sbt")) A$mc.sbt=read.table("sbt.mcmc", header=F)


fig1 <- function(A)
{	#Age-schedule information Length, Growth/maturity
	op=par(no.readonly=T)
	par(mfcol=c(2, 1))
	plot(A$age, A$la, xlab="Age", ylab="Length (cm)",type="o",  
	ylim=c(0,max(A$la)))
	plot(A$age, A$wa, xlab="Age", ylab="Weight (kg)/fecundity")
	lines(A$age, A$fa)
	legend("bottom",c("Weight-at-age", "Fecundity-at-age"), 
	lty=c(-1, 1), pch=c(1, -1), bty="n")
	
	fn=paste(prefix,"Fig1.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=4)
	par(op)
}

fig2 <- function(A)
{	#Total biomass and spawning biomass with reference points
	with(A,{
	matplot(yrs,cbind(bt, sbt),type="l", xlab="Year", 
		ylab="Biomass (t)", col=1, 
		ylim=c(0, 1.1*max(bt, sbt, bo)))
	
	abline(h=c(bo, bmsy), lty=3)
	legend("topright", c("total biomass", "spawning biomass"), 
		lty=c(1, 2), col=1, bty="n")
	})
	fn=paste(prefix,"Fig2.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
}

fig3 <- function(A)
{	#Observed & predicted relative abundance data
	with(A,{
		ii=iyr%in%yr	#match indices in case of retrospective
		plot(iyr[ii],it[ii], xlab="Year", ylab="Relative abundance", 
			ylim=c(0, max(it, pit)))
		q=exp(mean(log(it[ii])-log(pit[ii])))
		lines(iyr[ii], pit[ii]*q)
		legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), 
			pch=c(1, -1), bty="n")	
	})
	fn=paste(prefix,"Fig3.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
}

fig4 <- function(A)
{	#Estimated fishing mortality
	with(A, {
		mft=t(as.matrix(ft))
		mft[is.na(mft)]=0
		matplot(yr, mft, lty=1, xlab="Year",ylim=c(0, max(mft)),  
		 ylab="Fishing mortality", type="l")
		
	})
	fn=paste(prefix,"Fig4.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
}


fig5 <- function(A,...)
{	#Marginal distributions & priors for theta
	op=par(no.readonly=T)
	par(las=1,mar=c(5, 4, 1, 1), oma=c(1, 1, 0, 0))
	
	## Read control file to get bounds and priors for theta
	## ctrl=read.table(A$control.file, header=F, skip=13, nrow=6)
	
	with(A, {
		std=apply(mcmc[,1:6],2,sd)
		nr=length(std[std!=0])/2
		par(mfcol=c(nr, 2))
		for(i in 1:6){
			if(std[i]!=0){
				ps = mcmc[, i]  #posterior samples
				xl=range(ps)
				
				hist(ps,xlab=colnames(mcmc[i]),prob=T, 
					main="", ylab="", col="lightgrey",breaks=30, 
					xlim=xl, ...)
				
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
	fn=paste(prefix,"Fig5.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
	par(op)
	require(MCMCpack)
	mc=mcmc(A$post.samp[, 1:6], 1, thin=1)
	varnames(mc)=c("log(Ro)", "h", "log(M)", "log(Rbar)", "rho","vartheta")
	densityplot(mc, layout=c(2, 3),panel=function(x,...)
				{panel.densityplot(x,...)
				panel.abline(v=quantile(x,prob=c(0.025,0.5,0.975)),lty=2,col="grey")})
				
	xyplot(mc, layout=c(2, 3),col="grey", thin=1, 
		panel=function(x, ...)
		{	
			panel.lines(x, ...)
			panel.loess(x,span=2/3, ..., col.line="black")
			panel.loess(x,span=1/10, ..., col.line="red")
			
			#Cumulative quantile not working?
			#can't seem to get the proper y values.
			cquantile <- function(z,probs) {
		        cquant <- matrix(0, nrow = length(z), length(probs))
		        for (i in seq(along = z)) if (is.R()) {
		            cquant[i, ] <- quantile(z[1:i], probs = probs, names = FALSE)
		        }
		        else {
		            cquant[i, ] <- quantile(z[1:i], probs = probs)
		        }
		        cquant <- as.data.frame(cquant)
		        names(cquant) <- paste(formatC(100 * probs, format = "fg", 
		            wid = 1, digits = 7), "%", sep = "")
		        return(cquant)
			}
			mci=cquantile(x,c(0.025, 0.5, 0.975))
			#browser()
			#llines(mci,..., col.line="green")
		})
}

fig6 <- function(A, ...)
{	#Marginal distributions for MSY reference points.
	#plot Bo,  Bmsy,  MSY,  Fmsy
	op=par(no.readonly=T)
	with(A, {
		par(las=1,mar=c(5, 5, 1, 1), oma=c(1, 1, 0, 0))
		par(mfcol=c(2, 2))
		for(i in 7:10)
		{
			ps=mcmc[, i]
			xl=range(ps)
			hist(ps,xlab=colnames(mcmc[i]),prob=T, 
				main="", ylab="",
				xlim=xl, ...)
		}
		mtext(c("MSY reference points", "Probability density"), 
		c(1, 2),outer=T, line=-1, las=0)
	})
	fn=paste(prefix,"Fig6.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
	par(op)
}

fig7 <- function(A)
{	#KOBE plots
	## This routine needs some work to accomodate multiple
	## fleets. Also need to address the Fmsy calculation in iscam.
	with(A, {
		xx = sbt[1:length(yr)]/bmsy
		yy = ft/fmsy  #yy can be a matrix
		yy[yy==0]=NA; ii=!is.na(yy)

		matplot(xx, (yy[ii]), type="l", xlim=c(0,max(2,xx)), 
		ylim=c(0,max(2,yy[ii])),xlab="bstatus", ylab="fstatus")
		rect(0, 0, 1, 1, col="yellow", border=NA)
		rect(1, 1, max(2, xx),max(2, yy),col="yellow",border=NA)
		rect(1, 0, max(2, xx), 1, col="green", border=NA)
		rect(0, 1, 1, max(2, yy), col="red", border=NA)
		## add bayesian fried egg 
		## need to get marginal samples for ft and sbt
		## to correctly plot the fried egg uncertainty
		xxx=sbt[length(yr)]/mcmc$bmsy
		yyy=ft[length(yr)]/mcmc$fmsy
		fried.egg(xxx, yyy)
		lines(xx, yy[ii], type="l")
		text(xx, yy[ii], yr, cex=0.75)
	})
	fn=paste(prefix,"Fig7.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
}

fig8 <- function(A)
{	#DFO FMF plot (Fisheries Management Framework.)
	## This routine plots the spawning stock biomass/bmsy
	## versus the fishing mortality/fmsy using the 
	## DFO FMF as a plotting background.  Default
	## USR and LRP are 0.8Bmsy and 0.4 Bmsy
	op=par(no.readonly=T)
	par(las=1, mfcol=c(1, 1), mar=c(5, 5, 5, 4))
	with(A, {
		xx = sbt[1:length(yr)]/bmsy
		yy = ft/fmsy  #yy can be a matrix
		yy[yy==0]=NA; ii=!is.na(yy)
		
		yl=c(0,1.2*max(2,yy[ii]))
		matplot(xx, (yy[ii]), type="l", xlim=c(0,max(2,xx)),bty="l", 
		ylim=c(0,1.2*max(2,yy[ii])),xlab="Stock status", 
		ylab="Removal rate")
		rect(0, 0, 0.4, yl[2],col="salmon", border=NA)
		rect(0.4, 0, 0.8, yl[2],col="lightyellow", border=NA)
		rect(0.8, 0, max(2, xx), yl[2],col="lightgreen", border=NA)
		segments(0.4, 0, 0.4, yl[2], lwd=0.5)
		segments(0.8, 0, 0.8, yl[2], lwd=0.5)
		segments(0.8, 1.0, max(2, xx), 1.0)
		segments(0.4, 0,0.8, 1.0, lty=2)
		segments(0, 0,0.4,0, lty=2)
		
		text(0.4, 1.5, "Limit Reference\n Point", srt=90, cex=0.75)
		text(0.8, 1.5, "Upper Stock\n Reference", srt=90, cex=0.75)
		text(0.2, yl[2], "Critical\nZone", cex=0.75, pos=1)
		text(0.6, yl[2], "Cautious\nZone", cex=0.75, pos=1)
		text((max(2, xx)+0.8)/2, yl[2], "Healthy\nZone", cex=0.75, 
		pos=1)
		
		matlines(xx, yy[ii])
		ix=seq(1, length(xx), by=2)
		text(xx[ix], yy[ii][ix], yr[ix], cex=0.5)
	})
	fn=paste(prefix,"Fig8.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
	par(op)
}

fig9 <- function(A)
{	#3D Plots for time-varying selectivity curves
	## This routine uses a 3D perspective plot to 
	## show time-varying changes in age-specific 
	## selectivity.
	op = par(no.readonly=T)
	par(mgp=c(3, 3, 5))
	plot.sel<-function(x, y, z)
	{
		#z=exp(A$log_sel)*3
		x=A$yr
		y=A$age
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
			
		#wireframe(z, drap=TRUE, col=fill)
	}
	ix=1:length(A$yr)
	for(k in 1:A$ngear){
		plot.sel(A$yr, A$age, exp(A$log_sel[A$log_sel[,1]==k,-1]))
		file.name=paste(prefix, "Fig9",letters[k],".eps", sep="")
		if(savefigs) dev.copy2eps(file=file.name, height=8, width=8)
	}
	par(op)
	
}

fig10 <- function(A)
{	#Stock recruitment plots & St versus log(Rt/St)
	#Note the bias correction on Ro and rrt
	with(A, {
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
		plot(sbt[1:length(rt)],rt,xlim=c(0,max(sbt, st)),
			ylim=c(0,max(rt, rrt)),xlab="Spawning biomass (t)", 
			ylab="Age-1 recruits")
		lines(st, rrt, type="l")
		ro=ro*exp(-0.5*tau^2)
		points(bo, ro, pch="O", col=2)
		points(bo, ro, pch="+", col=2)
		fn=paste(prefix,"Fig10a.eps", sep="")
		if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
		
		plot(sbt[1:length(rt)], log(rt/sbt[1:length(rt)]), 
		    xlab="Spawning biomass (St)", ylab="Log(Rt/St)", 
			xlim=c(0,max(sbt, st)))
		lines(st, log(rrt/st))
		points(bo, log(ro/bo), pch="O", col=2)
		points(bo, log(ro/bo), pch="+", col=2)
		fn=paste(prefix,"Fig10b.eps", sep="")
		if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
	})
	
}

fig11 <- function(A)
{	#Residuals in abundance indices, catch, & surveys
	op=par(no.readonly=T)
	par(las=1, mfcol=c(2, 2))
	with(A, {
		#Catch residuals
		nyr=dim(ct)[2]
		matplot(yr, t(log(obs_ct[,1:nyr])-log(ct)),xlab= "Year", 
			ylab="Catch residuals", type="h")
		
		#Recruitment residuals
		plot(yr[-1], delta, type="h", xlab="Year", 
			ylab="Recruitment residual (delta)")
		
		#Relative abundance residuals
		ii=iyr%in%yr	#match indices in case of retrospective
		q=exp(mean(log(it[ii])-log(pit[ii])))
		plot(iyr[ii], epsilon[ii], xlab="Year", 
		ylab="Survey residual", type="h")
		
		#Observed vs. Predicted relative abundance
		xl=range(log(c(it[ii],pit[ii])))
		plot(log(it[ii]), log(pit[ii])+log(q),  
			xlab="Observed log(It)", ylab="Predicted log(It)", 
			xlim=xl, ylim=xl)
		abline(lm(log(pit[ii])+log(q)~log(it[ii])+0))
		abline(0, 1, col=2, lty=2)
		text(mean(xl),mean(xl), "1:1", srt=45, col=2, pos=3)
	})
	fn=paste(prefix,"Fig11.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
	par(op)
}

fig12 <- function(A)
{	#Spawning biomass depletion relative to reference points
	op=par(no.readonly=T)
	par(mfcol=c(1, 1))
	with(A, {
		dt.mci=t(apply(mc.sbt/mcmc$bo,2,quantile,probs=c(0.025,0.5,0.975)))
		bo.mci=quantile(mcmc$bo, probs=c(0.025,0.5,0.975))
		bmsy.mci=quantile(mcmc$bmsy, probs=c(0.025,0.5,0.975))
		plot(yrs, sbt/bo, type="n", xlab="Year", 
			ylab="Spawning biomass depletion", 
			ylim=c(0, max(sbt/bo, dt.mci, 1)))
		
		lrp=0.4*bmsy.mci[2]/bo.mci[2]
		usr=0.8*bmsy.mci[2]/bo.mci[2]
		
		rect(yrs[1],0,max(yrs), lrp, col="salmon", border=NA)
		rect(yrs[1],lrp,max(yrs), usr, col="lightyellow", border=NA)
		rect(yrs[1],usr,max(yrs), max(1,sbt/bo,dt.mci), col="lightgreen", border=NA)
		
		matlines(yrs, dt.mci, lty=c(2, 1, 2), col=1)
		text(median(yrs), lrp,"Limit reference point", cex=0.75)
		text(median(yrs), usr,"Upper stock reference", cex=0.75)
		text(yrs[1],mean(c(0, lrp)),"Critical Zone",pos=4)
		text(yrs[1],mean(c(lrp, usr)),"Cautious Zone", pos=4)
		text(yrs[1],mean(c(usr,max(1,sbt/bo,dt.mci))),"Healthy Zone", pos=4)
	})
	fn=paste(prefix,"Fig12.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
	par(op)
}


fig13 <- function(A)
{#Bubble plots for catch-at-age & residuals
	op=par(no.readonly=T)
	par(las=1, mfcol=c(2, 1), mar=c(0.1, 4, 0.1, 2), oma=c(5, 1, 1, 2))
	with(A, {
		igear=unique(A[, 2])
		for(i in igear)
		{
			O=subset(A, A[, 2]==i)
			P=subset(Ahat, Ahat[, 2]==i)
			ayr=O[,1]
			ii=ayr%in%yr
			bubble.plot(O[ii, 1],a_sage[i]:a_nage[i],O[ii, -(1:2)],
				scale=4, sample.size=FALSE, xlab="", xaxt="n")
			
			#Pearson residuals for multivariate logistic
			#log(o)-log(p)-mean(log(o)-log(p))
			o=0.001/(dim(O)[1]-2)
			r=log(O[ii,-(1:2)]+o)-log(P[,-(1:2)]+o)
			r=r-rowMeans(r)
			bubble.plot(O[ii, 1],a_sage[i]:a_nage[i],r,scale=4, 
				sample.size=FALSE, ylab="Pearson residuals", xlab="", 
				leg=TRUE)
			mtext("Year", 1, outer=T, line=3)
			
			file.name=paste(prefix, "Fig13",letters[i],".eps", sep="")
			if(savefigs) dev.copy2eps(file=file.name, height=8, width=8) 
		}
	})
	par(op)
}

fig14 <- function(A)
{	#Estimates of sage recruits
	with(A, {
		plot(yr[-1], rt, type="h", lwd=2, xlab="Year", 
			ylab=paste("Age", min(age), "recruits"), ylim=c(0, max(rt)))
		abline(h=median(rt), lty=2)
		abline(h=mean(rt), lty=3)
		legend("topright",c("Median","Average"),lty=c(2,3), bty="n")
	})
	fn=paste(prefix,"Fig14.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
}

fig15 <- function(A)
{	#Observed landings 
	with(A, {
		barplot(obs_ct,names.arg=yr, xlab="Year", ylab="Landings (t)")
	})
	fn=paste(prefix,"Fig15.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
}

fig16 <- function(A)
{	#Spawning biomass & depletion with credible intervals
	#Calls table 1 and jointly writes iscam.table1.tex
	table1(A)
	fn=paste(prefix,"Fig16.eps", sep="")
	if(savefigs) dev.copy2eps(file=fn, height=8, width=8)
}

fig17 <- function(A)
{	#Retrospective plots for spawning biomass
	rfiles = list.files(pattern="iscam.ret")
	sbt=matrix(NA,nrow=length(A$yrs),ncol=length(rfiles)+1)
	sbt[, 1]=A$sbt
	i=1
	for(fn in rfiles)
	{
		i=i+1
		R=read.rep(fn)
		sbt[1:length(R$yrs),i]=R$sbt
	}
	matplot(A$yrs,sbt[,1:6], type="l", lty=1, 
		ylim=c(0, max(sbt, na.rm=T)), xlab="Year", 
		ylab="Spawning biomass (t)")
}

## ________________________________________________________________ ##
## ________________________________________________________________ ##
## ________________________________________________________________ ##
## Tables

table1 <- function(A)
{
	op=par(no.readonly=T)
	par(mfcol=c(2, 1), las=1, mar=c(5, 4, 1, 1))
	with(A, {
		sb.mci=t(apply(mc.sbt,2,quantile,probs=c(0.025,0.5,0.975)))
		dt.mci=t(apply(mc.sbt/mcmc$bo,2,quantile,probs=c(0.025,0.5,0.975)))
		bo.mci=quantile(mcmc$bmsy, probs=c(0.025,0.5,0.975))
		
		t1 = cbind(yrs, sb.mci, dt.mci)
		t1 = apply(t1, 2, round, 3)
		
		
		matplot(t1[,1], t1[,2:4], type="l", lty=c(2, 1, 2),col=1, 
				xlab="Year", ylab="Spawning biomass (t)", 
				ylim=c(0, max(t1[,2:4])))
		#abline(h=bo.mci)
		matplot(t1[,1], t1[,5:7], type="l", lty=c(2, 1, 2),col=1, 
				xlab="Year", ylab="Spawning depletion (Bt/Bo)", 
				ylim=c(0, max(t1[,5:7])))
		
		filename<-paste(prefix,"iscam.Table1.tex",sep="")
			cap		<-"Recent trends in median estimate and 2.5\\% and 97.5\\% 
					credible intervals for spawning stock biomass, and
					spawning stock depletion. These estimates are based 
					on sampling the joint posterior distribution using MCMC."
			cgrp		<-c(" ","Spawning stock biomass","Depletion")
			ncgrp	<-c(1,3,3)
			colnames(t1)	<-c("Year",rep(c("2.5\\%","median","97.5\\%"),2))
		if(savefigs)
			tt <- latex(tail(t1, 10),file=filename,rowname=NULL,caption=cap,
				cgroup=cgrp,n.cgroup=ncgrp,label="iscam.T1")
	})
	par(op)
}

table2 <- function(A)
{	#MLE and Median values of the 6 leading parameters
	with(A, {
		mle = fit$est[1:6]
		std = fit$std[1:6]
		p=colnames(mcmc[,1:6])
		mci = t(apply(mcmc[,1:6],2, quantile, probs=c(0.5, 0.025, 0.975)))
		t1 = cbind(mle, std, mci)
		t1 = apply(t1,2, round, 3)
		
		cap <- "Maximum likelihood estimates (MLE) and standard deviations (SD)
				based on the inverse Hessian for the six leading parameters. Median
				values and the 95\\% credible interval based on posterior samples."
		colnames(t1) <- c("MLE","SD","Median", "2.5\\%","97.5\\%")
		rownames(t1) <- c("$\\ln(R_0)$","$h$","$\\ln(M)$","$\\ln(\\bar{R})$","$\\rho$","$\\vartheta$")
		filename<-paste(tprefix,"iscam.Table2.tex",sep="")
		if(savefigs)
			tt <- latex(t1,"",file=filename, caption=cap, label="iscam.T2")
		print(t1)
	})
}

