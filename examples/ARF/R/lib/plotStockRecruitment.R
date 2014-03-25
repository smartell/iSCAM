# Steven Martell
# Sept 2,  2012

.plotStockRecruit	<- function( M )
{
	n <- length(M)
	cat(".plotStockRecruit\n")
	
	mdf <- NULL
	for( i in 1:n )
	{
		SBt <- M[[i]]$sbt[1:(length(M[[i]]$yr)-M[[i]]$sage)]
		Rt <- M[[i]]$rep_rt[(1+M[[i]]$sage):length(M[[i]]$yr)]
		Year <- M[[i]]$yr[(1+M[[i]]$sage):(length(M[[i]]$yr))]
		kappa <- M[[i]]$kappa
		ro <- M[[i]]$ro
		sbo <-  M[[i]]$sbo
		tau <-  M[[i]]$tau
		rectype <- M[[i]]$rectype
		ssbt <- seq(0,sbo,length.out=length(SBt))
		
		if(M[[i]]$rectype==1) 
		{
			#Beverton-Holt
			Rt_hat=kappa*ro*ssbt/(sbo+(kappa-1)*ssbt)*exp(-0.5*tau^2)  
		}
		else if(M[[i]]$rectype==2)
		{
			#Ricker
			Rt_hat=kappa*ro*ssbt*exp(-log(kappa)*ssbt/sbo)/sbo *exp(-0.5*tau^2) 
		}


		df <- data.frame(Model=names(M[i]),Year=Year,SBt=SBt,Rt=Rt,Rt_hat=Rt_hat,ssbt=ssbt)
		
		mdf <- rbind(mdf,df)
	}
	print(head(mdf,3))

	p <- ggplot(mdf) + geom_point(aes(SBt,Rt,color=Year))
	p <- p + geom_line(aes(ssbt, Rt_hat),size=2)
	p <- p + labs(x="Spawning Stock Biomass",y="Recruitment")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}

#.plotOldStockRecruit	<- function( repObj )
#{
#	with(repObj, {
#		xx = sbt[1:(length(yr)-min(age))]
#		yy = rt
#		tt = yr[1:(length(yr)-min(age))]
#		
#		plot(xx, yy, type="n",ylim=c(0, max(yy, ro)),xlim=c(0, max(xx,sbo,bo)), 
#			xlab="Spawning biomass (t)", 
#			ylab=paste("Age-",min(age)," recruits\n", "(numbers)", sep=""), 
#			main=paste(stock))
#		grid()
#			
#		lines(xx, yy, type="l", col=colr(1, 0.7))
#		text(xx, yy, tt, cex=0.65)
#		#points(xx, yy)
#		
#		#points(xx[1],yy[1], pch=20, col="green")
#		#points(xx[length(xx)], yy[length(xx)], pch=20, col=2)
#		
#		
#		qn <- quantile(yy, probs=c(0.33, 0.66))
#		abline(h=quantile(yy, probs=c(0.33, 0.66)), col="grey")
#		text(0.9*max(xx,sbo,bo),qn, c("33 precentile", "66 percentile"), col="grey", pos=3, cex=0.75)
#		
#		plot.curve <- function(sbo, clr, ...)
#		{
#			st=seq(0, max(sbt, sbo), length=100)
#			if(rectype==1)
#			{
#				#Beverton-Holt
#				rrt=kappa*ro*st/(sbo+(kappa-1)*st)*exp(-0.5*tau^2)  
#			}
#			if(rectype==2)
#			{
#				#Ricker
#				rrt=kappa*ro*st*exp(-log(kappa)*st/sbo)/sbo *exp(-0.5*tau^2) 
#			}
#			lines(st, rrt, col=clr, ...)
#			ro=ro*exp(-0.5*tau^2)
#			points(sbo, ro, pch="O", col=clr, cex=1.25)
#			points(sbo, ro, pch="+", col=clr, cex=1.25)
#		}
#		plot.curve(sbo, clr=colr(1, 0.5), lwd=4);
#		plot.curve(bo, clr=colr(4, 0.5), lwd=2);
#		
#		
#		h <- round( kappa/(4+kappa), 2 )
#		legend("topleft", paste("Steepness =",h, sep=""), bty="n")
#		legend("topright", c("Historical S-R curve", "Projection S-R curve"), 
#		lty=1, col=c(colr(1, 0.5), colr(4, 0.5)), bty="n", lwd=c(4, 2))
#	})
#}

.plotRecruitsPerSpawner	<- function( M )
{
	require(ggplot2)

	n <- length(M)
	cat(".plotStockRecruit\n")
	
	mdf <- NULL
	
	for( i in 1:n )
	{
	SBt <- M[[i]]$sbt[1:(length(M[[i]]$yr)-M[[i]]$sage)]
	Rt <- M[[i]]$rep_rt[(1+M[[i]]$sage):length(M[[i]]$yr)]
	Year <- M[[i]]$yr[(1+M[[i]]$sage):(length(M[[i]]$yr))]
	kappa <- M[[i]]$kappa
	ro <- M[[i]]$ro
	sbo <-  M[[i]]$sbo
	bo <-  M[[i]]$bo
	tau <-  M[[i]]$tau
	rectype <- M[[i]]$rectype
	st <- seq(10,sbo,length.out=length(SBt))

	recruitment.model <-function(sbo, ...)
		{
			if(rectype==1)
			{
				rm  = kappa*ro*st/(sbo+(kappa-1)*st)*exp(-0.5*tau^2) 
			}
			if(rectype==2)
			{
				rm = kappa*ro*st*exp(-log(kappa)*st/sbo)/sbo*exp(-0.5*tau^2) 
			}
			return(rm)
		}
		
		rm = recruitment.model(sbo)
		rp = recruitment.model(bo)
		
		df  = data.frame(Model=names(M[i]),SSB=SBt,RS=log(Rt/SBt), R=Rt,rm=rm, rp=rp, st=st,  Year=Year)
		mdf <- rbind(mdf,df)
		}
		
		print(head(mdf,3))

		p <- ggplot(mdf) + geom_point(aes(x=SSB, y=RS, color=Year))
		p <- p + geom_line(aes(x=st, y=log(rm/st)), size=2)
		p <- p + labs(x="Spawning stock biomass (SSB)", y="Log(recruits/SSB)")
		p <- p + facet_wrap(~Model,scales="free")
		print(p + .THEME)
	}



.plotOldRecruitsPerSpawner	<- function( repObj )
{
	require(ggplot2)
	# Create data frame object from repObj
	with(repObj, {
		SSB  = sbt[1:(length(yr)-min(age))]
		R    = rt
		Year = yr[1:(length(yr)-min(age))]
		
		st   = seq(1e-5, max(sbt, sbo), length=length(SSB))
		recruitment.model <-function(sbo, ...)
		{
			
			if(rectype==1)
			{
				rm  = kappa*ro*st/(sbo+(kappa-1)*st)*exp(-0.5*tau^2) 
			}
			if(rectype==2)
			{
				rm = kappa*ro*st*exp(-log(kappa)*st/sbo)/sbo*exp(-0.5*tau^2) 
			}
			return(rm)
		}
		
		rm = recruitment.model(sbo)
		rp = recruitment.model(bo)
		obj  = data.frame(SSB=SSB,RS=log(R/SSB), R=R,rm=rm, rp=rp, st=st,  Year=Year)
		
		p <- ggplot(obj) + geom_point(aes(x=SSB, y=RS, col=Year))
		p <- p + geom_line(aes(x=st, y=log(rp/st)))
		p <- p + geom_line(aes(x=st, y=log(rm/st))) + labs(x="Spawning stock biomass (SSB)", y="Log(recruits/SSB)")
		p <- p + geom_point(x=sbo, y=log(ro*exp(-0.5*tau^2) /sbo), col=I("red"), shape=2)
		p <- p + geom_point(x=sbo, y=log(ro*exp(-0.5*tau^2) /bo), col=I("blue"), shape=3)
		print(p)
	})
}

# .plotoldStockRecruit	<- function( repObj )
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
# 		grid()
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
