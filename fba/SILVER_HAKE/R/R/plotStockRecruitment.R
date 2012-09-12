# Steven Martell
# Sept 2,  2012

.plotStockRecruit	<- function( repObj )
{
	with(repObj, {
		xx = sbt[1:(length(yr)-min(age))]
		yy = rt
		tt = yr[1:(length(yr)-min(age))]
		
		plot(xx, yy, type="n",ylim=c(0, max(yy, ro)),xlim=c(0, max(xx,sbo,bo)), 
			xlab="Spawning biomass (t)", 
			ylab=paste("Age-",min(age)," recruits\n", "(numbers)", sep=""), 
			main=paste(stock))
		grid()
			
		lines(xx, yy, type="l", col=colr(1, 0.7))
		text(xx, yy, tt, cex=0.65)
		#points(xx, yy)
		
		#points(xx[1],yy[1], pch=20, col="green")
		#points(xx[length(xx)], yy[length(xx)], pch=20, col=2)
		
		
		qn <- quantile(yy, probs=c(0.33, 0.66))
		abline(h=quantile(yy, probs=c(0.33, 0.66)), col="grey")
		text(0.9*max(xx,sbo,bo),qn, c("33 precentile", "66 percentile"), col="grey", pos=3, cex=0.75)
		
		plot.curve <- function(sbo, clr, ...)
		{
			st=seq(0, max(sbt, sbo), length=100)
			if(rectype==1)
			{
				#Beverton-Holt
				rrt=kappa*ro*st/(sbo+(kappa-1)*st)*exp(-0.5*tau^2)  
			}
			if(rectype==2)
			{
				#Ricker
				rrt=kappa*ro*st*exp(-log(kappa)*st/sbo)/sbo *exp(-0.5*tau^2) 
			}
			lines(st, rrt, col=clr, ...)
			ro=ro*exp(-0.5*tau^2)
			points(sbo, ro, pch="O", col=clr, cex=1.25)
			points(sbo, ro, pch="+", col=clr, cex=1.25)
		}
		plot.curve(sbo, clr=colr(1, 0.5), lwd=4);
		plot.curve(bo, clr=colr(4, 0.5), lwd=2);
		
		
		h <- round( kappa/(4+kappa), 2 )
		legend("topleft", paste("Steepness =",h, sep=""), bty="n")
		legend("topright", c("Historical S-R curve", "Projection S-R curve"), 
		lty=1, col=c(colr(1, 0.5), colr(4, 0.5)), bty="n", lwd=c(4, 2))
	})
}


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
