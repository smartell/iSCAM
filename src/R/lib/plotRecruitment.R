# Steven Martell
# Sept 2,  2011

.plotRecruitment	<- function( repObj )
{
	#plot age-a recruits.
	with(repObj, {
		xx = yr
		yy = exp(ln_rt)
		yy=yy
		yrange=c(0, max(yy, na.rm=T))
		
		plot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
			xlab="Year", main=paste(stock), 
			ylab=paste("Age-", min(age), " recruits", sep=""))
		grid()
		
		lines(xx, yy, type="h", lwd=3, col=colr(1, 0.7))
		
		
		#add 0.33 and 0.66 quantile lines
		qtl = quantile(yy, prob=c(1/3, 2/3))
		abline(h=qtl,lty=1, col=colr("darkgrey",0.75))
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		
	})
}