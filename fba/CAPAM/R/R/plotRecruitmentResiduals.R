# Steven Martell
# Aug 30,  2012

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
			xlab="Year", ylab="Recruitment residuals", main=paste(stock))
		
		points(xx, yy, cex=0.75)	
		lines(xx, yy, type="h", col="black")
		axis( side=1 )
		axis( side=2,  las=.VIEWLAS )
		box()
		grid()
	})
}
