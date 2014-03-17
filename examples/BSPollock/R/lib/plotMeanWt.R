# Steven Martell

.plotMeanwt	<- function( repObj )
{
	#plot mean weight-at-age by cohort
	with(repObj, {
		xx = yr		## xaxis labels
		yy = age	## yaxis labels
		nage=length(age)
		
		if(sum(par("mfcol"))==2)
		{
			xl = "Cohort year";xlm=""
			yl = "Weight-at-age (kg)";ylm=""
		}
		else
		{
			xlm = "Cohort year";xl=""
			ylm = "Weight-at-age (kg)";yl=""
		}
		
		plot(range(xx), range(wt_obs), type="n", axes=FALSE,
		xlab=xl, ylab=yl, main=paste(stock))
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		grid()
		
		for(i in 1:(dim(wt_obs)[1]-1))
		{
			yy = (diag(as.matrix(wt_obs[0:-i, ]))) 
			xx = 1:length(yy)+yr[i]-min(age)+1
			
			yy[yy==0]=NA;xx[yy==NA]=NA
			lines(xx,yy)
			
			points(xx[1],yy[1],pch=20,col="steelblue",cex=0.5)
			points(xx[nage],yy[nage],pch=20,col="salmon",cex=0.5)
		}
		for(i in 1:dim(wt_obs)[2]-1)
		{
			yy = diag(as.matrix(wt_obs[,-1:-i]))
			n = length(yy)
			xx = yr[1]:(yr[1]+n-1)
			lines(xx, yy)
			points(xx[n], yy[n], pch=20, col="salmon", cex=0.5)
		}
		
		mtext(xlm, side=1, outer=T, line=0)
		mtext(ylm, side=2, outer=T, line=0)
	})
}