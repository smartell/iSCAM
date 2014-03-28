# Steven Martell
# Aug 28,  2012

.plotSurveyResiduals	<- function( repObj, annotate=FALSE )
{
	#Plot residuals between observed and predicted relative abundance
	#indicies (epsilon)
	with(repObj, {
		if(is.matrix(epsilon)){
			xx = t(iyr)
			yy = t(epsilon)
			t1 = colSums(yy,na.rm=T)
			ng = length(t1[t1!=0])
		}else{
			xx = iyr
			yy = epsilon
			ng = 1
		}
		
		absmax = abs(max(yy, na.rm=TRUE))
		if(absmax< 1e-3) absmax=1
		yrange=c(-absmax, absmax)
		
		matplot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
			xlab="Year", ylab="Survey residual", main=paste(stock))
		
		matlines(xx, yy, type="h", col="black")
		matpoints(xx, yy,pch=1:ng, cex=0.75, col=1)
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		grid()
		if ( annotate )
		{
			n=dim(xx)[2]
			txt=paste("Survey",1:ng)
			
			mfg <- par( "mfg" )
			if ( mfg[1]==1 && mfg[2]==1 )
			legend( "top",legend=txt,
				bty='n',lty=1:ng,lwd=1,pch=-1,ncol=1 )
		}
	})
}
