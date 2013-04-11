# Rscript for plotting relative abundance data
# Steven Martell
# August 28,  2012

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
		n=nrow(t(as.matrix(yy)))
		yrange=c(0, max(yy, na.rm=TRUE))
		
		matplot(xx, yy, type="n", axes=FALSE,
			xlab="Year", ylab="Relative abundance", 
			ylim=yrange , main=paste(stock))
		
		matlines(xx, yy, col="black",type="o", pch=1:n)
		
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		grid()
		
		if( annotate )
		{
			txt=paste("Survey",1:n)
			legend("top", txt, lty=1:n, pch=1:n, bty="n")
		}
	})
}
