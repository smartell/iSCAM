# Steven Martell
# Aug 31,  2012

.plotBiomass	<- function( repObj, annotate=FALSE )
{
	#plot total biomass & spawning biomass 
	with(repObj, {
		xx=yr
		yy=cbind(bt[1:length(xx)], sbt[1:length(xx)])
		
		yrange=c(0, 1.2*max(yy, na.rm=TRUE))
		
		matplot(xx, yy, type="n",axes=FALSE,
				xlab="Year", ylab="Biomass (1000 t)",main=paste(stock), 
				ylim=yrange)
		
		matlines(xx,yy,
			type="l", col="black",
			ylim=c(0,max(yy,na.rm=T)))
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		grid()
		
		if ( annotate )
		{
			mfg <- par( "mfg" )
			if ( mfg[1]==1 && mfg[2]==1 )
			legend( "top",legend=c( "Pre-fishery biomass","Spawning biomass"),
				bty='n',lty=c(1,2),lwd=c(1,1),pch=c(-1,-1),ncol=1 )
		}
	})	
}
