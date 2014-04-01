# Steven Martell
# Aug 31,  2012

.plotSpawnBiomass <- function( M )
{
	n <- length(M)
	cat(".plotSpawnBiomass\n")

	mdf <- NULL
	for(i in 1:n)
	{
		fit = M[[i]]$fit
		sbt <- fit$est[fit$names=="sd_sbt"]
		std <- fit$std[fit$names=="sd_sbt"]
		bt <- data.frame(Model=names(M)[i],Year=M[[i]]$yrs,SBt=M[[i]]$sbt,Std=std)
		bt <- data.frame(bt,Bo=M[[i]]$bo)
		mdf <- rbind(mdf,bt)
	}

	p <- ggplot(mdf,aes(Year,SBt)) + geom_line(width=2)
	p <- p + geom_ribbon(aes(ymax=SBt+1.96*Std,ymin=SBt-1.96*Std),alpha=0.2)
	p <- p + geom_line(data=bt,aes(Year,Bo),col="blue")
	p <- p + labs(x="Year",y=paste("Spawning biomass",.UNITS))
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}


# .plotBiomass	<- function( repObj, annotate=FALSE )
# {
# 	#plot total biomass & spawning biomass 
# 	with(repObj, {
# 		xx=yr
# 		yy=cbind(bt[1:length(xx)], sbt[1:length(xx)])
		
# 		yrange=c(0, 1.2*max(yy, na.rm=TRUE))
		
# 		matplot(xx, yy, type="n",axes=FALSE,
# 				xlab="Year", ylab="Biomass (1000 t)",main=paste(stock), 
# 				ylim=yrange)
		
# 		matlines(xx,yy,
# 			type="l", col="black",
# 			ylim=c(0,max(yy,na.rm=T)))
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		grid()
		
# 		if ( annotate )
# 		{
# 			mfg <- par( "mfg" )
# 			if ( mfg[1]==1 && mfg[2]==1 )
# 			legend( "top",legend=c( "Pre-fishery biomass","Spawning biomass"),
# 				bty='n',lty=c(1,2),lwd=c(1,1),pch=c(-1,-1),ncol=1 )
# 		}
# 	})	
# }
