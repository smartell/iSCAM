# Steven Martell
# September 2,  2012

.plotDepletion <- function( M )
{
	n <- length( M )
	cat(".plotDepletion")
	mdf <- NULL
	for(i in 1:n)
	{
		df <- data.frame(Model=names(M)[i],Year=M[[i]]$yrs,Dt=M[[i]]$sbt/M[[i]]$sbo)
		mdf <- rbind(mdf,df)
	}

	p <- ggplot(mdf,aes(x=Year,y=Dt)) + geom_line()
	p <- p + labs(x="Year",y="Spawning biomass depletion")
	p <- p + ylim(0,max(mdf$Dt))
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}


.plotOldDepletion	<- function( repObj, annotate=FALSE )
{
	#plot the spawning biomass depletion level & reference points
	with(repObj, {
		xx=yr
		yy=sbt[1:length(xx)]/sbo
		yrange=c(0,1.1*max(yy, na.rm=TRUE))
		
		plot(xx, yy, type="n", axes=FALSE,
			xlab="Year", ylab="Spawning depletion",main=paste("stock"), 
			ylim=yrange)
		lines(xx, yy)
		rlvl=c(1.0, 0.8, 0.4)
		abline(h=rlvl*bmsy/sbo,lty=2,lwd=0.5)
		
		
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		grid()
		
		if ( annotate )
		{
			#mfg <- par( "mfg" )
			#if ( mfg[1]==1 && mfg[2]==1 )
			#legend( "top",legend=c( "Spawning biomass","MSY depletion level",
			#	"Upper stock reference","Limit reference point"),
			#	bty='n',lty=c(1,2,2,2),lwd=c(1,rlvl),pch=-1,ncol=2 )
			
			#Delinate critical zone,  cautious zone, healthy zone.
			rect(min(xx)-5,-0.5,max(xx)+5,0.4*bmsy/bo,col=colr("red",0.1), border=NA)
			rect(min(xx)-5,0.4*bmsy/bo,max(xx)+5,0.8*bmsy/bo,col=colr("yellow",0.1),border=NA)
			rect(min(xx)-5,0.8*bmsy/bo,max(xx)+5,1.5,col=colr("green",0.1), border=NA)
		}
	})
}
