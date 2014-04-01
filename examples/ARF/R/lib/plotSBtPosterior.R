# Steven Martell

.plotSbtPosterior	<- function( admbObj, depletion=FALSE, annotate=FALSE )
{
	## Median and 95% CI for spawning biomass or depletion
	print("	.plotSbtPosterior")
	with( admbObj, {
		xx <- yr
		yy <- t(apply( mcsbt, 2, quantile, probs=c(0.5, 0.025, 0.975) ))
		
		bo <- quantile( mcmc$bo, probs=c(0.5, 0.025, 0.975) )
		bmsy <- quantile( mcmc$bmsy, probs=c(0.5, 0.025, 0.975) )
		yl="Spawning biomass (1000 t)"
		
		if(depletion)
		{
			yy <- t(apply( mcsbt/mcmc$bo, 2, quantile, probs=c(0.5, 0.025, 0.975) ))
			yl = "Spawning biomass depletion"
		}
			

		matplot(xx, yy, type="n", xlab="Year", 
			ylab=yl, axes=FALSE, 
			ylim=c(0, 1.0*max(yy)), main=paste(stock) )
		
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		
		polygon(c(xx, rev(xx)), 
			c(yy[,2],rev(yy[,3])), 
			col=colr("red",0.25), border=NA)
		print(yy)
		lines(xx, yy[,1], lwd=1.5)
		#matlines(xx, yy)
		
		if(!depletion)
		{
			points(min(xx)-0.4,bo[1])
			arrows(min(xx)-0.4,bo[2], min(xx)-0.4, bo[3], length=0)
		}
		
		if ( annotate && depletion )
		{
			#mfg <- par( "mfg" )
			#if ( mfg[1]==1 && mfg[2]==1 )
			#legend( "top",legend=c( "Spawning biomass","MSY depletion level",
			#	"Upper stock reference","Limit reference point"),
			#	bty='n',lty=c(1,2,2,2),lwd=c(1,rlvl),pch=-1,ncol=2 )
			
			#Delinate critical zone,  cautious zone, healthy zone.
			rect(min(xx)-5,-0.5,max(xx)+5,0.4*bmsy[1]/bo[1],col=colr("red",0.1), border=NA)
			rect(min(xx)-5,0.4*bmsy[1]/bo[1],max(xx)+5,0.8*bmsy[1]/bo[1],col=colr("yellow",0.1),border=NA)
			rect(min(xx)-5,0.8*bmsy[1]/bo[1],max(xx)+5,1.5,col=colr("green",0.1), border=NA)
		}
		
		
	})
}
