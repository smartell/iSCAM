# Steven Martell
.plotCatchResiduals		<- function( repObj, annotate=FALSE )
{
	#Plot residuals between observed and predicted catches
	#residuals (epsilon=log(obs_ct)-log(ct))
	with(repObj, {
		epsilon=log(obs_ct)-log(ct)
		if(is.matrix(epsilon)){
			xx = yr
			yy = t(epsilon)
			t1 = colSums(yy,na.rm=T)
			ng = length(t1[t1!=0])
		}else{
			xx = yr
			yy = epsilon
			ng = 1
		}
		#browser()
		absmax = max(abs(yy), na.rm=TRUE)
		if(absmax<=1e-3)absmax=1
		yrange=c(-absmax, absmax)
		
		matplot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
			xlab="Year", ylab="Catch residual", main=paste(stock))
		
		matlines(xx, yy, type="h", col="black",lty=1)
		matpoints(xx, yy,pch=1:ng, cex=0.75, col=1)
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		grid()
		if ( annotate )
		{
			#n=dim(yy)[2]
			txt=paste("Gear",1:ng)
			
			mfg <- par( "mfg" )
			if ( mfg[1]==1 && mfg[2]==1 )
			legend( "top",legend=txt,cex=0.75, 
				bty='n',pch=1:ng,lwd=1,lty=-1,ncol=ng )
		}
	})
}