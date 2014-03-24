# R-script for catch residuals
# Sept 5, 2013
# Steven Martell

require(ggplot2)
require(reshape2)

.plotCatchResidual <- function(M)
{
	n <- length(M)

	mdf <- NULL
	for( i in 1:n )
	{
		ct <- data.frame( cbind(M[[i]]$dCatchData,M[[i]]$eta) )
		colnames(ct) <- c("Year","Gear","Area","Group","Sex","Type","Catch","Residual")
		ct <- data.frame(Model=names(M)[i],ct)
		mdf <- rbind(mdf,ct)
	}
	print(head(mdf,3))

	p <- ggplot(mdf,aes(x=factor(Year),Residual,fill=factor(Gear)))
	p <- p + geom_bar(width=0.75,position="dodge",stat='identity')
	p <- p + labs(x="Year",y="log residual")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}


# .plotCatchResiduals		<- function( repObj, annotate=FALSE )
# {
# 	#Plot residuals between observed and predicted catches
# 	#residuals (epsilon=log(obs_ct)-log(ct))
# 	with(repObj, {
# 		epsilon=log(obs_ct)-log(ct)
# 		if(is.matrix(epsilon)){
# 			xx = yr
# 			yy = t(epsilon)
# 			t1 = colSums(yy,na.rm=T)
# 			ng = length(t1[t1!=0])
# 		}else{
# 			xx = yr
# 			yy = epsilon
# 			ng = 1
# 		}
# 		#browser()
# 		absmax = max(abs(yy), na.rm=TRUE)
# 		if(absmax<=1e-3)absmax=1
# 		yrange=c(-absmax, absmax)
		
# 		matplot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
# 			xlab="Year", ylab="Catch residual", main=paste(stock))
		
# 		matlines(xx, yy, type="h", col="black",lty=1)
# 		matpoints(xx, yy,pch=1:ng, cex=0.75, col=1)
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		grid()
# 		if ( annotate )
# 		{
# 			#n=dim(yy)[2]
# 			txt=paste("Gear",1:ng)
			
# 			mfg <- par( "mfg" )
# 			if ( mfg[1]==1 && mfg[2]==1 )
# 			legend( "top",legend=txt,cex=0.75, 
# 				bty='n',pch=1:ng,lwd=1,lty=-1,ncol=ng )
# 		}
# 	})
# }