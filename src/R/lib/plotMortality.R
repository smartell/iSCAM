# Steven Martell
# Sept 6,  2012

.plotMortality <- function( M )
{
	n <- length(M)
	cat(".plotMortality\n")

	mdf <- NULL
	for(i in 1:n)
	{
		nn  <- length(M[[i]]$yr)
		if(M[[i]]$nsex==2) 
			sex = c(rep('Female',nn),rep('Male',nn))
		else
			sex = c(rep('Unsexed',nn))

		zt  <- data.frame(Model=names(M)[i],Year=M[[i]]$yr,Mt=rowMeans(M[[i]]$M))
		zt  <- data.frame(zt,Ft=rowMeans(M[[i]]$F),Sex=sex)
		mdf <- rbind(mdf,melt(zt,id.vars=c("Model","Year","Sex")))
	}
	print(head(mdf,3))

	p <- ggplot(mdf,aes(Year,value,fill=variable)) 
	p <- p + geom_area(position='stack',alpha=0.7)
	p <- p + labs(x="Year",y="Mean instantaneous mortality rate")
	p <- p + scale_fill_brewer(palette="Blues")
	p <- p + facet_wrap(Model~Sex,scales="free")
	print(p + .THEME)
  
}


# .plotOldMortality	<- function( repObj, annotate=FALSE )
# {
# 	# SJDM June 5, 2011 Changed to plot Average M and sum of average Fs by gear
# 	#plot average total mortality,  fishing mortality & natural mortality
# 	with(repObj, {
# 		xx=yr
# 		if(is.matrix(ft))
# 		{
# 			yy=t(as.matrix(ft))
# 		}
# 		else
# 		{
# 			yy=as.matrix(ft)
# 		}
# 		#n=nrow(t(as.matrix(yy)))
# 		icol=apply(yy,2,function(x){sum(cumsum(x))!=0.0})
# 		ng=length(icol[icol==T])
		
# 		yy = cbind( rowMeans(M_tot), yy[,icol] )
		
# 		csyy = t(apply(yy,1, cumsum))	#cumulative sum
			
# 		yrange=c(0, max(csyy, na.rm=TRUE))
# 		lw = c(1, rep(1,ng))
# 		lt = 1
# 		iclr = colr(1:(ng+1),0.5)
		
# 		matplot(xx, csyy, type="n", axes=FALSE, log="y",  
# 			xlab="Year", ylab="Mortality rate", main=paste(stock))
			
		
# 		#lines(xx, yy[,1], lwd=2, lty=1, col=1)
# 		matlines(xx, csyy, log="y", col=iclr, lwd=lw, lty=lt)
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		grid()
		
# 		ry=cbind(1.e-30,csyy)
# 		for(i in 1:dim(csyy)[2])
# 		{
# 			x2=c(xx, rev(xx))
# 			y2=c(ry[,i+1], rev(ry[,i]))
# 			#browser()
# 			polygon(x2, y2, border=NA,col=colr(i,0.2), log="y")
# 		}
		
		
		
# 		if ( annotate )
# 		{
# 			txt = c("Natural mortality", paste("Gear",1:ng))
# 			#mfg <- par( "mfg" )
# 			#if ( mfg[1]==1 && mfg[2]==1 )
# 			legend( "topright",legend=txt,col=iclr, 
# 				bty='n',lty=lt,lwd=5,pch=-1,ncol=2)
# 		}
# 	})
# }
