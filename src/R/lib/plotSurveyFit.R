# Steven Martell
# Sept 6,  2012

.plotSurveyfit	<- function( repObj, annotate=FALSE)
{
	with(repObj, {
		if(is.matrix(it)){
			xx = t(iyr)
			m = apply(it,1,max, na.rm=T)
			yy = t(pit/m)
			y2 = t(it/m)
		}else{
			xx = iyr
			yy = pit
			y2 = it
		}
		n=nrow(t(as.matrix(yy)))
		#n=dim(xx)[2]
		yrange=c(0, 1.15*max(yy, y2, na.rm=TRUE))
		
		matplot(xx, yy, type="n",axes=FALSE,ylim=yrange, 
			xlab="Year", ylab="Relative abundance", main=paste(stock))
		
		matlines(xx, yy, col=1:n, lty=1)
		matpoints(xx, y2, col=1:n, pch=1:n)
		
		axis( side=1 )
		axis( side=2, las=.VIEWLAS )
		box()
		grid()
		
		if ( annotate )
		{
			
			txt=paste("Survey ",1:n,", q=",round(q, 3), sep="")
			
			mfg <- par( "mfg" )
			#if ( mfg[1]==1 && mfg[2]==1 )
			legend( "top",legend=txt,
				bty='n',lty=1,lwd=1,pch=1:n,ncol=n, col=1:n )
				
			#print(q)
		}
	})
}



.plotSurveyfit <- function(M)
{
	n <- length(M)
	cat(".plotSurveyfit\n")
	mdf <- NULL
	for( i in 1:n )
	{
		it <- data.frame(cbind(M[[i]]$d3_survey_data,M[[i]]$it_hat))
		colnames(it) <- c("Year","Index","Gear","Area","Group","Sex","wt","timing","Index_hat")
		it <- data.frame(Model=names(M)[i],it)
		mdf <- rbind(mdf,it)
	}
	print(head(mdf,3))

	p <- ggplot(mdf) + geom_point(aes(Year,Index,colour = "Observed", name=" "))
	p <- p + geom_line(aes(Year,Index_hat,colour="Predicted"))
	p <- p + theme_set(theme_bw(12))
	p <- p + theme(legend.title=element_blank())
	p <- p + labs(x="Year",y="Relative abundance",linetype="Gear")
	p <- p + facet_wrap(~Model,scales="free")
	print(p)
	#print(p + .THEME)
}









