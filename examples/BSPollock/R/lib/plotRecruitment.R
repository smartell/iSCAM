# Steven Martell
# Sept 2,  2011

# Feb 26, 2014
.plotRecruitment <- function( M )
{
	n <- length(M)
	cat(".plotRecruitment\n")

	mdf <- NULL
	for(i in 1:n)
	{
		rt <- data.frame(Model=names(M)[i],Year=M[[i]]$yr,Rt=M[[i]]$rep_rt)
		mdf <- rbind(mdf,rt)
	}

	p <- ggplot(mdf,aes(factor(Year),Rt))
	p <- p + geom_bar(width=0.75,stat='identity')
	p <- p + labs(x="Year",y=paste("Recruitment"))
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}

# .plotOldRecruitment	<- function( repObj )
# {
# 	#plot age-a recruits.
# 	with(repObj, {
# 		xx = yr
# 		yy = exp(ln_rt)
# 		yy=yy
# 		yrange=c(0, max(yy, na.rm=T))
		
# 		plot(xx, yy, type="n", axes=FALSE, ylim=yrange, 
# 			xlab="Year", main=paste(stock), 
# 			ylab=paste("Age-", min(age), " recruits", sep=""))
# 		grid()
		
# 		lines(xx, yy, type="h", lwd=3, col=colr(1, 0.7))
		
		
# 		#add 0.33 and 0.66 quantile lines
# 		qtl = quantile(yy, prob=c(1/3, 2/3))
# 		abline(h=qtl,lty=1, col=colr("darkgrey",0.75))
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
		
# 	})
# }