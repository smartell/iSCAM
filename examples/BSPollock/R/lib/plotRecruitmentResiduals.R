# Steven Martell
# Aug 30,  2012

# R-script for a barplot for the observed catch.
# Steven Martell
# Sep 5, 2013

require(ggplot2)
require(reshape2)

.plotRecruitmentResidual <- function(M)
{
	n <- length(M)

	mdf <- NULL
	for( i in 1:n )
	{
		D  <- M[[i]]
		rt <- data.frame( t(rbind(D$yr[-(1:D$sage)],D$rt,D$delta)) )
		# ng <- (dim(rt)[2]-1)/2 # special case when ngroups > 1
		colnames(rt) <- c("Year","Recruitment","Residual")
		rt <- data.frame(Model=names(M)[i],rt)
		mdf <- rbind(mdf,rt)

		# ct <- data.frame(M[[i]]$catch_data)
		# colnames(ct) <- c("Year","Gear","Area","Group","Sex","Type","Catch")
		# ct <- data.frame(Model=names(M)[i],ct)
		# mdf <- rbind(mdf,ct)
	}
	print(head(mdf,3))
	print(tail(mdf,3))

	p <- ggplot(mdf,aes(x=factor(Year),Residual))
	p <- p + geom_bar(width=0.75,position="dodge",stat='identity')
	p <- p + labs(x="Year",y="Recruitment (log residual)")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}


# .plotRecruitmentResiduals	<- function( repObj )
# {
# 	#Plot the log residuals between the estimated recruits and
# 	#those obtained from the recruitment model (delta)
# 	with(repObj, {
# 		ii = 1:min(age)
# 		xx = yr[-ii]
# 		yy = delta
# 		absmax = abs(max(yy, na.rm=TRUE))
# 		if(absmax< 1e-3) absmax=1
# 		yrange=c(-absmax, absmax)
		
# 		plot(xx, yy, type="n",  axes=FALSE, ylim=yrange, 
# 			xlab="Year", ylab="Recruitment residuals", main=paste(stock))
		
# 		points(xx, yy, cex=0.75)	
# 		lines(xx, yy, type="h", col="black")
# 		axis( side=1 )
# 		axis( side=2,  las=.VIEWLAS )
# 		box()
# 		grid()
# 	})
# }
