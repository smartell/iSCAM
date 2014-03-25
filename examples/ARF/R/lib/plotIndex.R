# Rscript for plotting relative abundance data
# Steven Martell
# August 28,  2012

require(ggplot2)
require(reshape2)

.plotIndex <- function(M)
{
	n <- length(M)
	cat(".plotIndex\n")
	mdf <- NULL
	for( i in 1:n )
	{
		it <- data.frame(M[[i]]$d3_survey_data)
		colnames(it) <- c("Year","Index","Gear","Area","Group","Sex","wt","timing")
		it <- data.frame(Model=names(M)[i],it)
		mdf <- rbind(mdf,it)
	}
	print(head(mdf,3))

	p <- ggplot(mdf,aes(Year,Index,linetype=factor(Gear)))
	p <- p + geom_line()
	p <- p + labs(x="Year",y="Relative abundance",linetype="Gear")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}

.plotIndexResidual <- function(M)
{
	n <- length(M)
	cat(".plotIndexResidual\n")
	mdf <- NULL
	for( i in 1:n )
	{
		it <- data.frame( cbind(M[[i]]$d3_survey_data,M[[i]]$epsilon) )
		colnames(it) <- c("Year","Index","Gear","Area","Group","Sex","wt","timing","Residual")
		it <- data.frame(Model=names(M)[i],it)
		mdf <- rbind(mdf,it)
	}
	print(head(mdf,3))
	print(mdf)
	p <- ggplot(mdf,aes(x=factor(Year),Residual,fill=Gear))
	p <- p + geom_bar(width=0.75,position="dodge",stat='identity')
	p <- p + labs(x="Year",y="log residual",fill="Gear")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}


# .plotIndex	<- function( repObj, annotate=FALSE )
# {
# 	#line plot for relative abundance indices
# 	with(repObj, {
# 		if(is.matrix(it)){
# 			xx=t(iyr)
# 			yy=t(it)
# 		}else{
# 			xx=iyr
# 			yy=it
# 		}
# 		n=nrow(t(as.matrix(yy)))
# 		yrange=c(0, max(yy, na.rm=TRUE))
		
# 		matplot(xx, yy, type="n", axes=FALSE,
# 			xlab="Year", ylab="Relative abundance", 
# 			ylim=yrange , main=paste(stock))
		
# 		matlines(xx, yy, col="black",type="o", pch=1:n)
		
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		grid()
		
# 		if( annotate )
# 		{
# 			txt=paste("Survey",1:n)
# 			legend("top", txt, lty=1:n, pch=1:n, bty="n")
# 		}
# 	})
# }
