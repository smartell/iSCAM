require(reshape2)
require(ggplot2)

# |----------------------------------------------------------------------------------|
# | .plotWeightAtAge
# |----------------------------------------------------------------------------------|
# | - M[[1]]$wt_avg is the object being plotted.
# | - year area stock sex |age columns (sage, nage) of weight at age data |

.plotWeightAtAge <- function( M )
{
	n   <- length(M)
	mdf <- NULL

	for(i in 1:n)
	{
		age <- seq(M[[i]]$sage,M[[i]]$nage)
		wt  <- data.frame(Model=names(M[i]),M[[i]]$wt_avg)
		
		colnames(wt) <- c("Model","Year","Area","Stock","Sex",paste(age))
		mdf <- rbind(mdf,melt(wt,id.vars=c("Model","Year","Area","Stock","Sex")))
	}
	mdf <- cbind(mdf,Cohort=as.factor(mdf$Year-as.double(mdf$variable)))
	print(head(mdf,3))

	p <- ggplot(mdf,aes(x=Year,y=value,col=factor(Sex),linetype=variable))
	p <- p + stat_smooth(alpha=0.1,lineend="butt") + geom_point()
	p <- p + labs(x="Year",y="Weight",col="Sex",linetype="Age")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
	
}