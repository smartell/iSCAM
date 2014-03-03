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
	cat("Plot weight at age")
	for(i in 1:n)
	{
		age <- seq(M[[i]]$sage,M[[i]]$nage)
		wt  <- data.frame(Model=names(M[i]),M[[i]]$d3_wt_avg)
		
		colnames(wt) <- c("Model","Year","Area","Stock","Sex",paste(age))
		df  <- melt(wt,id.vars=c("Model","Year","Area","Stock","Sex"))
		mdf <- rbind(mdf,df)
	}

	mdf <- cbind(mdf,Cohort=as.factor(mdf$Year-as.double(mdf$variable)))
	if(diff(range(age)) > 10)
	{
		grp <- factor(mdf$variable,levels=pretty(unclass(mdf$variable)))
		mdf <- cbind(mdf,grp=grp)		
	}
	else
	{
		grp <- mdf$variable
		mdf <- cbind(mdf,grp=grp)
	}
	
	print(head(mdf,3))

	p <- ggplot(mdf,aes(x=Year,y=value,col=grp,linetype=factor(Sex),fill=factor(Sex)) )
	p <- p + stat_smooth(alpha=0.2,lineend="butt",method="loess") 
	p <- p + geom_point(data=mdf,aes(x=Year,y=value,col=grp))
	p <- p + labs(x="Year",y="Weight",col="Age",fill="Sex",linetype="Sex")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
	
}