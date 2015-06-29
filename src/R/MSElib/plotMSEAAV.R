#Catarina Wor
# Jun 29 2015

require(reshape)



.MSEplotAAV <- function(M, no=3, ci=95)
{
	df <- as.data.frame(M[[1]]$AAV.df)
	df2<- as.data.frame(M[[2]]$rawAAV.df)
	
	cat(".plotAAV\n")

	#if(ci==95){up=df$ct975;low=df$ct025}
	#if(ci==90){up=df$ct95;low=df$ct05}
	#if(ci==50){up=df$ct75;low=df$ct25}
	
	if(.OVERLAY)
	{		
		p <- ggplot(df2,aes(as.factor(Year),AAV,col=Scenario)) + geom_boxplot()
		p <- p + facet_grid(.~Procedure)
	}
	else
	{
		p <- ggplot(df2,aes(as.factor(Year),AAV,col=Scenario)) + geom_boxplot()
		p <- p + facet_grid(Scenario~Procedure,scales="free")
	}
	# p <- p + geom_line(data=bt,aes(Year,Bo),col="blue")
	p <- p + labs(x="Year",y=paste("Annual Absolute Variation in catch",.UNITS))
	print(p + .THEME)
}
