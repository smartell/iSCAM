#Catarina Wor
# Jun 29 2015

require(reshape)



.MSEplotCatch <- function(M, no=3, ci=95)
{
	df <- as.data.frame(M[[1]]$catch.df)
	df2<- as.data.frame(M[[2]]$rawcatch.df)
	
	cat(".plotCatch\n")

	if(ci==95){up=df$ct975;low=df$ct025}
	if(ci==90){up=df$ct95;low=df$ct05}
	if(ci==50){up=df$ct75;low=df$ct25}
	
	if(.OVERLAY)
	{
		
		p <- ggplot(df,aes(Year,ct50,col=gear)) + geom_line(width=2)
		p <- p +  geom_ribbon(aes(ymax=up, ymin=low,fill=Scenario),alpha=0.2)
		p <- p + facet_grid(.~Procedure)

		if(no>0)
			{
				itx <- sample(1:length(grep("ct",colnames(df2))),no)
				itxx <- grep("ct",colnames(df2))[itx]
				new.df <- df2[,c(1:4,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year","gear"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',linetype='variable', col='gear'))
			}

	}
	else
	{
		p <- ggplot(df,aes(Year,ct50),fill=gear) + geom_line(width=2)
		p <- p + geom_ribbon(aes(ymax=up, ymin=low,fill=gear),alpha=0.2)
		p <- p + facet_grid(Scenario~Procedure,scales="free")

		if(no>0)
			{
			
				itx <- sample(1:length(grep("ct",colnames(df2))),no)
				itxx <- grep("ct",colnames(df2))[itx]
				new.df <- df2[,c(1:4,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year","gear"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',linetype='variable',col="gear"))
			}
	}
	# p <- p + geom_line(data=bt,aes(Year,Bo),col="blue")
	p <- p + labs(x="Year",y=paste("Catch",.UNITS))
	print(p + .THEME)
}
