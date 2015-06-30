#Catarina Wor
# Jun 25 2015
#parei aqui
require(reshape)

.MSEplotSubLegal <- function(M, no=3, ci=95)
{
	df <- as.data.frame(M[[1]]$sublegal.df)
	df2<- as.data.frame(M[[2]]$rawsublegal.df)
	
	if(ci==95){up=dt975;low=dt025}
	if(ci==90){up=dt95;low=dt05}
	if(ci==50){up=dt75;low=dt25}

	cat(".plotSubLegalCatch\n")

	
	if(.OVERLAY)
	{
		
		p <- ggplot(df,aes(Year,dt50,col=gear)) + geom_line(width=2)
		p <- p +  geom_ribbon(aes(ymax=up,ymin=low,fill=Scenario),alpha=0.2)
		p <- p + facet_grid(.~Procedure)

		if(no>0)
			{

				itx <- sample(1:length(grep("dt",colnames(df2))),no)
				itxx <- grep("dt",colnames(df2))[itx]
				new.df <- df2[,c(1:4,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year", "gear"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',col='gear',linetype='variable'))

			}

	}
	else
	{
		p <- ggplot(df,aes(Year,p.Bt0.5),col=Scenario,fill=Scenario) + geom_line(width=2)
		p <- p + geom_ribbon(aes(ymax=up, ymin=low,col=Scenario,fill=Scenario),alpha=0.2)
		p <- p + facet_grid(Scenario~Procedure,scales="free")
	}
	# p <- p + geom_line(data=bt,aes(Year,Bo),col="blue")
	p <- p + labs(x="Year",y=paste("Spawning biomass",.UNITS))
	print(p + .THEME)
}


.MSEplotMilkaSpawnBiomass <- function(M, no=3, ci=95)
{
	df <- as.data.frame(M[[1]]$biomass.df)
	df2<- as.data.frame(M[[2]]$rawbiomass.df)
	
	cat(".plotSpawnBiomass\n")

	if(ci==95){up=df$t.Bt0.975;low=df$t.Bt0.025}
	if(ci==90){up=df$t.Bt0.95;low=df$t.Bt0.05}
	if(ci==50){up=df$t.Bt0.75;low=df$t.Bt0.25}
	
	if(.OVERLAY)
	{
		
		p <- ggplot(df,aes(Year,t.Bt0.5,col=Scenario)) + geom_line(width=2)
		p <- p +  geom_ribbon(aes(ymax=up, ymin=low,fill=Scenario),alpha=0.2)
		p <- p + facet_grid(.~Procedure)

		if(no>0)
			{
				itx <- sample(1:length(grep("m_sbt",colnames(df2))),no)
				itxx <- grep("m_sbt",colnames(df2))[itx]
				new.df <- df2[,c(1:3,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',col='variable'))
			}

	}
	else
	{
		p <- ggplot(df,aes(Year,m.Bt0.5),col=Scenario,fill=Scenario) + geom_line(width=2)
		p <- p + geom_ribbon(aes(ymax=up, ymin=low,col=Scenario,fill=Scenario),alpha=0.2)
		p <- p + facet_grid(Scenario~Procedure,scales="free")

		if(no>0)
			{
				itx <- sample(1:length(grep("m_sbt",colnames(df2))),no)
				itxx <- grep("m_sbt",colnames(df2))[itx]
				new.df <- df2[,c(1:3,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',col='variable'))
			}
	}
	# p <- p + geom_line(data=bt,aes(Year,Bo),col="blue")
	p <- p + labs(x="Year",y=paste("Spawning biomass",.UNITS))
	print(p + .THEME)
}
