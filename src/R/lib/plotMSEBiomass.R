#Catarina Wor
# Jun 25 2015



.MSEplotSpawnBiomass <- function( M )
{
	n <- length(M)
	cat(".plotSpawnBiomass\n")

	
	if(.OVERLAY)
	{
		p <- ggplot(biomass.df,aes(Year,p.Bt0.5,col=Scenario, sixe=Procedure)) + geom_line(width=2)
		p <- p +  geom_ribbon(aes(ymax=p.Bt0.975, ymin=p.Bt0.025,fill=Scenario),alpha=0.2)
	}
	else
	{
		p <- ggplot(biomass.df,aes(Year,p.Bt0.5),col=Scenario,fill=Scenario) + geom_line(width=2)
		p <- p + geom_ribbon(aes(ymax=p.Bt0.975, ymin=p.Bt0.025,col=Scenario,fill=Scenario),alpha=0.2)
		p <- p + facet_wrap(~Scenario +Procedure ,scales="free")
	}
	# p <- p + geom_line(data=bt,aes(Year,Bo),col="blue")
	p <- p + labs(x="Year",y=paste("Spawning biomass",.UNITS))
	print(p + .THEME)
}