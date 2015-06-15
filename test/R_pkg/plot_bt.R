## ------------------------------------------------------------------------- ##
## plot_bt                                                                   ##
## ------------------------------------------------------------------------- ##


# .plot_bt
# Purpose:    A generic time series plot using ggplot
# Parameters: M,  a list of model objects, 
# Parameters: SeriesType, a string corresponding to the type of series ("bt","sbt")
# Returns:    NULL
# Source:     
.plot_bt <-
function( M, which="bt", ... )
{
	# M is a list of model objects
	n	<- length(M)
	
	# Set up ylbls for plots depending on value of which
	case <- c("bt","sbt","ut")
	ylbl <- c("Biomass","Spawning Biomass", "Exploitation rate")
	iy   <- match(which, case)
	cat(ylbl[iy], "\n")
	
	# Construct data frame
	mm	<- NULL
	for( i in 1:n )
	{
		idx <- match(which,names(M[[i]]$rep))
		x	<- unlist(M[[i]]$rep$yrs)
		y	<- unlist(M[[i]]$rep[idx])
		
		if(length(x)!=length(y)) x	<- x[1:length(y)]
		df	<- data.frame(Year=x, Model=names(M)[i], y = y)
		mm	<- rbind(mm, melt(df, id.vars=c("Model", "Year"))) 
	}
	
	names(mm) <- c("Model", "Year", which, "Biomass")
	
	p	<- ggplot(mm, aes(x=Year, y=Biomass)) + geom_line()
	p	<- p + labs(y=ylbl[iy]) + theme_iscam()
	p	<- p + scale_y_continuous(limits=c(0, max(mm$Biomass)))
	p	<- p + facet_wrap(~Model, scales="free") 
	print(p)
	
}