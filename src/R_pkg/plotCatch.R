## ------------------------------------------------------------------------- ##
## plotCatch                                                                 ##
## ------------------------------------------------------------------------- ##


# .plotCatch
# Purpose:    Plot the catch data from an iscam model object    
# Parameters: M a list of model objects 
# Returns:    NULL
# Source:     
.plotCatch <-
function( M, ... )
{
	# M is a list of model objects
	n	<- length(M)
	
	# Construct data.frame for plotting
	mm	<- NULL
	for(i in 1:n)
	{
		print(names(M))
		ct	<- M[[i]]$rep$ct
		if(!is.matrix(ct)) ct = t(ct)
		df	<- data.frame(Year=M[[i]]$rep$yr, Model=names(M)[i], Catch=t(ct))
		mm	<- rbind(mm, melt(df, id.vars=c("Model", "Year")))
	}
	
	names(mm) <- c("Model", "Year", "Gear", "Catch")
	
	
	p	<- ggplot(mm,  aes(factor(Year),Catch, fill=Gear)) 
	p	<- p + geom_bar() + guides(fill=FALSE) + labs(x="Year", y="Catch (t)")
	p	<- p + theme_iscam(xlas=90)
	p	<- p + facet_wrap(~Model, scales="free")
	
	print(p)
}

