## plotIt.R  
## plot relative abundance indices


.plotIt	<-
function( M, ... )
{
	# M is a model object
	n	<- length(M)
	
	for(i in 1:n)
	{
		o	<- M[[i]]$rep
		x	<- o$iyr
		y	<- o$it
		if(is.vector(y))
		{
			y <- t(as.matrix(y))
			x <- t(as.matrix(x))
		}
		if(is.matrix(y))
		{
			mm				<- melt(t(y))
			mm[, 1]			<- as.vector(t(x))
			colnames(mm)	<- c("Year","Gear","It")
			mm	<- na.omit(mm)
			mm[, 2] <- factor(mm[, 2])
			p	<- ggplot(mm, aes(x=Year, y=It)) + labs(x="Year", y="Relative Abundance")
			p	<- p + geom_point() + theme_iscam()
			p	<- p + facet_wrap(~Gear, scale="free_y")
			print(p)
		}	
		
	}
}