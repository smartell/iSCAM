## ------------------------------------------------------------------------- ##
## plot_resid                                                                ##
## ------------------------------------------------------------------------- ##

.plot_resid <-
function( M,  which="eta", ... )
{
	# use which to determine with set of residuals to plot for each model.
	n	<- length(M)
	
	# Construct data frame for plotting
	mm	<- NULL
	for( i in 1:n )
	{
		O  <- M[[i]]$rep
		switch(which, 
			eta = {
				x <- O$yr
			}, 
			epsilon = {
				x <- O$iyr
			}, 
			delta = {
				x <- O$yr
			})
		print(x)
		idy <- match(which,names(O))
		y   <- O[[idy]]
		if(!is.matrix(y)) y <- t(y)
		
		df	<- data.frame(Year=x, Model=names(M)[i], Residual = t(y))
		print(df)
		mm	<- rbind(mm, melt(df, id.vars=c("Model", "Year")))
	}
	
	names(mm) <- c("Model", "Year", which, "Residual")
	class(mm) <- c(which, "data.frame")
	print(mm)
}

plot.eta <-
function( object, ... )
{
	p	<- ggplot(object, aes(x=Year, y=Residual, colour=eta))
	p	<- p + geom_point(alpha=0.8) + theme_iscam()
	p	<- p + facet_wrap(~Model, scales="free")
	print(p)
}

plot.epsilon <- 
function(object, ...)
{
	p	<- ggplot(object, aes(x=Year, y=Residual, colour=epsilon))
	p	<- p + geom_point(alpha=0.8, na.rm=TRUE) + theme_iscam()
	p	<- p + facet_wrap(~Model, scales="free")
	print(p)
}

plot.delta <-
function(object, ...)
{
	p	<- ggplot(object, aes(x=Year, y=Residual)) + geom_point()
	p	<- p + geom_hline(yintercept=0, col="grey50", size=0.5)
	
	p
}


deprecate.plot.resid <-
function(object, which="epsilon",  ...)
{
	#assign a class to all residual objects and plot
	# epsilon -> residuals in survey
	# delta -> residuals in stock-recruitment
	o	<- object
	
	switch(which,
		epsilon = {
			x	<- (o$iyr)
			y	<- (o$epsilon)
			if(is.vector(y)) y <- t(y)
			df		<- data.frame(melt(t(y)))
			df$X1	<- as.vector(t(x))
			df$X2	<- factor(df$X2, levels=unique(df$X2))
			df		<- na.omit(df)
			colnames(df)	<- c("Year","Gear","Residual")
		}, 
		delta = {
			x	<- (o$yr)
			y	<- (o$delta)
			n1	<- length(x) - length(y)+1
			n2	<- length(x)
			df	<- data.frame(Year=x[n1:n2], Residual=y)
		}
	)
	
	
	class(df) <- c(which, "data.frame")
	plot(df)

}
