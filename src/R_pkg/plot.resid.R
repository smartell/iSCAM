## plot residuals

plot.resid <-
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

plot.epsilon <- 
function(object, ...)
{
	p	<- ggplot(object, aes(x=Year, y=Residual, colour=Gear))
	p	<- p + geom_point() + geom_hline(yintercept=0, col="grey50", size=0.5) 
	
	p
}

plot.delta <-
function(object, ...)
{
	p	<- ggplot(object, aes(x=Year, y=Residual)) + geom_point()
	p	<- p + geom_hline(yintercept=0, col="grey50", size=0.5)
	
	p
}