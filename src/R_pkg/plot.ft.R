## S3 generics for plotting fishing mortality rates

## Gear specific fishing rates (average)
plot.ft <- 
function(object, ...)
{
	o	<- object
	x	<- o$yr
	y	<- o$ft
	if(is.vector(y))
	{
		y <- t(as.matrix(y))
	}
	
	if(is.matrix(y))
	{
		df				<- t(y)
		rownames(df)	<- x
		mm				<- melt(df)
		colnames(mm)	<- c("Year","Gear","Ft")
		
		p	<- ggplot(mm, aes(x=Year, y=Ft, group=Gear))
		p	+ geom_line(aes(color=Gear)) + ylim(0, max(y))
	}	
}

## Age-speicific fishing mortality rates
plot.F <- 
function(object, ...)
{
	o	<- object
	x	<- o$yr
	y	<- o$age
	z	<- o$F
	
	rownames(z)	<- x
	colnames(z)	<- y
	mm			<- melt(z)
	colnames(mm)<- c("Year", "Age", "Ft")
	
	p	<- ggplot(mm, aes(x=Year, y=Ft, group=Age))
	p + geom_line(aes(color=Age))
}