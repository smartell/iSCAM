## Plot Mean Weight at age over time.

.plotMeanWt	<-
function( M, ... )
{
	## M is a model object
	n	<- length(M)
	
	for(i in 1:n)
	{
		df	<- M[[i]]$rep$wt_obs
		colnames(df) <- paste("Age", M[[i]]$rep$age)
		df	<- data.frame(Year=M[[i]]$rep$yrs, df)
		m	<- melt(df, id=1)
		names(m) = c("Year", "Age", "Weight")
		p	<- ggplot(m, aes(Year, Weight, col=Age))
		print(p+stat_smooth()+geom_point())
	}
}


.ggplotMeanwt <- function( repObj )
{
	
	#Use ggplots to plot mean weight-at-age
	require(ggplot2)
	tmp <- repObj$wt_obs
	colnames(tmp)<-paste("Age", repObj$age)
	W <- data.frame(Year=repObj$yrs,tmp)
	M <- melt(W, id=1)
	names(M)=c("Year","Age","Weight")
	p <- ggplot(M,aes(Year,Weight,color=Age))
	p_title <- paste(repObj$stock)
	print(p+stat_smooth()+geom_point()+opts(title=p_title))
	
}
