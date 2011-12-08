## ------------------------------------------------------------------------- ##
## plotAgeComps                                                                 ##
## ------------------------------------------------------------------------- ##


.plotAgeComps <-
function( M, ... )
{
	## M is a list of model objects
	require(hacks)
	n	<- length(M)
	
	for( i in 1:n )
	{
		A	<- data.frame(M[[i]]$rep$A)
		Age	<- M[[i]]$rep$a_sage:M[[i]]$rep$a_nage
		colnames(A)	<- c("Year", "Gear", paste(Age))
		mm	<- melt(A, id=c("Year", "Gear"),variable_name="Age")
		print(head(mm))
		p	<- ggplot(mm, aes(x=Year, y=Age), ylab="Age") 
		p	<- p + geom_point(aes(size=value/max(value)), alpha=0.5)
		p	<- p + scale_colour_gradient(low="grey", high="black", alpha=0.5)
		p	<- p + scale_area() + facet_wrap(~Gear)
		print(p)
		
	}
}

.plotAgeHist <-
function( M, ... )
{
	## M is a list of model objects
	n	<- length(M)

	for( i in 1:n )
	{
		A	<- data.frame(M[[i]]$rep$A)
		Age	<- M[[i]]$rep$a_sage:M[[i]]$rep$a_nage
		colnames(A)	<- c("Year", "Gear", paste(Age))
		mm	<- melt(A, id=c("Year", "Gear"),variable_name="Age")
		print(head(mm))
		ig	<- unique(mm$Gear)
		for(j in ig)
		{
			sm	<- subset(mm, Gear==j)
			sm	<- cbind(sm, density=sm$value/median(sm$value))
			p	<- ggplot(sm, aes(x=Age, y=density)) 
			p	<- p + geom_bar() + facet_wrap(~Year)
			print(p)
		}
	}
}


.plotAgecomps_old	<- function(repObj, meanAge = FALSE )
{
	#Bubble plot of age-composition data
	#A is the observed age-comps
	#Ahat is the predicted age-comps (proportions)
	with( repObj, {
		if(!is.null(repObj$A)){
			nagear = unique(A[, 2])
			xrange = range(A[, 1])
			#par(mfcol=c(length(nagear), 1))
			for(i in nagear)
			{
				ac = subset(A, A[, 2]==i)
				xx = ac[, 1]
				zz = t(ac[, -1:-2])
				
				
				# plot proportions-at-age (cpro=TRUE)
				plotBubbles(zz, xval = xx, yval = age, cpro=TRUE, hide0=TRUE,  
					las=.VIEWLAS, xlab="Year", ylab="Age", frange=0.0, size=0.1, 
					bg=colr("steelblue", 0.5),main=paste(stock, "Gear", i), 
					xlim=xrange)
				
				if( meanAge )
				{
					tz = t(zz)
					p = t(tz/rowSums(tz))
					abar = colSums(t(tz/rowSums(tz))*age)
					sbar = sqrt(colSums(p*(1-p)*age))
					sbar = 1.96*colSums(sqrt(p*(1-p))/sqrt(age))
				
					lines( xx, abar, col=colr("steelblue", 0.75), lwd=2 )
				
					yy = c(exp(log(abar)+log(sbar)), rev(exp(log(abar)-log(sbar))))
					polygon(c(xx, rev(xx)),yy,border=NA,col=colr("steelblue",0.25))
				}
			}
		}
		else{print("There is no age-composition data")}
	})
}
