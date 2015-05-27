# Rscript for plotting age-comp bubble plots.
# Steven Martell
# Aug 28,  2012

.plotAgecomps	<- function(repObj, meanAge = FALSE )
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
					las=.VIEWLAS, xlab="Year", ylab="Age", frange=0.05, size=0.1, 
					bg=colr("steelblue", 0.5),main=paste(stock, "Gear", i), 
					xlim=xrange)
				
				grid()
				
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
