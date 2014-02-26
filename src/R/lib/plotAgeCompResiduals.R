# Steven Martell
# Aug 30,  2012

.plotAgeCompResiduals	<- function(repObj)
{
	#Bubble plot of age-composition data
	#A is the observed age-comps
	#Ahat is the predicted age-comps (proportions)
	with( repObj, {
		if(!is.null(repObj$d3_A)){
			nagear = unique(d3_A[, 2])
			xrange = range(d3_A[, 1])
			#par(mfcol=c(length(nagear), 1))
			j=0
			for(i in nagear)
			{
				j=j+1
				ac = subset(A_nu, A_nu[, 2]==i)
				xx = ac[, 1]
				zz = t(ac[, -1:-2])
			
				# plot residuals
				plotBubbles(zz, xval = xx, yval = age, rres=FALSE, hide0=TRUE,  
					las=.VIEWLAS, xlab="Year", ylab="Age", frange=0.05, size=0.5*age_tau2[j],
					bg=colr("white", 0.5), xlim=xrange,main=paste(stock, "Gear", i))
					
				grid()
				
				title(main=paste("Variance=",round(age_tau2[j],3)),line=-1,cex.main=0.75)
			}
			
		}
		else{print("There is no age-composition data")}
	})
}
