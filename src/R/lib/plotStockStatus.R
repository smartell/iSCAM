# Stock Status plot
# Steven Martell

.plotStockStatus	<- function( repObj )
{
	print("	.plotStockStatus")
	
	#KOBE plots
	## This routine needs some work to accomodate multiple
	## fleets. Also need to address the Fmsy calculation in iscam.
	with(repObj, {
		xx = Bstatus[1:length(yr)]
		yy = Fstatus
		n  = dim(yy)[1]
		for(i in 1:n)
		{
			if(sum(yy[i, ])!=0)
			{
				y = yy[i, ];y[y==0]=NA; ii=!is.na(y)
				plot(xx[ii], y[ii], type="l", xlim=c(0, max(2, xx)), 
				ylim=c(0, max(2, yy)), xlab="Spawning biomass/SBmsy", ylab="Ft/Fmsy")
				# Neopolitn colors
				rect(0, 0, 1, 1, col=colr("yellow", 0.5), border=NA)
				rect(1, 1, max(2, xx),max(2, yy),col=colr("yellow", 0.5),border=NA)
				rect(1, 0, max(2, xx), 1, col=colr("green", 0.5), border=NA)
				rect(0, 1, 1, max(2, yy), col=colr("red", 0.5), border=NA)
				print(ii)
				lines(xx[ii], y[ii], type="l", col=colr(1, 1))
				text(xx[ii], y[ii], yr, cex=0.5, col=colr(1, 1))
			}
		}
		#matplot(xx, t(yy), type="l", xlim=c(0,max(2,xx)), 
		#ylim=c(0,max(2,yy)),xlab="Spawning biomass/SBmsy", ylab="Ft/Fmsy")

		
		#text(xx, yy, paste(yr), cex=0.75, col=colr(1, 1))
	})
}
