# Stock Status plot
# Steven Martell

.plotStockStatus	<- function( M )
{
	require(ggplot2)

	n <- length(M)
	cat(".plotStockRecruit\n")
	
	mdf <- NULL
	
	for( i in 1:n )
	{
		Bbmsy <- M[[i]]$bt[1:length(M[[i]]$yr)]/M[[i]]$bmsy
		Ffmsy <- apply(M[[i]]$ft,2,sum)/sum(M[[i]]$fmsy)

		df <- data.frame(Model=names(M[i]),Year=M[[i]]$yr,Bbmsy=Bbmsy,Ffmsy=Ffmsy)
		df <- df[order(df$Year),]

		mdf <- rbind(mdf,df)

	}

	print(head(mdf,3))

	p <- ggplot(mdf) + geom_point(aes(Bbmsy,Ffmsy, color=Year))
	p <- p + ylim(0,max(Ffmsy)) + xlim(0,max(Bbmsy))
	p <- p + geom_vline(aes(xintercept=1))
	p <- p + geom_hline(aes(yintercept=1))
	p <- p + geom_path(aes(Bbmsy,Ffmsy,color=Year))
	p <- p + geom_point(aes(x=Bbmsy[length(Year)],y=Ffmsy[length(Year)]),size=4,color="Red")
	p <- p + geom_text(aes(label=Year[length(Year)], x=Bbmsy[length(Year)],y=Ffmsy[length(Year)]),hjust=0, vjust=0)
	p <- p + labs(x=expression(paste("B/","B"["MSY"])),y=expression(paste("F/","F"["MSY"])))
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)

}




.plotOldStockStatus	<- function( repObj )
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
