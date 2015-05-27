# Steven Martell
# Sept 7,  2012
.plotReferencePoints	<- function( admbObj )
{
	print("	.plotReferencePoints")
	op=par(no.readonly=T)
	with(admbObj, {
		par(las=1,mar=c(5, 5, 1, 1), oma=c(1, 1, 1, 0))
		par(mfcol=c(2, 2))
		for(i in 8:11)
		{
			ps=mcmc[, i]
			xl=range(ps)
			hist(ps,xlab=colnames(mcmc[i]),prob=T, 
				main="", ylab="",
				xlim=xl)#, ...)
		}
		mtext(c("MSY reference points", "Probability density",paste(stock)), 
		c(1, 2, 3),outer=T, line=-1, las=0)
	})
	
	par(op)
	
}
