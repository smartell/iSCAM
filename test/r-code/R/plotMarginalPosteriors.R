# Steven Martell
# Sept 7,  2012

.plotMarginalPosteriors	<- function( admbObj )
{
	print("	.plotMarginalPosteriors")
	#Marginal distributions & priors for theta
	op=par(no.readonly=T)
	par(las=1,mar=c(5, 4, 1, 1), oma=c(1, 1, 1, 0))

	## Read control file to get bounds and priors for theta
	## ctrl=read.table(A$control.file, header=F, skip=13, nrow=6)

	with(admbObj, {
		cat("Parameter controls\n", ctrl)
		std=apply(mcmc[,1:7],2,sd)
		nr=length(std[std!=0])
		if(nr > 6)
		{
			nRow=3; nCol=3
		}
		else if(nr>4)
		{
			nRow=3; nCol=2
		}
		else
		{
			nRow=2; nCol=2
		}
		par(mfcol=c(nRow, nCol))
		for(i in 1:7){
			if(std[i]!=0){
				ps = mcmc[, i]  #posterior samples
				xl=range(ps)

				hist(ps,xlab=colnames(mcmc[i]),prob=T, 
					main="", ylab="", col="lightgrey",breaks=30, 
					xlim=xl, border=NA)#, ...)
				
				
				## Add priors
				nfn=c("dunif","dnorm","dlnorm","dbeta","dgamma")
				pt = ctrl[i, 5]+1
				fn=match.fun(nfn[pt])
				p1=ctrl[i, 6]; p2=ctrl[i, 7]
				#browser()
				if(pt!=4)
					curve(unlist(lapply(x,fn,p1,p2)),
						xl[1],xl[2],add=T, col=colr(4, 0.7), lwd=2)
				else
					curve(unlist(lapply((x-ctrl[i,2])/
						 (ctrl[i,3]-ctrl[i,2])
						,fn,p1,p2)),xl[1],xl[2],add=T, col=colr(4, 0.7), lwd=2)
			}
		}
		mtext(c("Parameter", "Probability density",paste(stock)), c(1, 2, 3), 
			outer=T, line=-1, las=0)
	})
	par(op)
}
