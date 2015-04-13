# plotQ.R
# Steve Martell
# April 10, 2015


.plotQ <- function ( M )
{
	n <- length(M)
	cat(".plotQ\n")
	mdf <- NULL
	for( i in 1:n )
	{
		fit = M[[i]]$fit
		yr  = M[[i]]$d3_survey_data[,1]
		kk  = M[[i]]$d3_survey_data[,3]

		qt  = M[[i]]$qt
		qt  = na.omit(as.vector(t(qt)))
		ln.qt.mu = fit$est[grep("log_q_devs",fit$names)]
		ln.qt.sd = fit$std[grep("log_q_devs",fit$names)]
		ub  = ln.qt.mu + 1.96 * ln.qt.sd
		lb  = ln.qt.mu - 1.96 * ln.qt.sd

		
		df  <- data.frame(Model=names(M)[i],year=yr,gear=kk,qt=qt)
		names(df) <- c("Model","Year","Gear","Q")
		df <- subset(df,Q>0)

		df <- cbind(df,"delta"=ln.qt.mu,"ub"=ub,"lb"=lb)

		mdf <- rbind(mdf,subset(df,Q>0));

	}

	print(head(mdf))
	p <- ggplot(mdf) + geom_point(aes(Year,Q))
	p <- p + geom_point(aes(Year,delta)) + geom_pointrange(aes(Year,delta,ymax=ub,ymin=lb))
	p <- p + facet_grid(Gear~Model,scales="free")
	print(p + .THEME)
	
}