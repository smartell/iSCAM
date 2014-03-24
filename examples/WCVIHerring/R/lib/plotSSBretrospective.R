
.plotSSBretrospective <- function( M )
{
	cat(".plotSSBretrospective")
	n   <- length(M)
	mdf <- NULL
	for(i in 1:n)
	{
		df <- data.frame(Model=names(M)[i],TYear=max(M[[i]]$yr),Year=M[[i]]$yr,
		                 Sbt=M[[i]]$sbt[1:length(M[[i]]$yr)])
		nn <- length(M[[i]]$retSbt)
		ddf<- NULL
		for(j in 1:nn){
			ndf<- data.frame(Model=names(M)[i],
			            TYear=max(M[[i]]$retYrs[[j]]),
			            Year=M[[i]]$retYrs[[j]],
			            Sbt=M[[i]]$retSbt[[j]])
			ddf <- rbind(ddf,ndf)
		}
		df <- rbind(df,ddf)
		mdf <- rbind(mdf,df)
	}
	p <- ggplot(mdf,aes(Year,Sbt,color=factor(TYear)),size=2)+geom_line()
	p <- p + labs(x="Year",y=paste("Spawning biomass ",.UNITS),color="Terminal Year")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}

.plotSSBsquid <- function( M )
{
	cat(".plotSSBsquid \n")
	n   <- length(M)
	mdf <- NULL
	for(i in 1:n)
	{
		df <- data.frame(Model=names(M)[i],TYear=max(M[[i]]$yr),Year=M[[i]]$yr,
		                 Sbt=M[[i]]$sbt[1:length(M[[i]]$yr)])
		SB <- df$Sbt
		df$Sbt = 0
		nn <- length(M[[i]]$retSbt)
		ddf<- NULL
		xc <- NULL
		for(j in 1:nn){
			T  <- length(M[[i]]$retYrs[[j]])
			y  <- SB[1:T]
			x  <- M[[i]]$retSbt[[j]]
			z  <- 100*(x-y)/y
			ndf<- data.frame(Model=names(M)[i],
			            TYear=max(M[[i]]$retYrs[[j]]),
			            Year=M[[i]]$retYrs[[j]],
			            Sbt=z)
			ddf <- rbind(ddf,ndf)

			# Mohns rho
			xc<- c(xc,(x[T]-y[T])/y[T])

		}
		cat("Mohns Rho = ",sum(xc),"\n")

		df <- rbind(df,ddf)
		mdf <- rbind(mdf,df)
	}
	p <- ggplot(mdf,aes(Year,Sbt,color=factor(TYear)),size=2)+geom_line()
	p <- p + labs(x="Year",y=paste("Spawning biomass deviance ",.UNITS),color="Terminal Year")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}

.plotOldSSBretrospective <- function( repObj )
{
	plot.sb <- function(fn, ...)
	{
		if(file.exists(fn))
		{
			Obj = read.rep(fn)
			lines(Obj$yr, Obj$sbt[1:length(Obj$yr)],...)
			# points(max(Obj$yr), Obj$sbo, pch=20, ...)
			cat(fn, " ", Obj$bo, "\n")
		}
	}
	
	with(repObj, {
		plot( yr, sbt[1:length(yr)], type="l"
			, ylim=c(0,1.3*max(sbt)), lwd=2
			, xlab="Year", ylab="Spawnig biomass (1000 t)", 
			main=stock )
	
		fn = paste(Control.File,".ret",1:30, sep="")
		ix=0
		for(ifn in fn) 
		{
			ix = ix +1
			plot.sb(ifn, col=ix)
		}
		grid()
	})
}
