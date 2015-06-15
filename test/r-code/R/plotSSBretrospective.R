.plotSSBretrospective <- function( repObj )
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
