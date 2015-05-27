## generics for S3 classes for summary.admb summary.iscam


summary.admb <- function(object, ...)
{
	A 	<- object
	print(A)
}

print.admb <- function(x, digits = max(3, getOption("digits")-3), ...)
{
	cat("------------------------------------------------")
	cat("\n", paste("Number of estimated parameters", x$nopar, "of", x$npar, "\n"))
	cat(paste(" Objective function value ", x$nlogl, "\n"))
	cat(paste(" Max gradient componenet", x$maxgrad, "\n"))
	
	df <- data.frame(parameter=x$names,
		mode=format(x$est,digits=digits), 
		sd=format(x$std, digits=digits))
	cat("------------------------------------------------")
	cat("\n","Parameter estimates:\n")
	print(df)
	cat("------------------------------------------------") 
}

print.iscam <- function(object, digits = max(3, getOption("digits")-3), ...)
{
	A	<- object
	pname <- c(	"Unfished biomass (bo)", 
				"Steepness (h)", 
				"Natural mortality (M)", 
				"Maximum sustainable yield (MSY)", 
				"Spawnig biomass at MSY (Bmsy)", 
				"Average fishing rate at MSY (Fmsy)")
	h	<- A$kappa/(4+A$kappa)
	if(A$rectype==2) h<-0.2*A$kappa^(0.8)
	x	<- format(c(A$bo, h, A$m, A$msy, A$bmsy, A$fmsy), digits=digits)
	df	<- data.frame(variable=pname, value=x)
	print(df)
}