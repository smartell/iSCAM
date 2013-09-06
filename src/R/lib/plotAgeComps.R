# Rscript for plotting age-comp bubble plots.
# Steven Martell
# Aug 28,  2012
require(ggplot2)
require(reshape2)

.plotAgeComps <- function( M )
{
	n <- length(M)
	cat(".plotAgeComps\n")
	mdf <- NULL
	for( i in 1:n )
	{
		A   <- data.frame(M[[i]]$A)
		age <- seq(min(M[[i]]$a_sage),max(M[[i]]$a_nage))
		# year gear area group sex
		A   <- data.frame(Model=names(M)[i],A)
		colnames(A) <- c("Model","Year","Gear","Area","Group","Sex",paste(age))
		mdf <- rbind(mdf,A)

	}
	mdf <- melt(mdf,id.vars=c("Model","Year","Gear","Area","Group","Sex"))
	BroodYear <- mdf$Year-as.double(mdf$variable)
	mdf <- cbind(mdf,BroodYear)
	print(head(mdf,3))

	p <- ggplot(mdf,aes(factor(Year),variable,size=value))
	p <- p + geom_point(alpha=0.75,aes(colour=factor(BroodYear))) 
	p <- p + scale_area(range = c(0,20))
	p <- p + labs(x="Year",y="Age",size="Count")
	p <- p + facet_wrap(~Model+Sex,scales="free")
	p <- p + scale_colour_discrete(legend=FALSE)
	print(p + .THEME)
}

.plotAgeCompResiduals <- function( M )
{
	n <- length(M)
	cat(".plotAgeCompResiduals\n")
	mdf <- NULL
	for( i in 1:n )
	{
		A   <- cbind(M[[i]]$A[,1:5],M[[i]]$A_nu)
		A   <- data.frame(A)
		age <- seq(M[[i]]$sage,M[[i]]$nage)
		# year gear area group sex
		A   <- data.frame(Model=names(M)[i],A)
		colnames(A) <- c("Model","Year","Gear","Area","Group","Sex",paste(age))
		mdf <- rbind(mdf,A)
	}
	mdf <- melt(mdf,id.vars=c("Model","Year","Gear","Area","Group","Sex"))
	print(head(mdf,3))

	p <- ggplot(mdf,aes(factor(Year),variable,col=factor(sign(value)),size=abs(value)))
	p <- p + geom_point(alpha=0.75)
	p <- p + labs(x="Year",y="Age",size="Residual",col="")
	p <- p + facet_wrap(~Model+Sex,scales="free")
	print(p + .THEME)
}

# .plotAgecomps	<- function(repObj, meanAge = FALSE )
# {
# 	#Bubble plot of age-composition data
# 	#A is the observed age-comps
# 	#Ahat is the predicted age-comps (proportions)
# 	with( repObj, {
# 		if(!is.null(repObj$A)){
# 			nagear = unique(A[, 2])
# 			xrange = range(A[, 1])
# 			#par(mfcol=c(length(nagear), 1))
# 			for(i in nagear)
# 			{
# 				ac = subset(A, A[, 2]==i)
# 				xx = ac[, 1]
# 				zz = t(ac[, -1:-2])
				
				
# 				# plot proportions-at-age (cpro=TRUE)
# 				plotBubbles(zz, xval = xx, yval = age, cpro=TRUE, hide0=TRUE,  
# 					las=.VIEWLAS, xlab="Year", ylab="Age", frange=0.05, size=0.1, 
# 					bg=colr("steelblue", 0.5),main=paste(stock, "Gear", i), 
# 					xlim=xrange)
				
# 				grid()
				
# 				if( meanAge )
# 				{
# 					tz = t(zz)
# 					p = t(tz/rowSums(tz))
# 					abar = colSums(t(tz/rowSums(tz))*age)
# 					sbar = sqrt(colSums(p*(1-p)*age))
# 					sbar = 1.96*colSums(sqrt(p*(1-p))/sqrt(age))
				
# 					lines( xx, abar, col=colr("steelblue", 0.75), lwd=2 )
				
# 					yy = c(exp(log(abar)+log(sbar)), rev(exp(log(abar)-log(sbar))))
# 					polygon(c(xx, rev(xx)),yy,border=NA,col=colr("steelblue",0.25))
# 				}
				
# 			}
# 		}
# 		else{print("There is no age-composition data")}
# 	})
# }
