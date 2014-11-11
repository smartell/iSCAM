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
		A   <- data.frame(M[[i]]$d3_A)
		# Ensure proportions are being plotted.
		A[,-1:-5] <- A[,-1:-5]/rowSums(A[,-1:-5],na.rm=TRUE)
		age <- seq(min(M[[i]]$n_A_sage),max(M[[i]]$n_A_nage))
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
	p <- p + scale_area(range = c(0,10))
	p <- p + labs(x="Year",y="Age",size="Count")
	p <- p + facet_wrap(~Model+Sex+Gear,scales="free")
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
		A   <- cbind(M[[i]]$d3_A[,1:5],M[[i]]$A_nu)
		A   <- data.frame(A)
		age <- seq(min(M[[i]]$n_A_sage),max(M[[i]]$n_A_nage))
		# year gear area group sex
		A   <- data.frame(Model=names(M)[i],A)
		colnames(A) <- c("Model","Year","Gear","Area","Group","Sex",paste(age))
		mdf <- rbind(mdf,A)
	}
	mdf <- melt(mdf,id.vars=c("Model","Year","Gear","Area","Group","Sex"))
	print(head(mdf,3))

	p <- ggplot(mdf,aes(factor(Year),variable,col=factor(sign(value)),size=abs(value)))
	p <- p + geom_point(alpha=0.75)
	# p <- p + scale_size_area(max_size=6)
	p <- p + labs(x="Year",y="Age",size="Residual",colour="Sign")
	p <- p + facet_wrap(~Model+Sex+Gear,scales="free")
	print(p + .THEME)
}


.plotAgeSummary <- function( M )
{
	n <- length(M)
	cat(".plotAgeSummary\n")
	mdf <- NULL
	for( i in 1:n )
	{
		age <- seq(min(M[[i]]$n_A_sage),max(M[[i]]$n_A_nage))

		# Predicted data
		A   <- cbind(M[[i]]$d3_A[,1:5],M[[i]]$A_hat)
		A   <- data.frame(A)
		A[,-1:-5] <- A[,-1:-5]/rowSums(A[,-1:-5],na.rm=TRUE)
		agA <- aggregate(A[,-1:-5],by=list(A[,2],A[,3],A[,4],A[,5]),FUN=mean,na.rm=TRUE)
		colnames(agA) = c("Gear","Area","Group","Sex",paste(age))

		# Observed data
		O   <- data.frame(M[[i]]$d3_A)
		O[,-1:-5] <- O[,-1:-5]/rowSums(O[,-1:-5],na.rm=TRUE)
		agO <- aggregate(O[,-1:-5],by=list(O[,2],O[,3],A[,4],O[,5]),FUN=mean,na.rm=TRUE)
		colnames(agO) = c("Gear","Area","Group","Sex",paste(age))

		# Create data frame
		df  <- rbind(cbind(type="Predicted",agA),cbind(type="Observed",agO))
		# year gear area group sex
		df   <- data.frame(Model=names(M)[i],df)
		colnames(df) <- c("Model","Type","Gear","Area","Group","Sex",paste(age))
		mdf <- rbind(mdf,df)
	}
	mdf <- melt(mdf,id.vars=c("Model","Type","Gear","Area","Group","Sex"))
	print(head(mdf,3))

	p <- ggplot( mdf,aes(variable,value,col=Type) )
	p <- p + geom_point(alpha=0.75)
	# p <- p + scale_area(range = c(0,10))
	p <- p + labs(x="Age",y="Mean proportion",colour="Type")
	p <- p + facet_wrap(~Model+Sex+Gear,scales="free")
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
