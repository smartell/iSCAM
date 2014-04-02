# Steven Martell
# Sept 6,  2012

.plotSelectivity <- function( M )
{
	cat(".plotSelectivity")
	n <- length(M)
	mdf <- NULL
	for(i in 1:n)
	{
		df <- data.frame(Model=names(M)[i],logSel=M[[i]]$log_sel)
		colnames(df)<-c("Model","Gear","Sex","Year",M[[i]]$age)

		mdf <- rbind(mdf,melt(df,id=c("Model","Gear","Sex","Year")))
	}

	p <- ggplot(mdf,aes(x=Year,y=as.double(variable),z=exp(value)/max(exp(value))))
	p <- p + stat_contour(aes(colour = ..level..))
	p <- p + labs(x="Year",y="Age",colour="Selectivity")
	# p <- p + stat_contour(geom="polygon", aes(fill=exp(value)))
	p <- p + facet_wrap(~Model+Gear+Sex,scale="free")
	print(p + .THEME)
}

.plotSelex <- function( M )
{
	cat(".plotSelex\n")
	n   <- length(M)
	mdf <- NULL
	for(i in 1:n)
	{
		df <- data.frame(Model=names(M)[i],logSel=M[[i]]$log_sel)
		colnames(df)<-c("Model","Gear","Sex","Year",M[[i]]$age)

		mdf <- rbind(mdf,melt(df,id=c("Model","Gear","Sex","Year")))
	}

	p <- ggplot(mdf,aes(x=variable,y=exp(value),colour=factor(Sex))) + geom_point()
	p <- p + geom_line(aes(group=Sex))
	p <- p + labs(x="Age",y="Selectivity",color="Sex")
	p <- p + facet_grid(Year ~ Gear)
	print(p + .THEME)
}

# .plotOldSelectivity	<- function( repObj )
# {
# 	#plot the selectivity curves (3d plots)
# 	with(repObj, {
# 		#par(mgp=c(3, 3, 5))
# 		plot.sel<-function(x, y, z, ...)
# 		{
# 			#z=exp(A$log_sel)*3
# 			#x=A$yr
# 			#y=A$age
# 			z <- z/max(z)
# 			z0 <- 0#min(z) - 20
# 			z <- rbind(z0, cbind(z0, z, z0), z0)
# 			x <- c(min(x) - 1e-10, x, max(x) + 1e-10)
# 			y <- c(min(y) - 1e-10, y, max(y) + 1e-10)
# 			clr=colorRampPalette(c("honeydew","lawngreen"))
# 			nbcol=50
# 			iclr=clr(nbcol)
# 			nrz <- nrow(z)
# 			ncz <- ncol(z)
# 			zfacet <- z[-1, -1]+z[-1, -ncz]+z[-nrz, -1]+z[-nrz, -ncz]
# 			facetcol <- cut(zfacet, nbcol)
# 			fill <- matrix(iclr[facetcol],nr=nrow(z)-1,nc=ncol(z)-1)
# 			fill[ , i2 <- c(1,ncol(fill))] <- "white"
# 			fill[i1 <- c(1,nrow(fill)) , ] <- "white"

# 			par(bg = "transparent")
# 			persp(x, y, z, theta = 35, phi = 25, col = fill, expand=5, 
# 				shade=0.75,ltheta=45 , scale = FALSE, axes = TRUE, d=1,  
# 				xlab="Year",ylab="Age",zlab="Selectivity", 
# 				ticktype="simple", ...)
			
# 			#require(lattice)
# 			#wireframe(z, drap=TRUE, col=fill)
# 		}
# 		ix=1:length(yr)
# 		for(k in 1:ngear){
# 			plot.sel(yr, age, exp(log_sel[log_sel[,1]==k,-1]), 
# 			main=paste(stock, "Gear", k))
# 			#file.name=paste(prefix, "Fig9",letters[k],".eps", sep="")
# 			#if(savefigs) dev.copy2eps(file=file.name, height=8, width=8)
# 		}
		
# 	})
# }
