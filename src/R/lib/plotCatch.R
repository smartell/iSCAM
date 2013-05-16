# R-script for a barplot for the observed catch.
# Steven Martell
# Aug 28,  2012

require(ggplot2)
require(reshape2)

.plotCatch <- function(M)
{
	n <- length(M)

	mdf <- NULL
	for( i in 1:n )
	{
		ct <- data.frame(M[[i]]$catch_data)
		colnames(ct) <- c("Year","Gear","Area","Group","Sex","Type","Catch")
		ct <- data.frame(Model=names(M)[i],ct)
		mdf <- rbind(mdf,ct)
	}
	print(head(mdf,3))

	p <- ggplot(mdf,aes(x=factor(Year),Catch,fill=Gear))
	p <- p + geom_bar(width=0.75)
	p <- p + labs(x="Year",y="Catch (t)")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}

# .plotCatch	<- function( repObj, legend.txt=NULL )
# {
# 	#barplot of the observed catch
# 	with(repObj, {
# 		tmp = obs_ct
# 		if(!is.null(dim(tmp)))
# 		{
# 			iRows <- apply( tmp,1,function(x) { sum(diff(x))!=0.0 } )
# 			tmp   <- tmp[iRows, ]
# 		}

# 		barplot( tmp, names.arg=yr,axes=FALSE,  
# 			xlab="Year", ylab="Catch (1000 t)",main=paste(stock),  
# 			legend.text = legend.txt )
# 		axis( side=2, las=.VIEWLAS )
# 	})
# }
