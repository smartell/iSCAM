# R-script for a barplot for the observed catch.
# Steven Martell
# Aug 28,  2012

.plotCatch	<- function( repObj, legend.txt=NULL )
{
	#barplot of the observed catch
	with(repObj, {
		tmp = obs_ct
		if(!is.null(dim(tmp)))
		{
			iRows <- apply( tmp,1,function(x) { sum(diff(x))!=0.0 } )
			tmp   <- tmp[iRows, ]
		}

		barplot( tmp, names.arg=yr,axes=FALSE,  
			xlab="Year", ylab="Catch (1000 t)",main=paste(stock),  
			legend.text = legend.txt )
		axis( side=2, las=.VIEWLAS )
	})
}
