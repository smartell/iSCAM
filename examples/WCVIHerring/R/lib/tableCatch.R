.tableCatch	<- function( hdr )
{
	#hdr contains thie checked files for extracting the catch
	nRuns <- nrow(hdr)
	
	#Now loop over each file and assemble a catch table using cbind
	ctTable=NULL
	for(i in 1:nRuns)
	{
		repObj	<- read.rep(paste(hdr$Control.File[i],".rep", sep=""))
		
		#report only non-zero catches in the table.
		ct = t(as.matrix(repObj$ct))
		iCols <- apply( ct,2,function(x) { sum(diff(x))!=0.0 } )
		nCols <- length(iCols[iCols])
		ct[ct==0]=NA
		if(i==1)
		{
			tmp <- cbind(repObj$yr, signif(ct[, iCols], 3))
			colnames(tmp) = c("Year", paste("Gear",1:nCols))
			cgrp = c("Stock",hdr$Stock[i])
			ncgrp = c(1, nCols)
		}
		else
		{
			tmp <- signif(ct[, iCols], 3)
			colnames(tmp) = paste("Gear",1:nCols)
			cgrp = c(cgrp, hdr$Stock[i])
			ncgrp = c(ncgrp, nCols)
		}
		
		ctTable = cbind(ctTable, tmp)
		
		
	}
	print(ctTable)
	cap <- "Observed catch by gear type and year for each stock."
	
	fn=paste(.TABLEDIR, "CatchTable.tex", sep="")
	tmp <- latex(ctTable, file=fn, rowname=NULL, longtable=TRUE
		, landscape=TRUE, lines.page=60, cgroup=cgrp, n.cgroup=ncgrp
		, caption=cap, label="TableCatch", na.blank=TRUE, vbar=TRUE
		, size="footnotesize")
	
	cat("Latex table saved as", fn)
}