# Steven Martell
# R-script for doing the logistic regression

# R-script for plotting and creating the Decision table for Pacific Herring
# Author: Steven Martell
# Date  : Aug 11,  2012

.plotRiskTable<-function( hdr )
{
	print(".plotRiskTable")
	nRuns    <- nrow(hdr)
	TACprobs <- NULL
	#browser()
	plot.risk <- function(obj, ...)
	{
		#obj is a class object with x & y vectors
		with(obj, {
			tmp.x       <- jitter(x)
			tmp.y       <- y
			tmp.y[y<=1] <- 0
			tmp.y[y> 1] <- 1
			glm.fit     <- glm(tmp.y~tmp.x, family=binomial(logit))
			
			plot(tmp.x, tmp.y, pch=".", xlab="Catch option", ...)
			lines(sort(tmp.x), sort(glm.fit$fitted.values), col=2, lwd=2)
			
			# Now append to the obj with probabilities for each tac
			tac <- unique(x)
			ab  <- coef(glm.fit)
			xx  <- ab[1]+ab[2]*tac
			py  <- exp(xx)/(1+exp(xx))
			obj$py  <- py
			obj$tac <- tac
			return(obj)
		})
	}
	
	
	for(i in 1:nRuns)
	{
		repObj <- read.admb(hdr$Control.File[i])
		mcfile=paste(hdr$Control.File[i],".mcmc", sep="")
		if(file.exists(mcfile))
			repObj$mcmc = read.table( mcfile, header=TRUE )
		else
			cat("NB. MCMC file missing.")
		pcfile=paste(hdr$Control.File[i],".proj", sep="")
		if(file.exists(pcfile))
			repObj$proj <- read.table( pcfile, header=TRUE )
			
		print(head(repObj$proj))
		with(repObj, {
			ylbls <- c(
				"P(SB decline)", 
				"P(SB < 0.25 B0)", 
				"P(SB < 0.75 B0)",     #
				"P(SB < 0.40 Bmsy)",   #
				"P(SB < 0.80 Bmsy)",   #
				"P(U > Umsy)",         #
				"P(U > 1/2 Umsy)",     #
				"P(U > 2/3 Umsy)",     #
				"P(U3+ > 20%)",        
				"P(declining trend)")
			
			par(mfcol=c(3, 2), las=2, mar=c(4, 4, 1, 1))
			
			# Decision table for Report
			tac  <- unique(proj$tac)
			TACprobs = c(tac, TACprobs)
			icol = c(2, 3, 4, 5, 6, 10)
			for(j in icol)
			{
				O   <- list()
				O$x <- proj$tac
				O$y <- proj[[j]]
				class(O) <- "risk"
				O   <- plot(O, ylab=ylbls[j-1])
				print(cbind(O$tac, O$py))
				TACprobs = cbind(TACprobs, round(O$py, 3))
			}
			colnames(TACprobs)=c("TAC", ylbls[icol])
			print(TACprobs)			
			
			fn=paste(.TABLEDIR, "table:RiskTable",hdr$Stock[i],".tex", sep="")
			cap <- paste("Decision table for", hdr$Stock[i], "where tac is the total
			allowable catch (row),  and each of the columns referst to the probability
			of an undesirable outcome (e.g.,  P(SB decline) corresponds to the probability
			of the spawning stock biomass declining from 2012 to 2013).")
			tmp <- latex(TACprobs, file=fn,title="Risk",  longtable=FALSE
				, landscape=FALSE, cgroup=NULL, n.cgroup=NULL
				, caption=cap, label=paste("Table:HerringRisk",hdr$Stock[i], sep="")
				, na.blank=TRUE, vbar=FALSE)
			 
		})
	}
}
