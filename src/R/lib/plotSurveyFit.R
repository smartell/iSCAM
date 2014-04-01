require(ggplot2)
require(reshape2)
# fixed a plotting bug when NA present.
.plotSurveyFit <- function(M)
{
	#SJDM.  For the 95% Confidence interval for the observations.
	#  Assume ln(I_t) = ln(\hat(I_t)) + 0.5*\sigma^2
	# The 95% CI is then ci = ln(I_t)± 1.96*\sigma/it_wt
	# Plot (exp(ci))

	n <- length(M)
	cat(".plotSurveyfit\n")
	mdf <- NULL
	for( i in 1:n )
	{
		it_hat <- na.omit(as.vector(t(M[[i]]$it_hat)))
		it     <- M[[i]]$d3_survey_data[,2]
		sig<- M[[i]]$sig
		norm_it_wt <- M[[i]]$it_wt
		#fit <- M[[i]]$fit
		#rho <- fit$est[fit$names=="theta[6]"]
		#vartheta <- fit$est[fit$names=="theta[7]"]
		#sig <- sqrt(rho/vartheta)
		Ub <- exp(log(it)+1.96*sig/norm_it_wt)
		Ub <- Ub[1:length(Ub)]
		Lb <- exp(log(it)-1.96*sig/norm_it_wt)
		Lb <- Lb[1:length(Ub)]
		df <- data.frame(Model=names(M)[i],factor(M[[i]]$d3_survey_data[,1]),M[[i]]$d3_survey_data[,2:8],
						norm_it_wt,it,it_hat,Ub,Lb)
		colnames(df) <- c("Model","Year","Index","Gear","Area","Group","Sex","wt","timing", "norm_wt"
		                  ,"It","It_hat","low_bnd","high_bnd")
		mdf <- rbind(mdf,df)
	}
	print(head(mdf,3))

	limits <- aes(ymax=high_bnd, ymin=low_bnd)

	p <- ggplot(mdf) + geom_point(aes(Year,Index,color=factor(Gear)),size=2)
	p <- p + geom_line(aes(Year,It_hat,color=factor(Gear)),size=1)
	p <- p + geom_pointrange(aes(Year,It,ymax=high_bnd, ymin=low_bnd,color=factor(Gear)))
	# p <- p + geom_ribbon(aes(Year,ymax=high_bnd, ymin=low_bnd,fill=factor(Gear)),alpha=.3)
	p <- p + labs(x="Year",y="Relative abundance",linetype="Gear")
	p <- p + facet_wrap(~Model,scales="free")
	print(p + .THEME)
}


# Steven Martell
# Sept 6,  2012

# .plotOldSurveyfit	<- function( repObj, annotate=FALSE)
# {
# 	with(repObj, {
# 		if(is.matrix(it)){
# 			xx = t(iyr)
# 			m = apply(it,1,max, na.rm=T)
# 			yy = t(pit/m)
# 			y2 = t(it/m)
# 		}else{
# 			xx = iyr
# 			yy = pit
# 			y2 = it
# 		}
# 		n=nrow(t(as.matrix(yy)))
# 		#n=dim(xx)[2]
# 		yrange=c(0, 1.15*max(yy, y2, na.rm=TRUE))
		
# 		matplot(xx, yy, type="n",axes=FALSE,ylim=yrange, 
# 			xlab="Year", ylab="Relative abundance", main=paste(stock))
		
# 		matlines(xx, yy, col=1:n, lty=1)
# 		matpoints(xx, y2, col=1:n, pch=1:n)
		
# 		axis( side=1 )
# 		axis( side=2, las=.VIEWLAS )
# 		box()
# 		grid()
		
# 		if ( annotate )
# 		{
			
# 			txt=paste("Survey ",1:n,", q=",round(q, 3), sep="")
			
# 			mfg <- par( "mfg" )
# 			#if ( mfg[1]==1 && mfg[2]==1 )
# 			legend( "top",legend=txt,
# 				bty='n',lty=1,lwd=1,pch=1:n,ncol=n, col=1:n )
				
# 			#print(q)
# 		}
# 	})
# }










