# plotAgeBars.R
library(reshape2)
library(ggplot2)

.plotAgeBars <- function( M )
{
	n <- length( M )
	cat(".plogAgeBars\n")
	id <- grep("d3_A[1-9]",names(M[[1]]))
	iu <- grep("A_hat[1-9]",names(M[[1]]))
	hr <- c("Year","Gear","Area","Group","Sex","AgeErr")
	mdf <- NULL
	for(i in 1:n)
	{
		getDF <- function(x)
		{
			ix <- id[x]
			jx <- iu[x]

			age <- seq(M[[i]]$n_A_sage[x],M[[i]]$n_A_nage[x])
			df  <- data.frame(V="Observed",M[[i]][[ix]])
			df[,-1:-6] <- df[,-1:-6]/rowSums(df[,-1:-6],na.rm=TRUE)
			colnames(df) = c("V",hr,paste(age))

			dp  <- data.frame(V="Predicted",M[[i]][[ix]][,1:6],M[[i]][[jx]])
			dp[,-1:-6] <- dp[,-1:-6]/rowSums(dp[,-1:-6],na.rm=TRUE)
			colnames(dp) = c("V",hr,paste(age))
			
			return(rbind(df,dp))
		}

		B  <- lapply(1:length(id),getDF)
	}
	
	mB <- melt(B,id.vars=c("V",hr))
	mB$BroodYear <- mB$Year - as.integer(mB$variable)

	O <- subset(mB,V=="Observed")
	I <- subset(mB,V=="Predicted")

	p <- ggplot(O,aes(as.numeric(variable),as.double(value),fill=BroodYear))
	p <- p + geom_bar(stat="identity")
	p <- p + geom_line(data=I,aes(as.numeric(variable),value),color="red")
	p <- p + labs(x="Age",y="Proportion")
	p <- p + guides(title.position = "top")
	p <- p + facet_wrap(~Year,nrow=4)+coord_flip()
	print(p + .THEME + theme(legend.position="top") + theme_bw(8))
}