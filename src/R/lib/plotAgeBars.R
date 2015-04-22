# plotAgeBars.R
library(reshape2)
library(ggplot2)


.plotAgeBars <- function( M )
{
	n <- length( M )
	cat(".plotAgeBars\n")
	id <- grep("d3_A[1-9]",names(M[[1]]))
	iu <- grep("A_hat[1-9]",names(M[[1]]))
	hr <- c("Year","Gear","Area","Group","Sex","AgeErr")
	mdf <- NULL
	mcol<- -1:-7
	for(i in 1:n)
	{
		getDF <- function(x)
		{
			ix <- id[x]
			jx <- iu[x]

			age <- seq(M[[i]]$n_A_sage[x],M[[i]]$n_A_nage[x])
			df  <- data.frame(V="Observed",M[[i]][[ix]])
			df[,mcol] <- df[,mcol]/rowSums(df[,mcol],na.rm=TRUE)
			# cat("\nSum of observed = ",rowSums(df[,-1:-6]))
			colnames(df) = c("V",hr,paste(age))

			dp  <- data.frame(V="Predicted",M[[i]][[ix]][,1:6],M[[i]][[jx]])
			dp[,mcol] <- dp[,mcol]/rowSums(dp[,mcol],na.rm=TRUE)
			# cat("\nSum of predicted = ",rowSums(dp[,-1:-6]))
			colnames(dp) = c("V",hr,paste(age))
			
			dfp <- rbind(df,dp)
			dfp <- subset(dfp,Year %in% seq(M[[i]]$syr,M[[i]]$nyr))
			return(dfp)
		}

		B  <- lapply(1:length(id),getDF)
	}
	


	barplot <- function(B)
	{		
		mB <- melt(B,id.vars=c("V",hr))
		mB$BroodYear <- mB$Year - as.integer(mB$variable)
		

		O <- subset(mB,V=="Observed")
		I <- subset(mB,V=="Predicted")

		klbl <- .GEAR[O$Gear]
		albl <- .AREA[O$Area]
		glbl <- .GROUP[O$Group]
		slbl <- .SEX[O$Sex+1]

		p <- ggplot(O,aes(variable,as.double(value),fill=BroodYear))
		p <- p + geom_bar(stat="identity",alpha=0.7)
		p <- p + scale_fill_distiller(type = "div")
		p <- p + geom_line(data=I,aes(as.numeric(variable),value),color="red")
		p <- p + labs(x="Age",y="Proportion")
		p <- p + ggtitle(paste(klbl,":Area",albl,":Group",glbl,":Sex",slbl))
		p <- p + guides(title.position = "top")
		p <- p + facet_wrap(~Year)+coord_flip()
		return (p + .THEME + theme(legend.position="top") + theme_bw(8))
	}

	P <- lapply(B,barplot)
	print(P)
}