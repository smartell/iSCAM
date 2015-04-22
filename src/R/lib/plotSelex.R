# plotSelex.R
library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)

.plotSlx <- function( M )
{
	n <- length(M)
	cat("plotSlx\n")
	mdf <- NULL
	decade <- seq(1900,2050,by=10)
	for(i in 1:n)
	{
		df  <- data.frame(Model=names(M)[i],logSel=M[[i]]$log_sel)
		colnames(df) = c("Model","Gear","Area","Group","Sex","Year",M[[i]]$age)
		df$Decade <- decade[findInterval(df$Year,decade)] 
		mdf <- rbind(mdf,melt(df,id=c("Model","Gear","Area","Group","Sex","Year","Decade")))
	}

	f <- unique(mdf$Area)
	g <- unique(mdf$Group)
	h <- unique(mdf$Sex)
	k <- unique(mdf$Gear)
	d <- unique(mdf$Decade)
	ags <- as.matrix(expand.grid(k,f,g))

	fnp <- function(ff)
	{
		dfx <- mdf %>% filter(Gear==ff[1]&Area==ff[2]&Group==ff[3])

		p  <- ggplot(dfx,aes(as.numeric(variable),exp(value),col=factor(Year)))
		p  <- p + geom_line()
		p  <- p + labs(x="Age",y="Relative Selectivity")
		p  <- p + ggtitle(paste("Gear",.GEAR[ff[1]]))
		p  <- p + facet_wrap(~Decade+Sex,scales="free_y") 
		return(p + .THEME + theme_bw(8) + guides(col=FALSE))
	}
	
	P <- apply(ags,1,fnp)
	print(P)


}
