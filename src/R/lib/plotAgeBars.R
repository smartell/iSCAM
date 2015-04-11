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
			df   <- data.frame(M[[i]][[ix]])
			colnames(df) = c(hr,paste(age))
			
			
			return(df)
		}

		B  <- lapply(1:length(id),getDF)
	}
	
	mB <- melt(B,id.vars=hr)
	mB$BroodYear <- mB$Year - as.integer(mB$variable)

	p <- ggplot(subset(mB,L1==5),aes(as.integer(variable),value,fill=factor(BroodYear)))
	p <- p + geom_bar(stat="identity")
	p <- p + labs(x="Age",y="Proportion",fill=NULL)
	p <- p + guides(fill = guide_legend(nrow = 4,keyheight=0.5))
	p <- p + facet_wrap(~Year)+coord_flip()
	print(p + .THEME + theme(legend.position="top"))
}