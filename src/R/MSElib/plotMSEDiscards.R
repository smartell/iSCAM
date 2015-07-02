#Catarina Wor
# Jun 25 2015
#parei aqui
require(reshape)

.MSEplotWastage <- function(M, no=3, ci=95)
{
	df <- as.data.frame(M[[1]]$sublegal.df)
	df2<- as.data.frame(M[[2]]$rawsublegal.df)
	
	if(ci==95){up=df$wt975;low=df$wt025}
	if(ci==90){up=df$wt95;low=df$wt05}
	if(ci==50){up=df$wt75;low=df$wt25}

	cat(".plotWastage\n")

	df<-data.frame(cbind(df,up,low))


	
	if(.OVERLAY)
	{
		
		p <- ggplot(df,aes(Year,wt50,col=gear)) + geom_line(width=2)
		p <- p +  geom_ribbon(aes(ymax=up,ymin=low,fill=Scenario),alpha=0.2)
		p <- p + facet_grid(.~Procedure)

		if(no>0)
			{

				itx <- sample(1:length(grep("wt",colnames(df2))),no)
				itxx <- grep("wt",colnames(df2))[itx]
				new.df <- df2[,c(1:4,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year", "gear"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',col='gear',linetype='variable'))
				p <- p +   scale_linetype_manual(values=c("dashed","dotted","dotdash","longdash","twodash","1F","F1","4C88C488"))

			}
	}
	else
	{
		p <- ggplot(df,aes(Year,wt50),col=gear,fill=gear) + geom_line(width=2)
		p <- p + geom_ribbon(aes(ymax=up, ymin=low,col=gear,fill=gear),alpha=0.2)
		p <- p + facet_grid(Scenario~Procedure,scales="free")

		if(no>0)
			{

				itx <- sample(1:length(grep("wt",colnames(df2))),no)
				itxx <- grep("wt",colnames(df2))[itx]
				new.df <- df2[,c(1:4,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year", "gear"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',col='gear',linetype='variable'))
				p <- p +   scale_linetype_manual(values=c("dashed","dotted","dotdash","longdash","twodash","1F","F1","4C88C488"))

			}
	
	p <- p + labs(x="Year",y=paste("Wastage",.UNITS))
	print(p + .THEME)
	}
}

.MSEplotSubLegal <- function(M, no=3, ci=95)
{
	df <- as.data.frame(M[[1]]$sublegal.df)
	df2<- as.data.frame(M[[2]]$rawsublegal.df)
	
	if(ci==95){up=df$dt975;low=df$dt025}
	if(ci==90){up=df$dt95;low=df$dt05}
	if(ci==50){up=df$dt75;low=df$dt25}

	df<-data.frame(cbind(df,up,low))


	cat(".plotSubLegalCatch\n")

	
	if(.OVERLAY)
	{
		
		p <- ggplot(df,aes(Year,dt50,col=gear)) + geom_line(width=2)
		p <- p +  geom_ribbon(aes(ymax=up,ymin=low,fill=Scenario),alpha=0.2)
		p <- p + facet_grid(.~Procedure)

		if(no>0)
			{

				itx <- sample(1:length(grep("dt",colnames(df2))),no)
				itxx <- grep("dt",colnames(df2))[itx]
				new.df <- df2[,c(1:4,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year", "gear"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',col='gear',linetype='variable'))
				p <- p +   scale_linetype_manual(values=c("dashed","dotted","dotdash","longdash","twodash","1F","F1","4C88C488"))

			}
	}
	else
	{
		p <- ggplot(df,aes(Year,dt50),col=gear,fill=gear) + geom_line(width=2)
		p <- p + geom_ribbon(aes(ymax=up, ymin=low,col=gear,fill=gear),alpha=0.2)
		p <- p + facet_grid(Scenario~Procedure,scales="free")

		if(no>0)
			{

				itx <- sample(1:length(grep("dt",colnames(df2))),no)
				itxx <- grep("dt",colnames(df2))[itx]
				new.df <- df2[,c(1:4,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year", "gear"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',col='gear',linetype='variable'))
				p <- p +   scale_linetype_manual(values=c("dashed","dotted","dotdash","longdash","twodash","1F","F1","4C88C488"))

			}
	
	p <- p + labs(x="Year",y=paste("SubLegal Discards",.UNITS))
	print(p + .THEME)
	}
}

.MSEplotEfficiency <- function(M, no=3, ci=95)
{
	df <- as.data.frame(M[[1]]$sublegal.df)
	df2<- as.data.frame(M[[2]]$rawsublegal.df)
	
	if(ci==95){up=df$ef975;low=df$ef025}
	if(ci==90){up=df$ef95;low=df$ef05}
	if(ci==50){up=df$ef75;low=df$ef25}

	df<-data.frame(cbind(df,up,low))


	cat(".plotEfficiency\n")

	
	if(.OVERLAY)
	{
		
		p <- ggplot(df,aes(Year,ef50,col=gear)) + geom_line(width=2)
		p <- p +  geom_ribbon(aes(ymax=up,ymin=low,fill=Scenario),alpha=0.2)
		p <- p + facet_grid(.~Procedure)

		if(no>0)
			{

				itx <- sample(1:length(grep("ef",colnames(df2))),no)
				itxx <- grep("ef",colnames(df2))[itx]
				new.df <- df2[,c(1:4,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year", "gear"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',col='gear',linetype='variable'))
				p <- p +   scale_linetype_manual(values=c("dashed","dotted","dotdash","longdash","twodash","1F","F1","4C88C488"))

			}
	}
	else
	{
		p <- ggplot(df,aes(Year,ef50),col=gear,fill=gear) + geom_line(width=2)
		p <- p + geom_ribbon(aes(ymax=up, ymin=low,col=gear,fill=gear),alpha=0.2)
		p <- p + facet_grid(Scenario~Procedure,scales="free")

		if(no>0)
			{

				itx <- sample(1:length(grep("ef",colnames(df2))),no)
				itxx <- grep("ef",colnames(df2))[itx]
				new.df <- df2[,c(1:4,itxx)]
				new.df <- melt(new.df,id=c("Scenario","Procedure","Year", "gear"))
				
				p <- p + geom_line(data=new.df,aes_string(x="Year",y='value',col='gear',linetype='variable'))
				p <- p +  scale_linetype_manual(values=c("dashed","dotted","dotdash","longdash","twodash","1F","F1","4C88C488"))

			}
	
	p <- p + labs(x="Year",y=paste("Efficiency",.UNITS))
	print(p + .THEME)
	}
}

