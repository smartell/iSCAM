#helpers.R
require(ggplot2)
require(reshape2)
require(googleVis)
require(plyr)
load("data/MSE.Rdata")  #Creates a mse.data object
mse.DF <- mse.data$biomass.df


funnel.plot <- function(df,input)
{
	S <- input$scenario
	P <- input$procedure
	
	if( input$plotType=='Spawning biomass' )
	{
		ci95 <- aes(ymin=t.Bt0.025,ymax=t.Bt0.975)
		ci25 <- aes(ymin=t.Bt0.25,ymax=t.Bt0.75)
		p  <- ggplot(df,aes(Year,t.Bt0.5))+geom_line()
		p  <- p + geom_ribbon(ci95,alpha=0.15)
		p  <- p + geom_ribbon(ci25,alpha=0.25)
		p  <- p + labs(x="Year",y="Spawning biomass")
		
	}

	if( input$plotType=='Depletion' )
	{
		ci95 <- aes(ymin=t.Dt0.025,ymax=t.Dt0.975)
		ci25 <- aes(ymin=t.Dt0.25,ymax=t.Dt0.75)
		p  <- ggplot(df,aes(Year,t.Dt0.5))+geom_line()
		p  <- p + geom_ribbon(ci95,alpha=0.15)
		p  <- p + geom_ribbon(ci25,alpha=0.25)
		p  <- p + labs(x="Year",y="Spawning biomass depletion")
		
	}

	if( input$plotType=='Catch' )
	{
		ci <- aes(ymin=ct025,ymax=ct975)
		p  <- ggplot(df,aes(Year,ct50))+geom_line()
		p  <- p + geom_ribbon(ci,alpha=0.15)
		p  <- p + labs(x="Year",y="Catch")
		
	}

	if( input$plotType=='Sub-legal Catch' )
	{
		ci <- aes(ymin=dt025,ymax=dt975)
		p  <- ggplot(df,aes(Year,dt50))+geom_line()
		p  <- p + geom_ribbon(ci,alpha=0.15)
		p  <- p + labs(x="Year",y="Sub-legal Catch")
	}

	if( input$plotType=="AAV in Catch")
	{
		print(colnames(df))
		ci <- aes(ymin=AAV025,ymax=AAV975)
		p  <- ggplot(df,aes(Year,AAV50,color=factor(gear),fill=factor(gear)))+geom_line()
		p  <- p + geom_ribbon(ci,alpha=0.15)
		p  <- p + labs(x="Year",y="AAV in Catch")
		p  <- p + theme_bw(11) 
	}
	
	facets <- paste("Procedure",'~',"Scenario")
	# if( length(P) > 1 )
	# {
	# 	p <- p + aes_string(color="Procedure",fill="Procedure")
	# }
	
	# if( length(P)==1 )
	# {
	# 	p <- p + aes_string(color="Scenario",fill="Scenario")
	# }

	.LEGPOS <- 'top'
	print( p + theme( legend.position = .LEGPOS ) )
}


motionChart <- function(df,input)
{
	print(df)
	M1 <- gvisMotionChart(df,idvar="Procedure",timevar="Year")
	return(M1)
}