#helpers.R
require(ggplot2)
require(reshape2)
require(googleVis)
load("data/MSE.Rdata")  #Creates a mse.data object
mse.DF <- mse.data$biomass.df


funnel.plot <- function(df,input)
{
	S <- input$scenario
	P <- input$procedure
	
	if( input$plotType=='Spawning biomass' )
	{
		ci <- aes(ymin=t.Bt0.025,ymax=t.Bt0.975)
		p  <- ggplot(df,aes(Year,t.Bt0.5))+geom_line()
		p  <- p + geom_ribbon(ci,alpha=0.15)
		p  <- p + labs(x="Year",y="Spawning biomass")
		
	}
	
	if( input$plotType=='Catch' )
	{
		ci <- aes(ymin=ct025,ymax=ct975)
		p  <- ggplot(df,aes(Year,ct50))+geom_line()
		p  <- p + geom_ribbon(ci,alpha=0.15)
		p  <- p + labs(x="Year",y="Catch")
		
	}


	
	facets <- paste("Procedure",'~',"Scenario")
	if( length(P) > 1 )
	{
		p <- p + aes_string(color="Procedure",fill="Procedure")
	}
	
	if( length(P)==1 )
	{
		p <- p + aes_string(color="Scenario",fill="Scenario")
	}

	.LEGPOS <- 'top'
	print( p + theme( legend.position = .LEGPOS ) )
}


motionChart <- function(df,input)
{
	print(df)
	M1 <- gvisMotionChart(df,idvar="Procedure",timevar="Year")
	return(M1)
}