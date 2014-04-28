#helpers.R
require(ggplot2)
require(reshape2)
load("data/MSE.Rdata")  #Creates a mse.DF object


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
