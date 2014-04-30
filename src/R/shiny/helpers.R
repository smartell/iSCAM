#helpers.R
# LIBRARIES
require(ggplot2)
require(reshape2)
require(googleVis)
require(plyr)

# LOAD A mse.data DATA OBJECT
load("data/MSE.Rdata")  

# ——————————————————————————————————————————————————————————————————————————— #
# NOTES: The mse.data object is a list of data.frames that comes from the
#        saveMSEdataframe.R routine.  The following code disassbles this list
#        into several dataframes to be used in the Shiny application.
# DATA FRAMES:
#   - BIO.DF -> spawning biomass 
#   - CAT.DF -> catch related variables
#   - SUB.DF -> sublegal and wastage related variables
#   - MOT.DF -> data frame to be used with gvisMotionChart
# ——————————————————————————————————————————————————————————————————————————— #
BIO.DF <- mse.data$biomass.df
CAT.DF <- mse.data$catch.df
SUB.DF <- mse.data$sublegal.df

MRG.DF <- merge(BIO.DF,CAT.DF,by=c("Scenario","Procedure","Year"))
# Restricted data frame for gvisMotionChart for increased speed & less clutter.
hdr <- c("Scenario","Procedure","Year","t.Bt0.5","t.Dt0.5","ct50")
MOT.DF <- MRG.DF[,which(names(MRG.DF) %in% hdr)]


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

	if( input$plotType=='Wasteage' )
	{
		ci <- aes(ymin=wt025,ymax=wt975)
		p  <- ggplot(df,aes(Year,wt50))+geom_line()
		p  <- p + geom_ribbon(ci,alpha=0.15)
		p  <- p + labs(x="Year",y="Wasteage")
	}

	
	# facets <- paste("Procedure",'~',"Scenario")
	facets <- paste(input$facet_row,"~",input$facet_col)
	p <- p + facet_grid(facets)

	if( input$icolor != "." )
	{
		p <- p + aes_string(color=input$icolor,fill=input$icolor)
	}
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
	M1 <- gvisMotionChart(df,idvar="Procedure",timevar="Year",
	      options=list(showChartButtons = TRUE))
	return(M1)
}