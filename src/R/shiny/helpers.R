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
AAV.DF <- mse.data$AAV.df

MRG.DF <- merge(BIO.DF,CAT.DF,by=c("Scenario","Procedure","Year"))
MRG.DF <- merge(MRG.DF,AAV.DF,by=c("Scenario","Procedure","Year","gear","area","group"))
MSE.DF <- merge(MRG.DF,SUB.DF,by=c("Scenario","Procedure","Year","gear","area","sex","group"))
# Restricted data frame for gvisMotionChart for increased speed & less clutter.
hdr <- c("Scenario","Procedure","Year","t.Bt0.5","t.Dt0.5","ct50","AAV50")
MOT.DF <- MRG.DF[,which(names(MRG.DF) %in% hdr)]
MOT.DF$idvar <- paste0(MOT.DF$Scenario,MOT.DF$Procedure)
names(MOT.DF) <- c("Scenario","Procedure","Year","Median biomass","Depletion",
                   "Median catch","Median AAV","idvar")


tulip.plot <- function(df,input)
{
	.LEGPOS <- 'bottom'

	print(input$plotType)
	icol <- c("Scenario","Procedure","Year","gear","area","sex","group")
	if(input$plotType=='Spawning biomass')
	{
		icol <- c(icol,"t.Bt0.5","t.Bt0.025","t.Bt0.975")
	}
	if(input$plotType=='Depletion')
	{
		icol <- c(icol,"t.Dt0.5","t.Dt0.025","t.Dt0.975")	
	}
	if(input$plotType=='Catch')
	{
		icol <- c(icol,"ct50","ct025","ct975")
	}
	if(input$plotType=='Sub-legal Catch')
	{
		icol <- c(icol,"dt50","dt025","dt975")	
	}
	if(input$plotType=='Sub-legal Catch')
	{
		icol <- c(icol,"dt50","dt025","dt975")	
	}
	if(input$plotType=='AAV in Catch')
	{
		icol <- c(icol,"AAV50","AAV025","AAV975")
	}
	if(input$plotType=='Wastage')
	{
		icol <- c(icol,"wt50","wt025","wt975")	
	}
	if(input$plotType=='Efficiency')
	{
		icol <- c(icol,"ef50","ef025","ef975")	
	}
	if(input$plotType=='Fishing mortality')
	{
		icol <- c(icol,"Ft0.5","Ft0.025","Ft0.975")
	}
	sdf  <- df[,which(names(df) %in% icol)]
	n    <- dim(sdf)[2] - 2
	colnames(sdf)[n:(n+2)] <- c("lci","Median","uci")

	print(tail(sdf))

	# GGPLOT.
	ci <- aes(ymin=lci,ymax=uci)
	p <- ggplot(sdf,aes(x=Year,y=Median)) + geom_line()
	p <- p + geom_ribbon(ci,alpha=0.25)
	p <- p + labs(x="Year",y=input$plotType)

	if( input$icolor != "." )
	{
		p <- p + aes_string(fill=input$icolor,color=input$icolor)
	}

	facets <- paste(input$facet_row,"~",input$facet_col)
	p      <- p + facet_grid(facets)


	print(p + theme( legend.position = .LEGPOS ) )
	
}

# funnel.plot <- function(df,input)
# {
# 	S <- input$scenario
# 	P <- input$procedure
	
# 	if( input$plotType=='Spawning biomass' )
# 	{
# 		ci95 <- aes(ymin=t.Bt0.025,ymax=t.Bt0.975)
# 		ci25 <- aes(ymin=t.Bt0.25,ymax=t.Bt0.75)
# 		p  <- ggplot(df,aes(Year,t.Bt0.5))+geom_line()
# 		p  <- p + geom_ribbon(ci95,alpha=0.15)
# 		p  <- p + geom_ribbon(ci25,alpha=0.25)
# 		p  <- p + labs(x="Year",y="Spawning biomass")
		
# 	}

# 	if( input$plotType=='Depletion' )
# 	{
# 		ci95 <- aes(ymin=t.Dt0.025,ymax=t.Dt0.975)
# 		ci25 <- aes(ymin=t.Dt0.25,ymax=t.Dt0.75)
# 		p  <- ggplot(df,aes(Year,t.Dt0.5))+geom_line()
# 		p  <- p + geom_ribbon(ci95,alpha=0.15)
# 		p  <- p + geom_ribbon(ci25,alpha=0.25)
# 		p  <- p + labs(x="Year",y="Spawning biomass depletion")
		
# 	}

# 	if( input$plotType=='Catch' )
# 	{
# 		ci <- aes(ymin=ct025,ymax=ct975)
# 		p  <- ggplot(df,aes(Year,ct50,color=factor(gear),fill=factor(gear)))+geom_line()
# 		p  <- p + geom_ribbon(ci,alpha=0.15)
# 		p  <- p + labs(x="Year",y="Catch")
		
# 	}

# 	if( input$plotType=='Sub-legal Catch' )
# 	{
# 		ci <- aes(ymin=dt025,ymax=dt975)
# 		p  <- ggplot(df,aes(Year,dt50))+geom_line()
# 		p  <- p + geom_ribbon(ci,alpha=0.15)
# 		p  <- p + labs(x="Year",y="Sub-legal Catch")
# 	}

# 	if( input$plotType=="AAV in Catch")	
# 	{
# 		ci <- aes(ymin=AAV025,ymax=AAV975)
# 		p  <- ggplot(df,aes(Year,AAV50,color=factor(gear),fill=factor(gear)))+geom_line()
# 		p  <- p + geom_ribbon(ci,alpha=0.15)
# 		p  <- p + labs(x="Year",y="AAV in Catch")
# 	}

# 	if( input$plotType=='Wastage' )
# 	{
# 		ci <- aes(ymin=wt025,ymax=wt975)
# 		p  <- ggplot(df,aes(Year,wt50))+geom_line()
# 		p  <- p + geom_ribbon(ci,alpha=0.15)
# 		p  <- p + labs(x="Year",y="Wastage")
# 	}
	
# 	# facets <- paste("Procedure",'~',"Scenario")
# 	facets <- paste(input$facet_row,"~",input$facet_col)
# 	p <- p + facet_grid(facets)

# 	if( input$icolor != "." )
# 	{
# 		p <- p + aes_string(fill=factor(input$icolor))
# 	}
	
# 	# 95% CI
# 	# p <- p + geom_ribbon(ci,alpha=0.15)

	

# 	# if( length(P) > 1 )
# 	# {
# 	# 	p <- p + aes_string(color="Procedure",fill="Procedure")
# 	# }
	
# 	# if( length(P)==1 )
# 	# {
# 	# 	p <- p + aes_string(color="Scenario",fill="Scenario")
# 	# }

# 	.LEGPOS <- 'top'

# 	print( p + theme( legend.position = .LEGPOS ) )
# }


motionChart <- function(df,input)
{
	M1 <- gvisMotionChart(df,idvar="idvar",timevar="Year",
	                      xvar="Depletion",yvar="Median catch",
	                      sizevar="Median AAV",colorvar="Procedure",
	      options=list(showChartButtons = TRUE))
	return(M1)
}

