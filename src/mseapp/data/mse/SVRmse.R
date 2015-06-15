# SVRmse.R

# -------------------------------------------- #
# funnelPlot
# Pseudocode:
#   • MODEL
# 	- Using the input argument:
# 		- determine which data elements are needed
#		- subset the mse.data object as needed.
#
#   • VIEW
#	- Using the subsetDF :
#		- create ggplot object
#       - 
# 
# 	• CONTROLLER
# 	- Subset mse.data for appropriat list. DF
#	- subset DF for specific Year, Scenario, Procedures
# 	- determine which columns to put in ggplot2
# -------------------------------------------- #

.funnelPlot <- function(input)
{
	print(input$siScenario)
	print(input$siProcedure)

	# MODEL
	icol <- c("Scenario","Procedure","Recruitment","SizeAtAge","Year","gear","area","sex","group")
	# Subset mse.data by siPlotType
	if(input$siPlotType == "ssb")
	{
		df   <- mse.data[["biomass.df"]]
		icol <- c(icol,"t.Bt0.5","t.Bt0.025","t.Bt0.975")
		wrap <- facet_wrap(~Recruitment+SizeAtAge)
	}
	if(input$siPlotType == "fct")
	{
		df <- mse.data[["catch.df"]]
		icol <- c(icol,"ct50","ct025","ct975")
		wrap <- facet_wrap(~Recruitment+SizeAtAge+gear+area+sex)
	}
	if(input$siPlotType == "sbd")
	{
		df <- mse.data[["biomass.df"]]
		icol <- c(icol,"t.Dt0.5","t.Dt0.025","t.Dt0.975")
		wrap <- facet_wrap(~Recruitment+SizeAtAge)
	}



	# CONTROLLER sldr_year_range
	# Year      %in% input$years[1]:input$years[2] &
	# Filter df 
	sdf <- subset(df,
	                Year     %in% input$sldr_year_range[1]:input$sldr_year_range[2]
	              & Scenario %in% input$siScenario 
	              & Procedure %in% input$siProcedure
	              & Recruitment %in% input$siRecruitment
	              & SizeAtAge  %in% input$siGrowth)

	

	idf  <- sdf[,which(names(sdf) %in% icol)]
	print(head(idf))
	n    <- dim(idf)[2] - 2
	colnames(idf)[n:(n+2)] <- c("lci","Median","uci")

	print(tail(idf))
	# VIEW
	ci <- aes(ymin=lci,ymax=uci,fill=Procedure)

	p  <- ggplot(idf,aes(Year,Median,color=Procedure))
	p  <- p + geom_line()
	p  <- p + geom_ribbon(ci,alpha=0.25,linetype=0)
	p  <- p + labs(x="Year",y=input$siPlotType)
	print(p + .THEME + wrap)



}

.tablePeformanceMetric <- function(input,var = "t.Dt0.5")
{
	cat("Performance Metric table for ", var, "\n")
	if(var == "t.Dt0.5")     idf = "biomass.df"
	if(var == "t.Bt0.5")     idf = "biomass.df"
	if(var == "P.SSB.0.20.") idf = "biomass.df"
	if(var == "P.SSB.0.30.") idf = "biomass.df"
	if(var == "ct50"   )     idf = "catch.df" 
	if(var == "AAV50"  )     idf = "AAV.df"

	msedf <- mse.data[[idf]]
	df  <- subset(msedf,
            	  Year        %in% input$sldr_year_range[1]:input$sldr_year_range[2] 
                & Scenario    %in% input$siScenario                
                & Procedure   %in% input$siProcedure 
                & Recruitment %in% input$siRecruitment
	            & SizeAtAge   %in% input$siGrowth
             )
	mdf <- melt(df,id=c("Scenario","Procedure","Year","Recruitment","SizeAtAge"))
    tmp <- dcast(mdf,Procedure~Scenario,mean,na.rm=TRUE,margins="Scenario"
                 ,subset = .(variable==var))
    return(tmp)	
}


  