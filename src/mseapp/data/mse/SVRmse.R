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
	}
	if(input$siPlotType == "fct")
	{
		df <- mse.data[["catch.df"]]
		icol <- c(icol,"ct50","ct025","ct975")
	}



	# CONTROLLER
	# Filter df 
	sdf <- subset(df,Scenario %in% input$siScenario 
	              & Procedure %in% input$siProcedure
	              & Recruitment %in% input$siRecruitment)

	

	idf  <- sdf[,which(names(sdf) %in% icol)]
	print(head(idf))
	n    <- dim(idf)[2] - 2
	colnames(idf)[n:(n+2)] <- c("lci","Median","uci")

	print(tail(idf))
	# VIEW
	ci <- aes(ymin=lci,ymax=uci,fill=Procedure)
	# co <- aes(ymin=0.8*lci,ymax=1.2*uci,fill=Procedure)
	p  <- ggplot(idf,aes(Year,Median,color=Procedure))
	p  <- p + geom_line()
	p  <- p + geom_ribbon(ci,alpha=0.25,linetype=0)
	# p  <- p + geom_ribbon(co,alpha=0.25,linetype=0)
	p  <- p + labs(x="Year",y=input$siPlotType)
	print(p + theme_classic(14) + facet_wrap(~Recruitment))



}