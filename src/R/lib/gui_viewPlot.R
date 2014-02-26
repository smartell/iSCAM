# |----------------------------------------------------------------------------------|
# | gui_viewPlot.R
# |----------------------------------------------------------------------------------|
# | - main function is .viewPlot.R   
# | - .getRunObject is the main function that reads the iSCAM output files and 
# |   creates an list object from the particular model run.
# |
# | - gui_viewPlot.R list of functions:
# |
# | .viewPlot     : determine selected models and calls appropriate plotting routine.
# | .getRunObject : reads selected files and saves an .RData file of the object.
# | .readFile     : change the status of .BOOLREADFN based on gui hdr$Select
# |


# |----------------------------------------------------------------------------------|
# | .viewPlot
# |----------------------------------------------------------------------------------|
# | - determine which models have been selected for plotting and call the appropriate
# |   routine to do the plot.

.viewPlot <- function()
{
	guiInfo <- getWinVal(scope="L")
	cat("\n .viewPlot, plotType =",plotType,"\n")

	# | Determine which files have been selected
	hdr	  <- ifiles[ ifiles$Select, ]
	fn    <- hdr$Control.File
	nRuns <- nrow(hdr)

	# | Read objects & make global
	if(.BOOLREADFN)
	{
		M           <<- lapply(fn,.getRunObject)
		names(M)    <<- hdr$Label
		.BOOLREADFN <<- FALSE
		assign(hdr$Label,M,pos=1)		
	}

	if( plotType=="catch" )
	{
		.plotCatch( M )
	}
	else if( plotType=="survey" )
	{
		.plotIndex( M )
	}
	else if( plotType=="agecomps" )
	{
		.plotAgeComps( M )
	}
	else if( plotType=="weightatage" )
	{
		.plotWeightAtAge( M )
	}
	else if( plotType=="agecompsresid" )
	{
		.plotAgeCompResiduals( M )
	}
	else if( plotType=="agesummary" )
	{
		.plotAgeSummary( M )
	}
	else if( plotType=="catchresid" )
	{
		.plotCatchResidual( M )
	}
	else if( plotType=="surveyresid" )
	{
		.plotIndexResidual( M )
	}
	else if( plotType=="recresid" )
	{
		.plotRecruitmentResidual( M )
	}
	else if( plotType=="spawnbio" )
	{
		.plotSpawnBiomass( M )
	}
	else if( plotType=="mortality" )
	{
		.plotMortality( M )
	}

}

# |----------------------------------------------------------------------------------|
# | .getRunObject
# |----------------------------------------------------------------------------------|
# | - This function takes the control file name (fn) as an argument and reads the
# |   report file, par file, mcmc files, etc. and returns a single list object.

.getRunObject <- function(fn)
{
	cat("\n Reading",fn,"\n")
	tmp      <- read.admb(fn)

	# | simulation file
	simfile  <- paste(fn,".sim",sep="")
	if(file.exists(simfile))
	{
		tmp$sim      <- read.rep(simfile)
	}
	# | retrospective results
	retfile <- paste(fn,".ret",1:10,sep="")
	retSbt <- list()
	i <- 0
	for(ifn in retfile)
	{
		if(file.exists(ifn))
		{
			i <- i + 1
			sbt <- read.rep(ifn)$sbt
			nn  <- 1:(length(sbt)-1)
			retSbt[[i]] <- sbt[nn]
		}
	}
	tmp$retSbt <- retSbt

	save(tmp,file=paste(fn,".RData",sep=""))
	return(tmp)
}

# |----------------------------------------------------------------------------------|
# | .readFile
# |----------------------------------------------------------------------------------|
# | - function to set global variable .BOOLREADFN to true if user changes gui Select
.readFile <- function()
{
	.BOOLREADFN <<- TRUE
}
