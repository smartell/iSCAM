#|-----------------------------------------------------------------------------------|#
#| R4iSCAM.R                                                                         |#
#| AUTHOR: Steven Martell                                                            |#
#|-----------------------------------------------------------------------------------|#
#| NOTES, if using sublime text 3.0
#|    use ⌘ + \ to set working directory to the directory of the current file.
#|    use ⌘ + enter to source a single line or highligted lines.
#|    use ⌘ + B to source entire file.

# |----------------------------------------------------------------------------------|
# | DEFINITIONS
# |----------------------------------------------------------------------------------|
# | .PWD       <- Global Parent Working Directory for R-scripts
# | .FIGUREDIR <- Directory for saving figures.
# | .RFILES    <- List of R functions to source from the lib directory.
# .PWD        <- "/Users/stevenmartell1/Documents/iSCAM/examples/PacificHake14/R"
.PWD        <- "/Users/stevenmartell1/Documents/iSCAM-project/src/R"
.LIB        <- "../../../src/R/lib/"
.WIN        <- "../../../src/R/iScamWin2.txt"
setwd(.PWD)
.FIGUREDIR  <- "../FIGS/"
.RFILES     <- list.files(.LIB,pattern="\\.[Rr]$")
.VIEWTRCK   <- "iscamViewTracker.txt"
.BOOLREADFN <- TRUE
.OVERLAY    <- FALSE
require(ggplot2)
.THEME      <- theme_bw(11)
.UNITS      <- "(mlb)"

# | Labels for gear sex area and group index.
.GEAR  = c("Directed","Wasteage","Bycatch","Sport","Personal","Setline Survey")
.SEX   = c("F & M","Female","Male")
.AREA  = c("Coast wide")
.GROUP = c("Pacific Halibut")


.MODELDIRS   <- "../DATA"
.MODELNAME   <- list.files(.MODELDIRS,pattern="\\.RData",full.name=TRUE)
load(.MODELNAME)
names(M)     <- strsplit(basename(.MODELNAME),".RData")

for(nm in .RFILES) source(file.path(.LIB, nm), echo=FALSE)
.plotCatch( M )
.plotIndex( M )
.plotWeightAtAge( M )
.plotSpawnBiomass( M )
.plotDepletion( M )
.plotMortality( M )
.plotRecruitment( M )
.plotStockRecruit( M )
.plotRecruitsPerSpawner( M )
.plotSurveyFit( M )
.plotQ( M )
.plotCatchResidual( M )
.plotIndexResidual( M )
.plotRecruitmentResidual( M )
.plotAgeComps( M )
.plotAgeCompResiduals( M )
.plotAgeSummary( M )
.plotAgeBars( M )
.plotSlx( M )


# |----------------------------------------------------------------------------------|
# | guiView: Main function for starting the iSCAM Gui
# |----------------------------------------------------------------------------------|
# |
guiView  <- function()
 {
	setwd(.PWD)
	for(nm in .RFILES) source(file.path(.LIB, nm), echo=FALSE)
	.gui2("iSCAMView") 	
 }


# |----------------------------------------------------------------------------------|
# | .gui2: creates gui interface and reads .VIEWTRCK file
# |----------------------------------------------------------------------------------|
# | win    <- name of X11 window to close if one already exists.
# | [ ]    <- TODO, if .VIEWTRCK does not exist write header for file.
.gui2	<- function(win)
{
	require(PBSmodelling)
	
	trckExists <- file.exists( .VIEWTRCK )
	ifiles     <- NULL
	if (trckExists)
	{
		cat( "MSG (.hCamViewSetup): Viewer tracking file ",.VIEWTRCK, " found.\n" )
		ifiles <- read.table( file = .VIEWTRCK,  as.is=TRUE,  header=TRUE,  sep="," )
		# ifiles <- tmp
		print(ifiles)
		# browser()
	}
	else
	{
		cat( .VIEWTRCK," does not exist, please check file name. ")
	}
	if(exists(".PBSmod")) closeWin(win)
		
	createWin(.WIN)
	
	setWinVal(list(txtFigDir=.FIGUREDIR))
}


# |----------------------------------------------------------------------------------|
# | .saveMSEdataframe: Loads the rda files and creates a data.frame object for shiny.
# |----------------------------------------------------------------------------------|
# | 
.saveMSEdataframe <- function()
{
	.MSELIB   <- "../DATA/"
	.rdaFILES <- paste0(.MSELIB,list.files(.MSELIB,pattern=".Rdata",recursive=TRUE))
	.bn       <- basename(.rdaFILES)

	loadObj   <- function(fn){load(fn);return(sims)}
	S         <- lapply(.rdaFILES,loadObj)
	names(S)  <- substr(.bn,1,nchar(.bn)-6)
	lbl       <- strsplit(names(S),"_")
	quant     <- c(0.025,0.05,0.25,0.5,0.75,0.90,0.95,0.975)
	n         <- length(S)
	mse.DF    <- NULL
	for( i in 1:n)
	{
		#Scenario & Procedure from hdr Label
		cat("Scenario = ",lbl[[i]][2],"\t")
		cat("Procedure = ",lbl[[i]][3],"\n")

		Year  <- S[[i]][[1]]$yr
		nyr   <- length(Year)

		# Predicted Spawning biomass from iSCAM
		fn    <- function(X) return(X$sbt[1:nyr])
		sbt   <- sapply(S[[i]],fn)
		p_sbt <- as.data.frame(t(apply(sbt,1,quantile,probs=quant)))
		colnames(p_sbt) <- paste0("p.Bt",quant)

		# Reference Model spawning biomass from Milka
		fn    <- function(X) return(X$m_sbt[1:nyr])
		sbt   <- sapply(S[[i]],fn)
		m_sbt <- as.data.frame(t(apply(sbt,1,quantile,probs=quant)))
		colnames(m_sbt) <- paste0("t.Bt",quant)

		# Predicted catch from iSCAM
		fn    <- function(X) return(X$ct)
		ct    <- sapply(S[[i]],fn)
		p_ct  <- as.data.frame(t(apply(ct,1,quantile,probs=quant)))
		colnames(p_ct) <- paste0("Ct",quant)



		# DATA FRAME to concatenate
		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 Year      = Year,
		                 p_sbt, p_ct,
		                 m_sbt)

		mse.DF <- rbind(mse.DF,df)

	}
	save(mse.DF,file=paste0(.MSELIB,"MSE.Rdata"))
}

cat("Type: \n guiView()\n to start the iSCAM gui\n")