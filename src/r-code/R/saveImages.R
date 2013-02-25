# Steven Martell

.saveImages <- function()
{
	# Get list of files need for function argument.
	guiInfo <- getWinVal(scope="L")
	hdr     <- ifiles[ifiles$Select, ]
	nSelect <- nrow(hdr)
	
	# Figure Directory
	.FIGUREDIR <<-  txtFigDir
	
	for( i in 1:nSelect )
	{
		# Obtain the report Object for plotting
		repObj        <- getRepObj(hdr$Control.File[i])
		repObj$stock  <- hdr$Stock[i]
		if(chkPNG)
		{
			cat("Saving PNG files to ", .FIGUREDIR, "for stock", repObj$stock, "\n	")
			.pngPlots( repObj )
		}
	}
	
	graphics.off()
}

.pngPlots <- function( repObj )
{
	# This function quickly builds png files for all plots and saves them 
	# in the .FIGUREDIR folder.
	# rattle off the plots here.
	# FILE NAME FORMAT:
	#                 iSCAMfig:Stock:Catch.png
	fileprefix <- paste(.FIGUREDIR,"iSCAMfig:",repObj$stock, ":", sep="")
	H = 1980
	W = 1980
	Res = 200
	#gfn     <- paste(.FIGUREDIR, "fig:",repObj$stock, "iSCAM%d.png", sep="")
	
	#Input data
	gfn  <- paste(fileprefix, "Catch.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotCatch    ( repObj, legend.txt=NULL )
	dev.off()
	
	gfn  <- paste(fileprefix, "Index.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotIndex    ( repObj, annotate=TRUE )
	dev.off()
	
	gfn  <- paste(fileprefix, "AgeComps%d.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotAgecomps    ( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "MeanWt.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotMeanwt    ( repObj )
	dev.off()
	
	# Residuals
	gfn  <- paste(fileprefix, "CatchResidual.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotCatchResiduals ( repObj, annotate=TRUE )
	dev.off()
	
	gfn  <- paste(fileprefix, "SurveyResidual.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotSurveyResiduals ( repObj, annotate=TRUE )
	dev.off()
	
	gfn  <- paste(fileprefix, "RecruitmentResidual.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotRecruitmentResiduals ( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "AgeCompResidual%d.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotAgeCompResiduals ( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "Biomass.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotBiomass ( repObj, annotate=TRUE )
	dev.off()
	

	gfn  <- paste(fileprefix, "RetrospectiveBiomass.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotSSBretrospective ( repObj )
	dev.off()

	gfn  <- paste(fileprefix, "Depletion.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotDepletion ( repObj, annotate=TRUE )
	dev.off()
	
	gfn  <- paste(fileprefix, "Recruitment.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotRecruitment ( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "StockRecruitment.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotStockRecruit ( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "Mortality.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotMortality ( repObj, annotate=TRUE )
	dev.off()
	
	gfn  <- paste(fileprefix, "SurveyFit.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotSurveyfit ( repObj, annotate=TRUE )
	dev.off()
	
	gfn  <- paste(fileprefix, "Selectivity%d.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotSelectivity ( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "StockStatus.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotStockStatus ( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "MarginalPosteriors.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotMarginalPosteriors( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "ReferencePoints.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotReferencePoints( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "SbPosterior.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotSbtPosterior( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "SbDepletionPosterior.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotSbtPosterior( repObj, TRUE,  annotate=TRUE )
	dev.off()
	
	gfn  <- paste(fileprefix, "MCMCtrace.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotMCMCtrace( repObj )
	dev.off()
	
	gfn  <- paste(fileprefix, "MCMCpairs.png", sep="")
	png(gfn, width=W, height=H, res=Res)
	.plotMCMCpairs ( repObj )
	dev.off()
	
	
}