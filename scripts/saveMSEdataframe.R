# |----------------------------------------------------------------------------------|
# | .saveMSEdataframe: Loads the rda files and creates a data.frame object for shiny.
# |----------------------------------------------------------------------------------|
# | 
require(plyr)
.saveMSEdataframe <- function()
{
	.MSELIB   <- "../DATA/"
	.rdaFILES <- paste0(.MSELIB,list.files(.MSELIB,pattern="^mse_",recursive=TRUE))
	.bn       <- basename(.rdaFILES)

	loadObj    <- function(fn){load(fn);return(sims)}
	S          <- lapply(.rdaFILES,loadObj)
	names(S)   <- substr(.bn,1,nchar(.bn)-6)
	lbl        <- strsplit(names(S),"_")
	quant      <- c(0.025,0.05,0.25,0.5,0.75,0.90,0.95,0.975)
	n          <- length(S)
	biomass.df <- NULL
	catch.df   <- NULL
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

		# Reference model spawning biomass depletion from Milka
		fn    <- function(X) return(X$m_sbt[1:nyr]/X$m_dBo)
		sdt   <- sapply(S[[i]],fn)
		m_sdt <- as.data.frame(t(apply(sdt,1,quantile,probs=quant)))
		colnames(m_sdt) <- paste0("t.Dt",quant)

		# DATA FRAME to concatenate
		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 Year      = Year,
		                 p_sbt, m_sbt, m_sdt)

		biomass.df <- rbind(biomass.df,df)

		# Catch by gear.
		fn    <- function(X)
		{
			df <- as.data.frame(X$m_dCatchData)
			colnames(df) <- c("Year","gear","area","group","sex","type","value")
			return(df)
		}
		ct    <- lapply(S[[i]],fn)
		ctdf  <- ldply(ct,data.frame)
		m_ct  <- ddply(ctdf,.(Year,gear,area,sex,group),summarize,
		               ct025=quantile(value,0.025),
		               ct05=quantile(value,0.05),
		               ct25=quantile(value,0.25),
		               ct50=quantile(value,0.50),
		               ct75=quantile(value,0.75),
		               ct95=quantile(value,0.95),
		               ct975=quantile(value,0.975))

		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 m_ct)

		catch.df <- rbind(catch.df,df)
	}
	mse.data <- list(biomass.df = biomass.df, catch.df=catch.df)
	save(mse.data,file=paste0(.MSELIB,"MSE.Rdata"))
}
.saveMSEdataframe()