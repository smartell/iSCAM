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
	sublegal.df<- NULL
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

		# P(Bt < 0.20)
		fn <- function(X)(return(X$m_sbt[1:nyr]/X$m_dBo <= 0.2))
		tmp <- sapply(S[[i]],fn)
		pb20 <- apply(tmp,1,FUN=function(tmp){sum(tmp==TRUE)})
		pb20 <- pb20/dim(tmp)[2]

		# P(Bt < 0.30)
		fn <- function(X)(return(X$m_sbt[1:nyr]/X$m_dBo <= 0.3))
		tmp <- sapply(S[[i]],fn)
		pb30 <- apply(tmp,1,FUN=function(tmp){sum(tmp==TRUE)})
		pb30 <- pb30/dim(tmp)[2]

		# DATA FRAME to concatenate
		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 Year      = Year,
		                 p_sbt, m_sbt, m_sdt,
		                 "P(SSB<0.20)" = pb20,
		                 "P(SSB<0.30)" = pb30)

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
		               ct05 =quantile(value,0.05),
		               ct25 =quantile(value,0.25),
		               ct50 =quantile(value,0.50),
		               ct75 =quantile(value,0.75),
		               ct95 =quantile(value,0.95),
		               ct975=quantile(value,0.975))


		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 m_ct)

		catch.df <- rbind(catch.df,df)

		# Sublegal by gear
		fn  <- function(X)
		{
			df <- as.data.frame(X$m_dSubLegalData)
			colnames(df) <- c("Year","gear","area","group","sex","sublegal","waste")
			return(df)
		}
		dt    <- lapply(S[[i]],fn)
		dtdf  <- ldply(dt,data.frame)
		m_dt  <- ddply(dtdf,.(Year,gear,area,sex,group),summarize,
		               dt025=quantile(sublegal,0.025),
		               dt05 =quantile(sublegal,0.05),
		               dt25 =quantile(sublegal,0.25),
		               dt50 =quantile(sublegal,0.50),
		               dt75 =quantile(sublegal,0.75),
		               dt95 =quantile(sublegal,0.95),
		               dt975=quantile(sublegal,0.975))

		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 m_dt)

		sublegal.df <- rbind(sublegal.df,df)


	}
	mse.data <- list(biomass.df = biomass.df, catch.df=catch.df, sublegal.df=sublegal.df)
	save(mse.data,file=paste0(.MSELIB,"MSE.Rdata"))
}
.saveMSEdataframe()