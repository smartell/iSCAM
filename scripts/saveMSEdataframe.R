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

	#data frames for summarized performance metrics
	biomass.df <- NULL
	catch.df   <- NULL
	sublegal.df<- NULL
	AAV.df   <- NULL


	#calculated performance metrics for each MSE iteration
	rawbiomass.df <- NULL
	rawcatch.df   <- NULL
	rawsublegal.df<- NULL
	rawAAV.df   <- NULL


	for( i in 1:n)
	{
		#Scenario & Procedure from hdr Label
		cat("Scenario = ",lbl[[i]][2],"\t")
		cat("Procedure = ",lbl[[i]][3],"\n")

		Year  <- S[[i]][[1]]$yr
		nyr   <- length(Year)
		m_Year <- S[[i]][[1]]$m_dCatchData[,1]
		m_nyr <- length(m_Year)

		# Predicted Spawning biomass from iSCAM
		fn    <- function(X) return(X$sbt[1:nyr])
		sbt   <- sapply(S[[i]],fn)
		p_sbt <- t(apply(sbt,1,quantile,probs=quant))
		p_sbt <- as.data.frame(rbind(p_sbt, matrix(NA,ncol=ncol(p_sbt),nrow=m_nyr-nyr)))
		colnames(p_sbt) <- paste0("p.Bt",quant)

		p_sbt2 <- as.data.frame(rbind(sbt,matrix(NA,ncol=ncol(sbt),nrow=m_nyr-nyr)))
		colnames(p_sbt2) <- paste0("p.Bt",1:length(S[[i]]))
		
		

		# Reference Model spawning biomass from Milka
		fn    <- function(X) return(X$m_sbt[1:m_nyr])
		sbt   <- sapply(S[[i]],fn)
		m_sbt <- as.data.frame(t(apply(sbt,1,quantile,probs=quant)))
		colnames(m_sbt) <- paste0("t.Bt",quant)

		m_sbt2 <- as.data.frame(sbt)
		colnames(m_sbt2) <- paste0("m_sbt",1:length(S[[i]]))
		

		# Reference model spawning biomass depletion from Milka
		fn    <- function(X) return(X$m_sbt[1:m_nyr]/X$m_dBo)
		sdt   <- sapply(S[[i]],fn)
		m_sdt <- as.data.frame(t(apply(sdt,1,quantile,probs=quant)))
		colnames(m_sdt) <- paste0("t.Dt",quant)

		m_sdt2 <- as.data.frame(sdt)
		colnames(m_sdt2) <- paste0("t.Dt",1:length(S[[i]]))
		

		# P(Bt < 0.20)
		fn <- function(X)(return(X$m_sbt[1:m_nyr]/X$m_dBo <= 0.2))
		tmp <- sapply(S[[i]],fn)
		pb20 <- apply(tmp,1,FUN=function(tmp){sum(tmp==TRUE)})
		pb20 <- pb20/dim(tmp)[2]

		pb202 <- apply(tmp,2,FUN=function(tmp){as.numeric(tmp)})
		colnames(pb202) <- 1:length(S[[1]])
		
		# P(Bt < 0.30)
		fn <- function(X)(return(X$m_sbt[1:m_nyr]/X$m_dBo <= 0.3))
		tmp <- sapply(S[[i]],fn)
		pb30 <- apply(tmp,1,FUN=function(tmp){sum(tmp==TRUE)})
		pb30 <- pb30/dim(tmp)[2]

		pb302 <- apply(tmp,2,FUN=function(tmp){as.numeric(tmp)})

		# DATA FRAME to concatenate
		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 Year      = m_Year,
		                 p_sbt, m_sbt, m_sdt,
		                 "P(SSB<0.20)" = pb20,
		                 "P(SSB<0.30)" = pb30)

		biomass.df <- rbind(biomass.df,df)

		df2 <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 Year      = m_Year,
		                 p_sbt2, m_sbt2, m_sdt2,
		                 "P(SSB<0.20)" = pb202,
		                 "P(SSB<0.30)" = pb302)

		rawbiomass.df <- rbind(rawbiomass.df,df2)

		

		# Catch by gear.
		fn    <- function(X)
		{
			df <- as.data.frame(X$m_dCatchData)
			colnames(df) <- c("Year","gear","area","group","sex","type","value","sd")
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

		iteration <-rep(1:length(S[[i]]),each=nrow(ct[[1]]))
		m_ct2  <- cbind(ldply(ct,data.frame),iteration)
		


		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 m_ct)

		df2 <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 m_ct2)

		catch.df <- rbind(catch.df,df)

		rawcatch.df <- rbind(rawcatch.df,df)

		# Sublegal by gear
		fn  <- function(X)
		{
			df <- as.data.frame(X$m_dSubLegalData)
			colnames(df) <- c("Year","gear","area","group","sex","sublegal","waste","efficiency")
			return(df)
		}
		dt    <- lapply(S[[i]],fn)
		dtdf  <- ldply(dt,data.frame)
		m_dt  <- ddply(dtdf,.(Year,gear,area,group,sex),summarize,
		               dt025=quantile(sublegal,0.025),
		               dt05 =quantile(sublegal,0.05),
		               dt25 =quantile(sublegal,0.25),
		               dt50 =quantile(sublegal,0.50),
		               dt75 =quantile(sublegal,0.75),
		               dt95 =quantile(sublegal,0.95),
		               dt975=quantile(sublegal,0.975))

		iteration <-rep(1:length(S[[i]]),each=nrow(dt[[1]]))
		m_dt2  <- cbind(ldply(dt,data.frame)[,1:6],iteration)
		
		m_wt  <- ddply(dtdf,.(Year,gear,area,sex,group),summarize,
		               wt025=quantile(waste,0.025),
		               wt05 =quantile(waste,0.05),
		               wt25 =quantile(waste,0.25),
		               wt50 =quantile(waste,0.50),
		               wt75 =quantile(waste,0.75),
		               wt95 =quantile(waste,0.95),
		               wt975=quantile(waste,0.975))
		
		m_wt2 <- cbind(ldply(dt,data.frame)[,c(1:5,7)],iteration)


		m_ef  <- ddply(dtdf,.(Year,gear,area,sex,group),summarize,
		               ef025=quantile(efficiency,0.025),
		               ef05 =quantile(efficiency,0.05),
		               ef25 =quantile(efficiency,0.25),
		               ef50 =quantile(efficiency,0.50),
		               ef75 =quantile(efficiency,0.75),
		               ef95 =quantile(efficiency,0.95),
		               ef975=quantile(efficiency,0.975))

		m_ef2 <- cbind(ldply(dt,data.frame)[,c(1:5,8)],iteration)


		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 m_dt,m_wt,m_ef)

		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 m_dt2,m_wt2,m_ef2)

		sublegal.df <- rbind(sublegal.df,df)

		rawsublegal.df <- rbind(rawsublegal.df,df2)


		#AAV
		#calculate moving sum
		runsum<-function(x, int=5)
		{
			rsum=numeric(length=length(x))
			rsum[1:(int-1)]=0
	
			for(i in (int+1):length(x))
			{			
				rsum[i]<-sum(x[(i-int+1):i])
			} 
			return(rsum)
		}

		for(rr in 1:length(ct))
		{
			AAV=runsum(c(0,abs(diff(ct[[rr]]$value))))/runsum(ct[[rr]]$value)
			ct[[rr]]= cbind(ct[[rr]],AAV)
		}
		
		AAVdf  <- ldply(ct,data.frame)
		m_AAV  <- ddply(AAVdf,.(Year,gear,area,group),summarize,
		               AAV025=quantile(AAV,0.025,na.rm=T),
		               AAV05=quantile(AAV,0.05,na.rm=T),
		               AAV25=quantile(AAV,0.25,na.rm=T),
		               AAV50=quantile(AAV,0.50,na.rm=T),
		               AAV75=quantile(AAV,0.75,na.rm=T),
		               AAV95=quantile(AAV,0.95,na.rm=T),
		               AAV975=quantile(AAV,0.975,na.rm=T))

		m_AAV <- m_AAV[(nyr+1):m_nyr,]

		iteration <-rep(1:length(S[[i]]),each=nrow(ct[[1]]))
		m_AAVtmp  <- cbind(ldply(ct,data.frame),iteration)

		m_AAV2<-NULL
		for(a in 1:(length(S[[i]])))
		{
			tmp <- m_AAVtmp[m_AAVtmp$iteration==a,]
			tmp <- tmp[((nyr+1):m_nyr),]
			m_AAV2 <- rbind(m_AAV2,tmp)

		}

		
		df <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 m_AAV)
		df2 <- data.frame(Scenario  = lbl[[i]][2],
		                 Procedure = lbl[[i]][3],
		                 m_AAV2)                


		AAV.df <- rbind(AAV.df,df)

		rawAAV.df <- rbind(rawAAV.df,df2)


	}
	rawmse.data <- list(rawbiomass.df = rawbiomass.df, rawcatch.df=rawcatch.df, rawsublegal.df=rawsublegal.df,rawAAV.df=rawAAV.df)
	mse.data <- list(biomass.df = biomass.df, catch.df=catch.df, sublegal.df=sublegal.df,AAV.df=AAV.df)
	save(mse.data,file=paste0(.MSELIB,"MSE.Rdata"))
	save(rawmse.data,file=paste0(.MSELIB,"RAWMSE.Rdata"))
}
.saveMSEdataframe()