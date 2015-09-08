# --------------------------------------------------------------------------- #
# TOTAL MORTALITY ALLOCATION
# AUTHOR: Steven Martell
# This is an equilibrium model for allocating the total mortality rate among 
# one or more fisheries.

# TODO
# -[X] Change algorithm to estimate fbar only.
# --------------------------------------------------------------------------- #
library(plyr)
library(ggplot2)

# 
# STOCK PARAMETERS
# 
A     <- 30					# Plus group Age
age   <- 1:A 				# vector of ages
H     <- 2                  # number of sexes
sex   <- 1:H 				# vector of sex indexes (F=1,M=2)
bo    <- 520                # unfished female spawning stock biomass (520 from 2015 assessment)
# ro    <- 1.0            	# unfished equilibrium recruitment
h     <- 0.95				# steepness of the B-H SRR
kappa <- 4.0 * h / (1.0 - h)# recruitmetn compensation
m     <- c(0.15,0.16)		# sex-specific natural mortality rate
winf  <- c( 79.1688, 22.9836)
linf  <- c(150.8458,102.9795)
vbk	  <- c(  0.0795,  0.0975)
to	  <- c(  0.5970,  1.2430)
b	  <- 3.24
a	  <- 0.00000692
cm     <- c(0,0)				# power parameter for age-dependent M

ahat  <- c(11.589,1000)
ghat  <- c(1.7732,0.01)

theta <- list(A=A,bo=bo,h=h,kappa=kappa,m=m,age=age,linf=linf,
              winf=winf,vbk=vbk,ahat=ahat,ghat=ghat)


# 
# SELECTIVITY PARAMETERS - need to change these to the data given by ian
# slx1 -> Length at 50% selectivity.
# slx2 -> STD in length-at-50% selectivity.
# slx3 -> Shape parameter for exponential logistic (0-1, where 0=asymptotic)
# slx4 -> asymptote age.
# slim -> minimum size limit for each gear.
# dmr  -> Discard mortality rates for each gear.
glbl <- c("IFQ","PSC","SPT","PER")

#sel from Ian's assessment
slx1 <- c(68.326,38.409,69.838,69.838)
slx2 <- c(3.338,4.345,5.133,5.133)
slx3 <- c(0.000,0.072,0.134,0.134)
slx4 <- c(30,16,16,16)
slim <- c(82,00,82,00)
dmr  <- c(0.16,0.80,0.20,0.00)
slx  <- data.frame(sector=glbl,slx1=slx1,slx2=slx2,slx3=slx3,slx4=slx4)

sel1 <- slx[1,]
# aYPR -> Yield per recruit allocations.
aYPR <- c(0.66,0.18,0.14,0.02)
# aYPR <- c(0.666541,0.177337,0.134677,0.021445)
# aMPR -> Mortality per recruit allocations.
# aMPR <- c(0.43443443,0.43043043,0.12012012,0.01501502)
aMPR <- c(0.43856475,0.42604559,0.11420931,0.02118035)

# fixed PSC limit for status quo
pscLimit  = c(NA,7.75,NA,NA)

# MANAGEMENT PROCEDURES
fstar <- 0.107413
sprTarget <- 0.45
MP0   <- list(	fstar     = fstar,
				slx       = slx,
				pYPR      = aYPR/sum(aYPR),
				pMPR      = aMPR/sum(aMPR),
				pscLimit  = pscLimit,
				slim      = slim,
				dmr       = dmr,
				sprTarget = sprTarget,
				type      = "YPR")

MP1 <- MP2 <- AB0 <- JRG0 <- JRG1 <- JRG2 <- MP0
MP1$slx$slx3[2] = 0.1
MP2$slx$slx3[2] = 0.2
JRG1$slx$slx3[2] = 0.1
JRG2$slx$slx3[2] = 0.2

# Abundance based PSC
AB0$pscLimit = c(NA,NA,NA,NA)
AB1 <- AB2 <- BB0 <- AB0
AB1$slx$slx3[2] = 0.1
AB2$slx$slx3[2] = 0.2

# Mortality based footprint
BB0$type="MPR"
BB1 <- BB2 <- BB4 <- BB0
BB1$slx$slx3[2] = 0.1
BB2$slx$slx3[2] = 0.2
BB4$type="YPR"



# 
# AGE SCHEDULE INFORMATION
# 
.getAgeSchedules <- function(theta)
{
	with(theta,{
		# Length-at-age, weight-at-age, mature weight-at-age
		vonb <- function(age,linf,k,to) {return( linf*(1-exp(-k*(age+to))) )}
		la <- sapply(age,vonb,linf=linf,k=vbk,to=to)
		wa <- a*la^b
		ma <- sapply(age,plogis,location=ahat,scale=ghat);  ma[,1:7] <- 0
		fa <- ma * wa

		# Age-dependent natural mortality (Lorenzen)
		getM <- function(age,vbk,m,cm){
			t1 <- exp(vbk*(age+2))-1
			t2 <- exp(vbk*(age+1))-1
			sa <- (t1/t2)^(-m/vbk)
			mx=m*((log(t1)-log(t2))/vbk)^cm
			return(mx)
		}
		mx <- sapply(age,getM,m=m,vbk=vbk,cm=cm)

		# Survivorship at unfished conditions
		lx <- matrix(0,H,A)
		for(i in age)
		{
			if(i == min(age))
			{
				lx[,i] <- 1.0/H
			}
			else
			{
				lx[,i] <- lx[,i-1] * exp(-mx[,i-1])
			}
			if(i == A)
			{
				lx[,i] <- lx[,i] / (1-exp(-mx[,i]))
			}
		}
		phi.E <- sum(lx*fa)
		ro    <- bo/phi.E

		# List objects to return
		ageSc <- list(ro=ro,la=la,wa=wa,ma=ma,fa=fa,lx=lx,mx=mx,phi.E=phi.E)
		theta <- c(theta,ageSc)
		return(theta)
	})
}

# 
# GET SELECTIVITIES 
# mp - a list object representing the procedure
.getSelectivities <- function(mp,la)
{
	
	with (mp,{
		ngear <- dim(slx)[1]
		va <- array(0,dim=c(H,A,ngear))
		for(k in 1:ngear)
		for(h in sex)
		{

			sc <- xplogis(la[h,],slx$slx1[k],slx$slx2[k],slx$slx3[k])
			sc[slx$slx4[k]:A] <- sc[slx$slx4[k]-1]

			ra <- plogis(la[h,],slim[k],0.1*la[h,])
			da <- (1-ra)*dmr[k]
			va[h,,k] <- sc*(ra+da)
		}

		return(list(ngear=ngear,fstar=fstar,va=va))
	})
}


# 
# EXPONENTIAL LOGISTIC FUNCTION FOR SELECTIVITY
# 
xplogis <- function(x,mu,sd,g)
{
	s <- (1/(1-g))*((1-g)/g)^g
	p <- plogis(x,mu,sd)*exp(g*(mu-x)/sd)
	return(s*p)
}

# 
# EQUILIBRIUM MODEL
# 
eqModel <- function(theta,selex,type="YPR")
{

	with(c(theta,selex),{
		lz  <- matrix(1/H,nrow=H,ncol=A)
		za  <- matrix(0,nrow=H,ncol=A)
		qa  <- array(0,dim=c(H,A,ngear))
		pa  <- array(0,dim=c(H,A,ngear))
		dlz <- array(0,dim=c(H,A,ngear))

		# Survivorship under fished conditions at fstar
		fbar <- fstar
		lambda <- rep(1.0,length=ngear)
		for(iter in 1:4)
		{
			# total mortality and survival rates
			fe <- fbar * lambda
			# browser()
			for(h in sex)
			{
				# print(fe)
				if(dim(va)[3] > 1){
					fage   <- rowSums(fe*va[h,,])
				}
				else if(dim(va)[3] == 1){
					fage   <- fe * va[h,,]
				}
				za[h,] <- mx[h,] + fage
			}
			sa <- exp(-za)
			oa <- 1.0 - sa

			# per recruit yield & numbers for gear k
			for(k in 1:ngear)
			{
				qa[,,k] <- va[,,k] * wa * oa / za
				pa[,,k] <- va[,,k] * oa / za
			}

			#  survivorship
			for(j in 2:A)
			{
				lz[,j] <- lz[,j-1] * sa[,j-1]
				if(j == A)
				{
					lz[,j] <- lz[,j] / oa[,j]
				}

				# partial derivatives
				# for(k in 1:ngear)
				# {
				# 	dlz[,j,k] <- sa[,j-1]*(dlz[,j-1,k]-lz[,j-1]*va[,j-1,k])
				# 	if(j == A)
				# 	{
				# 		dlz[,j,k] <- dlz[,j,k]/oa[,j] - lz[,j-1]*sa[,j-1]*va[,j,k]*sa[,j]/oa[,j]^2
				# 	}
				# }
			}


			# Fmultipliers for fstar based on allocations
			qp    <- switch(type,YPR=qa,MPR=pa)
			ak    <- switch(type,YPR=pYPR,MPR=pMPR)
			phi.t <- 0
			for(h in sex)
			{
				phi.t <- phi.t + as.vector(lz[h,] %*% qp[h,,])
			}
			
			
			lam.t  <- ak / (phi.t/sum(phi.t))
			lambda <- lam.t / sum(lam.t)
			# cat(iter," lambda = ",lambda,"\n")
		}

		# incidence functions
		phi.e  <- sum(lz*fa)
		phi.q  <- phi.m <- dphi.e <- dre <- 0
		# dphi.q <- matrix(0,ngear,ngear) 
		for(h in sex)
		{
			# dphi.e <- dphi.e + as.vector(fa[h,] %*% dlz[h,,])
			phi.q  <- phi.q + as.vector(lz[h,] %*% qa[h,,])
			phi.m  <- phi.m + as.vector(lz[h,] %*% pa[h,,])
			# derivatives for yield equivalence
			# for(k in 1:ngear)
			# {
			# 	for(kk in 1:ngear)
			# 	{
			# 		va2 <- va[h,,k] * va[h,,kk]
			# 		dqa <- va2*wa[h,]*sa[h,]/za[h,] - va2*oa[h,]*wa[h,]/za[h,]^2
			# 		dpq  <- as.double(qa[h,,k] %*% dlz[h,,k] + lz[h,] %*% dqa)
			# 		dphi.q[k,kk] <- dphi.q[k,kk] + dpq
			# 	}
			# }
		}
		spr   <- (phi.e/phi.E)
		ispr  <- (phi.E/phi.e)

		# equilibrium recruitment & sensitivity
		re    <- max(0,ro*(kappa-ispr)/(kappa-1))
		be    <- re * phi.e
		# dre   <- ro * phi.E * dphi.e / (phi.e^2 * (kappa-1))
		
		# yield per recuit and yield
		ypr   <- fe * phi.q
		ye    <- re * ypr

		# mortality per recruit and mortality
		mpr   <- fe * phi.m
		me    <- re * mpr
		
		# Jacobian for yield
		# dye <- matrix(0,nrow=ngear,ncol=ngear)
		# for(k in 1:ngear)
		# {
		# 	for(kk in 1:ngear)
		# 	{
		# 		dye[k,kk] = fe[k]*re*dphi.q[k,kk] + fe[k]*phi.q[k]*dre[kk]
		# 		if(k == kk)
		# 		{
		# 			dye[k,kk] = dye[k,kk] + re*phi.q[k]
		# 		}
		# 	}
		# }
		# dye   <- re * as.vector(phi.q * diag(1,ngear)) + fe * phi.q * dre + fe * re * dphi.q
		# dye   <- re * (phi.q) + fe * phi.q * dre + fe * re * dphi.q
		
		# Yield equivalence
		# v  <- sqrt(diag(dye))
		# M  <- dye / (v %o% v)
		# print(v %o% v)
		# print(t(matrix(v,4,4)))

		# print(M)
		# cat("ye\n",ye,"\n")


		# Equilibrium Model output
		out <- list(
		            "fe"  = fe,
		            "ye"  = ye,
		            "me"  = me,
		            "be"  = be,
		            "re"  = re,
		            "spr" = spr,
		            "ypr" = ypr,
		            "mpr" = mpr,
		            # "dre" = dre,
		            # "dye" = as.vector(diag(dye)),
		            "fstar" = fstar,
		            "gear" = slx$sector
		            )

		return(out)
	})
}

run <- function(MP)
{
	theta <- .getAgeSchedules(theta)
	selex <- .getSelectivities(MP,theta$la)
	selex$pYPR <- MP$pYPR
	selex$pMPR <- MP$pMPR
	
	EM    <- eqModel(theta,selex,type=MP$type)
	return(EM)
}

plotSelex <- function()
{
	theta <- .getAgeSchedules(theta)
	S0    <- .getSelectivities(MP0,theta$la)$va[1,,2]
	S1    <- .getSelectivities(MP1,theta$la)$va[1,,2]
	S2    <- .getSelectivities(MP2,theta$la)$va[1,,2]
	df <- cbind(1:A,S0,S1,S2)
}


getFspr <- function(MP)
{
	fn <- function(log.fspr)
	{
		MP$fstar <- exp(log.fspr)
		spr  	 <- run(MP)$spr
		ofn   	 <- (spr-MP$sprTarget)^2
		return(ofn)
	}	
	fit <- optim(log(MP$fstar),fn,method="BFGS")

	MP$fstar = exp(fit$par)
	# print(fit)
	return(MP)
}

getFsprPSC <- function(MP)
{
	# add swtich statement here for type of allocation
	ak    <- switch(MP$type,YPR=MP$pYPR,MPR=MP$pMPR)
	bGear <- !is.na(MP$pscLimit)
	iGear <- which(!is.na(MP$pscLimit))
	pk    <- ak[!bGear]/sum(ak[!bGear])
	

	getFootPrint <- function(phi)
	{
		tmp        <- ak
		tmp[bGear] <- phi[-1]
		tmp[!bGear]<- (1-sum(tmp[bGear]))*pk
		return(tmp)
	}

	fn <- function(phi)
	{
		MP$fstar <- exp(phi[1])
		
		# MP$pYPR  <- tmp
		MP$pYPR  <- getFootPrint(phi)

		EM       <- run(MP)
		spr  	 <- EM$spr
		psc      <- EM$ye
		
		ofn   	 <- (spr-MP$sprTarget)^2 + 10*sum(na.omit((psc-MP$pscLimit)^2))
		
		return(ofn)
	}
	parms <- c(log(MP$fstar),ak[iGear])
	fit   <- optim(parms,fn,method="BFGS",control=list(maxit=1000))
	
	MP$fstar <- exp(fit$par[1])
	MP$pYPR  <- getFootPrint(fit$par)
	return(MP)
}

getFstar <- function(MP)
{
	# check if there is a fixed PSC limit
	if ( any(!is.na(MP$pscLimit)) ){
		print("PSC LIMIT")
		MP <- getFsprPSC(MP)
	}	
	else{
		print("NO PSC LIMIT")
		MP <- getFspr(MP)
	}

	return(MP)

}

runProfile <- function(MP)
{
	fbar <- seq(0,0.45,length=100)
	fn   <- function(log.fbar)
	{
		MP$fstar <- exp(log.fbar)
		em       <- run(MP)
		return(em)
	}
	runs <- lapply(log(fbar),fn)
	df   <- ldply(runs,data.frame)
	p <- ggplot(df,aes(fstar,ye))
	p <- p + geom_line(aes(col=gear),size=1.3)
	p <- p + labs(x="Fishing Rate (F*)",y="Yield",col="Sector")
	# p  <- p + geom_vline(aes(xintercept=fe[which.max(ye)],col=gear))
	# p  <- p + geom_line(aes(fe,dye/30,col=gear),data=df)
	# print(p+facet_wrap(~gear,scales="free_x"))
	print(p + theme_minimal(22) + ylim(c(0,1)))

	q <- ggplot(df,aes(fstar,spr)) 
	q <- q + geom_line(aes(col=gear),size=1.3)
	q <- q + labs(x="Fishing Rate (F*)",y="Spawning Potential Ratio (SPR)",col="Sector") 
	print(q + theme_minimal(22) + ylim(c(0,1)))

	r <- ggplot(df,aes(spr,1-exp(-fe))) 
	r <- r + geom_line(aes(col=gear),size=1.3)
	r <- r + labs(y="Harvest Rate (f)",x="Spawning Potential Ratio (SPR)",col="Sector") 
	print(r + theme_minimal(22) + ylim(c(0,0.3)))

	return(df)
}

yieldEquivalence <- function(MP)
{
	# base   <- run(MP)
	G    <- length(MP$pYPR)
	x    <- rep(1,length=G)
	D    <- rbind(rep(1,length=G),1 - diag(1,G))

	fn   <- function(x){
		# print(x)
		if(MP$type=="YPR")
		{
			ak      <- x * MP$pYPR
			MP$pYPR = ak / sum(ak)
		}
		if(MP$type=="MPR")
		{
			ak      <- x * MP$pMPR
			MP$pMPR = ak / sum(ak)	
		}
		ftmp     <- exp(getFspr(MP)$par)
		MP$fstar <- ftmp
		rtmp     <- run(MP)
		rtmp$x   <- x
		return(rtmp)
	}
	XX <- apply(D,1,fn)
	df <- ldply(XX,data.frame)
	Y  <- matrix(df$ye,ncol=G,byrow=TRUE)
	y  <- Y[1,]
	M  <- Y[-1,]
	E  <- t(t(M)/y)

}



main <- {
	# fspr <- exp(getFspr(MP0)$par)
	# MP0$fstar = fspr
	# M0 <- run(MP0)

	# status quo scenario (Fixed PSC limit)
	MP0 <- getFstar(MP0)
	M0  <- run(MP0)

	# Fixed PSC limit with excluder 1
	MP1 <- getFstar(MP1)
	M1  <- run(MP1)

	# Fixed PSC limit with excluder 2
	MP2 <- getFstar(MP2)
	M2  <- run(MP2)

	# Index-based PSC limit STQ
	AB0 <- getFstar(AB0)
	A0  <- run(AB0)
	AB1 <- getFstar(AB1)
	A1  <- run(AB1)
	AB2 <- getFstar(AB2)
	A2  <- run(AB2)

	# Index-based PSC limit MPR
	BB0 <- getFstar(BB0)
	B0  <- run(BB0)
	BB1 <- getFstar(BB1)
	B1  <- run(BB1)
	BB2 <- getFstar(BB2)
	B2  <- run(BB2)	
	BB4 <- getFstar(BB4)
	B4  <- run(BB4)	

	theta$c = c(0.8,0.8)
	JRG0 <- getFstar(JRG0)
	JG0  <- run(JRG0)

	JRG1 <- getFstar(JRG1)
	JG1  <- run(JRG1)

	JRG2 <- getFstar(JRG2)
	JG2  <- run(JRG2)


	# PROFILE OVER FSTAR
	# df <- runProfile(MP0)

	# COMPUTE YIELD EQUIVALENCE
	# E  <- yieldEquivalence(MP0)
}

# FIXED PSC LIMITS
cat("\nEquilibrium yield\n")
Ye <- print(data.frame(round(cbind(STQ=M0$ye,EX1=M1$ye,EX2=M2$ye)/2.2046,3)))

cat("\nEquilibrium yield Gauvin\n")
Yeg <- print(data.frame(round(cbind(STQ=JG0$ye,EX1=JG1$ye,EX2=JG2$ye)/2.2046,3)))



cat("\nYield Per Recruit \n")
YPR <- print(data.frame(round(cbind(STQ=M0$ypr,EX1=M1$ypr,EX2=M2$ypr),3)))

cat("\nYield Per Recruit Gauvin\n")
YPRg <- print(data.frame(round(cbind(STQ=JG0$ypr,EX1=JG1$ypr,EX2=JG2$ypr),3)))


cat("\nMortality Per Recruit\n")
MPR <- print(data.frame(round(cbind(STQ=M0$mpr,EX1=M1$mpr,EX2=M2$mpr),3)))

cat("\nMortality Per Recruit Gauvin\n")
MPRg <- print(data.frame(round(cbind(STQ=JG0$mpr,EX1=JG1$mpr,EX2=JG2$mpr),5)))


cat("\nYield Per Recruit Footprint\n")
YPRfp <- print(data.frame(round(data.frame(t(t(YPR)/colSums(YPR))),3))*100)

cat("\nMortality Per Recruit Footprint\n")
MPRfp <- print(data.frame(round(data.frame(t(t(MPR)/colSums(MPR))),3))*100)

##
##cat("\nMortality Per Recruit Footprint Gauvin\n")
##MPRfp <- print(data.frame(round(data.frame(t(t(MPRg)/colSums(MPRg))),3))*100)
##

# INDEX-BASED PSC LIMITS WITH ALLOCATION BASED ON YIELD
# cat("\nEquilibrium yield\n")
# Ye <- print(data.frame(round(cbind(STQ=A0$ye,EX1=A1$ye,EX2=A2$ye)/2.2046,3)))

# cat("\nYield Per Recruit\n")
# YPR <- print(data.frame(round(cbind(STQ=A0$ypr,EX1=A1$ypr,EX2=A2$ypr),3)))

# cat("\nMortality Per Recruit\n")
# MPR <- print(data.frame(round(cbind(STQ=A0$mpr,EX1=A1$mpr,EX2=A2$mpr),3)))

# cat("\nYield Per Recruit Footprint\n")
# YPRfp <- print(data.frame(round(data.frame(t(t(YPR)/colSums(YPR))),3))*100)

# cat("\nMortality Per Recruit Footprint\n")
# MPRfp <- print(data.frame(round(data.frame(t(t(MPR)/colSums(MPR))),3))*100)



# INDEX-BASED PSC LIMITS WITH ALLOCATION BASED ON MPR
# cat("\nEquilibrium yield\n")
# Ye <- print(data.frame(round(cbind(YPR=A0$ye,STQ=B0$ye,EX1=B1$ye,EX2=B2$ye)/2.2046,3)))

# cat("\nYield Per Recruit\n")
# YPR <- print(data.frame(round(cbind(STQ=B0$ypr,EX1=B1$ypr,EX2=B2$ypr),3)))

# cat("\nMortality Per Recruit\n")
# MPR <- print(data.frame(round(cbind(STQ=B0$mpr,EX1=B1$mpr,EX2=B2$mpr),3)))

# cat("\nYield Per Recruit Footprint\n")
# YPRfp <- print(data.frame(round(data.frame(t(t(YPR)/colSums(YPR))),3))*100)

# cat("\nMortality Per Recruit Footprint\n")
# MPRfp <- print(data.frame(round(data.frame(t(t(MPR)/colSums(MPR))),3))*100)

# CPUE <-print(round(Ye[,-1]/cbind(B0$fe,B1$fe,B2$fe),3))












