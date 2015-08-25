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
ro    <- 1.0            	# unfished equilibrium recruitment
h     <- 0.75				# steepness of the B-H SRR
kappa <- 4.0 * h / (1.0 - h)# recruitmetn compensation
m     <- c(0.15,0.16)		# sex-specific natural mortality rate
winf  <- c( 79.1688, 22.9836)
linf  <- c(150.8458,102.9795)
vbk	  <- c(  0.0795,  0.0975)
to	  <- c(  0.5970,  1.2430)
b	  <- 3.24
a	  <- 0.00000692
c     <- c(0,0)				# power parameter for age-dependent M

ahat  <- c(11.589,1000)
ghat  <- c(1.7732,0.01)

theta <- list(A=A,ro=ro,h=h,kappa=kappa,m=m,age=age,linf=linf,
              winf=winf,vbk=vbk,ahat=ahat,ghat=ghat)


# 
# SELECTIVITY PARAMETERS
# slx1 -> Length at 50% selectivity.
# slx2 -> STD in length-at-50% selectivity.
# slx3 -> Shape parameter for exponential logistic (0-1, where 0=asymptotic)
# slim -> minimum size limit for each gear.
# dmr  -> Discard mortality rates for each gear.
glbl <- c("IFQ","PSC","SPT","PER")
slx1 <- c(68.326,38.409,69.838,69.838)
slx2 <- c(3.338,4.345,5.133,5.133)
slx3 <- c(0.000,0.072,0.134,0.134)
slim <- c(82,00,82,82)
dmr  <- c(0.16,0.80,0.20,0.20)
slx  <- data.frame(sector=glbl,slx1=slx1,slx2=slx2,slx3=slx3)

sel1 <- slx[1,]
# aYPR -> Yield per recruit allocations.
aYPR <- c(0.50,0.40,0.05,0.05)
# aMPR -> Mortality per recruit allocations.
aMPR <- c(0.25,0.25,0.25,0.25)

# MANAGEMENT PROCEDURES
fstar <- 0.05
target_spr <- 0.40

MP0   <- list(fstar=fstar,
              slx=slx,
              pYPR=aYPR,
              pMPR=aMPR,
              slim=slim,
              dmr=dmr,
              type="YPR",
              target_spr = target_spr)

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
		getM <- function(age,vbk,m,c){
			t1 <- exp(vbk*(age+2))-1
			t2 <- exp(vbk*(age+1))-1
			sa <- (t1/t2)^(-m/vbk)
			mx=m*((log(t1)-log(t2))/vbk)^c
			return(mx)
		}
		mx <- sapply(age,getM,m=m,vbk=vbk,c=c)

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

		# List objects to return
		ageSc <- list(la=la,wa=wa,ma=ma,fa=fa,lx=lx,mx=mx,phi.E=phi.E)
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
			sc <- xplogis(la[h,],slx1[k],slx2[k],slx3[k])
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
		for(iter in 1:3)
		{
			# total mortality and survival rates
			fe <- fbar * lambda
			# browser()
			for(h in sex)
			{
				#print(fe)
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
				for(k in 1:ngear)
				{
					dlz[,j,k] <- sa[,j-1]*(dlz[,j-1,k]-lz[,j-1]*va[,j-1,k])
					if(j == A)
					{
						dlz[,j,k] <- dlz[,j,k]/oa[,j] - lz[,j-1]*sa[,j-1]*va[,j,k]*sa[,j]/oa[,j]^2
					}
				}
			}


			# Fmultipliers for fstar based on allocations
			qp    <- switch(type,YPR=qa,MPR=pa)
			ak    <- switch(type,YPR=pYPR,MPR=pQMPR)
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
		dphi.q <- matrix(0,ngear,ngear) 
		for(h in sex)
		{
			dphi.e <- dphi.e + as.vector(fa[h,] %*% dlz[h,,])
			phi.q <- phi.q + as.vector(lz[h,] %*% qa[h,,])
			phi.m <- phi.m + as.vector(lz[h,] %*% pa[h,,])
			# derivatives for yield equivalence
			for(k in 1:ngear)
			{
				for(kk in 1:ngear)
				{
					va2 <- va[h,,k] * va[h,,kk]
					dqa <- va2*wa[h,]*sa[h,]/za[h,] - va2*oa[h,]*wa[h,]/za[h,]^2
					dpq  <- as.double(qa[h,,k] %*% dlz[h,,k] + lz[h,] %*% dqa)
					dphi.q[k,kk] <- dphi.q[k,kk] + dpq
				}
			}
		}
		spr   <- (phi.e/phi.E)
		ispr  <- (phi.E/phi.e)

		# equilibrium recruitment & sensitivity
		re    <- max(0,ro*(kappa-ispr)/(kappa-1))
		dre   <- ro * phi.E * dphi.e / (phi.e^2 * (kappa-1))
		
		# yield per recuit and yield
		ypr   <- fe * phi.q
		ye    <- re * ypr

		# mortality per recruit and mortality
		mpr   <- fe * phi.m
		me    <- re * mpr
		
		# Jacobian for yield
		dye <- matrix(0,nrow=ngear,ncol=ngear)
		for(k in 1:ngear)
		{
			for(kk in 1:ngear)
			{
				dye[k,kk] = fe[k]*re*dphi.q[k,kk] + fe[k]*phi.q[k]*dre[kk]
				if(k == kk)
				{
					dye[k,kk] = dye[k,kk] + re*phi.q[k]
				}
			}
		}
		# dye   <- re * as.vector(phi.q * diag(1,ngear)) + fe * phi.q * dre + fe * re * dphi.q
		# dye   <- re * (phi.q) + fe * phi.q * dre + fe * re * dphi.q
		
		# Yield equivalence
		v  <- sqrt(diag(dye))
		M  <- dye / (v %o% v)
		# print(v %o% v)
		# print(t(matrix(v,4,4)))

		# print(M)
		# cat("ye\n",ye,"\n")


		# Equilibrium Model output
		out <- list(
		            "fe"  = fe,
		            "ye"  = ye,
		            "me"  = me,
		            "re"  = re,
		            "spr" = spr,
		            "ypr" = ypr,
		            "mpr" = mpr,
		            "dre" = dre,
		            "dye" = as.vector(diag(dye)),
		            "fstar" = fstar,
		            "gear" = slx$sector,
		            "dlz"  = dlz,
		            "lz"   = lz,
		            "ak"   = ak

		            )

		return(out)
	})
}


runModel <- function(MP)
{
	# Initialize model
	theta <- .getAgeSchedules(theta)
	selex <- .getSelectivities(MP,theta$la)
	selex$pYPR <- MP$pYPR
	selex$pMPR <- MP$pMPR
	
	
	
	fn    <- function(log.fbar)
	{
		MP$fstar = exp(log.fbar)
		theta <- .getAgeSchedules(theta)
		selex <- .getSelectivities(MP,theta$la)
		selex$pYPR <- MP$pYPR
		selex$pMPR <- MP$pMPR
		em <- eqModel(theta,selex,type=MP$type)
		print(em$spr)
		return((em$spr - MP$target_spr)^2)
	}
	# Perform nonlinear search to find Fspr = spr_target
	
	fit <- optim(log(MP$fstar),fn,method="BFGS")
	MP$fstar <- exp(fit$par)
	theta <- .getAgeSchedules(theta)
	selex <- .getSelectivities(MP,theta$la)
	selex$pYPR <- MP$pYPR
	selex$pMPR <- MP$pMPR
	
	EM    <- eqModel(theta,selex,type=MP$type)
	return(EM)
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

runProfile <- function(MP)
{
	fbar <- seq(0,0.32,length=100)
	fn   <- function(log.fbar)
	{
		MP$fstar <- exp(log.fbar)
		em       <- run(MP)
		return(em)
	}
	runs <- lapply(log(fbar),fn)
	df   <- ldply(runs,data.frame)
	p  <- ggplot(df,aes(fe,ye))
	p  <- p + geom_line(aes(col=gear),size=1.3)
	p  <- p + geom_vline(aes(xintercept=fe[which.max(ye)],col=gear))
	p  <- p + geom_line(aes(fe,dye/30,col=gear),data=df)
	print(p+facet_wrap(~gear,scales="free_x"))
	return(df)
}

runProfile2 <- function(MP)
{
	PSCaks<-seq(0.4,0.01,by=-0.01)
	

  	if(MP$type=="YPR"){
  		fn   <- function(PSCak)
		{
  			othak<- MP$pYPR[-2]/sum(MP$pYPR[-2])
			MP$pYPR[2]<- PSCak
			MP$pYPR[-2]<- (1-PSCak)*othak
			
			
			
			em       <- runModel(MP)
			return(em)
		}
	}else{
		fn   <- function(PSCak)
		{
			othak<- MP$pMPR[-2]/sum(MP$pMPR[-2])
			MP$pMPR[2]<- PSCak
			MP$pMPR[-2]<- (1-PSCak)*othak
			em       <- runModel(MP)
			return(em)
		}
	}	

  
	

	
	
	runs <- lapply(PSCaks,fn)
	df   <- ldply(runs,data.frame)

	AK_PSC <- rep(PSCaks, each=length(MP$pYPR))
	ye_PSC <- rep(df$ye[df$gear=="PSC"], each=length(MP$pYPR))
	
	gr<-unique(df$gear)
	tt_ye<-df$ye[df$gear==gr[1]]+df$ye[df$gear==gr[2]]+df$ye[df$gear==gr[3]]+df$ye[df$gear==gr[4]]
	to_ye<-df$ye[df$gear==gr[1]]+df$ye[df$gear==gr[3]]+df$ye[df$gear==gr[4]]
	ye_tot <- rep(tt_ye, each=length(MP$pYPR))
	ye_oth <- rep(to_ye, each=length(MP$pYPR))

	#mylm<-lm(df$ye[df$gear=="IFQ"]~df$ak[df$gear=="IFQ"])
	#plot(mylm)

	df1 <- cbind(df,AK_PSC,ye_PSC,ye_tot,ye_oth)
	df2 <- cbind(AK_PSC,ye_PSC,ye_tot,ye_oth)
	mdf<-melt(df1)

	p  <- ggplot(df1,aes(ye_PSC,ye_tot))
	p  <- p + geom_line(aes(col=gear),size=1.3)
	print(p)
	
	p  <- ggplot(mdf,aes(variable,value))
	p  <- p + geom_line(aes(col=gear),size=1.3)
	print(p)

	return(df1)
}

df1<-runProfile2(MP0)

p  <- ggplot(df1,aes(AK_PSC,round(spr,6)))
p  <- p + geom_line(aes(col=gear),size=1.3)
print(p)

p  <- ggplot(df1,aes(ye_PSC,ye/ye_PSC))
p  <- p + geom_line(aes(col=gear),size=1.3)
print(p)

p  <- ggplot(df1,aes(ye_PSC,ye_oth))
p  <- p + geom_line(aes(col=gear),size=1.3)
print(p)

p  <- ggplot(df1,aes(AK_PSC,ypr))
p  <- p + geom_line(aes(col=gear),size=1.3)

print(p+facet_wrap(~gear,scales="free_x"))
	
p  <- ggplot(df1,aes(ye_PSC,ye))
p  <- p + geom_line(aes(col=gear),size=1.3) 
print(p)
print(p+facet_wrap(~gear,scales="fixed"))
		



names(df1)
checkDerivatives <- function(MP)
{
	hh <- 1e-5
	mph <- mmh <- MP
	mph$fstar <- MP$fstar + hh 
	mmh$fstar <- MP$fstar - hh 
	ph <- run(mph)
	mh <- run(mmh)
	mp <- run(MP)

	ndlz <- (ph$lz-mh$lz)/(2*hh)
	matplot(t(ndlz))
	matlines(t(mp$dlz[,,1]))
}

main <- {
	df <- runProfile(MP0)
	# checkDerivatives(MP0)
}



