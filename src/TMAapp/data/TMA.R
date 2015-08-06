# --------------------------------------------------------------------------- #
# TOTAL MORTALITY ALLOCATION
# AUTHOR: Steven Martell
# This is an equilibrium model for allocating the total mortality rate among 
# one or more fisheries.

# TODO
# -[ ] Change algorithm to estimate fbar only.
# --------------------------------------------------------------------------- #
print("loading TMA.r")

library(ggplot2)

A     <- 25				# Plus group Age
ro    <- 1.0             # unfished equilibrium recruitment
h     <- 0.75			# steepness of the B-H SRR
kappa <- 4.0 * h / (1.0 - h)
m     <- 0.15			# instantaneous natural mortality rate
age   <- 1:A
mx    <- rep(m,length=A)
linf  <- 120
winf  <- 50
vbk   <- 1.5 * m
ahat  <- 11.5
ghat  <- 1.5

# 
# Selectivity parameters for each gear.
# 
s1    <- c( 85.0, 25.0, 45.0)
s2    <- c( 12.0, 5.50, 11.2)
s3    <- c( 00.0, 0.10, 0.05)

ng    <- length(s1) 		# number of fleets
gear  <- 1:ng


# 
# Age-dependent natural mortality (Lorenzen)
# 
t1 <- exp(vbk*(age+2))-1
t2 <- exp(vbk*(age+1))-1
sa <- (t1/t2)^(-m/vbk)
mx=m/vbk*(log(t1)-log(t2))

# 
# MANAGEMENT CONTROLS
# 
target_spr     <- 0.40
fspr           <- 0.10
# fe             <- rep(0.01,length=ng)
tma_allocation <- rep(1/ng,length=ng)
slim           <- c(NA,NA,NA)
dmr   		   <- c(0.16,0.90,0.2)
ak             <- tma_allocation


# 
# AGE-SCHEDULE INFORMATION
# 
la    <- linf*(1-exp(-vbk*age))
wa    <- winf / (linf^3)*la^3
ma    <- plogis(age,ahat,ghat)
fa    <- wa * ma
va    <- matrix(nrow=ng,ncol=A)
ra    <- matrix(nrow=ng,ncol=A)
qa    <- matrix(nrow=ng,ncol=A)
pa    <- matrix(nrow=ng,ncol=A)
dlz   <- matrix(nrow=ng,ncol=A)

lx    <- rep(1,length=A)
for(i in 2:A) lx[i] <- lx[i-1] * exp(-mx[i-1])
lx[A] <- lx[A]/(1-exp(-mx[A]))
phi.E <- as.double(lx %*% fa)




# 
# FISHERIES SELECTIVITY
#  
# exponential logistic curve from Grant Thompson.
# 
gplogis <- function(x,a,b,g)
{
	s <- (1/(1-g))*((1-g)/g)^g
	s <- s * exp(1/a*g*(b-x)) / (1+exp(1/a*(b-x)))
	return(s)
}

getSlx <- function(s1,s2,s3)
{
	for(k in gear)
	{
		# size-based selectivity
		sc 	   <- gplogis(la,a=s2[k],b=s1[k],g=s3[k])
		if(!is.na(slim[k]))
			ra[k,] <- plogis(la,slim[k],1.0)
		if(is.na(slim[k]))
			ra[k,] <- 0
		va[k,] <- sc*(ra[k,]+(1.0-ra[k,])*dmr[k])
	}
	return(va)
}

# 
# Survivorship under fished conditions.
# 
# va <- getSlx(s1,s2,s3)
# za <- mx + colSums(fe*va)
# sa <- exp(-za)
# oa <- 1.0 - sa
# lz <- rep(1,length=A) 
# for(j in 2:A)
# {
# 	lz[j] <- lz[j-1] * sa[j-1]	
# 	if(j == A)
# 	{
# 		lz[A] <- lz[A] / (1.0 - sa[A])	
# 	} 
# }


# p0 <- mx/za*oa
# p1 <- fe[1]*va[1,]/za*oa
# p2 <- fe[2]*va[2,]/za*oa

# ta <- oa - p0

# l1 <- sum(p1)/sum(ta)
# l2 <- sum(p2)/sum(ta)






#
# EQUILIBRIUM MODEL
# 
equilibriumModel <- function(fbar, type="YPR")
{

	lambda <- rep(1.0,length=length(gear))



	# 
	# Survivorship under fished conditions.
	#
	for(iter in 1:2)
	{
		fe <- fbar * lambda
		va <- getSlx(s1,s2,s3)
		za <- mx + colSums(fe*va)
		sa <- exp(-za)
		oa <- 1.0 - sa
		lz <- rep(1,length=A) 

		for(k in gear)
		{
			qa[k,] <- va[k,] * wa * oa / za		
			pa[k,] <- va[k,] * oa / za
			dlz[k,1] <- 0
		}
		
		# 
		# F multipliers (lambda) based on YPR or MPR
		# 
		qp     <- switch(type,YPR = qa, MPR = pa)
		phi.t  <- as.vector(lz %*% t(qp))
		lam.t  <- ak / (phi.t/sum(phi.t))
		lambda <- lam.t / mean(lam.t)

		cat(iter," lambda = ",lambda,"\n")
	}
	print("cheguei aqui")

	for(j in 2:A)
	{
		lz[j] <- lz[j-1] * sa[j-1]	
		if(j == A)
		{
			lz[A] <- lz[A] / (1.0 - sa[A])	
		} 

		for(k in gear)
		{
			dlz[k,j] <- sa[j-1]*(dlz[k,j-1]-lz[j-1]*va[k,j-1])
			if(j == A)
			{
				dlz[k,j] <- dlz[k,j]/oa[j] 
				- lz[j-1] * sa[j-1] * va[k,j] * sa[j] / (oa[j])^2
			}
		}
	}
	
	

	# incidence functions
	phi.m  <- as.vector(lz %*% t(pa))
	phi.e  <- as.double(lz %*% fa)
	phi.q  <- as.vector(lz %*% t(qa))
	dphi.e <- as.vector(fa %*% t(dlz))
	spr    <- phi.e / phi.E
	ispr   <- phi.E / phi.e
	dphi.q <- matrix(nrow=ng,ncol=ng)
	dre    <- rep(0,length=ng)
	for (k in gear) 
	{
		for(kk in gear)
		{
			va2 <- (va[k,] * va[kk,])
			t0  <- oa / za
			t1  <- lz*wa*va2/za
			t3  <- sa - t0

			dphi.q[k,kk] <- as.double(qa[k,] %*% dlz[k,] + t1 %*% t3)
		}
		dre[k]  <- ro * phi.E * dphi.e[k] / (phi.e^2 *(kappa-1))
	}


	# equilibrium recruitment
	re    <- max(0,ro * (kappa - ispr) / (kappa - 1.0))

	# mortality per recruit
	mpr   <- fe * phi.m

	# equilibrium catch
	ypr   <- fe * phi.q
	ye    <- re * ypr

	# Jacobian for yield
	dye   <- re * phi.q + fe * phi.q * dre + fe * re * dphi.q

	# Mitigation
	v     <- sqrt(diag(dye))
	M     <- dye / (v %o% v)
	
	# cat("SPR = ",   round(spr,3))
	# cat("\n Re = ", round(re,3))
	# cat("\n ye = ", round(ye,3))

	out <- list("spr"   = spr,
	         "fe"       = fe,
	         "ypr"      = as.double(ypr),
	         "yield"    = as.double(ye),
	         "mpr"      = as.vector(mpr),
	         "M"        = as.matrix(M)
	         )
	return(out)
}
equilibriumModel(fbar=0.1)






fnB <- function(log.fbar)
{
	fe  <- exp(log.fbar)
	em  <- as.list(equilibriumModel(fe,type="YPR"))
	spr <- em$spr

	f1  <- (spr - target_spr)^2
	return(f1)
}

fnC <- function(log.fbar)
{
	fe  <- exp(log.fbar)
	em  <- as.list(equilibriumModel(fe,type="MPR"))
	spr <- em$spr
	
	f1  <- (spr - target_spr)^2
	return(f1)
}


getFs <- function(TMAParams)
{
	fitB <- optim(log(fspr),fnB,method="BFGS",hessian=TRUE)
	fitC <- optim(log(fspr),fnC,method="BFGS",hessian=TRUE)
	
	fb   <- exp(fitB$par)
	fc   <- exp(fitC$par)

	runB <- equilibriumModel(fb,type="YPR")
	runC <- equilibriumModel(fc,type="MPR")

	dfB  <- with(runB,data.frame(method="B",fe,ypr,yield,mpr))
	dfC  <- with(runC,data.frame(method="C",fe,ypr,yield,mpr))
	df   <- rbind(dfB,dfC)
	return(df)
}


# equilibriumModel(fe)
# fd <- seq(0,0.2,by=0.01)
# fb <- seq(0,0.4,by=0.01)
# fs <- expand.grid(fd,fb)
# names(fs) <- c("F.d","F.b")

# s1    <- c( 8.0, 3.0)
# s2    <- c( 1.2, 1.2)
# s3    <- c( 0.0, 0.1)
# EQM <- as.data.frame(cbind(fs,t(apply(fs,1,equilibriumModel))))

# s1    <- c( 9.0, 3.0)
# s2    <- c( 0.2, 1.2)
# s3    <- c( 0.0, 0.1)
# EQF <- as.data.frame(cbind(fs,t(apply(fs,1,equilibriumModel))))

# s1    <- c( 8.0, 2.0)
# s2    <- c( 1.2, 3.2)
# s3    <- c( 0.0, 0.1)
# EQB <- as.data.frame(cbind(fs,t(apply(fs,1,equilibriumModel))))





# # 
# # GRAPHICS
# # 
# spr_brk <- seq(0.05,0.50,by=0.05)
# df <- data.frame(age=age,lx=lx,lz=lz,fa=fa)

# p <- ggplot(df)
# p <- p + geom_area(aes(x=age,y=lx*fa),alpha=0.3,fill="black")
# p <- p + geom_area(aes(x=age,y=lz*fa),alpha=0.3,fill="red")
# p <- p + labs(x="Age",y="Spawning biomass per recruit")
# print(p)

# # p <- ggplot(EQM,aes(F.d,F.b,z=ypr1/ypr2)) 
# # p <- p + stat_contour(alpha=0.5,aes(color=..level..))
# # p <- p + stat_contour(data=EQM,aes(F.d,F.b,z=ypr1/ypr2),breaks=1,color="black")
# # p <- p + stat_contour(data=EQM,aes(F.d,F.b,z=spr),breaks=0.35,color="red")
# # print(p)

# p <- ggplot(EQM,aes(F.d,F.b,z=spr)) 
# p <- p + stat_contour(breaks=spr_brk,alpha=0.5,aes(color=..level..))
# p <- p + stat_contour(breaks=0.35,colour="red",size=1.5,alpha=0.5)
# p <- p + stat_contour(data=EQF,aes(F.d,F.b,z=spr),breaks=0.35,colour="green",size=1.5,alpha=0.5)
# p <- p + stat_contour(data=EQB,aes(F.d,F.b,z=spr),breaks=0.35,colour="orange",size=1.5,alpha=0.5)
# p <- p + labs(x="Directed fishery F",y="Bycatch F",col="SPR")
# print(p)

# # not sure what this is good for.
# br <- seq(2,10,by=1)
# p <- ggplot(EQM) 
# p <- p + stat_contour(aes(F.d,F.b,z=yield1,color =..level..),breaks=br,alpha=0.1)
# p <- p + stat_contour(aes(F.d,F.b,z=yield2,size =..level..),breaks=br,alpha=0.1)
# p <- p + geom_contour(aes(F.d,F.b,z=spr),breaks=0.35,col="red",alpha=0.5,size=1.5)
# p <- p + stat_contour(data=EQB,aes(F.d,F.b,z=yield1,color=..level..),breaks=br,size=1.5)
# p <- p + stat_contour(data=EQB,aes(F.d,F.b,z=yield2,size=..level..),breaks=br,alpha=0.5)
# p <- p + stat_contour(data=EQF,aes(F.d,F.b,z=spr),breaks=0.35,colour="green",size=1.5,alpha=0.5)
# p <- p + stat_contour(data=EQB,aes(F.d,F.b,z=spr),breaks=0.35,colour="orange",size=1.5,alpha=0.5)
# # p <- p + stat_contour(breaks=0.35,colour="red")
# p <- p + labs(x="Directed fishery F",y="Bycatch F",col="Removals",size="Bycatch")
# print(p)

# p <- ggplot(EQM)
# p <- p + stat_contour(aes(F.d,F.b,z=yield1+yield2,col=..level..))
# p <- p + geom_contour(aes(F.d,F.b,z=spr),breaks=0.35,col="red",alpha=0.5,size=1.5)
# p <- p + labs(x="Directed fishery F",y="Bycatch F",col="Yield")
# print(p)

# # p <- ggplot(EQM,aes(spr,yield1,color=factor(yield2))) + geom_line()

# # v <- ggplot(SPR,aes(Halibut.Fishery,Bycatch.Fishery,z=SPR)) 
# # v <- v + stat_contour(breaks=seq(0.05,0.50,by=0.05))

# plot.va <- function()
# {
# 	for(k in gear)
# 	{
# 		va[k,] = gplogis(age,a=s2[k],b=s1[k],g=s3[k])
# 	}
# 	df  <- data.frame(Age=age,"Directed"=va[1,],"Bycatch"=va[2,])
# 	mdf <- melt(df,id.vars="Age")
# 	p <- ggplot(mdf) + geom_line(aes(Age,value,col=variable))
# 	p <- p + labs(x="Age",y="Selectivity-at-age")
# 	print(p)
# }
# m











