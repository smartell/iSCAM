# --------------------------------------------------------------------------- #
# TOTAL MORTALITY ALLOCATION
# AUTHOR: Steven Martell
# This is an equilibrium model for allocating the total mortality rate among 
# one or more fisheries.

# TODO
# -[ ] Change algorithm to estimate fbar only.
# --------------------------------------------------------------------------- #

library(plyr)
library(ggplot2)

A     <- 25				# Plus group Age
ro    <- 1.0             # unfished equilibrium recruitment
h     <- 0.75			# steepness of the B-H SRR
kappa <- 4.0 * h / (1.0 - h)
m     <- 0.15			# instantaneous natural mortality rate
age   <- 1:A
mx    <- rep(m,length=A)
linf  <- 198.929
vbk   <- 0.078
to    <- 0.164
a     <- 0.00000692
b     <- 3.24
winf  <- a*linf^b
ahat  <- 11.5
ghat  <- 1.5

theta <- list(A=A,ro=ro,h=h,kappa=kappa,m=m,age=age,linf=linf,
              winf=winf,vbk=vbk,ahat=ahat,ghat=ghat)

# 
# Selectivity parameters for each gear.
# 
s1    <- c( 85.0, 25.0, 45.0, 85.0 )
s2    <- c( 12.0, 5.50, 11.2, 12.0 )
s3    <- c( 00.0, 0.03, 0.03, 00.0 )
slx   <- list(s1=s1,s2=s2,s3=s3)
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
target_spr <- 0.40
fspr       <- 0.10
slim       <- c(82,NA,82,00)
dmr        <- c(0.16,0.80,0.2,0.00)
ak         <- c(0.727,0.165,0.104,0.013)
# ak         <- rep(1/ng,length=ng)

# List object for the default harvest policy inputs.
hpDefaultInput <- list(allocation = ak,
                       target_spr = target_spr,
                       fspr       = fspr,
                       dmr        = dmr,
                       sizeLimit  = slim,
                       type       = "YPR",
                       slx        = slx)

hpSTQ  <- hpDefaultInput;
hpMSY  <- hpSTQ; hpMSY$allocation <- c(1,rep(0,length=ng-1))


# 
# AGE-SCHEDULE INFORMATION
# 
la    <- linf*(1-exp(-vbk*age))
wa    <- winf / (linf^3)*la^3
ma    <- plogis(age,ahat,ghat)
fa    <- wa * ma

lx    <- rep(1,length=A)
for(i in 2:A) lx[i] <- lx[i-1] * exp(-mx[i-1])
lx[A] <- lx[A]/(1-exp(-mx[A]))
phi.E <- as.double(lx %*% fa)

va    <- matrix(nrow=ng,ncol=A)
ra    <- matrix(nrow=ng,ncol=A)
qa    <- matrix(nrow=ng,ncol=A)
pa    <- matrix(nrow=ng,ncol=A)
dlz   <- matrix(nrow=ng,ncol=A)

getAgeSchedules <- function(theta)
{
	with(theta,{
		# 
		# AGE-SCHEDULE INFORMATION
		# 
		la    <- linf*(1-exp(-vbk*(age+to)))
		wa    <- a*la^b
		ma    <- plogis(age,ahat,ghat)
		fa    <- wa * ma

		# 
		# Age-dependent natural mortality (Lorenzen)
		# 
		t1 <- exp(vbk*(age+2))-1
		t2 <- exp(vbk*(age+1))-1
		sa <- (t1/t2)^(-m/vbk)
		mx=m/vbk*(log(t1)-log(t2))


		lx    <- rep(1,length=A)
		for(i in 2:A) lx[i] <- lx[i-1] * exp(-mx[i-1])
		lx[A] <- lx[A]/(1-exp(-mx[A]))
		phi.E <- as.double(lx %*% fa)
		
		ageSc <- list(la=la,wa=wa,ma=ma,fa=fa,lx=lx,mx=mx,phi.E=phi.E)
		theta <- c(theta,ageSc)

		return(theta)
	})
	
}




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

getSelex <- function(slxPars)
{
	# slxPars = s1,s2,s3,slim,dmr
	with(slxPars,{
		ngear <- length(s1)
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
	})
}


#
# EQUILIBRIUM MODEL
# 
equilibriumModel <- function(theta, type="YPR")
{



	with(theta,{
		# 
		# Survivorship under fished conditions.
		#
		fbar   <- fspr
		lambda <- rep(1.0,length=length(gear))
		for(iter in 1:5)
		{
			fe <- fbar * lambda
			za <- mx + colSums(fe*va)
			sa <- exp(-za)
			oa <- 1.0 - sa
			

			for(k in gear)
			{
				qa[k,] <- va[k,] * wa * oa / za		
				pa[k,] <- va[k,] * oa / za
				dlz[k,1] <- 0
			}
			
			lz    <- rep(1,length=A)
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

			# 
			# F multipliers (lambda) based on YPR or MPR
			# 
			qp     <- switch(type,YPR = qa, MPR = pa)
			phi.t  <- as.vector(lz %*% t(qp))
			lam.t  <- allocation / (phi.t/sum(phi.t))
			lambda <- lam.t / mean(lam.t)

			# cat(iter," lambda = ",lambda,"\n")
		

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
		Id    <- diag(1,ng)
		dye   <- re * (phi.q * Id) + fe * phi.q * dre + fe * re * dphi.q
		# browser()
		# Yield Equivalence
		v     <- sqrt(diag(dye))
		M     <- dye / (v %o% v)
		
		# cat("SPR = ",   round(spr,3))
		# cat("\n Re = ", round(re,3))
		# cat("\n ye = ", round(ye,3))

		out <- list("spr"   = spr,
		         "fe"       = fe,
		         "gear"     = gear,
		         "fspr"     = as.double(fbar),
		         "lambda"   = lambda,
		         "ypr"      = as.double(ypr),
		         "yield"    = as.double(ye),
		         "dye"      = as.vector(diag(dye)),
		         "mpr"      = as.vector(mpr),
		         "M"        = as.matrix(M)
		         )
		return(out)
	})
}





runModel <- function(theta,hp)
{
	# Initialize model
	theta <- getAgeSchedules(theta)
	theta <- c(theta,list(va=getSelex(hp)),hp)
	type  <- hp$type
	
	fn    <- function(log.fbar)
	{
		theta$fspr = exp(log.fbar)
		em <- as.list(equilibriumModel(theta,type))
		print(em$spr)
		return((em$spr - target_spr)^2)
	}
	# Perform nonlinear search to find Fspr = spr_target
	fit <- optim(log(theta$fspr),fn,method="BFGS")
	theta$fspr = exp(fit$par)
	EM  <- as.list(equilibriumModel(theta,type))
	return(EM)
}

runProfile <- function(theta,hp)
{
	# Initialize model
	theta <- getAgeSchedules(theta)
	theta <- c(theta,list(va=getSelex(hp)),hp)
	type  <- hp$type
	print(type)
	fn    <- function(log.fbar)
	{
		theta$fspr = exp(log.fbar)
		em <- as.list(equilibriumModel(theta,type))		
		return(em)
	}
	fbar <- seq(0,0.2,length=50)
	runs <- lapply(log(fbar),fn)
	df   <- ldply(runs,data.frame)
	
	return(df)
}

# .THEME <- theme(
#     panel.grid.minor = element_blank(), 
#     panel.grid.major = element_blank(),
#     panel.background = element_rect(fill="cyan",colour=NA),
#     plot.background = element_blank()
#    )

.THEME <- function (base_size = 12, base_family = "") 
{
    theme_grey(base_size = base_size, base_family = base_family)%+replace% 
    theme(axis.text = element_text(size = rel(0.8),colour="white"), 
          axis.ticks = element_line(colour = "black"), 
          strip.text = element_text(size = rel(0.9),colour="white"),
          legend.key = element_rect(colour = NA), 
          panel.background = element_rect(fill = "black",colour = NA), 
          panel.border = element_rect(fill = NA,colour = "grey90"), 
          panel.grid.major = element_line(colour = "grey50",size = 0.2), 
          panel.grid.minor = element_line(colour = "grey80",size = 0.5), 
          strip.background = element_rect(fill = "grey80",colour = "grey50", size = 0.2),
          plot.background  = element_rect(fill = "black")
          )
}

.GEAR <- c("Commercial","Bycatch","Sport","Other")
main <- function()
{
	rm <- runModel(theta,hpSTQ)
	df <- runProfile(theta,hpSTQ)

	p  <- ggplot(df,aes(fspr,yield)) + geom_line(aes(col=factor(gear)))
	print(p)
	p  <- ggplot(df,aes(spr,fe)) + geom_line(aes(col=factor(gear)))
	p  <- p + geom_vline(xintercept=target_spr,col="grey")
	print(p)

	p  <- ggplot(df,aes(spr,yield)) 
	p  <- p + geom_line(aes(colour=factor(.GEAR[gear])))
	p  <- p + labs(x="SPR",y="Yield")
	print(p + .THEME(28))
	return (p)
}

#