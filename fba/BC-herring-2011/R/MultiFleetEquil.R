# R-script for developing an altorithm for multifleet MSY calculation
# Author: Steven Martell
# Date:   June 8,  2012
set.seed(2313)
# Model parameters
A   <- 10
K   <- 3			# No. of fleets
age <- 1:A
#pk  <- c(0.3428571,  0.3800000,  0.2771429)
pk  <- c(0.1, 0.45, 0.45)  # Allocation among fleets
ro  <- 1.0
h   <- 0.75
m   <- 0.55
wa  <- (1-exp(-m/1.5*age))^3
fa  <- wa*plogis(age, 4, 0.3)
ah  <- c(1.5, 3, 4)#rnorm(K, 3.0, 1.5)
gh  <- rnorm(K, 0.5, 0.2)
V   <- sapply(age,plogis,location=ah,scale=gh)
matplot(t(V), type="o")

theta <- c(ro=ro, h=h, m=m, wa=wa, fa=fa, V=V)


# FUNCTIONS
equil  <- function(fe=0)
{
	with(list(theta), {
		lx     <- exp(-m*(age-1))
		lx[A]  <- lx[A]/(1.-exp(-m))
		
		phi.E  <- sum(lx*fa)
		kap    <- 4*h/(1-h)
		so     <- kap/phi.E
		beta   <- (kap-1)/(ro*phi.E)
		
		# survivorship under fished conditions.
		lambda <- pk/mean(pk)
		for(iter in 1:15)
		{
			za     <- m + colSums(fe*lambda*V)
			sa     <- exp(-za)
			oa     <- 1-sa
			qa     <- t(t(V)/za*oa)
		
			lz     <- rep(1, length=A)
			for(i in age)
			{
				if(i>min(age))
				{
					lz[i]  <- lz[i-1]*sa[i-1]
				}
			
				if(i==A)
				{
					lz[A]  <- lz[A]/oa[A]
				}
			}
			phi.e  <- sum(lz*fa)
			phi.q  <- colSums(wa*lz*t(qa))
			re     <- ro*(kap-phi.E/phi.e)/(kap-1.)
			ye     <- (fe*lambda)*re*phi.q
			
			# Update f-multiplier
			ak     <- ye/(sum(ye)+1.e-30)
			lambda <- lambda *(pk/(ak+1e-30))
			lambda <- lambda/mean(lambda)
			if(sum(abs(pk-ak)) <=1e-8) 
			{
				cat("Iteration ", iter, "\n")
				break
			}
			#lambda <- lambda*(ak/pk)/sum(lambda*(ak/pk)+1.e-30)
			#//fm=fm*(ak/pk)/sum(fm*(ak/pk))
			#print(c(iter, re, fe, sum(fe*lambda), ak/pk, ye))
		}
		return(list(Ye=sum(ye), fe=fe, ye=ye, lambda=lambda, re=re, be=re*phi.e, bo=ro*phi.E))
	})
}

.list2vector <- function(A, arg="Ye")
{
	ii <- match( arg, names(A[[1]]) )
	fn <- function(i) A[[i]][[ii]]
	i  <- 1:length(A)
	xx <- sapply(i, fn)
	return(xx)
}

fe <- seq(0,2.5,by=0.02)
AA <- sapply(fe,equil,simplify=F)
Ye <- .list2vector(AA, "Ye")
ye <- t(apply(.list2vector(AA,"ye"),2,cumsum))
plot(fe,Ye,type="n", ylim=c(0, max(ye)))
matlines(fe, ye, lty=1, lwd=2)
ix <- which.max(ye[,3])
MSY <- ye[ix, 3]
Fmsy <- fe[ix]
cat("MSY  = ", MSY, "\nFmsy =", Fmsy)
