# July 25,  2012 
# In response to issue 1 in the Git repository
# A simple example in R to develope an algorithm for
# solving for Fmsy for multiple fleets with allocation.
set.seed(123)
# POPULATION PARAMETERS
A     <- 20
age   <- 1:A
ro    <- 1
h     <- 0.75
vbk   <- 0.2
m     <- 1.5*vbk
theta <- c(ro=ro, h=h, vbk=vbk, m=m)

# FISHING FLEET PARAMETERS
ngear  <- 2
lambda <- rep(1/ngear, length=ngear)		# allocation to each fleet
ah     <- rep(3,ngear)
#ah     <- rlnorm(ngear, log(log(3)/m), 0.3)
gh     <- 0.5*ah
V      <- sapply(age, plogis, location=ah, scale=gh)
matplot(age, t(V), type="o")

equil <- function(theta, fe=0)
{
	with(as.list(theta), {
		# age-schedules
		wa   <- (1-exp(-vbk*age))^3
		fa   <- wa; fa[1:4] <- 0
		lx   <- exp(-m)^(age-1); lx[A] <- lx[A]/(1-exp(-m))
		phie <- sum(lx*fa)
		reck <- 4*h/(1-h)
		alfa <- reck/phie
		beta <- (reck-1)/(ro*phie)
		
		# fished conditions
		lz   <- rep(1, length=A)
		za   <- m + colSums(fe*V)
		sa   <- exp(-za)
		oa   <- (1-sa)
		qa   <-  t(t(V)*wa/za*oa)
		dlz.df    <- matrix(0, nrow=ngear, ncol=A)
		ddlz.df   <- matrix(0, nrow=ngear, ncol=A)
		dphif.df  <- rep(0, ngear)
		ddphif.df <- rep(0, ngear)
		dphiq.df  <- wa[1]*V[,1]/za[1]*(sa[1]-oa[1]/za[1])
		ddphiq.df <- -wa[1]*V[,1]^3*sa[1]/za[1] - 2*wa[1]*V[,1]^3*sa[1]/za[1]^2 + 2*wa[1]*V[,1]^3*oa[1]/za[1]^3
		for(i in 2:A)
		{
			lz[i]   <- lz[i-1]*sa[i-1]
			dlz.df[,i]  <- sa[i-1]*(dlz.df[,i-1] - lz[i-1]*V[,i-1])
			ddlz.df[,i] <- sa[i-1]*(ddlz.df[,i-1] + lz[i-1]*V[,i-1]^2)
			
			if(i==A)
			{
				lz[i]   <- lz[i]/oa[i]
				dlz.df[,i]  <- dlz.df[,i]/oa[i] - lz[i-1]*sa[i-1]*V[,i]*sa[i]/oa[i]^2
				
				ddlz.df[,i] <- ddlz.df[,i]/oa[i] + 2*lz[i-1]*V[,i-1]^2*sa[i-1]/oa[i]
				               + 2*lz[i-1]*V[,i-1]*sa[i-1]*V[,i]*sa[i]/oa[i]^2
				               + 2*lz[i-1]*sa[i-1]*V[,i]^2*sa[i]/oa[i]^3
				               + lz[i-1]*sa[i-1]*V[,i]^2*sa[i]/oa[i]^2
			}
			dphif.df  <- dphif.df  + fa[i]*dlz.df[,i]
			ddphif.df <- ddphif.df + fa[i]*ddlz.df[,i]
			
			dphiq.df  <- dphiq.df  + qa[,i]*dlz.df[,i] + lz[i]*wa[i]*V[,i]/za[i]*(sa[i]-oa[i]/za[i])
			
			ddphiq.df <- ddphiq.df + ddlz.df[,i]*qa[,i] + 2*dlz.df[,i]*wa[i]*V[,i]^2*sa[i]/za[i]
			             - 2*dlz.df[,i]*wa[i]*V[,i]^2*oa[i]/za[i]^2 - lz[i]*wa[i]*V[,i]^3*sa[i]/za[i]
			             - 2*lz[i]*wa[i]*V[,i]^3*sa[i]/za[i]^2 + 2*lz[i]*wa[i]*V[,i]^3*oa[i]/za[i]^3
		}
		
		
		phif    <- sum(lz*fa)
		phiq    <- colSums(lz*t(qa))
		dre.df  <- ro/(reck-1)*phie/phif^2*dphif.df
		ddre.df <- -2*ro*phie*dphif.df^2/(phif^3*(reck-1)) + ro*phie*ddphif.df/(phif^2*(reck-1))
		re      <- ro*(reck-phie/phif)/(reck-1)
		ye      <- fe*re*phiq
		dye.df  <- re*phiq + fe*phiq*dre.df + fe*re*dphiq.df
		ddye.df <- 2*dre.df*phiq + 2*re*dphiq.df + fe*ddre.df*phiq + 2*fe*dre.df*dphiq.df + fe*re*ddphiq.df
		#ddye.dj <- fe*ddre.df*phiq + 2*fe*dre.df*dphiq.df + fe*re*ddphiq.df
		
		U		<- phiq*outer(ddre.df,fe) + outer(dre.df*dphiq.df,2*fe) + re*outer(ddphiq.df,fe)
		
		# Jacobian matrix
		J       <- U
		diag(J) <- ddye.df
		
		invJ    <- -solve(J)
		delta   <- dye.df %*% invJ
		#cat(dye.df, "\t")
		#browser()
		print(sqrt(sum(dye.df^2)))
		cat("dye\n", dye.df, "\n")
		#cat("Jacobian\n", J, "\n")
		#browser()
		return(fe+delta)
	})
}

f<-rep(0.2, ngear)
equil(theta, f)

f <- rep(0.6*m, ngear)
for(iter in 1:100)
{
	f <- as.vector(equil(theta, f))
}
