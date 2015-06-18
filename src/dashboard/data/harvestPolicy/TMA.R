# --------------------------------------------------------------------------- #
# TOTAL MORTALITY ALLOCATION
# AUTHOR: Steven Martell
# This is an equilibrium model for allocating the total mortality rate among 
# one or more fisheries.
# --------------------------------------------------------------------------- #


library(ggplot2)


A     <- 50				 	# Plus group Age
ro    <- 1.0             	# unfished equilibrium recruitment
h     <- 0.75			 	# steepness of the B-H SRR
kappa <- 4.0 * h / (1.0 - h)
m     <- 0.15			 	# instantaneous natural mortality rate
age   <- 1:A 				# age sequence
winf  <- 100 				# maximum avetage weigth
k     <- 1.5 * m 			# vb k
ahat  <- 11.5  			
ghat  <- 1.5
s1    <- c( 8.0, 3.0)
s2    <- c( 1.2, 1.2)
s3    <- c( 0.0, 0.1)
ng    <- length(s1) 		# number of fleets
gear  <- 1:ng


# 
# MANAGEMENT CONTROLS
# 
target_spr     <- 0.40
fe             <- c(0.15,0.05)
tma_allocation <- c(0.50,0.50)
slim           <- c(32,580)


# 
# AGE-SCHEDULE INFORMATION
# 
wa    <- winf*(1-exp(-k*age))^3
ma    <- plogis(age,ahat,ghat)
fa    <- wa * ma
lx    <- exp(-m)^(age-min(age)); lx[A] <- lx[A]/(1-exp(-m))
va    <- matrix(nrow=ng,ncol=A)
qa    <- matrix(nrow=ng,ncol=A)
dlz   <- matrix(nrow=ng,ncol=A)
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

# 
# Survivorship under fished conditions.
# 
for(k in gear)
{
	va[k,] = gplogis(age,a=s2[k],b=s1[k],g=s3[k])
}
za <- m + colSums(fe*va)
sa <- exp(-za)
oa <- 1.0 - sa
lz <- rep(1,length=A) 


for(j in 2:A)
{
	lz[j] <- lz[j-1] * sa[j-1]	
	if(j == A)
	{
		lz[A] <- lz[A] / (1.0 - sa[A])	
	} 
}










#
# EQUILIBRIUM MODEL
# 
equilibriumModel <- function(fe)
{


	# 
	# Survivorship under fished conditions.
	# 
	for(k in gear)
	{
		va[k,] = gplogis(age,a=s2[k],b=s1[k],g=s3[k])
	}
	za <- m + colSums(fe*va)
	sa <- exp(-za)
	oa <- 1.0 - sa
	lz <- rep(1,length=A) 

	for(k in gear)
	{
		qa[k,] <- va[k,] * wa * oa / za

		dlz[k,1] <- 0
	}

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

	# Jacobian for SPR
	dspr <- dphi.e
	print(dspr)
	# dspr <- matrix(nrow=ng,ncol=ng)
	# for (k in gear) 
	# {
	# 	# for (kk in gear) 
	# 	# {
	# 	# 	if(k == kk)
	# 	# 	{
	# 	# 		dspr[k,kk] = dphi.e[k] * phi.e
	# 	# 	}
	# 	# 	else
	# 	# 	{
	# 	# 		dspr[k,kk] = dphi.e[kk]
	# 	# 	}
	# 	# }
	# }
	# dspr   <-   dspr / phi.E
	# invJ   <- - solve(dspr)
	# fstp   <-   (spr*tma_allocation-target_spr) %*% invJ
	
	# print(jacobi)

	# equilibrium recruitment
	re    <- max(0,ro * (kappa - ispr) / (kappa - 1.0))


	# equilibrium catch
	ypr   <- fe * phi.q
	ye    <- re * ypr

	# cat("SPR = ",   round(spr,3))
	# cat("\n Re = ", round(re,3))
	# cat("\n ye = ", round(ye,3))

	out <- c("spr"   = spr,
	         "dspr"  = dspr,
	         "ypr"   = as.double(ypr),
	         "yield" = as.double(ye)
	         )
	return(out)
}


equilibriumModel(fe)
fd <- seq(0,0.2,by=0.01)
fb <- seq(0,0.4,by=0.01)
fs <- expand.grid(fd,fb)
names(fs) <- c("F.d","F.b")

EQM <- as.data.frame(cbind(fs,t(apply(fs,1,equilibriumModel))))

s1    <- c( 9.0, 3.0)
s2    <- c( 0.2, 1.2)
s3    <- c( 0.0, 0.1)
EQF <- as.data.frame(cbind(fs,t(apply(fs,1,equilibriumModel))))

s1    <- c( 8.0, 9.0)
s2    <- c( 1.2, 1.2)
s3    <- c( 0.0, 0.1)
EQB <- as.data.frame(cbind(fs,t(apply(fs,1,equilibriumModel))))





# 
# GRAPHICS
# 
spr_brk <- seq(0.05,0.50,by=0.05)
df <- data.frame(age=age,lx=lx,lz=lz,fa=fa)

p <- ggplot(df)
p <- p + geom_area(aes(x=age,y=lx*fa),alpha=0.3,fill="black")
p <- p + geom_area(aes(x=age,y=lz*fa),alpha=0.3,fill="red")
p <- p + labs(x="Age",y="Spawning biomass per recruit")
print(p)


p <- ggplot(EQM,aes(F.d,F.b,z=spr)) 
p <- p + stat_contour(breaks=spr_brk,alpha=0.5,aes(color=..level..))
p <- p + stat_contour(breaks=0.35,colour="red",size=1.5,alpha=0.5)
p <- p + stat_contour(data=EQF,aes(F.d,F.b,z=spr),breaks=0.35,colour="green",size=1.5,alpha=0.5)
p <- p + labs(x="Directed fishery F",y="Bycatch F",col="SPR")
print(p)

# not sure what this is good for.
br <- seq(2,10,by=1)
p <- ggplot(EQM) 
p <- p + stat_contour(aes(F.d,F.b,z=yield1,color =..level..),breaks=br,alpha=0.1)
p <- p + stat_contour(aes(F.d,F.b,z=yield2,size =..level..),breaks=br,alpha=0.1)
p <- p + geom_contour(aes(F.d,F.b,z=spr),breaks=0.35,col="red",alpha=0.5,size=1.5)
p <- p + stat_contour(data=EQB,aes(F.d,F.b,z=yield1,color=..level..),breaks=br,size=1.5)
p <- p + stat_contour(data=EQB,aes(F.d,F.b,z=yield2,size=..level..),breaks=br,alpha=0.5)
p <- p + stat_contour(data=EQF,aes(F.d,F.b,z=spr),breaks=0.35,colour="green",size=1.5,alpha=0.5)
p <- p + stat_contour(data=EQB,aes(F.d,F.b,z=spr),breaks=0.35,colour="orange",size=1.5,alpha=0.5)
# p <- p + stat_contour(breaks=0.35,colour="red")
p <- p + labs(x="Directed fishery F",y="Bycatch F",col="Removals",size="Bycatch")
print(p)

p <- ggplot(EQM)
p <- p + stat_contour(aes(F.d,F.b,z=yield1+yield2,col=..level..))
p <- p + geom_contour(aes(F.d,F.b,z=spr),breaks=0.35,col="red",alpha=0.5,size=1.5)
p <- p + labs(x="Directed fishery F",y="Bycatch F",col="Yield")
print(p)

# p <- ggplot(EQM,aes(spr,yield1,color=factor(yield2))) + geom_line()

# v <- ggplot(SPR,aes(Halibut.Fishery,Bycatch.Fishery,z=SPR)) 
# v <- v + stat_contour(breaks=seq(0.05,0.50,by=0.05))














