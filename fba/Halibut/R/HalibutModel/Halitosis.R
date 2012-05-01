# -------------------------------------------------------------------------- ##
# Halitosis.R
# Written by Steve Martell,  UBC
# Date: Jan 8, 2012
# 
# Feb 7, 2011.  Added multiple growth groups to allow for cumulative effects
# of size selective fishing.  Need to show how mean weigth-at-age decreases
# with increasing fishing mortality rates.
# -------------------------------------------------------------------------- ##

# -------------------------------------------------------------------------- ##
# Libraries
require(Hmisc)
# -------------------------------------------------------------------------- ##

# Data and other constants
A	<- 30	# maximum age.
G	<- 11	# number of growth groups
S	<- 2	# number of sexes
dim	<- c(A, G, S)
age	<- 1:A	# vector of ages
pg	<- dnorm(seq(-3, 3, length=G), 0, 1); pg <- pg/sum(pg)

# Population parameters (female, male)
bo		<- 100.0			# unfished female spawning biomass
h		<- 0.75				# steepness
m		<- c(0.15, 0.18)	# natural mortality rate
dm		<- 0.17				# discard mortality rate
a50		<- 10.91			# age at 50% maturity
k50		<- 1.406			# std at 50% maturity
a		<- 6.92e-6			# length-weight allometry
b		<- 3.24				# length-weight allometry
linf	<- c(145, 110)		# Range female 145-190, male 110-155 (cm)
k		<- c(0.1, 0.12)		# eyeballed growth pars from Clark & Hare 2002.

# Selectivity parameters (cm)
lhat	<- 97.132
ghat	<- 1/0.1667
slim	<- 81.28
cvlm	<- 0.1

# Halibut prices (10-20) (20-40) (40+)
# $6.75  $7.30  $7.50  In Homer Alaska.


tsasm	<- 
function(fe=0, slim=0, dm=0.17)
{
	# A two sex age structured model. #
	
	
	# Age-schedule information
	# Survivorship array lx(A, G, S)
	lx	<- array(0, dim)
	la	<- array(0, dim)
	wa	<- array(0, dim)
	fa	<- array(0, dim)
	pa	<- array(0, dim)	# price per pound 
	ma	<- plogis(age, a50, k50)
	for(i in 1:S)
	{
		lx[,,i]  <- exp(-m[i])^(age-1)
		lx[A,,i] <- lx[A,,i]/(1-exp(-m[i]))
		
		# growth
		'vonb'  <- function(linf,k) len <- linf*(1-exp(-k*age))
		dev     <- linf[i]*0.3
		linf.g  <- seq(linf[i]-dev, linf[i]+dev, length=G)
		la[,,i] <- sapply(linf.g, vonb,k=k[i])
		wa[,,i] <- a*la[,,i]^b
		
		# maturity (this assumes maturity at a fixed age)
		fa[,,i] <- ma*wa[,,i]
	}
	
	# price premiums based on fish weight
	pa[wa<10]  <- 0
	pa[wa>=10] <- 6.75
	pa[wa>=20] <- 7.30
	pa[wa>=40] <- 7.50
	
	# lx	<- sapply(age, function(age) exp(-m)^(age-1) )
	# lx[, A] <- lx[, A]/(1-exp(-m))
	# la	<- sapply(age, function(age) linf*(1-exp(-k*age)))
	# wa	<- a*la^b
	# ma	<- plogis(age, a50, k50)
	# fa	<- t(t(wa)*ma)
	
	# Derivation of S-R parameters (Ricker model)
	# Re = Ro*[(log(kap)-log(phi.E/phi.e))/log(kap)]
	# Assume 50:50 sex ratio at birth
	kap		<- (5*h)^(5/4)
	phi.E	<- sum(t(t(lx[,,1]*fa[,,1])*pg))
	#phi.E	<- as.double(lx[1, ] %*% fa[1, ])
	ro		<- bo/phi.E
	
	# Length-based selectivity (length-based -> age-based)
	sc	<- array(0, dim)
	sr	<- array(0, dim)
	sd	<- array(0, dim)
	va	<- array(0, dim)
	std	<- cvlm*slim+1.e-30
	for(i in 1:S)
	{
		sc[,,i]  <- plogis(la[,,i],location=lhat, scale=ghat)
		sr[,,i]  <- plogis(la[,,i],location=slim, scale=std)
		sd[,,i]  <- 1-sr[,,i]
		va[,,i]  <- sc[,,i]*(sr[,,i]+sd[,,i]*dm)
	}
	# sc	<- t(apply(la,1,plogis,location=lhat,scale=ghat))
	# std	<- cvlm*slim+1.e-30
	# sr	<- t(apply(la,1,plogis,location=slim,scale=std))
	# sd	<- 1-sr
	# va	<- sc*(sr+sd*dm)		# age-specific probability of dying due to F
	
	# Age-specific total mortality, survival, retention, and discard rate.
	za	<- array(0, dim)
	sa	<- array(0, dim)
	qa	<- array(0, dim)
	da	<- array(0, dim)
	for(i in 1:S)
	{
		za[,,i]  <- m[i] + fe*va[,,i]
		sa[,,i]  <- exp(-za[,,i])
		qa[,,i]  <- (sc[,,i]*sr[,,i]) * (1-sa[,,i])/za[,,i]
		da[,,i]  <- (sc[,,i]*sd[,,i]) * (1-sa[,,i])/za[,,i]
	}
	# za	<- m+fe*va
	# sa	<- exp(-za)
	# qa	<- (sc*sr)*(1.-sa)/za	# fraction retained
	# da	<- (sc*sd)*(1.-sa)/za	# fraction discarded
	
	# Survivorship under fished conditions lz(A, G, S)
	lz	<- array(1, dim)
	for(i in 1:S)
	{
		for(j in 2:A)
		{
			lz[j,,i] <- lz[j-1,,i]*exp(-za[j-1,,i])
		}
		lz[A,,i] <- lz[A,,i]/(1-exp(-za[A,,i]))
	}
	
	# lz	<- matrix(1, 2, A)
	# for(i in 2:A)
	# {
	# 	lz[, i] <- lz[, i-1]*exp(-za[, i-1])
	# 	if(i==A)
	# 	{
	# 		lz[, A] <- lz[, A]/(1-exp(-za[, i-1]))
	# 	}
	# }
	
	# Incidence functions
	phi.e	<- sum( t(lz[,,1]*fa[,,1])*pg )
	#phi.e	<- as.double(lz[1, ] %*% fa[1, ])
	
	# Equilibrium calculations
	t1		<- log(phi.E/(kap*phi.e))
	t2		<- (log(kap)*phi.e)
	re		<- max(0, -(t1*ro*phi.E)/t2)
	be		<- re * phi.e
	ye		<- 0
	de		<- 0
	ypr		<- 0
	wbar	<- rep(0, S)
	for(i in 1:S)
	{
		ye  <- ye + sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
		de	<- de + sum( re * fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
		ypr <- ypr + sum( fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
		
		# Average weigth of a 10-year old fish (female & male)
		tmp		<- t(lz[,,i]*wa[,,i])*pg
		tmpn	<- t(lz[,,i])*pg
		wbar[i] <- weighted.mean(wa[10,,i], tmpn[,10])
	}
	spr		<- phi.e/phi.E
	
	# Need to calculate landed value per recruit versus slim and fe
	# using the price and size categories from above.
	# How many millions of dollars in Halibut are discarded each year?
	landed.value	<- 0
	discard.value	<- 0
	for(i in 1:S)
	{
		t1 <- sum(re * fe * t(lz[,,i]*wa[,,i]*qa[,,i]*pa[,,i])*pg)
		landed.value  <- landed.value + t1
		t2 <- sum(re * fe * dm * t(lz[,,i]*wa[,,i]*da[,,i]*pa[,,i])*pg)
		discard.value <- discard.value + t2
	}
	
	# b<- seq(0, 2*bo, length=100)
	# r<- kap*b*exp(-log(kap)*b/bo)/phi.E
	# plot(b, r)
	# points(bo, ro, pch=20, col=2)
	# points(be, re, pch=20, col=3)
	return(list(re=re, be=be, ye=ye, 
		de=de, spr=spr, ypr=ypr, 
		dep=be/bo, wbar.f=wbar[1], wbar.m=wbar[2], 
		landed.value=landed.value, 
		discard.value=discard.value))
}

.equil	<-
function(arg="ye", dm=0.17)
{
	fn	<- function(fe=0, slim=0, dm=dm)
	{
		tmp <- tsasm(fe, slim, dm)
		idx <- which(names(tmp)==arg)
		return(as.double(tmp[idx]))
	}
	V	<- Vectorize(fn, c("fe", "slim"))
	fe	<- seq(0, 0.45, length=20)
	sl	<- seq(60, 100, length=20)
	Z	<- outer(fe, sl, V, dm)
	obj	<- list(x=fe, y=sl, Z=Z)
	class(obj) <- "isopleth"
	return(obj)
}

plot.isopleth <- 
function(obj, ...)
{
	# Plots the contour plot for the contour class
	x <- obj$x
	y <- obj$y
	z <- obj$Z
	
	contour(x, y, z, ...)
}

# A key question is for every pound of bycatch what is the corresponding
# yield loss to the directed fishery. This is computed by IPHC as the 
# yield loss ration = (Wt. of future yield loss)/(Wt. of bycatch).  The wt. of 
# the bycatch is straight forward. The yield loss is the difference between
# the yield obtained with discard mortality =0 and discard mortality =0.17


SPR <- .equil("spr", dm=dm)
SPR0<- .equil("spr", dm=0)
YE  <- .equil("ye", dm=dm)
YE0 <- .equil("ye", dm=0)
BE  <- .equil("be", dm=dm)
BE0 <- .equil("be", dm=0)
DE	<- .equil("de", dm=dm)
W.F <- .equil("wbar.f", dm=dm)
W.M <- .equil("wbar.m", dm=dm)
LV	<- .equil("landed.value", dm=dm)
DV	<- .equil("discard.value", dm=dm)

# REPORT SECTION
par(mfcol=c(1, 1), las=1)
isolvl <- c(0.35, seq(0, 1, by=0.1))
isolwd <- c(2, rep(1, 11))
xl     <- "Fishing mortality"
yl     <- "Size limit (cm)"
plot(SPR,xlab=xl,ylab=yl,levels=isolvl,lwd=isolwd,main="Spawn potential ratio")
plot(YE ,xlab=xl,ylab=yl,main="Equilibrium yield")
plot(DE ,xlab=xl,ylab=yl,main="Discarded yield")
X = DE
X$Z = (YE0$Z-YE$Z)/(DE$Z)
plot(X, xlab=xl,ylab=yl,main="Yield loss ratio")

#The following in the spawning biomass per recruit lost per 
#unit of discard. This should be the spawning biomass,  not SPR
SE = DE
SE$Z = (BE0$Z-BE$Z)/(DE$Z)
plot(SE, add=TRUE, col="blue", levels=seq(0, 10, by=.25))

E=DE
E$Z = YE$Z/(YE$Z+DE$Z)
plot(E, add=TRUE, col="red")

par(mfcol=c(2, 2))
plot(W.F, xlab=xl, ylab=yl, main="Mean weight of age-10 females")
plot(W.M, xlab=xl, ylab=yl, main="Mean weight of age-10 males")
plot(LV, xlab=xl, ylab=yl, main="Landed Value ($$)")
X = LV
X$Z = DV$Z/LV$Z
plot(X, xlab=xl, ylab=yl, main="Discard Value/Landed Value")
