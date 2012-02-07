# -------------------------------------------------------------------------- ##
# Halitosis.R
# Written by Steve Martell,  UBC
# Date: Jan 8, 2012
# -------------------------------------------------------------------------- ##

# -------------------------------------------------------------------------- ##
# Libraries
require(Hmisc)
# -------------------------------------------------------------------------- ##

# Data and other constants
A	<- 30	# maximum age.
age	<- 1:A	# vector of ages

# Population parameters (female, male)
bo		<- 1.0				# unfished female spawning biomass
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

tsasm	<- 
function(fe=0, slim=0, dm=0.17)
{
	# A two sex age structured model. #
	
	# Age-schedule information
	# Survivorship matrix lx(sex, age)
	lx	<- sapply(age, function(age) exp(-m)^(age-1) )
	lx[, A] <- lx[, A]/(1-exp(-m))
	la	<- sapply(age, function(age) linf*(1-exp(-k*age)))
	wa	<- a*la^b
	ma	<- plogis(age, a50, k50)
	fa	<- t(t(wa)*ma)
	
	# Derivation of S-R parameters (Ricker model)
	# Re = Ro*[(log(kap)-log(phi.E/phi.e))/log(kap)]
	# Assume 50:50 sex ratio at birth
	kap		<- (5*h)^(5/4)
	phi.E	<- as.double(lx[1, ] %*% fa[1, ])
	ro		<- bo/phi.E
	
	# Length-based selectivity (convert weight-at-age to length-at-age)
	sc	<- t(apply(la,1,plogis,location=lhat,scale=ghat))
	std	<- cvlm*slim+1.e-30
	sr	<- t(apply(la,1,plogis,location=slim,scale=std))
	sd	<- 1-sr
	va	<- sc*(sr+sd*dm)		# age-specific probability of dying due to F
	
	# Age-specific total mortality rate.
	za	<- m+fe*va
	sa	<- exp(-za)
	qa	<- (sc*sr)*(1.-sa)/za	# fraction retained
	da	<- (sc*sd)*(1.-sa)/za	# fraction discarded
	
	# Survivorship under fished conditions
	lz	<- matrix(1, 2, A)
	for(i in 2:A)
	{
		lz[, i] <- lz[, i-1]*exp(-za[, i-1])
		if(i==A)
		{
			lz[, A] <- lz[, A]/(1-exp(-za[, i-1]))
		}
	}
	
	# Incidence functions
	phi.e	<- as.double(lz[1, ] %*% fa[1, ])
	
	# Equilibrium calculations
	t1		<- log(phi.E/(kap*phi.e))
	t2		<- (log(kap)*phi.e)
	re		<- max(0, -(t1*ro*phi.E)/t2)
	be		<- re * phi.e
	ye		<- sum(re * fe * lz*wa*qa)
	de		<- sum(re * fe * dm * lz*wa*da)
	
	# Per recruit calculations
	ypr		<- sum(fe * lz*wa*qa)
	spr		<- phi.e/phi.E
	
	# b<- seq(0, 2*bo, length=100)
	# r<- kap*b*exp(-log(kap)*b/bo)/phi.E
	# plot(b, r)
	# points(bo, ro, pch=20, col=2)
	# points(be, re, pch=20, col=3)
	return(list(re=re, be=be, ye=ye, 
		de=de, spr=spr, ypr=ypr, 
		dep=be/bo))
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
	fe	<- seq(0, 0.4, length=100)
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
# yield loss ration = (Wt. of yield loss)/(Wt. of bycatch).  The wt. of 
# the bycatch is straight forward. The yield loss is the difference between
# the yield obtained with discard mortality =0 and discard mortality =0.17


SPR <- .equil("spr", dm=dm)
YE  <- .equil("ye", dm=dm)
YE0 <- .equil("ye", dm=0)
DE	<- .equil("de", dm=dm)


# REPORT SECTION
par(mfcol=c(2, 2), las=1)
isolvl <- c(0.35, seq(0, 1, by=0.1))
isolwd <- c(2, rep(1, 11))
xl     <- "Fishing mortality"
yl     <- "Size limit (cm)"
plot(SPR,xlab=xl,ylab=yl,levels=isolvl,lwd=isolwd,main="Spawn biomass per recruit")
plot(YE ,xlab=xl,ylab=yl,main="Equilibrium yield")
plot(DE ,xlab=xl,ylab=yl,main="Discarded yield")
X = DE
X$Z = (YE0$Z-YE$Z)/(DE$Z)
plot(X, xlab="Fishing Mortality",ylab="Size limit",main="Yield loss ratio")



