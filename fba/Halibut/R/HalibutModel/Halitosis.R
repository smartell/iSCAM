# -------------------------------------------------------------------------- ##
# Halitosis.R
# Written by Steve Martell,  UBC
# Date: Jan 8,  2012
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
a50		<- 10.91			# age at 50% maturity
k50		<- 1.406			# std at 50% maturity
a		<- 6.92e-6			# length-weight allometry
b		<- 3.24				# length-weight allometry
linf	<- c(145, 110)		# Range female 145-190, male 110-155 (cm)
k		<- c(0.1, 0.12)		# eyeballed growth pars from Clark & Hare 2002.

# Selectivity parameters (cm)
lhat	<- 80
ghat	<- 08

tsasm	<- function()
{}
	# A two sex age structured model. #
	
	# Age-schedule information
	# Survivorship matrix lx(sex, age)
	lx	<- sapply(age, function(age) exp(-m)^age)
	la	<- sapply(age, function(age) linf*(1-exp(-k*age)))
	wa	<- a*la^b
	ma	<- plogis(age, a50, k50)
	fa	<- t(t(wa)*ma)
	
	# Derivation of S-R parameters (Ricker model)
	# Re = Ro*[(log(kap)-log(phi.E/phi.e))/log(kap)]
	# Assume 50:50 sex ratio at birth
	kap		<- (5*h)^(5/4)
	phi.E	<- as.double(lx[1, ] %*% fa[1, ])
	ro		<- 2.0*bo/phi.E
	
	# Length-based selectivity (convert weight-at-age to length-at-age)
	sa	<- t(apply(la,1,plogis,location=lhat,scale=ghat))
	
#}










