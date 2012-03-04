# Feb 21, 2012
# OBJECTIVE:  to develop an numerical routine to calculate Fmsy 
#             for multiple fleets where there is an allocation
#             specified a priori for each of the participlating
#             fishing fleets.  This is key for developing reference
#             points in cases where there are two or more fleets
#             with different selectivity curves.  
# 
# 
# 
# PSEUDO CODE:
#            - Specify model dimensions & stock parameters
#            - Set up unfished conditions & stock-recruit parameters
#            - Write a function to calculate equilibrium yield for each gear

# DATA_SECTION:
A   = 20;			# Maximum age (min age =1)
G   = 3;			# Number of fishing gears
age = seq(1, A);	# vector of ages


# PARAMETER_SECTION
# Stock parameters
ro     = 1;			#Unfished recruitment
reck   = 12;		#Recruitment compensation
m      = 0.25;		#Instantaneous natural mortality rate
Winf   = 1;			#Asymptotic body weight (kg)
k      = 0.25;		#vonBert growth coefficient
am     = 7;			#age @ 50% maturity (assume logistic maturity ojive)
std.am = 1.5		#StDev in age at maturity

# Gear selectivity parameters & allocation (ak)
ah     = c(5, 9, 4);		#age @ 50% vulnerability in logistic selctivity
std.ah = c(2, 1, 3);		#StDev in logistic selectivity
ak     = c(0.5, 0.2, 0.3);	#Allocation (by weight) to each gear type.



# PROCEDURE_SECTION

# Unfished conditions
lx   = exp(-m)^(age-min(age))	#survivroship curve
lx[A]= lx[A]/(1.-exp(-m))
wa   = Winf*(1-exp(-k*age))^3	#mean weight at age
ma   = plogis(age,am,std.am)	#maturity at age schedule
phie = sum(lx*ma*wa)			#Eggs per recruit ~ or spawning stock biomass per recruit.
so   = reck/phie;				#maximum juvenile survival rate
b    = (reck-1)/(ro*phie);		#density dependent term in Beverton Holt model

# selectivity matrix va[1,G,1,A]
va   = sapply(age, plogis, location=ah, scale=std.ah)

# Equilibrium yield for a given fishing rate vector (fe)
"calc.ye" <-
function(fe=0)
{
	# The objective of this routine is to calculate 
	# the equilibrium yields for each fleet for a 
	# given vector of fishing rates (fe).  This can 
	# be used to numerically solve for the correct fe 
	# values to ensure the allocation constraints are 
	# satisfied.
	
	# Note that fe*phiq/sum(fe*phiq) is the per-recruit
	# yield for each gear.  At this point fe can be 
	# adjusted to satisfy thie allocation constraints.
	
	lz = rep(1, A)
	za = m + colSums(fe*va)		#age-specific total mortality
	sa = exp(-za)
	qa = t(va)*(1-sa)/za
	for(i in 2:A)
	{
		lz[i]  = lz[i-1]*sa[i-1]		
		if(i==A)
		{
			lz[i]  = lz[i]/(1-sa[i])
		}
	}
	
	
	phie = sum(lz*ma*wa)
	phiq = colSums(lz*wa*qa)
	re = max( c(0, so*phie-1/(b*phie)) )
	ye = re*fe*phiq
	return(ye)
}


# Use the following little algorithm to adjust 
# fishing mortality rate multipliers to satisfy allocation.
"equilibrium" <-
function(fe = 0, arg="ye")
{
	fm=rep(1/G, G)
	fk = fe*fm
	for(i in 1:10)
	{
		ye = calc.ye(fk)
		if(fe != 0)
		{
			pk = ye/sum(ye)
			fm=fm*(ak/pk)/sum(fm*(ak/pk))
			fk = fe*fm
		}
		if(sum(abs(ak-pk))<1e-6) break;
	}
	return(list(ye=ye, fk=fk))
}
