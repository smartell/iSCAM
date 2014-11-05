# runHalitosis_cpp.R
# library(inline)
# library(Rcpp)
# library(RcppArmadillo)
# # setwd("/Users/stevenmartell1/Documents/iSCAM-project/src/R/shiny")
sourceCpp("data/halitosis.cpp")



# |---------------------------------------------------------------------------|
# | MODEL DIMENSIONS
# |---------------------------------------------------------------------------|
# | Stock -> is a list object with all parameters and output.
   A	<- 30								# maximum age.
   G	<- 11								# number of growth groups
   S	<- 2								# number of sexes
   dim	<- c(A, G, S)				# array dimensions
   age	<- 1:A					 	# vector of ages
   pg	<- dnorm(seq(-1.96, 1.96, length=G), 0, 1); 
   pg  <- pg/sum(pg) 			# proportion assigned to each growth-type group.
   
   Stock <- list(A=A,G=G,S=S,dim=dim,age=age,pg=pg)


# |---------------------------------------------------------------------------|
# | Population parameters 476.8910
# |---------------------------------------------------------------------------|
   bo				<- 476.891							# unfished female spawning biomass
   h				<- 0.75				  				# steepness
   m        <- c(0.201826,0.169674) # natural mortality rate
   linf     <- c(158.86,101.525)		# asymptotic length (cm)
   vonk     <- c(0.0684,0.0842)			# vonk
   to       <- c(-5.2,-4.838)   		# time at zero length
   p        <- c(1.451,1.0424)  		# vonb Power parameter.
   cv       <- c(0.1,0.1)				# CV in length-at-age.
   a50      <- c(10.91,1000)     		# age at 50% maturity.
   k50      <- c(1.406,1)     		# std at 50% maturity.
	 a				<- rep(6.821e-6, 2) 		# length-weight allometry (Clark 1992)
	 b				<- rep(3.24, 2)					# length-weight allometry (CLark 1992)

   # dm		<- 0.16				# discard mortality rate
   cm		<- 0				# Size-dependent natural mortality rate (-0.5, 0.5)
   dev   <- seq(-1.96, 1.96, length=G)
   if(G==1) dev <- 0
   Stock <- c(Stock,list(bo=bo,h=h,m=m,linf=linf,vonk=vonk,to=to,p=p,cv=cv,a50=a50,k50=k50))
   Stock <- c(Stock,list(a=a,b=b,cm=cm,dev=dev))


# |---------------------------------------------------------------------------|
# | MODEL PROCEDURES
# |---------------------------------------------------------------------------|
# | mp1 -> is a list object with all parameters and output.

mp1 <- list(slim=82,
            ulim=1500,
            dmr=0.16,
            selex_50=55,
            selex_95=70,
            selex_bycatch50=45,
            selex_bycatch95=55,
            selex_asymptote=0.65,
            mid_points = seq(50,200, by = 2.5),
            fe = seq(0,0.5,by=0.01),
            bycatch = 5)



# |---------------------------------------------------------------------------|
# | Main Model Run
# |---------------------------------------------------------------------------|
# | 



equilibrium_model_cpp <- function(size_limit=c(32,100),
                              discard_mortality_rate=0.16,
                              spr_target=c(0.2,0.30),
                              selex_fishery=c(34,40),
                              selex_bycatch=c(24,40),
                              selex_asymptote=0.65,
                              num_bycatch=8,
                              five=0,
                              ten=5,
                              twenty=5,
                              forty=5,
                              linf_dev=0,
                              vonk_dev=0,
                              maternal_effect=1,
                              prefix="a")
{
	# -----------------------------------------
	# Inputs
	# -----------------------------------------

	# mimumum size limit
	size.limits = size_limit * 2.54 
	
	# discard mortality rate
	discard.mortality.rate = discard_mortality_rate

	# 50% selectivity at length
	selex.fishery  = selex_fishery * 2.54

	# 95% selectivity at length
	selex.bycatch  = selex_bycatch * 2.54

	# Bycatch mortality in Mlbs
	num.bycatch   = num_bycatch

	# Price vector
	price = c(five,ten,twenty,forty)

	Stock$linf = Stock$linf +(linf_dev/100*Stock$linf)
	Stock$vonk = Stock$vonk +(vonk_dev/100*Stock$vonk)
	Stock$mate = maternal_effect;

	# Procedure list
	mp1 <- list(slim=size_limit[1]*2.54,
            	ulim=size_limit[2]*2.54,
            	dmr=discard_mortality_rate,
            	selex_50=selex_fishery[1]*2.54,
	            selex_95=selex_fishery[2]*2.54,
	            selex_bycatch50=selex_bycatch[1]*2.54,
	            selex_bycatch95=selex_bycatch[2]*2.54,
	            selex_asymptote=selex_asymptote,
	            mid_points = seq(50,200, by = 2.5),
	            fe = seq(0,0.5,by=0.01),
	            bycatch = num_bycatch)
	print("Running cpp equilibrium model")

	mod <- 	new(Equilibrium,Stock)
	out <-	mod$calcLifeTable(Stock)
	df  <-  mod$runModel(mp1)
	df  <-  cbind(prefix=prefix,
	              ssb_limit=spr_target[1],
	              ssb_threshold=spr_target[2],
	              slim=size_limit[1]*2.54,
	              dm = discard_mortality_rate,
	              bycatch = num_bycatch,
	              df)

	return(df)
}




