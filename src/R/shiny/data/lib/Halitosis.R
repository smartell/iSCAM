#EquilibriumModel.R
# |---------------------------------------------------------------------------|
# | Libraries
# |---------------------------------------------------------------------------| 
   require(Hmisc)
   require(ggplot2)
   require(reshape2)
   # require(Riscam)
   require(grid)            #arrow function
   require(plyr)

   #source("Selex.R")
   #source("mytheme.R")

mytheme <- function (base_size = 12, base_family = "") 
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
        theme(axis.text = element_text(size = rel(0.8)), axis.ticks = element_line(colour = "black"), 
            legend.key = element_rect(colour = "grey80"), panel.background = element_rect(fill = "white", 
                colour = NA), panel.border = element_rect(fill = NA, 
                colour = "grey50"), panel.grid.major = element_line(colour = "grey90", 
                size = 0.2), panel.grid.minor = element_line(colour = "grey98", 
                size = 0.5), strip.background = element_rect(fill = "grey80", 
                colour = "grey50"), strip.background = element_rect(fill = "grey80", 
                colour = "grey50"),
                legend.position = 'bottom')
}


# |---------------------------------------------------------------------------|
# | Data and other constants
# |---------------------------------------------------------------------------|
# | Stock -> is a list object with all parameters and output.
   regArea <- "2B"
   A	<- 30					# maximum age.
   G	<- 3					# number of growth groups
   S	<- 2					# number of sexes
   dim	<- c(A, G, S)			# array dimensions
   age	<- 1:A					# vector of ages
   pg	<- dnorm(seq(-1.96, 1.96, length=G), 0, 1); 
   pg  <- pg/sum(pg) 			# proportion assigned to each growth-type group.
   
   Stock <- list(A=A,G=G,S=S,dim=dim,age=age,pg=pg)


# |---------------------------------------------------------------------------|
# | Population parameters 
# |---------------------------------------------------------------------------|
   bo		<- 476.891			# unfished female spawning biomass
   h		<- 0.75				# steepness
   m        <- c(0.201826,0.169674) # natural mortality rate
   linf     <- c(158.86,101.525)		# asymptotic length
   vonk     <- c(0.0684,0.0842)	# vonk
   to       <- c(-5.2,-4.838)   # time at zero length
   p        <- c(1.451,1.0424)  # vonb Power parameter.
   cv       <- c(0.1,0.1)		# CV in length-at-age.
   a50      <- rep(10.91,2)     # age at 50% maturity.
   k50      <- rep(1.406,2)     # std at 50% maturity.
	 a				<- rep(6.821e-6, 2) # length-weight allometry (Clark 1992)
	 b				<- rep(3.24, 2)		# length-weight allometry (CLark 1992)

   # dm		<- 0.16				# discard mortality rate
   cm		<- 0				# Size-dependent natural mortality rate (-0.5, 0.5)
   Stock <- c(Stock,list(bo=bo,h=h,m=m,linf=linf,vonk=vonk,to=to,p=p,cv=cv,a50=a50,k50=k50))
   Stock <- c(Stock,list(a=a,b=b,cm=cm))

# |---------------------------------------------------------------------------|
# | Commercial selectivities from Stewart 2012.                               
# |---------------------------------------------------------------------------|
   bin   <- seq(60, 130, by=10)
   CSelL <- matrix(
   (data=c(0,  0.0252252,  0.250685,  0.617268,  1,  1.36809,  1.74292,  2.12022, 
   		0,  0.0151914,  0.144236,  0.513552,  1,  1.48663,  1.97208,  2.45702)
   ), nrow=8, ncol=2)
   

   # THis is from the survey selectivity
   SSelL <- matrix(
	data=c(0,  0.285875,  0.579811,  0.822709,  1,  1.13862,  1.29717,  1.46621, 
		0,  0.246216,  0.454654,  0.68445,  1,  1.29683,  1.58379,  1.86742)
	, nrow=8, ncol=2)
	
   #CSelL <- t(t(CSelL)/apply(CSelL,2,max))
   CSelL <- t(t(SSelL)/apply(SSelL,2,max))
   #CSelL <- CSelL/max(CSelL)
   #CSelL <- SSelL/max(SSelL)
   slim	<- 81.28
   ulim	<- 1500
   cvlm	<- 0.1

   fe  <- seq(0, 0.50, by=0.01)#sequence of fishing mortality rates.
   sizelimit <- 2.54 * c(32)
   discardmortality <- c(0.16)
   bycatch  <- c(10)

   bycatchSel <- c(0, 0, 0.379083, 0.923116, 1, 0.748264, rep(0.650509,length=29))
# |---------------------------------------------------------------------------|


# |---------------------------------------------------------------------------|
# | Selex.R
# |---------------------------------------------------------------------------|
# | Compute the age-specific capture probability based on size selectivity.
# | 



# |---------------------------------------------------------------------------|
# | Calculate vulnerability-at-age given P(l|a) and P(l)
# |---------------------------------------------------------------------------|
# | Args: la -> vector of mean length-at-age
# |       sa -> vector of std in length-at-age
# |       pl -> probabilty of capturing an individual of a given length l.
# |       xl -> vector of mid points for length frequency distributions.
# | Pseudocode:
# |  1. Construct an age-length key given la and sa
# |  2. Compute va = P(l)*P(l|a)
.calcPage <- function(la, sa, pl, xl)
{
	A   <- length(la)
	L   <- length(xl)
	hbw <- 0.5*(xl[2] - xl[1])
	ALK <- matrix(0, nrow=L, ncol=A)
	
	ALK <-  sapply(xl+hbw,pnorm,mean=la,sd=sa) - sapply(xl-hbw,pnorm,mean=la,sd=sa)
	
	va  <-  as.vector(pl %*% t(ALK))
	return(va)
}

.va2vl <- function(la, sa, va, xl,linf,vbk,to,p)
{
	A   <- length(la)
	L   <- length(xl)

	 x = c(0,0.25,0.5,1.0)
	 c = c(1,2,2,2)
	wk = exp(-x^2/2)
	wk = wk/sum(c*wk)

	invfn <- function(mu) {
		- (-vbk*to+log(-exp(log(mu/linf)/p))+1)/vbk
	}
	k = -3:3

	
	# hbw <- 0.5*(xl[2] - xl[1])
	# ALK <- matrix(0, nrow=L, ncol=A)
	
	# ALK <-  sapply(xl+hbw,pnorm,mean=la,sd=sa) - sapply(xl-hbw,pnorm,mean=la,sd=sa)
	
	# vl  <-  as.vector(pa %*% ALK)
	# return(vl)
}


.calcALK <- function(la, sa, xl)
{
	A   <- length(la)
	L   <- length(xl)
	hbw <- 0.5*(xl[2] - xl[1])
	ALK <- matrix(0, nrow=L, ncol=A)
	
	ALK <-  sapply(xl+hbw,pnorm,mean=la,sd=sa) - sapply(xl-hbw,pnorm,mean=la,sd=sa)
	
	return(ALK)
}

.eplogis <- function (x, a, b, g) 
{
    sx = (1/(1 - g)) * ((1 - g)/g)^g * (exp(a * g * (b - x))/(1 + 
        exp(a * (b - x))))
}



# |---------------------------------------------------------------------------|
# | Age Schedule Information.
# |---------------------------------------------------------------------------|
.calcLifeTable <- function( Stock ) 
{

	gvonb <- function(t, linf, vbk, to, p=1)
	{
		l <- linf*(1.0-exp(-vbk*(t - to)))^p
		
		return(l)
	}

	with(Stock, {
		lx	   <- array(1, dim)
		la	   <- array(0, dim)
		sd_la  <- array(0, dim)
		wa	   <- array(0, dim)
		fa	   <- array(0, dim)
		pa     <- array(0, dim)
		M      <- array(0, dim)
		ma     <- plogis(age, a50, k50)

		for(i in 1:S)
		{
			# Lenght-at-age
			mu    <- gvonb(age,linf[i],vonk[i],to[i],p[i])
			sigma <- cv[i] * mu
			dev   <- seq(-1.96, 1.96, length=G)
			if(G==1) dev <- 0

			la[,,i]    <- sapply(dev,fn<-function(dev){la=mu+dev*sigma})
			sd_la[,,i] <- sqrt(1/G*(cv[i]*mu)^2)
			wa[,,i]    <- a[i]*la[,,i]^b[i]
			fa[,,i]    <- ma*wa[,,i]

			# Size dependent natural mortality rate 
			# M_l = M (l_a/l_r)^c
			l_r     <- 100
			delta   <- (la[,,i]/l_r)^cm / mean((la[,,i]/l_r)^cm)
			M[,,i]  <- m[i] * delta
			
			# Survivorship
			for(j in 2:A)
			{
				lx[j,,i] <- lx[j-1,,i]*exp(-M[j-1,,i])
			}
			lx[A,,i] <- lx[A,,i]/(1-exp(-M[A,,i]))
		}

		# price premiums based on fish weight
		pa[wa<10]  <- 0.00
		pa[wa>=10] <- 6.75
		pa[wa>=20] <- 7.30
		pa[wa>=40] <- 7.50

		Stock$lx	   = lx	
		Stock$la	   = la	
		Stock$sd_la    = sd_la
		Stock$wa	   = wa
		Stock$pa       = pa
		Stock$fa	   = fa	
		Stock$M        = M
		Stock$ma       = ma

		return(Stock)
	})

}

# |---------------------------------------------------------------------------|
# | Calculate size-based selectivities and joint capture probability
# |---------------------------------------------------------------------------|
# | This function calls .calcPage(la, sa, pl, xl) in Selex.R
.calcSelectivities <- function(Stock,slim=0,ulim=1500,cvlm=0.1,dm=0)
{
	# TODO Fix selectivities such that va is roughly the same with G=1 or G=11 groups.
	# Calculate capture probabilities
	with(Stock, {
		# Length-interval midpoints for integration
		xl  <- seq(5,200,by=2.5)
		# cat("S50 \t",selex_50)
		# Length-based selectivity (length-based -> age-based)
		sc	<- array(0, dim)  #size capture
		sr	<- array(0, dim)  #size retention
		sd	<- array(0, dim)  #size discard
		va	<- array(0, dim)  #retained
		vd	<- array(0, dim)  #bycatch fishery
		std	<- cvlm*slim+1.e-30	
		for(i in 1:S)
		{
			# pl       <- approx(bin, CSelL[,i], xl, yright=1, yleft=0)$y
			# pl       <- 1.0/(1.0+exp(-log(19)*(xl-selex_50)/(selex_90-selex_50)))
			pl       <- .plogis95(xl,selex_50,selex_90)
			sc[,,i]  <- .calcPage(la[,,i],sd_la[,,i],pl,xl)
			#sc[,,i]  <- approx(bin, CSelL[,i], la[,,i], yleft=0, yright=1)$y
			
			# retention proability
			pr       <-  plogis(xl, slim, 0.1) - plogis(xl, ulim, 0.1)
			sr[,,i]  <- .calcPage(la[,,i],sd_la[,,i],pr,xl)
			
			#sr[,,i]  <-  plogis(la[,,i],location=slim, scale=std) - plogis(la[,,i],location=ulim, scale=std)
			sd[,,i]  <- 1-sr[,,i]
			va[,,i]  <- sc[,,i]*(sr[,,i]+sd[,,i]*dm)
			
			# bycatch fishery selecitvity
			pd       <- plogis(xl, 40,  5) - (1-0.65)*plogis(xl, 100.3, 5.0)
			vd[,,i]  <- .calcPage(la[,,i],sd_la[,,i], pd, xl)
			

		}
		
		
		Stock$sc <- sc	# Length-based commercial selectivity.
		Stock$sr <- sr	# Age-specific retention probability.
		Stock$sd <- sd	# Age-specific discard probability.
		Stock$va <- va	# Joint capture probability.
		Stock$vd <- vd	# Discard probability in trawl fishery.
		Stock$dm <- dm  # Discard mortality rate
		
		return(Stock)
	})
}

.plogis95 <- function(x,s50,s90)
{
	1.0/(1.0+exp(-log(19)*(x-s50)/(s90-s50)))
}

# |---------------------------------------------------------------------------|
# | Calculate stock recruitment relationship.                       
# |---------------------------------------------------------------------------|
# |
.calcSRR <- function(Stock)
{
	with(Stock, {
		# Unfished SPR  (phi.E)
		phi.E	<- sum(t(t(lx[,,1]*fa[,,1])*pg))
		
		# Unfished recruitment (ro)
		ro		<- bo/phi.E
		
		# Beverton-Holt model
		kap <- 4*h/(1-h)

		# Ricker Model
		# kap <- (5*h)^(5/4)
		
		Stock$phi.E <- phi.E
		Stock$ro    <- ro
		Stock$kap   <- kap
		
		return(Stock)
	})
}

# |---------------------------------------------------------------------------|
# | Age-structure equilibrium model asem                            
# |---------------------------------------------------------------------------|
# | fe is the equilibrium fishing mortality rate.
.asem <- function(fe=0, Stock, ct=0)
{
	# | Psuedocode:
	# | 1. Calculate age-specific total mortality, retention, and discard rates
	# | 2. Calculate survivorship with fe>0.
	# | 3. Calculate equilibrium recruitment (re) and biomass (be)
	# | 4. Calculate yield per recruit,  spawning biomass per recruit,  yield, discards.
	# | 5. Calculate average weight-at-age.
	with(Stock, {
		# 1. Age-specific total mortality, survival, retention, and discard rate.
		za	<- array(0, dim)
		sa	<- array(0, dim)
		qa	<- array(0, dim)
		da	<- array(0, dim)
		ta  <- array(0, dim)	# yield per recruit in trawl discard fishery.
		
		bycatch <- ct
		# bapprox <- bo * 0.15/(0.15+fe)
		# fd      <- bycatch/bapprox
		fd      <- 0;
		#if(ct>0)
		#cat("fe = ", fe, " fd = ", fd, "\n")
		for(iter in 1:25)
		{


			for(i in 1:S)
			{
				za[,,i]  <- M[,,i] + fe*va[,,i] + fd*vd[,,i]
				sa[,,i]  <- exp(-za[,,i])
				qa[,,i]  <- (sc[,,i]*sr[,,i]) * (1-sa[,,i])/za[,,i]
				da[,,i]  <- (sc[,,i]*sd[,,i]) * (1-sa[,,i])/za[,,i]
				ta[,,i]  <- (vd[,,i])         * (1-sa[,,i])/za[,,i]
			}
			
			# 2. Survivorship under fished conditions lz(A, G, S)
			lz	<- array(1, dim)
			for(i in 1:S)
			{
				for(j in 2:A)
				{
					lz[j,,i] <- lz[j-1,,i]*exp(-za[j-1,,i])
				}
				lz[A,,i] <- lz[A,,i]/(1-exp(-za[A,,i]))
			}
			
			# 3. Calculate equilibrium recruitment and biomass
			phi.e	<- sum( t(lz[,,1]*fa[,,1])*pg )
			# Beverton-Holt model
			t1      <- phi.E/phi.e
			t2      <- (kap-t1)
			re      <- max(0, ro*t2/(kap-1))
			# Ricker model
			#t1		<- log(phi.E/(kap*phi.e))
			#t2		<- (log(kap)*phi.e)
			#re		<- max(0, -(t1*ro*phi.E)/t2)
			
			be		<- re * phi.e
			if(re ==0 ) break();
			
			# 3b. Calculate bycatch per recruit to determine F in bycatch fisheries.
			de      <- 0
			for(i in 1:S)
			{
				de	<- de + re * sum( t(lz[,,i]*wa[,,i]*ta[,,i])*pg )
			}
			if(bycatch < de)
				fd <- -log(1.0-bycatch/de)
			else
				fd <- -log(0.01) 
			# cat("fd = ",fd,"\t de = ",de,"re = ",re," \n")
		}
		# 4. Calculate yield per recruit, spawning biomass per recruit,  yield, discards.
		ye		<- 0
		ne    <- 0
		de		<- 0
		we    <- 0
		yev   <- 0
		dev   <- 0  # Value of wastage.
		byv   <- 0 	# Value of discarded bycatch.
		ypr		<- 0
		dpr   <- 0
		wpr   <- 0
		bpr   <- 0
		spr		<- phi.e/phi.E
		for(i in 1:S)
		{

			ye  <- ye + sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			ne  <- ne + sum( re * fe * t(lz[,,i]*qa[,,i])*pg )
			we	<- we + sum( re * fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
			de	<- de + sum( re * fe * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
			bpr <- bpr + sum( t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			ypr <- ypr + sum( fe * t(lz[,,i]*wa[,,i]*qa[,,i])*pg )
			dpr <- dpr + sum( fe * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
			wpr <- wpr + sum( fe * dm * t(lz[,,i]*wa[,,i]*da[,,i])*pg )
			n18 <- re * fe * t(lz[18,,i]*qa[,,i])*pg
			
			#
			# landed value of retained fish.
			#
			yev <- yev + sum( re * fe * t(lz[,,i]*wa[,,i]*qa[,,i]*pa[,,i])*pg )
			# 
			# value of discarded fish.
			# 
			dev <- dev + sum( re * fe * t(lz[,,i]*wa[,,i]*da[,,i]*pa[,,i])*pg )
			#
			# value of bycaught fish
			#
			byv <- byv + sum( re * fd * t(lz[,,i]*wa[,,i]*ta[,,i]*pa[,,i])*pg )
		}
		cbar <- ye / ne
		
		# 5. Calculate average weight-at-age
		wbar <- matrix(0, nrow=S, ncol=A)
		for(i in 1:S)
		{
			tmp      <- lz[,,i]/rowSums(lz[,,i])
			wbar[i,] <- rowSums(wa[,,i]*tmp)
		}

		
		# 6. Calculate average weight of the landed catch.

		Stock$lz   <- lz
		Stock$re   <- re
		Stock$be   <- be
		Stock$ye   <- ye
		Stock$de   <- de
		Stock$we   <- we
		Stock$fd   <- fd
		Stock$bpr  <- bpr
		Stock$ypr  <- ypr
		Stock$spr  <- spr
		Stock$dpr  <- dpr
		Stock$wpr  <- wpr
		Stock$wbar <- wbar
		Stock$cbar <- cbar
		Stock$yev  <- yev 
		Stock$dev  <- dev
		Stock$byv  <- byv

		
		return(Stock)
	})
}

# |---------------------------------------------------------------------------|
# | Run Model and Construct Data Frame for fe, slim, dm, bycatch combinations.
# |---------------------------------------------------------------------------|
# | This function creates a large data frame object with id.vars for fe, slim
# | dm and bycatch levels, and colums for each performance measure (ye, spr)
# | The order of operations is important here:
# | 2. .calcLifeTable
# | 3. .calcSelectivities
# | 4. .calcSRR (Stock recruitment relationship)
# | 5. .calcEquilibrium (data frame of values versus fe)

.runModel <- function(df)
{
	# print(df)
	with(as.list(df), {
		Stock$selex_50 = selex_50
		Stock$selex_90 = selex_90
		Stock$selex_bycatch50 = selex_bycatch50
		Stock$selex_bycatch90 = selex_bycatch90
		M1 <- .calcLifeTable(Stock)
		M1 <- .calcSelectivities(M1,slim=slim,ulim=ulim,cvlm=0.1,dm=dm)
		M1 <- .calcSRR(M1)
		M1 <- .asem(fe,M1,bycatch)
		
		out <- c(fe=fe,slim=slim,dm=dm,bycatch=bycatch,
		         Ye=M1$ye,De=M1$de,We=M1$we,Be=M1$be,Re=M1$re,
		         SPR=M1$spr,YPR=M1$ypr,DPR=M1$dpr,WPR=M1$wpr,BPR=M1$bpr,
		         Cbar=M1$cbar,EE=M1$ye/(M1$ye+M1$de),Fd = M1$fd,
		         YEv=M1$yev,DEv=M1$dev,BYv=M1$byv,
		         wbar_f=M1$wbar[1,18],wbar_m=M1$wbar[2,18])
		return(out)
	})
}




equilibrium_model <- function(size_limit=c(32,100),
                              discard_mortality_rate=0.16,
                              selex_fishery=c(34,40),
                              selex_bycatch=c(24,40),
                              num_bycatch=8,
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

	print("Running equilibrium model")
	
	df <- expand.grid(fe=fe,
	                  slim=size.limits[1],
	                  ulim=size.limits[2],
	                  dm=discard.mortality.rate,
	                  bycatch=num.bycatch,
	                  selex_50=selex.fishery[1],
	                  selex_90=selex.fishery[2],
	                  selex_bycatch50 = selex.bycatch[1],
	                  selex_bycatch90 = selex.bycatch[2])
	
	SA <- cbind(prefix=prefix,as.data.frame(t(apply(df,1,.runModel))))
	return(SA)
}

# Was going to add logistic 5095 curve to calcSelectivities.
# template<class T, class T2>
# 	const T plogis95(const T &x, const T2 &s50, const T2 &s95)
# 	{
# 		dvar_vector selex	= T2(1.0)/(T2(1.0)+(exp(-log(19)*((x-s50)/(s95-s50)))));
#     selex /= selex(selex.indexmax());	
# 		return selex;
# 	}


