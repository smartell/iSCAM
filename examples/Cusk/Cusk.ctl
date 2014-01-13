## ____________________________________________________________________________ ##
##                            CUSK CONTROLS
## ____________________ CONTROLS FOR ESTIMATED PARAMETERS _____________________ ##
##  Prior descriptions:
##   -0 uniform (0,0)
##   -1 normal (p1=mu,p2=sig)
##   -2 lognormal (p1=log(mu),p2=sig)
##   -3 beta (p1=alpha,p2=beta)
##   -4 gamma(p1=alpha,p2=beta)
## ____________________________________________________________________________ ##
7   ## npar
##  ival        lb      ub      phz   prior     p1      p2      parameter name
## ____________________________________________________________________________ ##
    2.5         -5.0    15       4     0         -5.0    15     #log_ro/msy
    0.60        0.2     1.0      4     3         30.0    20.0     #steepness/fmsy
    -1.89712   -5.0    0.0      -2     2         -1.89712  0.05 #log.m
    2.1         -5.0    15       1     0         -5.0    15     #log_avgrec
    7.20        -5.0    15      -1     0         -5.0    15     #log_recinit
    0.50        0.001   0.999    3     3         125.0   125.0  #rho
    5.00        0.01    500      3     4         0.01    0.01   #kappa (precision)
## ____________________________________________________________________________ ##
## ------------------------------------------------------------------------- ##
## CONTROL PARAMETERS FOR AGE/SIZE COMPOSITION DATA FOR na_gears             ##
## ------------------------------------------------------------------------- ##
## Likelihood type for each gear:
##     -1 : multivariate logistic (dmvlogistic)
##     -2 : multinomial, sample size based on input data
##     -3 : logistic_normal, 3 flavors, no autocorrelation, AR1, AR2.
## ------------------------------------------------------------------------- ##
## Number of columns == na_gears.
   1                        ## Gear Index
   3                        ## Likelihood type
   0.00                     ## Minimum proportion for aggregation
   0.00                     ## Small constant to add to comps & renormalize
   -2                       ## phase for phi1 estimation: bounded (-1,1) AR1
   -2                       ## phase for phi2 estimation: bounded (0,1)  AR2 
   -12345                   ## int check (-12345)
## ------------------------------------------------------------------------- ##
##                                                                           ##

## _________________________SELECTIVITY PARAMETERS_____________________________ ##
## OPTIONS FOR SELECTIVITY:
##      1) logistic selectivity parameters
##      2) selectivity coefficients
##      3) a constant cubic spline with age-nodes
##      4) a time varying cubic spline with age-nodes
##      5) a time varying bicubic spline with age & year nodes.
##      6) fixed logistic (set isel_type=1, and estimation phase to -1)
## Gear 1 fishery
## isel_type
    1		1		1		1
## Age at 50% selectivity (logistic)
    8.5		8.5		4.5		4.5
## STD at 50% selectivity (logistic)
    0.50	0.25	0.25	0.25
## No. of age nodes for each gear (0 to ignore).
    5		5		5		5
## No. of year nodes for each gear (0 to ignore).
    17      5		5		5
## Estimation phase
    -1      -1		-1		-1
## Penalty weight for 2nd differences w=1/(2*sig^2)
    12.5    12.5	12.5	12.5
## Penalty weight for dome-shaped selectivity 1=1/(2*sig^2)
    3.125	3.125	3.125	3.125
## ____________________________________________________________________________ ##

## ____________________________________________________________________________ ##
##                             Priors for Survey q                              ##
## ____________________________________________________________________________ ##
## nits  #number of surveys
	3
## priors	0=uniform density
## 			1=normal density
##			2=random walk in q
	2		2		2
## prior log(mean)
	0		0		0
## prior sd
	1		1		0
## ____________________________________________________________________________ ##

## _______________________OTHER MISCELLANEOUS CONTROLS_________________________ ##
0           ## 1  verbose ADMB output (0=off, 1=on)
1           ## 2  recruitment model (1=beverton-holt, 2=ricker)
0.05        ## 3  std in observed catches in first phase.
0.025       ## 4  std in observed catches in last phase.
1           ## 5  Assume unfished in first year (0=FALSE, 1=TRUE)
0.01		## 6  Minimum proportion to consider in age-proportions for dmvlogistic
0.20		## 7  Mean fishing mortality for regularizing the estimates of Ft
0.01		## 8  std in mean fishing mortality in first phase
5.00		## 9  std in mean fishing mortality in last phase
-1			## 10 phase for estimating m_deviations (use -1 to turn off mdevs)
0.1			## 11 std in deviations for natural mortality
3			## 12 number of estimated nodes for deviations in natural mortality
0.0	        ## 13 fraction of total mortality that takes place prior to spawning
1           ## 14 switch for age-composition likelihood (1=dmvlogistic,2=dmultinom)
## ____________________________________________________________________________ ##
## eofc
999