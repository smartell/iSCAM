## ____________________________________________________________________________ ##
##                            PACIFIC HAKE CONTROLS
## ___________________CONTROLS FOR ESTIMATED PARAMETERS________________________ ##
##  Prior descriptions:
##                      -0 uniform (0,0)
##                      -1 normal (p1=mu,p2=sig)
##                      -2 lognormal (p1=log(mu),p2=sig)
##                      -3 beta (p1=alpha,p2=beta)
##                      -4 gamma(p1=alpha,p2=beta)
## ____________________________________________________________________________ ##
6   ## npar
##  ival        lb      ub      phz     prior    p1      p2      parameter name
## ____________________________________________________________________________ ##
    6.0         -5.0    15       4       0       -5.0    15      #log_ro/msy 
    0.65        0.2     1.0      4       3       3       2       #steepness/fmsy
    -1.096614   -5.0    0.0      2       1       -1.0966 0.05    #log.m
    6.0         -5.0    15       1       0       -5.0    15      #log_avgrec
    0.2         0.001   0.999    3       3       25.0    100     #rho
    1.25        0.01    500      3       4       1.01    1.01    #kappa (precision)
## ____________________________________________________________________________ ##


## _________________________SELECTIVITY PARAMETERS_____________________________ ##
## OPTIONS FOR SELECTIVITY:
##      1) logistic selectivity parameters
##      2) selectivity coefficients
##      3) a constant cubic spline with age-nodes
##      4) a time varying cubic spline with age-nodes
##      5) a time varying bicubic spline with age & year nodes.
##      6) fixed logistic (set isel_type=1, and estimation phase to -1)
## Gear 1:3 fishery:  Gear 4 survey
## isel_type
    1        1			1		1	1
## Age at 50% selectivity (logistic)
    1.5      2.5		2.5		2.05	2.05
## STD at 50% selectivity (logistic)
    1.0      0.5		0.2		0.05	0.05
## No. of age nodes for each gear (0 to ignore).
    5        5			5		0		0
## No. of year nodes for each gear (0 to ignore).
    12       3			10		0		0
## Estimation phase
    2        2			2		-2		-2
## Penalty weight for 2nd differences w=1/(2*sig^2)
    12.5     12.5		12.5	12.5	12.5
## Penalty weight for dome-shaped selectivity 1=1/(2*sig^2)
    3.125    200.0		200.0	200.0	200.0
## ____________________________________________________________________________ ##

## _______________________OTHER MISCELLANEOUS CONTROLS_________________________ ##
0			## verbose ADMB output (0=off, 1=on)
1			## recruitment model (1=beverton-holt, 2=ricker)
0.05		## std in observed catches in first phase.
0.025		## std in observed catches in last phase.
0			## Assume unfished in first year (0=FALSE, 1=TRUE)
## ____________________________________________________________________________ ##


## eofc
999