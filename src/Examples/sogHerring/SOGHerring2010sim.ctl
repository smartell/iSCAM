## ____________________________________________________________________________ ##
##                            SOG HERRING CONTROLS
## ___________________CONTROLS FOR ESTIMATED PARAMETERS________________________ ##
##  Prior descriptions:
##                      -0 uniform (0,0)
##                      -1 normal (p1=mu,p2=sig)
##                      -2 lognormal (p1=log(mu),p2=sig)
##                      -3 beta (p1=alpha,p2=beta)
##                      -4 gamma(p1=alpha,p2=beta)
## ____________________________________________________________________________ ##
7   ## npar
##  ival        lb      ub      phz     prior    p1      p2      parameter name
## ____________________________________________________________________________ ##
    7.60        -5.0    15       4       0       -5.0    15     #log_ro/msy 
    0.80        0.2     1.0      4       0       1.1     1.1    #steepness/fmsy
    -0.7985077  -5.0    0.0     -3       0       -0.7985 0.2    #log.m
    7.60        -5.0    15       3       0       -5.0    15     #log_avgrec
    7.60        -5.0    15       3       0       -5.0    15     #log_recinit
    0.05        0.001   0.999   -3       3       1.01    1.01   #rho
    4999        0.01    5000    -3       4       1.01    1.01   #vartheta
## ____________________________________________________________________________ ##

## ____________________________________________________________________________ ##
## _________________________SELECTIVITY PARAMETERS_____________________________ ##
## OPTIONS FOR SELECTIVITY:
##      1) logistic selectivity parameters
##      2) selectivity coefficients
##      3) a constant cubic spline with age-nodes
##      4) a time varying cubic spline with age-nodes
##      5) a time varying bicubic spline with age & year nodes.
##      6) fixed logistic (set isel_type=1, and estimation phase to -1)
## Gear 1:3 fishery:  Gear 4-5 survey
## isel_type
    1        1          1       1       1
## Age at 50% selectivity (logistic)
    1.5      2.5        3.5     2.05    2.05
## STD at 50% selectivity (logistic)
    0.5      0.5       0.5     0.05    0.05
## No. of age nodes for each gear (0 to ignore).
    5        5          5       0       0
## No. of year nodes for each gear (0 to ignore).
    12       3          10      0       0
## Estimation phase
    -2      -2         -2      -2      -2
## Penalty weight for 2nd differences w=1/(2*sig^2)
    12.5     12.5       12.5    12.5    12.5
## Penalty weight for dome-shaped selectivity 1=1/(2*sig^2)
    3.125    200.0      200.0   200.0   200.0
## ____________________________________________________________________________ ##

## ____________________________________________________________________________ ##
##                             Priors for Survey q                              ##
## ____________________________________________________________________________ ##
## nits  #number of surveys
    2
## priors 0=uniform density     1=normal density
    0       0
## prior log(mean)
    0       0
## prior sd
    1       1
## ____________________________________________________________________________ ##

## _______________________OTHER MISCELLANEOUS CONTROLS_________________________ ##
0           ##   1 verbose ADMB output (0=off, 1=on)
1           ##   2 recruitment model (1=beverton-holt, 2=ricker)
0.015       ##   3 std in observed catches in first phase.
0.00001     ##   4 std in observed catches in last phase.
0           ##   5 Assume unfished in first year (0=FALSE, 1=TRUE)
0.00        ##   6 Minimum proportion to consider in age-proportions for dmvlogistic
0.05        ##   7 Mean fishing mortality for regularizing the estimates of Ft
0.01        ##   8 std in mean fishing mortality in first phase
5.00        ##   9 std in mean fishing mortality in last phase
-3          ##   10 phase for estimating m_deviations (use -1 to turn off mdevs)
0.01        ##   11 std in deviations for natural mortality
0.99        ##   12 fraction of total mortality that takes place prior to spawning
1           ##   13 switch for age-composition likelihood (1=dmvlogistic,2=dmultinom)
## ____________________________________________________________________________ ##


## eofc
999