## ____________________________________________________________________________ ##
##                            4WVX Herring CONTROLS
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
    1.6         -5.0    15       4       0       -5.0     15     #log_ro/msy
    0.750        0.2     1.0      4       3       5       5       #steepness/fmsy
##    -1.469      -5.0    0.0      2       1       -1.469  0.05    #log.m
    -1.60944      -5.0    0.0      -1       1       -1.469  0.05    #log.m
    6.0         -5.0    15       1       0       -5.0    15      #log_avgrec
    0.1         0.001   0.999   -3       3       25.0    100      #rho
    1.25        0.01    500     -3       4       1.01    1.01    #kappa (precision)
## ____________________________________________________________________________ ##
#
#
## _________________________SELECTIVITY PARAMETERS_____________________________ ##
## OPTIONS FOR SELECTIVITY:
##      1) logistic selectivity parameters
##      2) selectivity coefficients
##      3) a constant cubic spline with age-nodes
##      4) a time varying cubic spline with age-nodes
##      5) a time varying bicubic spline with age & year nodes.
##      6) fixed logistic (set isel_type=1, and estimation phase to -1)
## Gear 1 fishery:  Gear 2 survey
## isel_type
    5        1
## Age at 50% selectivity (logistic)
    3.5      2.05
## STD at 50% selectivity (logistic)
   1.0      0.5
## No. of age nodes for each gear (0 to ignore).
    5        0
## No. of year nodes for each gear (0 to ignore).
    20       0
## Estimation phase
    2        -2
## Penalty weight for 2nd differences w=1/(2*sig^2)
    12.5     12.5
## Penalty weight for dome-shaped selectivity 1=1/(2*sig^2)
    3.125    200.0
## ____________________________________________________________________________ ##
## ____________________________________________________________________________ ##
##                             Priors for Survey q                              ##
## ____________________________________________________________________________ ##
## nits  #number of surveys
    1
## priors 0=uniform density     1=normal density
    1
## prior log(mean)
    0
## prior sd
    1
## ____________________________________________________________________________ ##
#
#
## _______________________OTHER MISCELLANEOUS CONTROLS_________________________ ##
0           ## verbose ADMB output (0=off, 1=on)
2           ## recruitment model (1=beverton-holt, 2=ricker)
0.05        ## std in observed catches in first phase.
0.025       ## std in observed catches in last phase.
0           ## Assume unfished in first year (0=FALSE, 1=TRUE)
0.02        ## Minimum proportion to consider in age-proportions for dmvlogistic
1.30        ## Mean fishing mortality for regularizing the estimates of Ft
0.01        ## std in mean fishing mortality in first phase
5.00        ## std in mean fishing mortality in last phase
-3           ## phase for estimating m_deviations (use -1 to turn off mdevs)
0.01        ## std in deviations for natural mortality
12			## 11 number of estimated nodes for deviations in natural mortality
0.99        ## fraction of total mortality that takes place prior to spawning
1           ## 13 switch for age-composition likelihood (1=dmvlogistic,2=dmultinom)
## ____________________________________________________________________________ ##
#
## eofc
999
