## ------------------------------------------------------------------------- ##
## CONTROL FILE TEMPLATE                                                     ##
## ------------------------------------------------------------------------- ##
##
##
## ------------------------------------------------------------------------- ##
## CONTROLS FOR LEADING PARAMETERS                                           ##
##  Prior descriptions:                                                      ##
##                      -0 uniform      (0,0)                                ##
##                      -1 normal       (p1=mu,p2=sig)                       ##
##                      -2 lognormal    (p1=log(mu),p2=sig)                  ##
##                      -3 beta         (p1=alpha,p2=beta)                   ##
##                      -4 gamma        (p1=alpha,p2=beta)                   ##
## ------------------------------------------------------------------------- ##
## npar
7
## ival         lb      ub      phz     prior   p1      p2      #parameter   ##
	7.0          0.0     10      4       0       0.0     12.0 #log_ro      ##
	0.75         0.2     1.0     4       3       3.00    2.00 #steepness   ##
	-1.89712    -3.0     2.0    -4       1       -1.203  0.15 #log_m g&b   ##
	6.5          0.0     10      1       0       0.0     10   #log_avgrec  ##
	7.5          0.0     10      1       0       0.0     10   #log_recinit ##
	0.00         0.00    1.00   -3       3       3.00    5.00 #rho         ##
	0.6          0.01    5.0    -3       4       1.01    1.01 #sigma_r     ##

## ------------------------------------------------------------------------- ##
##
## ———————————————————————————————————————————————————————————————————————————————————— ##
## CONTROL PARAMETERS FOR AGE/SIZE COMPOSITION DATA FOR na_gears                        ##
## ———————————————————————————————————————————————————————————————————————————————————— ##
## Likelihood type for each gear:
##     -1 : multivariate logistic (dmvlogistic)
##     -2 : multinomial, sample size based on input data
##     -3 : logistic_normal, no autocorrelation, AR1, AR2.
##     -4 : logistic_normal, AR1
##     -5 : logistic_normal, AR2
##     -6 : multinomial with estimated effective sample size.
## ———————————————————————————————————————————————————————————————————————————————————— ##
## Number of columns == na_gears.
	1        1       6       6       6       6      ## : Gear Index
	2        2       2       2       2       2      ## : Likelihood type
	0.000    0.000   0.000   0.000   0.000   0.000  ## : Minimum proportion for aggregation & compression
	0.0000   0.0000  0.0000  0.0000  0.0000  0.0000 ## : Small constant to add to comps & renormalize
	-1       -1      -1      -1      -1      -1     ## : phase for log_age_tau2 estimation.
	-2       -2      -2      -2      -2      -2     ## : phase for phi1 estimation : bounded (-1,1) AR1
	-2       -2      -2      -2      -2      -2     ## : phase for phi2 estimation : bounded (0,1)  AR2
	-2       -2      -2      -2      -2      -2      ## : phase for degrees of freedom for student T.
	-12345                 ## : int check (-12345)
## ———————————————————————————————————————————————————————————————————————————————————— ##



# ## ———————————————————————————————————————————————————————————————————————————————————— ##
# ## SELECTIVITY CONTROLS                                                                 ##
# ## ———————————————————————————————————————————————————————————————————————————————————— ##
# ## - Each gear must have a selectivity and retention curve.
# ## • Index       = gear index for selectivity curve.
# ## • sel_type    = type of selectivity function (see Legend).
# ## • sel_mu      = mean age/length 50% selectivity.
# ## • sel_sd      = std in 50% selectivity
# ## • sex_dep     = 0 -> no;  1 -> offset for sex 2.
# ## • size_nodes  = # of nodes for age/size cubic spline.
# ## • year_nodes  = # of nodes for time varying bicubic spline.
# ## • phz_mirror  = phase of estimation (-ve phase to mirror selextivity index)
# ## • lam1        = penalty weight for 2nd differences (w = 1/(2•sig^2)).
# ## • lam2        = penalty weight for dome-shaped selectivity.
# ## • lam3        = penalty weight for time-varying selectivity.
# ## • start_block = year index for first year of selectivity curve.
# ## ———————————————————————————————————————————————————————————————————————————————————— ##
# ## sel_nBlocks    ret_nBlocks    ## Gear Description.
#    1              1              ## Commercial retained
#    1              1              ## Commercial discards
#    1              1              ## Bycatch in non-directed fisheries.
#    1              1              ## Sport
#    1              1              ## Personal use
#    1              1              ## Setline survey.
# ## ———————————————————————————————————————————————————————————————————————————————————— ##
# ## slx_dControls
# ##        sel   sel  sel  sex  size   year  phz                       start  end        ##
# ## Index  type  mu   sd   dep  nodes  nodes mirror lam1  lam2  lam3 | block  block      ##
# ## ———————————————————————————————————————————————————————————————————————————————————— ##
# ## Selectivity P(capture of all size/age)
#    1      3     10   2.0  1    5      0     2      0.0   0.0   0.0    1888
#    2      1     3.0  1.5  0    5      0    -2      0.0   0.0   0.0    1888
#    3      1     4.0  2.5  0    5      0    -3      0.0   0.0   0.0    1888
#    4      1     4.0  2.0  0    5      0    -4      0.0   0.0   0.0    1888
#    5      3     10   2.0  1    5      0    -1      0.0   0.0   0.0    1888
#    6      5     3.0  2.0  1    5      0     2      0.0   0.0   0.0    1888
# ## ———————————————————————————————————————————————————————————————————————————————————— ##
# ## ret_dControls
# ##        sel   sel  sel  sex  size   year  phz                       start  end        ##
# ## Index  type  mu   sd   dep  nodes  nodes mirror lam1  lam2  lam3 | block  block      ##
# ## ———————————————————————————————————————————————————————————————————————————————————— ##
# ## Retention P(retaining size/age)
#   -1      3     10   2.0  1    5      0     2      0.0   0.0   0.0    1888
#   -2      1     3.0  1.5  0    5      0    -2      0.0   0.0   0.0    1888
#   -3      1     4.0  2.5  0    5      0    -3      0.0   0.0   0.0    1888
#   -4      1     4.0  2.0  0    5      0    -4      0.0   0.0   0.0    1888
#   -5      3     10   2.0  1    5      0    -1      0.0   0.0   0.0    1888
#   -6      5     3.0  2.0  1    5      0     2      0.0   0.0   0.0    1888










##  TO BE DEPRECATED
## ------------------------------------------------------------------------- ##
## SELECTIVITY PARAMETERS Columns for gear                                   ##
## OPTIONS FOR SELECTIVITY (isel_type):                                      ##
##      1) logistic selectivity parameters                                   ##
##      2) selectivity coefficients                                          ##
##      3) a constant cubic spline with age-nodes                            ##
##      4) a time varying cubic spline with age-nodes                        ##
##      5) a time varying bicubic spline with age & year nodes.              ##
##      6) fixed logistic (set isel_type=6, and estimation phase to -1)      ##
##      7) logistic function of body weight.                                 ##
##      8) logistic with weight deviations (3 parameters)                    ##
##      11)logistic selectivity with 2 parameters based on mean length       ##
##      12)length-based selectivity coefficients with spline interpolation   ##
##      sig=0.05 0.10 0.15 0.20 0.30 0.40 0.50                               ##
##      wt =200. 50.0 22.2 12.5 5.56 3.12 2.00                               ##
## ------------------------------------------------------------------------- ##
  3      1     1     1     3     5     # 1  -selectivity type ivector(isel_type) for gear
  82     3.0   4.0   4.0   3.0   65.0  # 2  -Age/length at 50% selectivity (logistic)
  5.5    1.5   2.5   2.0   2.5   5.5   # 3  -STD at 50% selectivity (logistic)
  5      0     0     0     5     5     # 4  -No. of age nodes for each gear (0=ignore)
  0      0     0     0     0     5     # 5  -No. of year nodes for 2d spline(0=ignore)
  2     -2    -3    -4    -1     2     # 6  -Phase of estimation (-1 for fixed)
  0.0    0.0   0.0   0.0   0.0   12.5  # 7  -Penalty wt for 2nd differences w=1/(2*sig^2)
  200.0  200.0 200.0 200.0 200.0 200.0 # 8  -Penalty wt for dome-shaped w=1/(2*sig^2)
  50.5   12.5  12.5  12.5  50.5  12.5  # 9  -Penalty wt for time-varying selectivity
  1      1     1     1     1     1     # 10 -n_sel_blocks (number of selex blocks)
## ------------------------------------------------------------------------- ##
## Start year of each time block: 1 row for each gear
1888 
1888
1888  
1888  
1888 
1970
## 
## ———————————————————————————————————————————————————————————————————————————————————— ##
## TIME VARYING NATURAL MORTALIIY RATES                                                 ##
## ———————————————————————————————————————————————————————————————————————————————————— ##
## TYPE: 
##      0 = constant natural mortality
##      1 = Random walk (deviates constrained by variance in M)
##      2 = Cubic Spline (deviates constrined by nodes & node-placement)
## ———————————————————————————————————————————————————————————————————————————————————— ##
  0
## Phase of estimation
  3
## STDEV in m_dev for Random walk
  0.01
## Number of nodes for cubic spline
  6
## Year position of the knots (vector must be equal to the number of nodes)
  1996 1998 1999 2001 2002 2013
## ———————————————————————————————————————————————————————————————————————————————————— ##
##
##
## ———————————————————————————————————————————————————————————————————————————————————— ##
## ABUNDANCE OBSERVATION MODELS
## ———————————————————————————————————————————————————————————————————————————————————— ##
## QTYPE:
##		0 = FIXED SURVEY Q (specify log(mean) for prior log(mean))
##		1 = CONSTANT Q     (use MLE for q and optional informative prior)
##      2 = RANDOM WALK Q  (use prior mean & sd for penalized random walk)
## ———————————————————————————————————————————————————————————————————————————————————— ##
2	     				# -number of surveys (n_it_nobs) 
 2   2					# -QTYPE (see legend above)
 0   0					# -prior log(mean)
 0   0					# -prior sd (set to 0 for uniformative prior)
 1   1                  # -Estimation phase
## ———————————————————————————————————————————————————————————————————————————————————— ##
##

## ------------------------------------------------------------------------- ##
## OTHER MISCELANEOUS CONTROLS                                               ##
## ------------------------------------------------------------------------- ##
0         # 1  -verbose ADMB output (0=off, 1=on)
1         # 2  -recruitment model (1=beverton-holt, 2=ricker)
1.0       # 3  -std in observed catches in first phase.
1.0       # 4  -std in observed catches in last phase.
0         # 5  -Assume unfished in first year (0=FALSE, 1=TRUE)
1.00      # 6  -Maternal effects power parameter (1.0 = no maternal effects)
0.30      # 7  -Mean fishing mortality for regularizing the estimates of Ft
0.10      # 8  -std in mean fishing mortality in first phase
2.00      # 9  -std in mean fishing mortality in last phase
-3        # 10 -DEPRECATED phase for estimating m_deviations (use -1 to turn off mdevs)
0.1       # 11 -DEPRECATED std in deviations for natural mortality 
12        # 12 -DEPRECATED number of estimated nodes for deviations in natural mortality
1.00      # 13 -fraction of total mortality that takes place prior to spawning
82        # 14 -number of perspective years to start assessment at.
0         # 15 -switch for IFD distribution in selectivity simulations
##
## ------------------------------------------------------------------------- ##
## MARKER FOR END OF CONTROL FILE (eofc)
## ------------------------------------------------------------------------- ##
999