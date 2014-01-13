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
## npar                                                                      ##
7                                                                            ##
## ival         lb      ub      phz     prior   p1      p2      #parameter   ##
   7.1         -5.0    15.0     4       0      -5.0    15.0     #log_ro      ##
   0.95         0.2     1.0     4       3       1.01    1.01    #steepness   ##
  -1.06421     -5.0     0.0    -2       2      -1.06421 0.05    #log_m       ##
   7.1         -5.0    15.0     1       0      -5.0    15.0     #log_avgrec  ##
   7.20        -5.0    15.0    -1       0      -5.0    15.0     #log_recinit ##
   0.50         0.001   0.999   3       3     125.0   125.0     #rho         ##
   5.00         0.01  500.0     3       4       0.01    0.01    #vartheta    ##
## ------------------------------------------------------------------------- ##
##
## ------------------------------------------------------------------------- ##
## CONTROL PARAMETERS FOR AGE/SIZE COMPOSITION DATA FOR na_gears             ##
## ------------------------------------------------------------------------- ##
## Likelihood type for each gear:
##     -1 : multivariate logistic (dmvlogistic)
##     -2 : multinomial, sample size based on input data
##     -3 : logistic_normal, 3 flavors, no autocorrelation, AR1, AR2.
## ------------------------------------------------------------------------- ##
## Number of columns == na_gears.
   1 						## Gear Index
   3                        ## Likelihood type
   0.00                     ## Minimum proportion for aggregation
   0.00                     ## Small constant to add to comps & renormalize
   -2                       ## phase for phi1 estimation: bounded (-1,1) AR1
   -2                       ## phase for phi2 estimation: bounded (0,1)  AR2 
   -12345                   ## int check (-12345)
## ------------------------------------------------------------------------- ##
##																			 ##
##
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
##      11) logistic selectivity with 2 parameters based on mean length      ##
##      12) length-based selectivity coefficients with spline interpolation  ##
##      sig=0.05 0.10 0.15 0.20 0.30 0.40 0.50                               ##
##      wt =200. 50.0 22.2 12.5 5.56 3.12 2.00                               ##
## ------------------------------------------------------------------------- ##
	1		        # 1  -selectivity type ivector(isel_type) for gear
	3.5		        # 2  -Age/length at 50% selectivity (logistic)
	0.25	        # 3  -STD at 50% selectivity (logistic)
	7				# 4  -No. of age nodes for each gear (0=ignore)
	12				# 5  -No. of year nodes for 2d spline(0=ignore)
	3				# 6  -Phase of estimation (-1 for fixed)
	12.5			# 7  -Penalty wt for 2nd differences w=1/(2*sig^2)
	12.5 			# 8  -Penalty wt for dome-shaped w=1/(2*sig^2)
	12.5			# 9  -Penalty wt for time-varying selectivity
	1               # 10 -n_sel_blocks (number of selex blocks)
## ------------------------------------------------------------------------- ##
1965
##
##
##
## ------------------------------------------------------------------------- ##
## PRIORS FOR SURVEY Q                                                       ##
## Prior type:                                                               ##
##			0 - uninformative prior                                          ##
##			1 - normal prior density for log(q)                              ##
##			2 - random walk in q                                             ##
## ------------------------------------------------------------------------- ##
1					# -number of surveys (nits)
0					# -prior type (see legend above)
0					# -prior log(mean)
0					# -prior sd
## ------------------------------------------------------------------------- ##
##

## ------------------------------------------------------------------------- ##
## OTHER MISCELANEOUS CONTROLS                                               ##
## ------------------------------------------------------------------------- ##
0     # 1  -verbose ADMB output (0=off, 1=on)                                ##
1     # 2  -recruitment model (1=beverton-holt, 2=ricker)                    ##
0.05  # 3  -std in observed catches in first phase.                          ##
0.025 # 4  -std in observed catches in last phase.                           ##
1     # 5  -Assume unfished in first year (0=FALSE, 1=TRUE)                  ##
0.01  # 6  -Minimum proportion to consider in age-comps for dmvlogistic      ##
0.20  # 7  -Mean fishing mortality for regularizing the estimates of Ft      ##
0.01  # 8  -std in mean fishing mortality in first phase                     ##
5.00  # 9  -std in mean fishing mortality in last phase                      ##
-1	  # 10 -phase for estimating m_deviations (use -1 to turn off mdevs)     ##
0.1	  # 11 -std in deviations for natural mortality                          ##
3	  # 12 -number of estimated nodes for deviations in natural mortality    ##
0.0	  # 13 -fraction of total mortality that takes place prior to spawning   ##
1     # 14 -switch for age-comp likelihood (1=dmvlogistic,2=dmultinom)       ##
0     # 15 -switch for IFD distribution in selectivity simulations			 ##
##                                                                           ##
## ------------------------------------------------------------------------- ##
## MARKER FOR END OF CONTROL FILE (eofc)
## ------------------------------------------------------------------------- ##
999