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
5
## ival         lb      ub      phz     prior   p1      p2      #parameter   ##
   15.9          -5.0    15      4       0       -5.0    15     #log_ro      ##
   0.75         0.2     1.0     4       3       1.01    1.01    #steepness   ##
#  -0.824       -3.0    2.0     -1      1       -1.74   0.1     #log_m_f     ##
#  -0-869       -3.0    2.0     -1      1       -1.74   0.1     #log_m_m     ##
   17.93703     -5.0    15      1       0       -5.0    15      #log_avgrec  ##
#  15.18892     -5.0    15      1       0       -5.0    15      #log_recinit ##
#  15.56507     -5.0    15      1       0       -5.0    15      #log_recinit ##
   0.5          0.01    0.99    -3      3       1.01    1.01    #rho         ##
   0.8          0.01    5.0     -3      4       1.01    1.01    #vartheta    ##
## ------------------------------------------------------------------------- ##
##
## ------------------------------------------------------------------------- ##
## CONTROLS FOR SEX BASED PARAMETERS (nsex arrays, 9 rows, 7 cols)           ##
## ------------------------------------------------------------------------- ##
## FEMALE                                                                    ##
## ival         lb      ub      phz     prior   p1      p2      #parameter   ##
## ------------------------------------------------------------------------- ##
   15.18892     -5.0    15      1       0       -5.0    15      #log_recinit ##
   -1.89712     -3.0    2.0     -1      1       -1.74   0.1     #log_m_f     ##
   145.0        0.0     200     -1      0       0.0     200     #linf        ##
   0.10         0.01    1.0     -1      0       0.01    1.0     #vonk        ##
   -1.2197      -2.0    0.0     -1      0       -2.0    0.0     #to          ##
   6.921e-6     0.0     1.0     -1      0       0.0     1.0     #a           ##
   3.24         2.0     3.5     -1      0       2.0     3.5     #b           ##
   11.49        0.0     30.     -1      0       0.0     30.     #ah          ##
   1.776        0.0     30.     -1      0       0.0     30.     #gh          ##   
## ------------------------------------------------------------------------- ##
## MALE                                                                      ##
## ival         lb      ub      phz     prior   p1      p2      #parameter   ##
## ------------------------------------------------------------------------- ##
   15.34718     -5.0    15      1       0       -5.0    15      #log_recinit ##
   -1.99897     -3.0    2.0     -1      1       -1.74   0.1     #log_m_f     ##
   110.0        0.0     200     -1      0       0.0     200     #linf        ##
   0.12         0.01    1.0     -1      0       0.01    1.0     #vonk        ##
   -1.2197      -2.0    0.0     -1      0       -2.0    0.0     #to          ##
   6.921e-6     0.0     1.0     -1      0       0.0     1.0     #a           ##
   3.24         2.0     3.5     -1      0       2.0     3.5     #b           ##
   11.49        0.0     30.     -1      0       0.0     30.     #ah          ##
   1.776        0.0     30.     -1      0       0.0     30.     #gh          ##   
## ------------------------------------------------------------------------- ##

##
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
##comm  sprt  pers  o32b  032w  u32b  u32w     
  11    11    11    11    11    11    11       # -selectivity type ivector(isel_type) for gear
  97.13 97.13 97.13 97.13 97.13 97.13 97.13    # -Age/length at 50% selectivity (logistic)
  6     6     6     6     6     6     6        # -STD at 50% selectivity (logistic)
  0     0     0     0     0     0     0        # -No. of age nodes for each gear (0=ignore)
  0     0     0     0     0     0     0        # -No. of year nodes for 2d spline(0=ignore)
  -2    -2    -2    -2    -2    -2    -2       # -Phase of estimation (-1 for fixed)
  2     2     2     2     2     2     2        # -Penalty wt for 2nd differences w=1/(2*sig^2)
  3     3     3     3     3     3     3        # -Penalty wt for dome-shaped w=1/(2*sig^2)
  81.3  0     0     0     0     0     0        # -Size limit (cm)
  0.17  0     0     0     0.17  0     0.17     # -Discard mortality rate
## ------------------------------------------------------------------------- ##
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
2					# -number of surveys (nits)
0	2   			# -prior type (see legend above)
0	0   			# -prior log(mean)
0	0.01   			# -prior sd
## ------------------------------------------------------------------------- ##
##

## ------------------------------------------------------------------------- ##
## OTHER MISCELANEOUS CONTROLS                                               ##
## ------------------------------------------------------------------------- ##
0           # 1  -verbose ADMB output (0=off, 1=on)
1           # 2  -recruitment model (1=beverton-holt, 2=ricker)
0.100       # 3  -std in observed catches in first phase.
0.0707      # 4  -std in observed catches in last phase.
0           # 5  -Assume unfished in first year (0=FALSE, 1=TRUE)
0.00        # 6  -Minimum proportion to consider in age-proportions for dmvlogistic
0.20        # 7  -Mean fishing mortality for regularizing the estimates of Ft
0.01        # 8  -std in mean fishing mortality in first phase
2.00        # 9  -std in mean fishing mortality in last phase
-3          # 10 -phase for estimating m_deviations (use -1 to turn off mdevs)
0.1         # 11 -std in deviations for natural mortality
12          # 12 -number of estimated nodes for deviations in natural mortality
0.50        # 13 -fraction of total mortality that takes place prior to spawning
1           # 14 -switch for age-composition likelihood (1=dmvlogistic,2=dmultinom)
#81.28       # 15 -Size limit (cm) for retention (logistic with 10% CV)
#0.17        # 16 -Base discard mortality rate (age-size independent)
##
## ------------------------------------------------------------------------- ##
## MARKER FOR END OF CONTROL FILE (eofc)
## ------------------------------------------------------------------------- ##
999