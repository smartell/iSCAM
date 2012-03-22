  //  This is the 2011 halibut assessment. Selectivity is entirely
  //  determined by mean length. Females and males are counted 
  //  separately. Ages are smeared. There are alternative formulations
  //  of length-specific selectivity. As of 2006, the model allows
  //  different catchability, fishing mortality, selectivity,
  //  and natural mortality for the two sexes. In 2007, a penalty
  //  on the variance of log(survey catchability) is added as a way
  //  of implementing annual catchability deviations.
  //
  // New in 2005 is the option of estimating bycatch selectivities.
  //
  // Naming conventions:
  //
  // Catch, tag C as in e.g. CF (commercial F): commercial catch.
  // Surv, tag S: setline survey
  // CPUE: setline catch per skate in TOTAL number, 
  // WPUE: setline catch per skate of LEGAL fish in weight.
  // Q: catchability (CQ, SQ).
  // SelL: estimated length-specific selectivity (CSelL, SSelL).
  // SelX: a vector of proportional multipliers applied top-down
  //       from the length of full selectivity to generate a SelL
  //       vector (CSelLX, SSelLX).
  // Sel: calculated age-specific selectivity (CSel, SSel).
  // Flag: a 0/1 switch to turn an option on or off.
  // RSS: a residual sum of squares, weighted inversely by the
  //      variance of each data point, e.g. Catch_RSS.
  // PSS: a penalty sum of squares.
  // SD: the standard deviation of a data point, supplied as data.
  //     Also a pseudo-SD used to penalize differences in random walks.
  // SDW: a weight applied to an SD in computing an RSS.
  // Lambda: a weight applied to an RSS in computing the total RSS.
  //         (In production fits all SDW's and Lambda's are 1.)
  // W: average weight at age, e.g CatchW.
  // P: tag for estimation phase, e.g. PC, PSJ.
  // 

GLOBALS_SECTION
  //
  #include <admodel.h>
  #include <ourlib.hpp>
  // Global container for simulation flag.
  int Simulation = 0;
  // Symbolic constants.
  # define NA -99.0 // Missing data or prediction.
  # define CheckValue -12345  // Place marker in .dat and .pin files.
  # define Sex1 1   // Identifies sex as an index in array declarations.
  # define Sex2 2
  # define Sex3 3
  //

TOP_OF_MAIN_SECTION
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(200000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  arrmblsize=2000000;
  //

DATA_SECTION
  //
  // Control parameters at head of .dat file.-------------------
  //
  // Simulation flag.
  //
  init_int	SimFlag
  //
  // Area (2A=21, 2B=22 etc.)
  //
  init_int Area
  //
  // First year in model fit; may be later than in data.
  //
  init_int	FirstYear
  //
  // Last year in model fit; may be earlier than in data.
  //
  init_int	LastYear
  //
  // Number of mean survey lengths (points) at 10 cm intervals
  // at which commercial and setline selectivity are estimated
  // (or calculated, or fixed, according to the selectivity option).
  //
  init_int	SLptN
  //
  // List of said lengths.
  //
  init_ivector	SLptList(1,SLptN)
  //
  // Element of SLptList at which selectivity is fixed at one.
  //
     init_int	SLptOne
  //
  // Fixed length-specific coastwide selectivities, used to compute EBioFixed.
  //
  init_vector	FSelL(1,SLptN)
  // 
  // Minimum and maximum allowed values for fishing mortality.
  //
  init_number	CFParMin
  init_number	CFParMax
  //
  // Number of years when commercial catchability is estimated.
  //
  init_int	CQEstN
  //
  // List of said years.
  //
  init_ivector	CQEstYears(1,CQEstN)
  //
  // Corresponding path types (flat=0 or ramped=1)
  //
  init_ivector CQEstPathTypes(1,CQEstN)
  //
  // Pseudo-SD for computing penalty on annual changes in CQ.
  //
  init_number	SteadyCQ_SD
  //
  // Maximum allowed scaler (C-hook Q)/(J-hook Q), like 2.5 or 3.
  // 
  init_number	HookCorrMax
  //
  // Number of years when commercial selectivity is estimated.
  //
  init_int	CSelLEstN
  //
  // List of said years.
  //
  init_ivector	CSelLEstYears(1,CSelLEstN)
  //
  // Corresponding path types
  //
  init_ivector	CSelLEstPathTypes(1,CSelLEstN)
  //
  // Form of commercial selectivity schedule:
  // 1 = free except zero at SLptList(1) and one at SLptOne
  // 2 = domed with peak at SLptOne
  // 3 = aymptotic with sel = 1 at & after SLptOne
  //
  init_int	CSelLForm
  // 
  // Number of years when survey catchability is estimated.
  //
  init_int	SQEstN
  //
  // List of said years.
  //
  init_ivector	SQEstYears(1,SQEstN)
  //
  // Corresponding path types (flat=0 or ramped=1)
  //
  init_ivector SQEstPathTypes(1,SQEstN)
  //
  // Penalty pseudo-sd on variance of log(SQEst).
  //
  init_int SteadySQFlag
  init_number SteadySQ_SD
  //
  // Flag to require zero slope of regression of SQ on year.
  // Operative only if SteadySQFlag is set.
  //
  init_int	TrendlessSQFlag
  //
  // Number of years when survey selectivity is estimated.
  //
  init_int	SSelLEstN
  //
  // List of said years.
  //
  init_ivector	SSelLEstYears(1,SSelLEstN)
  //
  // Corresponding path types
  //
  init_ivector	SSelLEstPathTypes(1,SSelLEstN)
  //
  // Form of survey selectivity schedule:
  // 1 = free except zero at SLptList(1) and one at SLptOne
  // 2 = domed with peak at SLptOne
  // 3 = aymptotic with sel = 1 at & after SLptOne
  //
  init_int	SSelLForm
  // 
  // Number of lengths (points) at 10 cm intervals at which
  // bycatch selectivity is estimated. These have to be the same 
  // as the intervals used in the bycatch data, which are
  // seq(0, 120, 10). The last length is the plus length
  // (i.e., 120+).
  //
  init_int	BLptN
  //
  // List of said lengths.
  //
  init_ivector	BLptList(1,BLptN)
  //
  // Length at which bycatch selectivity at length is set to 1.0.
  // Must be one of the interior elements of BLptList.
  //
  init_int	BLptOne
  //
  //
  // Number of years when bycatch selectivity is estimated.
  //
  init_int	BSelLEstN
  //
  // List of said years.
  //
  init_ivector	BSelLEstYears(1,BSelLEstN)
  //
  // Corresponding path types
  //
  init_ivector	BSelLEstPathTypes(1,BSelLEstN)
  //
  // Estimate bycatch selectivity by fitting bycatch in number at
  // length (flag!=0) or not (flag=0). In that case the supplied 
  // bycatch selectivity at length is used throughout and the 
  // computed bycatch in weight is fitted to the data value 
  // in BycatchWt.
  //
  init_int	FitBycatchNosFlag
  //
  // Limit on estimated year-class strength at age 1 (millions).
  //
  init_number	RMMax
  //
  // Smear ages (flag!=0) or not.
  //
  init_int	SmearAgesFlag
  //
  // Smooth CSelL and SSelL (flag!=0) or not.
  //
  init_int	SelSmoothFlag
  //
  // Pseudo-SD used to penalize second differences in selectivity
  // schedules if SelSmoothFlag is set.
  //
  init_number	SelSmooth_SD
  //
  // Variance scalers for total and sex-specific catch at age.
  // (This value actually scales standard deviation not variance.)
  //
  init_number	CatchTau_T
  init_number	CatchTau_F
  init_number	CatchTau_M
  //
  // Variance scalers for total and sex-specific commercial CPUE 
  // at age in number.
  //
  init_number	CCPUETau_T
  init_number	CCPUETau_F
  init_number	CCPUETau_M
  //
  // Variance scalers for total commercial CPUE in number & weight.
  //
  init_number	CCPUETotTau
  init_number	CWPUETotTau
  //
  // Variance scalers for total and sex-specific proportion at age
  // in survey catches.
  //
  init_number	SurvPropTau_T
  init_number	SurvPropTau_F
  init_number	SurvPropTau_M
  //
  // Variance scalers for total and sex-specific survey CPUE at age
  // in number.
  //
  init_number	SurvCPUETau_T
  init_number	SurvCPUETau_F
  init_number	SurvCPUETau_M
  //
  // Variance scalers for total survey CPUE in number and legal
  // survey CPUE in weight.
  //
  init_number	SurvCPUETotTau
  init_number	SurvWPUETotTau
  //
  // Variance scaler for bycatch in number at length.
  //
  init_number	BycatchTau
  //
  // Flag for robustifying residuals (in phase 2).
  //
  init_int	RobustifyFlag
  //
  // Weight for Catch at age/sex.
  //
  init_number	CatchLambda
  //
  // Weight for commercial CPUE at age/sex.
  //
  init_number	CCPUELambda
  //
  // Weight for total commercial CPUE in number.
  //
  init_number	CCPUETotLambda
  //
  // Weight for total commercial CPUE in weight.
  //
  init_number	CWPUETotLambda
  //
  // Weight for Survey proportion at age/sex.
  //
  init_number	SurvPropLambda
  //
  // Weight for Survey CCPUE at age/sex.
  //
  init_number	SurvCPUELambda
  //
  // Weight for total survey CPUE in number.
  //
  init_number	SurvCPUETotLambda
  //
  // Weight for total survey CPUE in weight.
  //
  init_number	SurvWPUETotLambda
  //
  // Weight for Bycatch.
  //
  init_number	BycatchLambda
  //
  // Estimate female M (flag!=0) or not (flag=0).
  //
  init_int	EstMFlag
  //
  // Estimate male M (flag!=0) or not (flag=0).
  //
  init_int 	EstMMFlag
  //
  // Setting this flag to zero makes selectivity a function of mean
  // length at age (not at age/sex) and leaves the sex-specific RSS
  // terms out of RSSTotal. (old)
  //
  init_int	SexFlag
  //
  // Setting this flag to zero removes sex-specific Survey and Commercial CPUE from the objective
  // function, as well as the dividing by 8 to account for double fitting
  //
  init_int	MultiFitFlag
  //
  // Last year of commercial and survey data to be included in objective 
  // function calculations. Normally LastYear but sometimes less.
  //
  init_int	CLastYearFit
  init_int	SLastYearFit
  //
  // Flag for fitting a specified value of F in LastYear, and said value.
  //
  init_int	FitLastFFlag
  init_number	LastF
  //
  // Flag for fitting to CatchN and including it in RSS
  init_int	FitCatchNFlag
  //
  // Flag for estimating separate parameters for males.
  // 0 = unisex; 1 = separate catchability and selectivity.
  //
  init_int	SplitMalesFlag
  //
  // Check value at end of settings and flags.
  //
  init_int 	DataCheck1
  //
  !!  if (DataCheck1 != CheckValue)
  !!     {
  !!     cerr << "Flags and settings are not in order." << endl;
  !!     exit(1);
  !!     } else
  !!     cout << "Flags and settings are in order." << endl;
  //
  // Commercial and survey data to be fitted.--------------------------------
  //
  // First and last years in the data files.
  //
  init_int	YearFirstDat
  init_int	YearLastDat
  //
  // Min and max observed ages in data.
  //
  init_int	AgeMinDat
  init_int	AgeMaxDat
  //
  // Observed plus ages.
  //
  init_int	AgePlusSurf
  init_int	AgePlusBurn
  //
  // Min and max true ages in "data" values of size
  // at true age.
  //
  init_int	AgeTrueMinDat
  init_int	AgeTrueMaxDat
  //
  // Number and list of 10 cm intervals in bycatch data.
  // Each length is the low end of a 10 cm interval.
  // The last is a plus length. Retained in this version
  // only because they appear in the consolidated
  // data file single.dat.
  //
  init_int	LintBN
  init_ivector	LintBList(1,LintBN)
  //
  //
  // Catch at age.
  //
  init_matrix	Catch_F_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	Catch_M_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	Catch_T_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Standard error thereof.
  //
  init_matrix	Catch_SE_F_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	Catch_SE_M_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	Catch_SE_T_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Commercial CPUE at age/sex in number.
  //
  init_matrix	CCPUE_F_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	CCPUE_M_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	CCPUE_T_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Standard error thereof.
  //
  init_matrix	CCPUE_SE_F_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	CCPUE_SE_M_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	CCPUE_SE_T_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Mean weight at observed age in the commercial catch
  //
  init_matrix	CatchW_F(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	CatchW_M(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	CatchW_T(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Mean weight at true age in the commercial catch
  //
  init_matrix	CatchWTrue_F(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	CatchWTrue_M(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	CatchWTrue_T(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  //
  // Total commercial CPUE in number and standard error.
  //
  init_vector	CCPUETot_Dat(YearFirstDat,YearLastDat)
  init_vector	CCPUETot_SE_Dat(YearFirstDat,YearLastDat)
  //
  // Total commercial CPUE in weight and standard error.
  //
  init_vector	CWPUETot_Dat(YearFirstDat,YearLastDat)
  init_vector	CWPUETot_SE_Dat(YearFirstDat,YearLastDat)
  //
  // Survey proportion at age/sex in number.
  //
  init_matrix	SurvProp_F_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvProp_M_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvProp_T_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Standard error thereof.
  //
  init_matrix	SurvProp_SE_F_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvProp_SE_M_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvProp_SE_T_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Survey CPUE at age/sex in number.
  //
  init_matrix	SurvCPUE_F_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvCPUE_M_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvCPUE_T_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Standard error thereof.
  //
  init_matrix	SurvCPUE_SE_F_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvCPUE_SE_M_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvCPUE_SE_T_Dat(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Total survey CPUE in number and standard error.
  //
  init_vector	SurvCPUETot_Dat(YearFirstDat,YearLastDat)
  init_vector	SurvCPUETot_SE_Dat(YearFirstDat,YearLastDat)
  //
  // Total survey CPUE in weight and standard error.
  //
  init_vector	SurvWPUETot_Dat(YearFirstDat,YearLastDat)
  init_vector	SurvWPUETot_SE_Dat(YearFirstDat,YearLastDat)
  //
  // Mean survey length at observed age.
  //
  init_matrix	SurvL_F(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvL_M(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvL_T(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Mean survey length at true age.
  //
  init_matrix	SurvLTrue_F(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	SurvLTrue_M(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	SurvLTrue_T(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  //
  // Mean survey weight at observed age.
  //
  init_matrix	SurvW_F(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvW_M(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvW_T(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Mean survey weight at true age.
  //
  init_matrix	SurvWTrue_F(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	SurvWTrue_M(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	SurvWTrue_T(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  //
  // Proportion of legal size at observed age.
  //
  init_matrix	SurvPropLegal_F(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvPropLegal_M(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvPropLegal_T(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  // Mean weight of legal sized fish at observed age.
  //
  init_matrix	SurvLegalW_F(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvLegalW_M(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  init_matrix	SurvLegalW_T(YearFirstDat,YearLastDat,AgeMinDat,AgeMaxDat)
  //
  //
  // Total bycatch in number by 10 cm length interval and standard error.
  //
  init_matrix	Bycatch_Dat(YearFirstDat,YearLastDat,1,BLptN)
  init_matrix	Bycatch_SE_Dat(YearFirstDat,YearLastDat,1,BLptN)
  //
  // Mean bycatch length at true age.
  //
  init_matrix	BycatchLTrue_F(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	BycatchLTrue_M(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	BycatchLTrue_T(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  //
  // Standard deviation of bycatch length at true age.
  //
  init_matrix	BycatchLSDTrue_F(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	BycatchLSDTrue_M(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	BycatchLSDTrue_T(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  //
  // 
  //
  init_matrix	BycatchWTrue_F(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	BycatchWTrue_M(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  init_matrix	BycatchWTrue_T(YearFirstDat,YearLastDat,AgeTrueMinDat,AgeTrueMaxDat)
  //
  // Total commercial, bycatch, sport catch, and personal use catch in weight (lb).
  //
  init_vector	CatchWt_Dat(YearFirstDat,YearLastDat)
  init_vector	DiscardWt_Dat(YearFirstDat,YearLastDat)
  init_vector	BycatchWt_Dat(YearFirstDat,YearLastDat)
  init_vector	SportCatchWt_Dat(YearFirstDat,YearLastDat)
  init_vector	PersUseWt_Dat(YearFirstDat,YearLastDat)
  //
  // End commercial and survey data.-----------------------------------------------
  //
  // Miscellaneous working arrays.
  //
  // Distribution of surface and burn ages at true age.
  //
  init_matrix	TrueToSurf(1,AgePlusSurf,1,AgeTrueMaxDat)
  init_matrix	TrueToBurn(1,AgePlusBurn,1,AgeTrueMaxDat)
  //
  // Proportion mature at age (true or observed).
  //
  init_vector	PropMat(AgeMinDat,AgeTrueMaxDat)
  //
  //
  // At this point the data file should have a line with
  // the single check value -12345.
  //
  init_int 	DataCheck2
  //
  !! if (DataCheck2 != -12345)
  !!    {cerr << "Data are not in order." << endl; exit(1);} else
  !!     cout << "Data are in order." << endl;
  //
  // Data files loaded.------------------------------------------
  //
  // Assign global Simulation.
  //
  !! Simulation = SimFlag;
  //
  // Working variables for selecting and organizing the observations.
  //
  // Internal model arrays always have ages 1 through LastAge = AgeTrueMaxDat.
  // This code also assumes that AgeTrueMinDat is 1.
  //
  int		LastAge
  !!		LastAge = AgeTrueMaxDat;
  //
  int		FirstAgeDat
  !!		FirstAgeDat = AgeMinDat;
  int		LastAgeDat
  !!		LastAgeDat = AgeMaxDat;
  //
  // SMinAge and FirstAgeDat are synonymous.
  // This code assumes SMinAge is the first age used in survey fits.
  // Likewise SMaxAge and LastAgeDat.
  int		SMinAge
  !!		SMinAge = FirstAgeDat;
  int		SMaxAge
  !!		SMaxAge = LastAgeDat;
  //
  // Only ages CMinAge and higher are used in commercial fits.
  int		CMinAge
  !!		CMinAge = 8;
  //
  int		AgePlus      // Holder used in calcs.
  int		AgePlusFirst // AgePlus in FirstYear.
  !!  if (FirstYear <= 2001) AgePlusFirst = AgePlusSurf; else
  !!                         AgePlusFirst = AgePlusBurn;
  //
  // Vectors to hold true and observed age distns.
  vector	TrueAgeFreqs(1,LastAge)
  vector	ObsAgeFreqs(1,LastAgeDat)
  //
  int		y	// Always indexes year.
  int		a	// Always indexes observed age.
  int		b	// Always indexes true age.
  int		s	// Always indexes sex (1=female, 2=male, 3=total).
  int		v	// Always indexes length point/interval.
  int		w	// Always indexes waypoint.
  int		istep	// Indexes steps in waypoint calculations.
  int		p	// Indexes selectivity param 1,...,SSelLParN.
  //
  // Lengths etc. used in computing predicted bycatch by length interval.
  number	LBot
  number	zBot
  number	LTop
  number	zTop
  number	LPorp
  //
  // Internal robustification flag, left off in phase 1 even if
  // RobustifyFlag is set.
  //
  int		Robust
  //
  // Load sex-specific matrices into 3darrays so that the same code
  // will be run on both sexes in a loop. This step also limits
  // the observations that are predicted to the years in the fit. 
  // Other items that are not predicted (like size at
  // age) are not limited to years in the fit so that e.g.
  // SurvW(LastYear+1,a,s) can be used for model outputs.
  //
  3darray	Catch(FirstYear,LastYear,SMinAge,SMaxAge,Sex1,Sex3)
  matrix	Catch_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	Catch_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	Catch_T(FirstYear,LastYear,SMinAge,SMaxAge)
  3darray	Catch_SE(FirstYear,LastYear,SMinAge,SMaxAge,Sex1,Sex3)
  matrix	Catch_SE_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	Catch_SE_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	Catch_SE_T(FirstYear,LastYear,SMinAge,SMaxAge)
  3darray	CCPUE(FirstYear,LastYear,SMinAge,SMaxAge,Sex1,Sex3)
  matrix	CCPUE_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	CCPUE_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	CCPUE_T(FirstYear,LastYear,SMinAge,SMaxAge)
  3darray	CCPUE_SE(FirstYear,LastYear,SMinAge,SMaxAge,Sex1,Sex3)
  matrix	CCPUE_SE_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	CCPUE_SE_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	CCPUE_SE_T(FirstYear,LastYear,SMinAge,SMaxAge)
  3darray	SurvProp(FirstYear,LastYear,SMinAge,SMaxAge,Sex1,Sex3)
  matrix	SurvProp_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvProp_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvProp_T(FirstYear,LastYear,SMinAge,SMaxAge)
  3darray	SurvProp_SE(FirstYear,LastYear,SMinAge,SMaxAge,Sex1,Sex3)
  matrix	SurvProp_SE_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvProp_SE_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvProp_SE_T(FirstYear,LastYear,SMinAge,SMaxAge)
  3darray	SurvCPUE(FirstYear,LastYear,SMinAge,SMaxAge,Sex1,Sex3)
  matrix	SurvCPUE_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvCPUE_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvCPUE_T(FirstYear,LastYear,SMinAge,SMaxAge)
  3darray	SurvCPUE_SE(FirstYear,LastYear,SMinAge,SMaxAge,Sex1,Sex3)
  matrix	SurvCPUE_SE_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvCPUE_SE_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvCPUE_SE_T(FirstYear,LastYear,SMinAge,SMaxAge)
  //
  3darray	CatchWTrue(YearFirstDat,YearLastDat,1,LastAge,Sex1,Sex3)
  3darray	SurvL(YearFirstDat,YearLastDat,SMinAge,SMaxAge,Sex1,Sex3)
  3darray	SurvLTrue(YearFirstDat,YearLastDat,1,LastAge,Sex1,Sex3)
  3darray	SurvW(YearFirstDat,YearLastDat,SMinAge,SMaxAge,Sex1,Sex3)
  3darray	SurvWTrue(YearFirstDat,YearLastDat,1,LastAge,Sex1,Sex3)
  3darray	SurvPropLegal(YearFirstDat,YearLastDat,SMinAge,SMaxAge,Sex1,Sex3)
  3darray	SurvLegalW(YearFirstDat,YearLastDat,SMinAge,SMaxAge,Sex1,Sex3)
  3darray	BycatchLTrue(YearFirstDat,YearLastDat,1,LastAge,Sex1,Sex3)
  3darray	BycatchLSDTrue(YearFirstDat,YearLastDat,1,LastAge,Sex1,Sex3)
  3darray	BycatchWTrue(YearFirstDat,YearLastDat,1,LastAge,Sex1,Sex3)
  //
  // Arrays of observations that are to be fitted.
  //
  !! for (y=FirstYear; y<=LastYear; y++)
  !! for (a=SMinAge; a<=SMaxAge; a++)
  !!   {
  !!   Catch(y,a,1) = Catch_F_Dat(y,a);
  !!   Catch(y,a,2) = Catch_M_Dat(y,a);
  !!   Catch(y,a,3) = Catch_T_Dat(y,a);
  !!   Catch_F(y,a) = Catch_F_Dat(y,a);
  !!   Catch_M(y,a) = Catch_M_Dat(y,a);
  !!   Catch_T(y,a) = Catch_T_Dat(y,a);
  //
  !!   Catch_SE(y,a,1) = Catch_SE_F_Dat(y,a);
  !!   Catch_SE(y,a,2) = Catch_SE_M_Dat(y,a);
  !!   Catch_SE(y,a,3) = Catch_SE_T_Dat(y,a);
  !!   Catch_SE_F(y,a) = Catch_SE_F_Dat(y,a);
  !!   Catch_SE_M(y,a) = Catch_SE_M_Dat(y,a);
  !!   Catch_SE_T(y,a) = Catch_SE_T_Dat(y,a);
  //
  !!   CCPUE(y,a,1) = CCPUE_F_Dat(y,a);
  !!   CCPUE(y,a,2) = CCPUE_M_Dat(y,a);
  !!   CCPUE(y,a,3) = CCPUE_T_Dat(y,a);
  !!   CCPUE_F(y,a) = CCPUE_F_Dat(y,a);
  !!   CCPUE_M(y,a) = CCPUE_M_Dat(y,a);
  !!   CCPUE_T(y,a) = CCPUE_T_Dat(y,a);
  //
  !!   CCPUE_SE(y,a,1) = CCPUE_SE_F_Dat(y,a);
  !!   CCPUE_SE(y,a,2) = CCPUE_SE_M_Dat(y,a);
  !!   CCPUE_SE(y,a,3) = CCPUE_SE_T_Dat(y,a);
  !!   CCPUE_SE_F(y,a) = CCPUE_SE_F_Dat(y,a);
  !!   CCPUE_SE_M(y,a) = CCPUE_SE_M_Dat(y,a);
  !!   CCPUE_SE_T(y,a) = CCPUE_SE_T_Dat(y,a);
  //
  !!   SurvProp(y,a,1) = SurvProp_F_Dat(y,a);
  !!   SurvProp(y,a,2) = SurvProp_M_Dat(y,a);
  !!   SurvProp(y,a,3) = SurvProp_T_Dat(y,a);
  !!   SurvProp_F(y,a) = SurvProp_F_Dat(y,a);
  !!   SurvProp_M(y,a) = SurvProp_M_Dat(y,a);
  !!   SurvProp_T(y,a) = SurvProp_T_Dat(y,a);
  //
  !!   SurvProp_SE(y,a,1) = SurvProp_SE_F_Dat(y,a);
  !!   SurvProp_SE(y,a,2) = SurvProp_SE_M_Dat(y,a);
  !!   SurvProp_SE(y,a,3) = SurvProp_SE_T_Dat(y,a);
  !!   SurvProp_SE_F(y,a) = SurvProp_SE_F_Dat(y,a);
  !!   SurvProp_SE_M(y,a) = SurvProp_SE_M_Dat(y,a);
  !!   SurvProp_SE_T(y,a) = SurvProp_SE_T_Dat(y,a);
  //
  !!   SurvCPUE(y,a,1) = SurvCPUE_F_Dat(y,a);
  !!   SurvCPUE(y,a,2) = SurvCPUE_M_Dat(y,a);
  !!   SurvCPUE(y,a,3) = SurvCPUE_T_Dat(y,a);
  !!   SurvCPUE_F(y,a) = SurvCPUE_F_Dat(y,a);
  !!   SurvCPUE_M(y,a) = SurvCPUE_M_Dat(y,a);
  !!   SurvCPUE_T(y,a) = SurvCPUE_T_Dat(y,a);
  //
  !!   SurvCPUE_SE(y,a,1) = SurvCPUE_SE_F_Dat(y,a);
  !!   SurvCPUE_SE(y,a,2) = SurvCPUE_SE_M_Dat(y,a);
  !!   SurvCPUE_SE(y,a,3) = SurvCPUE_SE_T_Dat(y,a);
  !!   SurvCPUE_SE_F(y,a) = SurvCPUE_SE_F_Dat(y,a);
  !!   SurvCPUE_SE_M(y,a) = SurvCPUE_SE_M_Dat(y,a);
  !!   SurvCPUE_SE_T(y,a) = SurvCPUE_SE_T_Dat(y,a);
  !!   }
  //
  // Arrays of observations that are not to be fitted.
  !! for (y=YearFirstDat; y<=YearLastDat; y++)
  !! for (a=SMinAge; a<=SMaxAge; a++)
  !!   {  
  !!   SurvL(y,a,1) = SurvL_F(y,a);
  !!   SurvL(y,a,2) = SurvL_M(y,a);
  !!   SurvL(y,a,3) = SurvL_T(y,a);
  !!   SurvW(y,a,1) = SurvW_F(y,a);
  !!   SurvW(y,a,2) = SurvW_M(y,a);
  !!   SurvW(y,a,3) = SurvW_T(y,a);
  !!   SurvPropLegal(y,a,1) = SurvPropLegal_F(y,a);
  !!   SurvPropLegal(y,a,2) = SurvPropLegal_M(y,a);
  !!   SurvPropLegal(y,a,3) = SurvPropLegal_T(y,a);
  !!   SurvLegalW(y,a,1) = SurvLegalW_F(y,a);
  !!   SurvLegalW(y,a,2) = SurvLegalW_M(y,a);
  !!   SurvLegalW(y,a,3) = SurvLegalW_T(y,a);
  !!   }
  //
  // Arrays of true values.
  //
  !!   for (y=YearFirstDat; y<=YearLastDat; y++)
  !!   for (b=1; b<=LastAge; b++)
  !!      {
  !!      CatchWTrue(y,b,1) = CatchWTrue_F(y,b);
  !!      CatchWTrue(y,b,2) = CatchWTrue_M(y,b);
  !!      CatchWTrue(y,b,3) = CatchWTrue_T(y,b);
  !!      SurvLTrue(y,b,1) = SurvLTrue_F(y,b);
  !!      SurvLTrue(y,b,2) = SurvLTrue_M(y,b);
  !!      SurvLTrue(y,b,3) = SurvLTrue_T(y,b);
  !!      SurvWTrue(y,b,1) = SurvWTrue_F(y,b);
  !!      SurvWTrue(y,b,2) = SurvWTrue_M(y,b);
  !!      SurvWTrue(y,b,3) = SurvWTrue_T(y,b);
  !!      BycatchLTrue(y,b,1) = BycatchLTrue_F(y,b);
  !!      BycatchLTrue(y,b,2) = BycatchLTrue_M(y,b);
  !!      BycatchLTrue(y,b,3) = BycatchLTrue_T(y,b);
  !!      BycatchLSDTrue(y,b,1) = BycatchLSDTrue_F(y,b);
  !!      BycatchLSDTrue(y,b,2) = BycatchLSDTrue_M(y,b);
  !!      BycatchLSDTrue(y,b,3) = BycatchLSDTrue_T(y,b);
  !!      BycatchWTrue(y,b,1) = BycatchWTrue_F(y,b);
  !!      BycatchWTrue(y,b,2) = BycatchWTrue_M(y,b);
  !!      BycatchWTrue(y,b,3) = BycatchWTrue_T(y,b);
  !!      }
  //
  // Observed bycatch at length (to be fitted maybe).
  //
  matrix	Bycatch(FirstYear,LastYear,1,BLptN)
  matrix	Bycatch_SE(FirstYear,LastYear,1,BLptN)
  //
  !!	for (y=FirstYear; y<=LastYear; y++)
  !!	for (v=1; v<=BLptN; v++)
  !!	   {
  !!	   Bycatch(y,v) = Bycatch_Dat(y,v);
  !!	   Bycatch_SE(y,v) = Bycatch_SE_Dat(y,v);
  !!	   }
  //
  // Reload observation vectors to match Year range.
  //
  vector	CCPUETot(FirstYear,LastYear)
  vector        CCPUETot_SE(FirstYear,LastYear)
  vector	CWPUETot(FirstYear,LastYear)
  vector        CWPUETot_SE(FirstYear,LastYear)
  vector	SurvCPUETot(FirstYear,LastYear)
  vector        SurvCPUETot_SE(FirstYear,LastYear)
  vector	SurvWPUETot(FirstYear,LastYear)
  vector        SurvWPUETot_SE(FirstYear,LastYear)
  vector	CatchWt(FirstYear,LastYear)
  vector	DiscardWt(FirstYear,LastYear)
  vector	BycatchWt(FirstYear,LastYear)
  vector	SportCatchWt(FirstYear,LastYear)
  vector	PersUseWt(FirstYear,LastYear)
  //
  !!   for (y=FirstYear; y<=LastYear; y++)
  !!      {
  !!      CCPUETot(y) = CCPUETot_Dat(y);
  !!      CCPUETot_SE(y) = CCPUETot_SE_Dat(y);
  !!      CWPUETot(y) = CWPUETot_Dat(y);
  !!      CWPUETot_SE(y) = CWPUETot_SE_Dat(y);
  !!      SurvCPUETot(y) = SurvCPUETot_Dat(y);
  !!      SurvCPUETot_SE(y) = SurvCPUETot_SE_Dat(y);
  !!      SurvWPUETot(y) = SurvWPUETot_Dat(y);
  !!      SurvWPUETot_SE(y) = SurvWPUETot_SE_Dat(y);
  !!      for (v=1; v<=BLptN; v++) Bycatch(y,v) = Bycatch_Dat(y,v);
  !!      CatchWt(y) = CatchWt_Dat(y);
  !!      DiscardWt(y) = DiscardWt_Dat(y);
  !!      BycatchWt(y) = BycatchWt_Dat(y);
  !!      SportCatchWt(y) = SportCatchWt_Dat(y);
  !!      PersUseWt(y) = PersUseWt_Dat(y);
  !!      }
  //
  // Compute total catch in number (CMinAge+) by year for computing CatchN_RSS.
  //
  vector	CatchN(FirstYear,LastYear)
  !!  for (y=FirstYear; y<=LastYear; y++)
  !!     {
  !!     if (y<=2001) AgePlus = AgePlusSurf; else AgePlus = AgePlusBurn;
  !!     CatchN(y) = 0.0;
  !!     for (a=CMinAge; a<=AgePlus; a++) CatchN(y) += Catch_T(y,a);
  !!     }
  //
  // Last year-class estimable from the data; depends on when last 
  // survey was done.
  //
  int		LastYearClass
  int		LastSurvYear
  !!  LastSurvYear = 0;
  !!  for (y=FirstYear; y<=SLastYearFit; y++)
  !!     if (SurvCPUETot(y) != NA) LastSurvYear = y;
  !!  LastYearClass = LastSurvYear - SMinAge;
  !!  if (CLastYearFit - CMinAge > LastSurvYear - SMinAge)
  !!      LastYearClass = CLastYearFit - CMinAge;
  //
  // Variables used for waypoint/path calculations.
  //
  int		StartYear
  int		StopYear
  //
  // Determine number of survey selectivity parameters to be estimated.
  // Set vOne to index of SLptOne in SLptList (same for SSelL and CSelL).
  //
  int 		SSelLParN
  int		vOne
  !!  vOne = 0;
  !!  for (v=1; v<=SLptN; v++)
  !!     if (SLptList(v)==SLptOne) vOne = v;
  !!  if (SSelLForm==1) SSelLParN = SLptN - 2;
  !!  if (SSelLForm==2) SSelLParN = SLptN - 1;
  !!  if (SSelLForm==3) SSelLParN = vOne - 1;
  //
  // Determine number of commercial selectivity parameters to be estimated.
  //
  int 		CSelLParN
  !!  if (CSelLForm==1) CSelLParN = SLptN - 2;
  !!  if (CSelLForm==2) CSelLParN = SLptN - 1;
  !!  if (CSelLForm==3) CSelLParN = vOne - 1;
  //
  // Bycatch is simpler.
  //
  int 		BSelLParN
  int		vOneB
  !!  BSelLParN = BLptN - 1;
  !!  vOneB = 0;
  !!  for (v=1; v<=BLptN; v++)
  !!     if (BLptList(v)==BLptOne) vOneB = v;
  //
  // If (SmearAgesFlag == 0), substitute simple transfer-and-accumulate
  // matrices for TrueToSurf and TrueToBurn.
  //
  !!  if (SmearAgesFlag==0)
  !!     {
  !!     TrueToSurf = TrueToSurf - TrueToSurf;
  !!     for (a=1; a<=AgePlusSurf; a++) TrueToSurf(a,a) = 1.0;
  !!     if (LastAge > AgePlusSurf)
  !!        for (b=AgePlusSurf+1; b<=LastAge; b++) 
  !!           TrueToSurf(AgePlusSurf,b) = 1.0;
         //
  !!     TrueToBurn = TrueToBurn - TrueToBurn;
  !!     for (a=1; a<=AgePlusBurn; a++) TrueToBurn(a,a) = 1.0;
  !!     if (LastAge > AgePlusBurn)
  !!        for (b=AgePlusBurn+1; b<=LastAge; b++) 
  !!           TrueToBurn(AgePlusBurn,b) = 1.0;
  !!     }
  //
  //
  // GOOD_COP_SECTION ----------------------------------------------------
  //
  // Various sanity and consistency checks.
  //
  // Make sure SLptOne is an element of SLPtList and that
  // there are more than zero selectivity params to estimate.
  // Similar for BLptOne.
  //
  !!  if (vOne==0)
  !!     {cout << "SLptOne is not one of SLptList" << endl; exit(1);}
  !!  if (SSelLParN<=0)
  !!     {cout << "SSelLParN <= 0" << endl; exit(1);}
  !!  if (vOneB==0)
  !!     {cout << "BLptOne is not one of BLptList" << endl; exit(1);}
  //
  // Make sure catchabilities and selectivities are estimated
  // in FirstYear. Technically it should also be checked that
  // PathType is zero for the last waypoint, but it is assumed 
  // to be zero in the calcs so it is not checked here.
  //
  !!  if (CQEstYears(1) != FirstYear)
  !!     {cerr << "CQ must be estimated in FirstYear" << endl; exit(1);}
  !!  if (CSelLEstYears(1) != FirstYear)
  !!     {cerr << "CSelL must be estimated in FirstYear" << endl; exit(1);}
  !!  if (SQEstYears(1) != FirstYear)
  !!     {cerr << "SQ must be estimated in FirstYear" << endl; exit(1);}
  !!  if (SSelLEstYears(1) != FirstYear)
  !!     {cerr << "SSelL must be estimated in FirstYear" << endl; exit(1);}
  !!  if (BSelLEstYears(1) != FirstYear)
  !!     {cerr << "BSelL must be estimated in FirstYear" << endl; exit(1);}
  //
  // Make sure CQ is not estimated in 1983.
  //
  !!  for (w=1; w<=CQEstN; w++)
  !!     if (CQEstYears(w)==1983)
  !!        {cerr << "CQ cannot be estimated in 1983" << endl; exit(1);}
  //
  // Make sure CLastYearFit and SLastYearFit are in range.
  //
  !!  if (CLastYearFit < FirstYear || CLastYearFit > LastYear)
  !!    {cerr << "CLastYearFit out of range" << endl; exit(1);}
  !!  if (SLastYearFit < FirstYear || SLastYearFit > LastYear)
  !!    {cerr << "SLastYearFit out of range" << endl; exit(1);}
  //
  // Make sure LastYearClass is in range.
  //
  !!  if (LastYearClass < FirstYear)
  !!    {cerr << "LastYearClass < FirstYear" << endl; exit(1);}
  //

  // RUBE_GOLDBERG_SECTION
  //
  // Parameters like BSelLX that may or may not be estimated are assigned
  // a positive or negative estimation phase here to accomplish that.
  //
  int 	P_BSelLX
  int 	FitBycatchWtFlag
  !!  if (FitBycatchNosFlag == 0) 
  !!        {
  !!        P_BSelLX = -1; 
  !!        FitBycatchWtFlag = 1;
  !!        } else
  !!        {
  !!        P_BSelLX = 1;
  !!        FitBycatchWtFlag = 0;
  !!        }
  //
  // Estimate female M?
  int 	PM
  !! if (EstMFlag != 0) PM = 1; else PM = -1;
  //
  // Estimate male M?
  int PMM
  !! if (EstMMFlag != 0) PMM = 1; else PMM = -1;
  //
  // Estimate separate catchability/selectivity
  // parameters for males?
  //
  // Phase for male catchability (& fishing mortality) scalers.
  int P_QMX
  !! if (SplitMalesFlag == 1) P_QMX = 1; else P_QMX = -1;
  //
  // Phase for separate male selectivity parameters.
  int P_SelL_M
  !! if (SplitMalesFlag == 1) P_SelL_M = 1; else P_SelL_M = -1;
  //
  // Robustify residuals in a phase 2?
  //
  int	P_ProForma
  !! if (RobustifyFlag != 0) P_ProForma = 2; else P_ProForma = 1;
  //

PARAMETER_SECTION
  //
  // Objective function.------------------------------------
  //
  // Total SS including penalties.
  // 
  objective_function_value MinusLogL
  number Deviance
  //
  // Total sum of squares.
  //
  number	SSTotal
  //
  // Just the deviations.
  //
  number	RSSTotal
  number	RSSNobs
  //
  // Just the penalties.
  //
  number	PSSTotal
  //
  // Model variables including model parameters.-------------
  //
  // Natural mortality. Normally not estimated.
  //
  matrix	M(FirstYear,LastYear,Sex1,Sex2)
  init_bounded_number	MPar(0.05,0.25,PM)
  init_bounded_number	MPar_M(0.05,0.25,PMM)
  //
  // Full-recuitment commercial fishing mortality CF(y,s).
  //
  matrix	CF(FirstYear,LastYear,Sex1,Sex2)
  init_bounded_vector	CFPar(FirstYear,LastYear,CFParMin,CFParMax)
  //
  // Commercial catchablity CQ. The parameter is CQEst x 1E6.
  //
  matrix	CQ(FirstYear,LastYear,Sex1,Sex2)
  init_bounded_vector	CQEstPar(1,CQEstN,1.0E-2,1.0E2)
  vector	CQEst(1,CQEstN)
  init_bounded_number	CQMX(0.1,10.0,P_QMX)
  !! if (SplitMalesFlag == 0) CQMX = 1.0; // GOOD_COP.
  number	StartVal  // Used in path & waypoint calcs.
  number 	StopVal
  number	StepVal
  number	PathLength
  //
  // Commercial selectivity at length CSelL.
  //
  3darray	CSelL(FirstYear,LastYear,1,SLptN,Sex1,Sex2) // Working array.
  matrix	CSelL_F(FirstYear,LastYear,1,SLptN) // For REPORT only.
  matrix	CSelL_M(FirstYear,LastYear,1,SLptN) // For REPORT only.
  matrix	CSelLy(Sex1,Sex2,1,SLptN)  // For computing CSel w approx(). 
  init_bounded_matrix	CSelLPar(1,CSelLEstN,1,CSelLParN,0.0,1.0)
  init_bounded_matrix	CSelLPar_M(1,CSelLEstN,1,CSelLParN,0.0,1.0,P_SelL_M)
  3darray	CSelLParSplit(1,CSelLEstN,1,CSelLParN,Sex1,Sex2)
  3darray	CSelLEst(1,CSelLEstN,1,SLptN,Sex1,Sex2)
  //
  // Estimated commercial selectivity at age CSel.
  //
  3darray	CSel(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  matrix	CSel_F(FirstYear,LastYear,1,LastAge) // For REPORT only.
  matrix	CSel_M(FirstYear,LastYear,1,LastAge) // For REPORT only.
  //
  // Commercial selectivity at age calculated from fixed coastwide
  // length-specific selectivities.
  //
  3darray	FSel(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  matrix	FSel_F(FirstYear,LastYear,1,LastAge)
  matrix	FSel_M(FirstYear,LastYear,1,LastAge)
  //
  // Fully selected sublegal discard mortality.
  //
  vector	DF(FirstYear,LastYear)
  init_bounded_vector	DFPar(FirstYear,LastYear,0.0,0.1)
  //
  // Sublegal discard selectivity at length; same for all years, both sexes.
  //
  vector	DSelL(1,SLptN)
  !! for (v=1; v<=SLptN; v++)
  !!    {if  (SLptList(v)==80) DSelL(v) = 1.0; else DSelL(v) = 0.0;}
  //
  // Sublegal discard selectivity at age.
  //
  3darray	DSel(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  matrix	DSel_F(FirstYear,LastYear,1,LastAge)
  matrix	DSel_M(FirstYear,LastYear,1,LastAge)
  //
  // Survey catchablity SQ. The param is SQEst x 1E6.
  //
  matrix	SQ(FirstYear,LastYear,Sex1,Sex2)
  init_bounded_vector	SQEstPar(1,SQEstN,1.0E-2,1.0E2)
  vector	SQEst(1,SQEstN)
  init_bounded_number	SQMX(0.1,10.0,P_QMX)
  !! if (SplitMalesFlag == 0) SQMX = 1.0; // GOOD_COP.
  number	SQ_bhat  // Slope of log(SQ) over time.
  //
  // Survey selectivity at length SSelL.
  //
  3darray	SSelL(FirstYear,LastYear,1,SLptN,Sex1,Sex2)
  matrix	SSelL_F(FirstYear,LastYear,1,SLptN) // For REPORTing.
  matrix	SSelL_M(FirstYear,LastYear,1,SLptN) // For REPORTing.
  matrix	SSelLy(Sex1,Sex2,1,SLptN) // For computing SSel w approx()
  init_bounded_matrix	SSelLPar(1,SSelLEstN,1,SSelLParN,0.0,1.0)
  init_bounded_matrix	SSelLPar_M(1,SSelLEstN,1,SSelLParN,0.0,1.0,P_SelL_M)
  3darray	SSelLParSplit(1,SSelLEstN,1,SSelLParN,Sex1,Sex2)
  3darray	SSelLEst(1,SSelLEstN,1,SLptN,Sex1,Sex2)
  vector	StartVec(1,SLptN)  // Used in path & waypoint calcs.
  vector	StopVec(1,SLptN)
  vector	StepVec(1,SLptN)
  vector	SelVec(1,SLptN) // Used in SelSmooth calcs.
  //
  // Survey selectivity at age SSel.
  //
  3darray	SSel(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  matrix	SSel_F(FirstYear,LastYear,1,LastAge) // For REPORTing.
  matrix	SSel_M(FirstYear,LastYear,1,LastAge) // For REPORTing.
  //
  // Full-recruitment bycatch fishing mortality BF(y).
  //
  init_bounded_vector	BF(FirstYear,LastYear,0.0,1.0)
  //
  // Bycatch selectivity at length BSelL.
  //
  3darray	BSelL(FirstYear,LastYear,1,BLptN,Sex1,Sex2)
  matrix	BSelL_F(FirstYear,LastYear,1,BLptN) // For REPORTing.
  matrix	BSelL_M(FirstYear,LastYear,1,BLptN) // For REPORTing.
  matrix	BSelLy(Sex1,Sex2,1,BLptN) // For computing BSel w approx()
  init_bounded_matrix	BSelLPar(1,BSelLEstN,1,BSelLParN,0.001,1000.0,P_BSelLX)
  matrix	BSelLEst(1,BSelLEstN,1,BLptN)
  //
  // Survey selectivity at age BSel.
  //
  3darray	BSel(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  matrix	BSel_F(FirstYear,LastYear,1,LastAge) // For REPORTing.
  matrix	BSel_M(FirstYear,LastYear,1,LastAge) // For REPORTing.
  //
  vector	StartVecB(1,BLptN)  // Used in path & waypoint calcs
  vector	StopVecB(1,BLptN)   // for bycatch.
  vector	StepVecB(1,BLptN)
  //
  // Recreational (sport) fishing mortality and selectivity.
  // Length-specific RSelL=SSelL, so age-specific RSel=SSel
  //
  init_bounded_vector	RFPar(FirstYear,LastYear,0.0,1.0)
  matrix	RF(FirstYear,LastYear,Sex1,Sex2)
  3darray	RSel(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  //
  // Personal use (subsistence) fishing mortality.
  // As above, PSel=SSel.
  //
  init_bounded_vector	PFPar(FirstYear,LastYear,0.0,1.0)
  matrix	PF(FirstYear,LastYear,Sex1,Sex2)
  3darray	PSel(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  //
  // Other mortality stuff.
  vector	OldZ(Sex1,Sex2) // Used in setting InitN.
  // Following arrays are year/age/sex-specific mortalities;
  // e.g. Cf(y,b,s) = CF(y,s) * CSel(y,b,s),
  //      Bf(y,b,s) = BF(y) * BSel(y,b,s).
  // Z(y,b,s) is sum of fishing mortalities and M(y,s).
  3darray	Cf(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  3darray	Df(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  3darray	Bf(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  3darray	Rf(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  3darray	Pf(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  3darray	 Z(FirstYear,LastYear,1,LastAge,Sex1,Sex2)
  matrix	Z_F(FirstYear,LastYear,1,LastAge) // For reporting.
  matrix	Z_M(FirstYear,LastYear,1,LastAge) // For reporting.
  //
  // Full population matrix (Jan. 1), recruitment and initial numbers.
  // Note that the working matrices are extended beyond the
  // last year to the year after to do projections routinely.
  // INitN and R parameters are numbers * 1E-6.
  //
  3darray	N(FirstYear,LastYear+1,1,LastAge,Sex1,Sex3)
  matrix	N_F(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
  matrix	N_M(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
  matrix	N_T(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
  3darray	NBar(FirstYear,LastYear,1,LastAge,Sex1,Sex3)
  matrix	NBar_F(FirstYear,LastYear,1,LastAge) // For REPORTing.
  matrix	NBar_M(FirstYear,LastYear,1,LastAge) // For REPORTing.
  matrix	NBar_T(FirstYear,LastYear,1,LastAge) // For REPORTing.
                // Commercial exploitable nos. at true age.
  3darray	CEN(FirstYear,LastYear+1,1,LastAge,Sex1,Sex3)
  matrix	CEN_F(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
  matrix	CEN_M(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
  matrix	CEN_T(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
                // Commercial exploitable nos. at true age calculated 
                // with fixed coastwide length-specific selectivities.
  3darray	FEN(FirstYear,LastYear+1,1,LastAge,Sex1,Sex3)
  matrix	FEN_F(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
  matrix	FEN_M(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
  matrix	FEN_T(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
                // Survey exploitable nos. at true age.
  3darray	SEN(FirstYear,LastYear+1,1,LastAge,Sex1,Sex3)
  matrix	SEN_F(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
  matrix	SEN_M(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
  matrix	SEN_T(FirstYear,LastYear+1,1,LastAge) // For REPORTing.
  //
  matrix	InitN(2,LastAge,Sex1,Sex3)
  init_bounded_matrix	InitNSeenPar(2,AgePlusFirst,Sex1,Sex2,0.0,RMMax)
  matrix	InitNSeen(2,AgePlusFirst,Sex1,Sex2)
  matrix	R(FirstYear,LastYear+1,Sex1,Sex3)
  init_bounded_vector	RSeenPar(FirstYear,LastYearClass+1,0.0,RMMax)
  init_bounded_number	PropFem(0.4,0.6,-1)
  vector	RSeen(FirstYear,LastYearClass+1)
  number	RMean
  //
  // Pro forma parameter to estimate in phase 2 so that robustification
  // can be switched on after the ordinary least-squares estimate has
  // been located. Its true value is 0; it is initialized to 0.
  // Phase P_ProForma = 1 or 2 depending on RobustifyFlag.
  //
  init_number	ProForma(P_ProForma)
  //
  // Check value at end of .pin file.
  //
  init_number	ParamCheck(-1)
  //
  !! if (ParamCheck != CheckValue)
  !!     {
  !!     cerr << "Parameters are not in order." << endl;
  !!     cerr << "Check value is " << ParamCheck << endl;
  !!     exit(1);
  !!     } else
  !!     cout << "Parameters are in order." << endl;
  //
  // End model parameters.
  //
  // PREDICTIONS
  //
  // Predictions of the observations and related quantities.
  //
  // The population matrix redistributed to smeared ages.
  //
  matrix	NS_F(FirstYear,LastYear+1,SMinAge,SMaxAge)
  matrix	NS_M(FirstYear,LastYear+1,SMinAge,SMaxAge)
  //
  // The commercial exploitable population redistributed to smeared ages.
  //
  matrix	CENS_F(FirstYear,LastYear+1,SMinAge,SMaxAge)
  matrix	CENS_M(FirstYear,LastYear+1,SMinAge,SMaxAge)
  //
  // The commercial exploitable population, calculated with the fixed
  // coastwide length-specific selectivities, redistributed to smeared ages.
  //
  matrix	FENS_F(FirstYear,LastYear+1,SMinAge,SMaxAge)
  matrix	FENS_M(FirstYear,LastYear+1,SMinAge,SMaxAge)
  //
  // The survey exploitable population redistributed to smeared ages.
  //
  matrix	SENS_F(FirstYear,LastYear+1,SMinAge,SMaxAge)
  matrix	SENS_M(FirstYear,LastYear+1,SMinAge,SMaxAge)
  //
  // Commercial catch at age/sex.
  //
  matrix	Catch_Pred_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	Catch_Pred_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	Catch_Pred_T(FirstYear,LastYear,SMinAge,SMaxAge)
  number	Catch_RSS
  number	Catch_RSS_F
  number	Catch_RSS_M
  number	Catch_RSS_T
  number	Catch_n
  number	Catch_n_F
  number	Catch_n_M
  number	Catch_n_T
  number	Catch_rms
  number	Catch_rms_F
  number	Catch_rms_M
  number	Catch_rms_T
  vector	RSS_list(1,3)	// CalcRRSS() returns a vector.
  //
  // Commercial CPUE at age/sex.
  //
  matrix	CCPUE_Pred_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	CCPUE_Pred_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	CCPUE_Pred_T(FirstYear,LastYear,SMinAge,SMaxAge)
  number	CCPUE_RSS
  number	CCPUE_RSS_F
  number	CCPUE_RSS_M
  number	CCPUE_RSS_T
  number	CCPUE_n
  number	CCPUE_n_F
  number	CCPUE_n_M
  number	CCPUE_n_T
  number	CCPUE_rms
  number	CCPUE_rms_F
  number	CCPUE_rms_M
  number	CCPUE_rms_T
  //
  // Commercial CPUE in total.
  //
  vector	CCPUETot_Pred(FirstYear,LastYear)
  number	CCPUETot_RSS
  number	CCPUETot_n
  number	CCPUETot_rms
  //
  // Commercial WPUE in total.
  //
  vector	CWPUETot_Pred(FirstYear,LastYear)
  number	CWPUETot_RSS
  number	CWPUETot_n
  number	CWPUETot_rms
  //
  // Survey proportion at age/sex.
  //
  matrix	SurvProp_Pred_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvProp_Pred_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvProp_Pred_T(FirstYear,LastYear,SMinAge,SMaxAge)
  number	SurvProp_RSS
  number	SurvProp_RSS_F
  number	SurvProp_RSS_M
  number	SurvProp_RSS_T
  number	SurvProp_n
  number	SurvProp_n_F
  number	SurvProp_n_M
  number	SurvProp_n_T
  number	SurvProp_rms
  number	SurvProp_rms_F
  number	SurvProp_rms_M
  number	SurvProp_rms_T
  //
  // Survey CPUE at age/sex.
  //
  matrix	SurvCPUE_Pred_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvCPUE_Pred_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SurvCPUE_Pred_T(FirstYear,LastYear,SMinAge,SMaxAge)
  number	SurvCPUE_RSS
  number	SurvCPUE_RSS_F
  number	SurvCPUE_RSS_M
  number	SurvCPUE_RSS_T
  number	SurvCPUE_n
  number	SurvCPUE_n_F
  number	SurvCPUE_n_M
  number	SurvCPUE_n_T
  number	SurvCPUE_rms
  number	SurvCPUE_rms_F
  number	SurvCPUE_rms_M
  number	SurvCPUE_rms_T
  //
  // Survey CPUE in total.
  //
  vector	SurvCPUETot_Pred(FirstYear,LastYear)
  number	SurvCPUETot_RSS
  number	SurvCPUETot_n
  number	SurvCPUETot_rms
  //
  // Survey WPUE in total.
  //
  vector	SurvWPUETot_Pred(FirstYear,LastYear)
  number	SurvWPUETot_RSS
  number	SurvWPUETot_n
  number	SurvWPUETot_rms
  //
  // Total sublegal discard mortality in weight by year.
  //
  vector	DiscardWt_Pred(FirstYear,LastYear)
  number	DiscardWt_RSS
  number	DiscardWt_n
  number	DiscardWt_rms
  //
  // Pro forma sd used for computing the RSS for all catch in weight totals.
  //
  vector	XCatchWt_SE(FirstYear,LastYear) // Pro forma for computing an RSS.
  !!		XCatchWt_SE.fill_seqadd(10000.0,0.0);
  //
  // Bycatch in number at length.
  //
  matrix	Bycatch_Pred(FirstYear,LastYear,1,BLptN)
  number	Bycatch_RSS
  number	Bycatch_n
  number	Bycatch_rms
  //
  // Total bycatch in weight by year.
  //
  vector	BycatchWt_Pred(FirstYear,LastYear)
  number	BycatchWt_RSS
  number	BycatchWt_n
  number	BycatchWt_rms
  //
  // Total recreational catch in weight by year.
  //
  vector	SportCatchWt_Pred(FirstYear,LastYear)
                // Smeared model catches computed with SSel:
  matrix	SportCatch_Pred_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SportCatch_Pred_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	SportCatch_Pred_T(FirstYear,LastYear,SMinAge,SMaxAge)
  number	SportCatchWt_RSS
  number	SportCatchWt_n
  number	SportCatchWt_rms
  //
  // Total personal use catch in weight by year.
  //
  vector	PersUseWt_Pred(FirstYear,LastYear)
                // Smeared model catches computed with SSel:
  matrix	PersUse_Pred_F(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	PersUse_Pred_M(FirstYear,LastYear,SMinAge,SMaxAge)
  matrix	PersUse_Pred_T(FirstYear,LastYear,SMinAge,SMaxAge)
  number	PersUseWt_RSS
  number	PersUseWt_n
  number	PersUseWt_rms
  //
  // Total catch in number each year. (Determines CF(y) if y > CLastYearFit.)
  //		
  vector	CatchN_Pred(FirstYear,LastYear)
  number	CatchN_RSS
  number	CatchN_n
  number	CatchN_rms
  vector	XCatchN_SE(FirstYear,LastYear) // Pro forma for computing an RSS.
  !!		XCatchN_SE.fill_seqadd(1000.0,0.0);
  //
  // Penalty on relative annual changes in commercial catchability CQ,
  // computed with SteadyCQ_SD.
  //
  number	SteadyCQ_PSS
  //
  // Penalty on unsmooth selectivities.
  //
  number	SelSmooth_PSS
  //
  // Penalty on deviation of CF(LastYear) from LastF if FitLastFFlag is set;
  // pseudo-sd for calculations.
  //
  number	LastF_PSS
  number	LastF_SD
  !!		LastF_SD = 0.0001;
  //
  // Penalty on uneven sex ratio in InitN estimates.
  // UnevenSexRatio_SD = 0.01 normally.
  //
  number	UnevenSexRatio_PSS
  number	UnevenSexRatio_SD
  !!		UnevenSexRatio_SD = 0.01;
  //
  // Penalty on large differences between successive InitN estimates.
  // Assessed on log differences.
  // WildYearClass_SD = 20.0/(LastYear - FirstYear) in closed-area fits.
  // WildYearClass_SD = 2.3 (order of magnitude) in coastwide fits.
  //
  number	WildYearClass_PSS
  number	WildYearClass_SD
  !!		WildYearClass_SD = 2.3;  // 2.3 = log(10)
  //
  // Penalty on large log differences between successive RSeen estimates.
  //
  number	WildR_PSS
  number	WildR_SD
  !!		WildR_SD = 2.3; // 2.3 = log(10) 
  //
  // Penalty on variance of annual estimates of survey catchability SQ,
  // computed with SteadySQ_SD.
  //
  number	SteadySQ_PSS
  //
  // Penalty on departure of annual estimates of SQ from zero slope,
  // computed with SQ_bhat_SD and computed along with SteadySQ_PSS 
  // (and not computed if SteadySQFlag==0).
  //
  number	TrendlessSQ_PSS
  number	SQ_bhat_SD
  !!		if (TrendlessSQFlag == 0) SQ_bhat_SD = 100.0; else SQ_bhat_SD = 0.0001;
  //
  // That's all the predictions of the observations. 
  //
  // Model catches in all fisheries in number at true age.
  //
  3darray	CCatch(FirstYear,LastYear,1,LastAge,Sex1,Sex3)
  matrix	CCatch_F(FirstYear,LastYear,1,LastAge)
  matrix	CCatch_M(FirstYear,LastYear,1,LastAge)
  matrix	CCatch_T(FirstYear,LastYear,1,LastAge)
  //
  3darray	DCatch(FirstYear,LastYear,1,LastAge,Sex1,Sex3)
  matrix	DCatch_F(FirstYear,LastYear,1,LastAge)
  matrix	DCatch_M(FirstYear,LastYear,1,LastAge)
  matrix	DCatch_T(FirstYear,LastYear,1,LastAge)
  //
  3darray	BCatch(FirstYear,LastYear,1,LastAge,Sex1,Sex3)
  matrix	BCatch_F(FirstYear,LastYear,1,LastAge)
  matrix	BCatch_M(FirstYear,LastYear,1,LastAge)
  matrix	BCatch_T(FirstYear,LastYear,1,LastAge)
  //
  3darray	RCatch(FirstYear,LastYear,1,LastAge,Sex1,Sex3)
  matrix	RCatch_F(FirstYear,LastYear,1,LastAge)
  matrix	RCatch_M(FirstYear,LastYear,1,LastAge)
  matrix	RCatch_T(FirstYear,LastYear,1,LastAge)
  //
  3darray	PCatch(FirstYear,LastYear,1,LastAge,Sex1,Sex3)
  matrix	PCatch_F(FirstYear,LastYear,1,LastAge)
  matrix	PCatch_M(FirstYear,LastYear,1,LastAge)
  matrix	PCatch_T(FirstYear,LastYear,1,LastAge)
  //
  3darray	MCatch(FirstYear,LastYear,1,LastAge,Sex1,Sex3)
  matrix	MCatch_F(FirstYear,LastYear,1,LastAge)
  matrix	MCatch_M(FirstYear,LastYear,1,LastAge)
  matrix	MCatch_T(FirstYear,LastYear,1,LastAge)
  //
  // Model CPUE in numbers at true age.
  //
  matrix	CCPUE_True_F(FirstYear,LastYear,1,LastAge)
  matrix	CCPUE_True_M(FirstYear,LastYear,1,LastAge)
  matrix	CCPUE_True_T(FirstYear,LastYear,1,LastAge)
  //
  matrix	SurvProp_True_F(FirstYear,LastYear,1,LastAge)
  matrix	SurvProp_True_M(FirstYear,LastYear,1,LastAge)
  matrix	SurvProp_True_T(FirstYear,LastYear,1,LastAge)
  //
  matrix	SurvCPUE_True_F(FirstYear,LastYear,1,LastAge)
  matrix	SurvCPUE_True_M(FirstYear,LastYear,1,LastAge)
  matrix	SurvCPUE_True_T(FirstYear,LastYear,1,LastAge)
  //
  // OUTPUTS
  //
  // Routine model outputs.
  //
  vector		R6(FirstYear,LastYear+1)
  vector		R8(FirstYear,LastYear+1)
  vector		N8Plus_F(FirstYear,LastYear+1)
  vector		N8Plus_M(FirstYear,LastYear+1)
  vector		TBio8Plus_F(FirstYear,LastYear+1)
  vector		TBio8Plus_M(FirstYear,LastYear+1)
  vector		LBio(FirstYear,LastYear+1)
  vector		SBio(FirstYear,LastYear+1)
  vector		EBio(FirstYear,LastYear+1)
  vector	 	EBioFixed(FirstYear,LastYear+1)
  sdreport_number	LastEBioFixed
  likeprof_number	LastEBioFixedToo
  likeprof_number	LastSBio
  vector		EBioSurv(FirstYear,LastYear+1)
  sdreport_vector	R8SD(FirstYear,LastYear+1)
  //

PROCEDURE_SECTION
  //
  // cout << "Hi, Mommy!" << endl;
  //
  // Initialize objective function components.
  //
  SSTotal = 0.0;
  RSSTotal = 0.0;
  RSSNobs = 0.0;
  PSSTotal = 0.0;
  //
  // Fill the working arrays of age- and sex-specfic values
  // M, CF, CQ, CSel, SQ, SSel, BF, BSel etc. according to their
  // various parameterizations.
  //
     Fill_M(); 
  //
     Fill_CF();
  //
     if (FitLastFFlag != 0) Compute_LastF_PSS();
  //
     Fill_CQ();
  //
     Compute_SteadyCQ_PSS();
  //
     Fill_CSelL();
  //
     Fill_DF();
  //
     Fill_SQ();
  //
     if (SteadySQFlag != 0) Compute_SteadySQ_PSS();
  //
     Fill_SSelL();
  //
     Fill_RF_and_PF();
  //
     Fill_BSelL();
  //
     if (SelSmoothFlag==1) Compute_SelSmooth_PSS();
  //
     Compute_selectivity_at_age();
  //
     Compute_mortality_rates_at_age();
  //
     Fill_InitN();
  //
     Fill_R();
  //
     Initialize_N();
  //
     Compute_removals_and_survivors();
  //
  // Maybe switch on the internal robustification flag.
     if (RobustifyFlag != 0 && current_phase() == 2) Robust = 1;
        else Robust = 0;
  //
     Compute_Catch_RSS();
  //
     Compute_CCPUE_RSS();
  //
     Compute_SurvProp_RSS();
  //
     Compute_SurvCPUE_RSS();
  //
     Compute_CPUE_Totals_RSS();
  //
     Compute_DiscardWt_RSS();
  //
     Compute_Bycatch_RSS();
     Compute_BycatchWt_RSS();
  //
     Compute_SportCatchWt_RSS();
  //
     Compute_PersUseWt_RSS();
  //
     Compute_CatchN_RSS();
  //
     SSTotal = RSSTotal + PSSTotal;
     SSTotal += sq(ProForma); // Initialized to zero and kept there.
     if (MultiFitFlag==0)
	     MinusLogL = SSTotal/2; else
	MinusLogL=SSTotal/8.0;
  //
     if (sd_phase()) Compute_outputs();
  //
  // cout << "Bye, Mommy!" << endl;
  //

FUNCTION Fill_M
  //
  // Fill in the natural mortality matrix.
  //
  for (y=FirstYear; y<=LastYear; y++)
     {
     M(y,1) = MPar;
     if (EstMMFlag != 0) M(y,2) = MPar_M; else M(y,2) = MPar;
     }

FUNCTION Fill_CF
  //
  // Fill in the commercial fishing mortality matrix
  //
  for (y=FirstYear; y<=LastYear; y++)
     {
     CF(y,1) = CFPar(y); // Females.
     CF(y,2) = CQMX * CFPar(y); // Males.
     }
   //

FUNCTION Compute_LastF_PSS
  //
  // Assess a heavy penalty on deviations from the specified LastF.
  //
  LastF_PSS = sq( (CFPar(LastYear) - LastF)/LastF_SD );
  PSSTotal += LastF_PSS;
  //


FUNCTION  Fill_CQ
  //
  // Fill in working matrix CQ from waypoint and path type params.
  // GOOD_COP has made sure that CQ is estimated in FirstYear.
  //
  // Scale the parameter.
  //
  for (w=1; w<=CQEstN; w++) CQEst(w) = CQEstPar(w) * 1.0E-6;
  //
  for (w=1; w<=CQEstN; w++)
     {
     // Starting value on this segment is CQEst(w).
     StartYear = CQEstYears(w);
     StartVal = CQEst(w);
     // Figure out StopYear and StopVal.
     if (w == CQEstN)
        {
        // Last (or only) waypoint. Use this value for remainder
        // of series regardless of value of CQEstPathTypes(w).
        StopYear = LastYear;
        StopVal = StartVal;
        } else
        {
        // First or intermediate waypoint. Determine StopYear, StopVal.
        StopYear = CQEstYears(w+1);
        if (CQEstPathTypes(w) == 0) 
           StopVal = StartVal; else StopVal = CQEst(w+1);
        }
     // Compute PathLength and StepVal.
     PathLength = StopYear - StartYear;
     StepVal = (StopVal - StartVal);
     if (PathLength > 0) StepVal = StepVal/PathLength;
     // Assign values point by point.
     for (istep=0; istep<=PathLength; istep++)
        {
        y = StartYear + istep;
        CQ(y,1) = StartVal + istep * StepVal; // Females.
        CQ(y,2) = CQMX * (StartVal + istep * StepVal); // Males.
        }
     // Done with this segment.
     } // End loop on segment index w.
  //

FUNCTION Compute_SteadyCQ_PSS
  //
  // Compute SteadyCQ_PSS
  //
  dvariable DeltaCQ;
  SteadyCQ_PSS = 0.0;
  for (s=1; s<=2; s++)
     for (y=FirstYear+1; y<=CLastYearFit; y++)
        {
        DeltaCQ = (CQ(y,s) - CQ(y-1,s)) / CQ(y-1,s); 
        if (y==1984 && DeltaCQ > HookCorrMax-1)
           SteadyCQ_PSS += sq( (DeltaCQ - (HookCorrMax-1)) / 0.01 );
        if (y!=1984)
           SteadyCQ_PSS += sq(DeltaCQ/SteadyCQ_SD);
        }
  PSSTotal += SteadyCQ_PSS;
  //

FUNCTION Fill_CSelL
  //
  // Fill in working array CSelL from path & waypoint params.
  // GOOD_COP has assured that an estimate is made in FirstYear.
  //
  // Compute implied estimates CSelLEst at waypoints.
  // Interpretation of params depends on CSelLForm.
  //
  // First load split-sex parameter array CSelLParSplit.
  //  If (SplitMalesFlag) load CSelLPar_M into male tier else CSelLPar.
  for (w=1; w<=CSelLEstN; w++)
     for (p=1; p<=CSelLParN; p++)
        {
        CSelLParSplit(w,p,1) = CSelLPar(w,p); // Females
        if (SplitMalesFlag == 0) CSelLParSplit(w,p,2) = CSelLPar(w,p);
                            else CSelLParSplit(w,p,2) = CSelLPar_M(w,p);
        } // End loops on w, p.
  //
  // Build split-sex selectivity array CSelLEst.
  //
  for (s=1; s<=2; s++)
  {
  if (CSelLForm==1)
  // Free except zero at SLptList(1) and one at SLptOne. 
  // Scale by 5 because params are all (0,1).
  for (w=1; w<=CSelLEstN; w++)
     {
     CSelLEst(w,1,s) = 0.0;
     CSelLEst(w,vOne,s) = 1.0;
     for (v=2; v<vOne; v++) CSelLEst(w,v,s) = 5.0 * CSelLParSplit(w,v-1,s);
     for (v=vOne+1; v<=SLptN; v++) CSelLEst(w,v,s) = 5.0 * CSelLParSplit(w,v-2,s);
     } // End loop on waypoint w.
  else if (CSelLForm==2)
  // Domed on either side of SLptOne
  for (w=1; w<=CSelLEstN; w++)
      {
      CSelLEst(w,vOne,s) = 1.0;
      for (v=vOne-1; v>=1; v--) CSelLEst(w,v,s) = CSelLEst(w,v+1,s) * CSelLParSplit(w,v,s);
      for (v=vOne+1; v<=SLptN; v++) CSelLEst(w,v,s) = CSelLEst(w,v-1,s) * CSelLParSplit(w,v-1,s);
      } // End loop on waypoint w.
  else if (CSelLForm==3)
  // Asymptotic.
  for (w=1; w<=CSelLEstN; w++)
      {
      CSelLEst(w,vOne,s) = 1.0;
      for (v=vOne-1; v>=1; v--) CSelLEst(w,v,s) = CSelLEst(w,v+1,s) * CSelLParSplit(w,v,s);
      for (v=vOne+1; v<=SLptN; v++) CSelLEst(w,v,s) = 1.0;
      } // End loop on waypoint w.
  //
  // Assign waypoints and compute paths.
  //
  for (w=1; w<=CSelLEstN; w++)
     {
     // Starting value on this segment is CSelLEst(w).
     StartYear = CSelLEstYears(w);
     for (v=1; v<=SLptN; v++)
        StartVec(v) = CSelLEst(w,v,s);
     // Figure out StopYear and StopVec.
     if (w == CSelLEstN)
        {
        // Last (or only) waypoint. Use this value for remainder
        // of series regardless of value of CSelLEstPathTypes(w).
        StopYear = LastYear;
        StopVec = StartVec;
        } else
        {
        // First or intermediate waypoint. Determine StopYear, StopVec.
        StopYear = CSelLEstYears(w+1);
        if (CSelLEstPathTypes(w) == 0) 
           StopVec = StartVec;
           else for(v=1; v<=SLptN; v++) StopVec(v) = CSelLEst(w+1,v,s);
        }
     // Compute PathLength and StepVec.
     PathLength = StopYear - StartYear;
     StepVec = (StopVec - StartVec);
     if (PathLength > 0) StepVec = StepVec/PathLength;
     // Assign values point by point.
     for (istep=0; istep<=PathLength; istep++)
     for (v=1; v<=SLptN; v++)
        CSelL(StartYear+istep,v,s) = StartVec(v) + istep * StepVec(v);
     // Done with this segment.
     } // End loop on waypoint.
  } // End loop on sex.
  //

FUNCTION Fill_DF
  //
  // Fill in the sublegal discard mortality matrix
  //
  for (y=FirstYear; y<=LastYear; y++)
     DF(y) = DFPar(y);
   //

FUNCTION  Fill_SQ
  //
  // Fill in working matrix SQ from waypoint and path type params.
  // GOOD_COP has made sure that SQ is estimated in FirstYear.
  //
  // Scale the parameter.
  //
  for(w=1; w<=SQEstN; w++) SQEst(w) = SQEstPar(w) * 1.0E-6;
  //
  for (w=1; w<=SQEstN; w++)
     {
     // Starting value on this segment is SQEst(w).
     StartYear = SQEstYears(w);
     StartVal = SQEst(w);
     // Figure out StopYear and StopVal.
     if (w == SQEstN)
        {
        // Last (or only) waypoint. Use this value for remainder
        // of series regardless of value of SQEstPathTypes(w).
        StopYear = LastYear;
        StopVal = StartVal;
        } else
        {
        // First or intermediate waypoint. Determine StopYear, StopVal.
        StopYear = SQEstYears(w+1);
        if (SQEstPathTypes(w) == 0) 
           StopVal = StartVal; else StopVal = SQEst(w+1);
        }
     // Compute PathLength and StepVal.
     PathLength = StopYear - StartYear;
     StepVal = (StopVal - StartVal);
     if (PathLength > 0) StepVal = StepVal/PathLength;
     // Assign values point by point.
     for (istep=0; istep<=PathLength; istep++)
        {
        y = StartYear + istep;
        SQ(y,1) = StartVal + istep * StepVal; // Females.
        SQ(y,2) = SQMX * (StartVal + istep * StepVal); // Males.
        } // End loop on istep.
     // Done with this segment.
     } // End loop on w.
  //

FUNCTION Compute_SteadySQ_PSS
  //
  // Compute penalty on variance in log(SQEst) -- post 1984 only.
  // Also compute penalty on any trend (departure from zero slope).
  //
  //sig_PSS = 0.0;
  SteadySQ_PSS = 0.0;
  dvariable SQSum;
  dvariable SQMean;
  dvariable YearSum;
  dvariable YearMean;
  dvariable xySum;
  dvariable xxSum;
  SQSum = 0.0;
  YearSum = 0.0;
  int nyrs = 0;
  for (y=FirstYear; y<=SLastYearFit; y++)
     {
     if (y < 1984) continue;
     nyrs += 1;
     SQSum += SQ(y,1);
     YearSum += y;
     }
  if (nyrs > 1)
     {
     SQMean = SQSum/nyrs;
     YearMean = YearSum/nyrs;
     xxSum = 0.0;
     xySum = 0.0;
     for (y=FirstYear; y<=SLastYearFit; y++)
        {
        if (y < 1984) continue;
        SteadySQ_PSS += sq( (log(SQ(y,1)) - log(SQMean) ) / SteadySQ_SD );
        xxSum += sq(y - YearMean);
        xySum += (y - YearMean) * (log(SQ(y,1)) - log(SQMean));
        } // End loop on y.
     } // End condition on nyrs.
  SQ_bhat = xySum / xxSum;
  TrendlessSQ_PSS = sq(SQ_bhat/SQ_bhat_SD);
  if (TrendlessSQFlag != 0) PSSTotal += TrendlessSQ_PSS;
  if (SteadySQFlag != 0) PSSTotal += SteadySQ_PSS;
  //

FUNCTION Fill_SSelL
  //
  // Fill in working array SSelL from path & waypoint params.
  // GOOD_COP has assured that an estimate is made in FirstYear.
  //
  // Compute implied estimates SSelLEst at waypoints.
  // Interpretation of params depends on SSelLForm.
  //
  // First load split-sex parameter array SSelLParSplit.
  //  If (SplitMalesFlag) load SSelLPar_M into male tier else SSelLPar.
  for (w=1; w<=SSelLEstN; w++)
     for (p=1; p<=SSelLParN; p++)
        {
        SSelLParSplit(w,p,1) = SSelLPar(w,p); // Females
        if (SplitMalesFlag == 0) SSelLParSplit(w,p,2) = SSelLPar(w,p);
                            else SSelLParSplit(w,p,2) = SSelLPar_M(w,p);
        } // End loops on w, p.
  //
  // Build split-sex selectivity array SSelLEst.
  //
  for (s=1; s<=2; s++)
  {
  if (SSelLForm==1)
  // Free except zero at SLptList(1) and one at SLptOne. 
  // Scale by 5 because params are all (0,1).
  for (w=1; w<=SSelLEstN; w++)
     {
     SSelLEst(w,1,s) = 0.0;
     SSelLEst(w,vOne,s) = 1.0;
     for (v=2; v<vOne; v++) SSelLEst(w,v,s) = 5.0 * SSelLParSplit(w,v-1,s);
     for (v=vOne+1; v<=SLptN; v++) SSelLEst(w,v,s) = 5.0 * SSelLParSplit(w,v-2,s);
     } // End loop on waypoint w.
  else if (SSelLForm==2)
  // Domed on either side of SLptOne
  for (w=1; w<=SSelLEstN; w++)
      {
      SSelLEst(w,vOne,s) = 1.0;
      for (v=vOne-1; v>=1; v--) SSelLEst(w,v,s) = SSelLEst(w,v+1,s) * SSelLParSplit(w,v,s);
      for (v=vOne+1; v<=SLptN; v++) SSelLEst(w,v,s) = SSelLEst(w,v-1,s) * SSelLParSplit(w,v-1,s);
      } // End loop on waypoint w.
  else if (SSelLForm==3)
  // Asymptotic.
  for (w=1; w<=SSelLEstN; w++)
      {
      SSelLEst(w,vOne,s) = 1.0;
      for (v=vOne-1; v>=1; v--) SSelLEst(w,v,s) = SSelLEst(w,v+1,s) * SSelLParSplit(w,v,s);
      for (v=vOne+1; v<=SLptN; v++) SSelLEst(w,v,s) = 1.0;
      } // End loop on waypoint w.
  //
  // Assign waypoints and compute paths.
  //
  for (w=1; w<=SSelLEstN; w++)
     {
     // Starting value on this segment is SSelLEst(w).
     StartYear = SSelLEstYears(w);
     for (v=1; v<=SLptN; v++)
        StartVec(v) = SSelLEst(w,v,s);
     // Figure out StopYear and StopVec.
     if (w == SSelLEstN)
        {
        // Last (or only) waypoint. Use this value for remainder
        // of series regardless of value of SSelLEstPathTypes(w).
        StopYear = LastYear;
        StopVec = StartVec;
        } else
        {
        // First or intermediate waypoint. Determine StopYear, StopVec.
        StopYear = SSelLEstYears(w+1);
        if (SSelLEstPathTypes(w) == 0) 
           StopVec = StartVec;
           else for(v=1; v<=SLptN; v++) StopVec(v) = SSelLEst(w+1,v,s);
        }
     // Compute PathLength and StepVec.
     PathLength = StopYear - StartYear;
     StepVec = (StopVec - StartVec);
     if (PathLength > 0) StepVec = StepVec/PathLength;
     // Assign values point by point.
     for (istep=0; istep<=PathLength; istep++)
     for (v=1; v<=SLptN; v++)
        SSelL(StartYear+istep,v,s) = StartVec(v) + istep * StepVec(v);
     // Done with this segment.
     } // End loop on waypoint.
  } // End loop on sex.
  //

FUNCTION Fill_RF_and_PF
  //
  // Fill in the recreational and psonal use fishing mortality matrices.
  //
  for (y=FirstYear; y<=LastYear; y++)
     {
     RF(y,1) = RFPar(y); // Females.
     RF(y,2) = SQMX * RFPar(y); // Males.
     PF(y,1) = PFPar(y); // Females.
     PF(y,2) = SQMX * PFPar(y); // Males.
     }
   //

FUNCTION Fill_BSelL
  //
  // Fill in working array BSelL from path and waypoint params.
  // GOOD_COP has assured that an estimate is made in FirstYear.
  //
  // Compute implied estimates BSelLEst at waypoints.
  // For bycatch selectivity, each row of the param matrix BSelLPar
  // is a set of multipliers to be applied on either side of BLptOne,
  // so a row has BSelLParN = BLptN-1 elements.
  // The limits placed on BSelLPar determine whether any of the
  // bycatch selectivities can exceed 1.0.
  //
  for (w=1; w<=BSelLEstN; w++)
     {
     BSelLEst(w,vOneB) = 1.0;
     for (v=vOneB-1; v>=1; v--) BSelLEst(w,v) = BSelLEst(w,v+1) * BSelLPar(w,v);
     for (v=vOneB+1; v<=BLptN; v++) BSelLEst(w,v) = BSelLEst(w,v-1) * BSelLPar(w,v-1);
     }
  //
  // Assign waypoints and compute paths.
  //
  for (w=1; w<=BSelLEstN; w++)
     {
     // Starting value on this segment is BSelLEst(w).
     StartYear = BSelLEstYears(w);
     StartVecB = BSelLEst(w);
     // Figure out StopYear and StopVecB.
     if (w == BSelLEstN)
        {
        // Last (or only) waypoint. Use this value for remainder
        // of series regardless of value of BSelLEstPathTypes(w).
        StopYear = LastYear;
        StopVecB = StartVecB;
        } else
        {
        // First or intermediate waypoint. Determine StopYear, StopVecB.
        StopYear = BSelLEstYears(w+1);
        if (BSelLEstPathTypes(w) == 0) 
           StopVecB = StartVecB; else StopVecB = BSelLEst(w+1);
        }
     // Compute PathLength and StepVecB.
     PathLength = StopYear - StartYear;
     StepVecB = (StopVecB - StartVecB);
     if (PathLength > 0) StepVecB = StepVecB/PathLength;
     // Assign values point by point.
     for (istep=0; istep<=PathLength; istep++)
     for (v=1; v<=BLptN; v++)
     for (s=1; s<=2; s++)
        BSelL(StartYear+istep,v,s) = StartVecB(v) + istep * StepVecB(v);
     // Done with this segment.
     }
  //
     
FUNCTION Compute_selectivity_at_age
  //
  // Compute all the age-specific selectivities from length-specific
  // schedules and mean (true) length in either the survey or the bycatch.
  //
  for (y=FirstYear; y<=LastYear; y++)
     {
     // Construct year-specific matrices CSelLy and SSelLy.
     // Each row is a sex-specific vector that can be used with approx().
     for (s=1; s<=2; s++)
     for (v=1; v<= SLptN; v++)
        {
        CSelLy(s,v) = CSelL(y,v,s);
        SSelLy(s,v) = SSelL(y,v,s);
        }
     // Same for BSelLy.
     for (s=1; s<=2; s++)
     for (v=1; v<=BLptN; v++)
        BSelLy(s,v) = BSelL(y,v,s);
     //
     // Run through ages and sexes doing the interpolation.
     for (b=1; b<=LastAge; b++)
     for (s=1; s<=2; s++)
        {
        // First hook-and-line fisheries.
        if (b < SMinAge)
           {
           CSel(y,b,s) = 0.0;
           DSel(y,b,s) = 0.0;
           FSel(y,b,s) = 0.0;
           SSel(y,b,s) = 0.0;
           } else
           {
           if (SexFlag==1)
              { // Use sex-specific mean length at age.
              CSel(y,b,s) = approx(SLptList, CSelLy(s), SurvLTrue(y,b,s), 0);
              DSel(y,b,s) = approx(SLptList,     DSelL, SurvLTrue(y,b,s), 0);
              FSel(y,b,s) = approx(SLptList,     	, SurvLTrue(y,b,s), 0);
              SSel(y,b,s) = approx(SLptList, SSelLy(s), SurvLTrue(y,b,s), 0);
              } else
              { // Use overall mean length at age for both sexes.
              CSel(y,b,s) = approx(SLptList, CSelLy(s), SurvLTrue(y,b,3), 0);
              DSel(y,b,s) = approx(SLptList,     DSelL, SurvLTrue(y,b,3), 0);
              FSel(y,b,s) = approx(SLptList,     FSelL, SurvLTrue(y,b,3), 0);
              SSel(y,b,s) = approx(SLptList, SSelLy(s), SurvLTrue(y,b,3), 0);
              }
           }
        RSel(y,b,s) = SSel(y,b,s);
        PSel(y,b,s) = SSel(y,b,s);
        // Bycatch.
        if (SexFlag==1)
           BSel(y,b,s) = approx(BLptList, BSelLy(s), BycatchLTrue(y,b,s), 0);
           else BSel(y,b,s) = approx(BLptList, BSelLy(s), BycatchLTrue(y,b,3), 0);
        } // End age/sex loop.
     } // End year loop.
  //

FUNCTION Compute_mortality_rates_at_age
  //
  // Compute all the year/age/sex-specific mortalities.
  //
  for (y=FirstYear; y<=LastYear; y++)
  for (b=1; b<=LastAge; b++)
  for (s=1; s<=2; s++)
     {
     Cf(y,b,s) = CSel(y,b,s) * CF(y,s);
     Df(y,b,s) = DSel(y,b,s) * DF(y);
     Bf(y,b,s) = BSel(y,b,s) * BF(y);
     Rf(y,b,s) = RSel(y,b,s) * RF(y,s);
     Pf(y,b,s) = PSel(y,b,s) * PF(y,s);
      Z(y,b,s) = Cf(y,b,s) + Df(y,b,s) + Bf(y,b,s) + Rf(y,b,s) + Pf(y,b,s) + M(y,s);
     }

FUNCTION Fill_InitN
  //
  // Fill in the first row of the population arrays InitN. To avoid some 
  // poorly determined estimates of initial abundance of the oldest (true) 
  // age groups, they are not estimated but instead computed from the 
  // initial abundance at true age AgePlusFirst (which is really an 
  // observed age), like 20. As a result the total abundance estimate for
  // the fish that were (true age) 20+ in 1974 should be OK, but their 
  // distn among ages will be assigned.
  //
  // Scale the parameter.
  //
  for (b=2; b<=AgePlusFirst; b++)
     for(s=1; s<=2; s++)
        InitNSeen(b,s) = InitNSeenPar(b,s) * 1.0E6;
  //
  // Seen well enough.
  for (b=2; b<=AgePlusFirst; b++)
     {
     for (s=1; s<=2; s++) InitN(b,s) = InitNSeen(b,s);
     InitN(b,3) = InitN(b,1) + InitN(b,2);
     }
  // Not seen well enough.
  // Get Z by sex in FirstYear.
  for (s=1; s<=2; s++)
     OldZ(s) = CF(FirstYear,s) + M(FirstYear,s);
  //
  for (b=AgePlusFirst+1; b<=LastAge; b++)
     {
     for (s=1; s<=2; s++)
        InitN(b,s) = InitN(AgePlusFirst,s) * 
                     exp(-(b - AgePlusFirst) * OldZ(s));
     InitN(b,3) = InitN(b,1) + InitN(b,2);
     }
  // Jack up the plus group.
  for (s=1; s<=2; s++) InitN(LastAge,s) /= (1.0 - exp(-OldZ(s)));
  InitN(LastAge,3) = InitN(LastAge,1) + InitN(LastAge,2);
  //
  // Compute UnevenSexRatio_PSS.
  dvariable pf; // Prop. female at age.
  dvariable pflast = 0.5;
  UnevenSexRatio_PSS = 0.0;
  // Penalize deviations from 0.5 through age 8.
  // Penalize first differences thereafter.
  for (b=2; b<=AgePlusFirst; b++)
     {
     pf = InitN(b,1)/(InitN(b,1) + InitN(b,2));
     if (b <= 8)
        UnevenSexRatio_PSS += sq( (pf - 0.5)/UnevenSexRatio_SD );
        else
        UnevenSexRatio_PSS += sq( (pf - pflast)/UnevenSexRatio_SD );
     pflast = pf;
     }
  PSSTotal += UnevenSexRatio_PSS;
  //
  // Compute WildYearClass_PSS.
  dvariable nt;
  dvariable ntlast = InitN(2,3);
  WildYearClass_PSS = 0.0;
  for (b=3; b<=AgePlusFirst; b++)
     {
     nt = InitN(b,3);
     WildYearClass_PSS += sq( log(nt/ntlast)/WildYearClass_SD );
     }
  PSSTotal += WildYearClass_PSS;
  //

FUNCTION Fill_R
  //
  // Fill in the recruitment matrix R.
  // Set unseen recruitments to mean of seen ones.
  //
  // Scale the parameter.
  //
  for(y=FirstYear; y<=LastYearClass+1; y++) RSeen(y) = RSeenPar(y) * 1.0E6;
  //
  // Seen.
  WildR_PSS = 0.0;
  for (y=FirstYear; y<=LastYearClass+1; y++)
     {
     R(y,1) = PropFem * RSeen(y);
     R(y,2) = (1.0 - PropFem) * RSeen(y);
     R(y,3) = RSeen(y);
     if (y > FirstYear) WildR_PSS += sq( log(RSeen(y)/RSeen(y-1)) / WildR_SD );
     }
  PSSTotal += WildR_PSS;
  // Unseen.
  RMean = 0.5 * mean(RSeen);
  for (y=LastYearClass+2; y<=LastYear+1; y++)
     {
     for (s=1; s<=2; s++) R(y,s) = RMean;
     R(y,3) = 2.0 * RMean;
     }
  //

FUNCTION Initialize_N
  //
  // Drop InitN and R into N.
  //
  for (b=2; b<=LastAge; b++)
     for (s=1; s<=3; s++)
        N(FirstYear,b,s) = InitN(b,s);
  //
  for (y=FirstYear; y<=LastYear+1; y++)
     for (s=1; s<=3; s++)
        N(y,1,s) = R(y,s);
  //

FUNCTION  Compute_removals_and_survivors
  //
  // Run the numbers: compute removals and survivors for all
  // sources of mortality, years, ages, and sexes.
  //
  // These are the core calculations of the model.
  //
  for (y=FirstYear; y<=LastYear; y++)
  for (b=1; b<=LastAge; b++)
  for (s=1; s<=2; s++)
     {
     // Average number present during the year.
     NBar(y,b,s) = N(y,b,s) * (1 - exp(-Z(y,b,s))) / Z(y,b,s);
     //
     // All removals.
     CCatch(y,b,s) = Cf(y,b,s) * NBar(y,b,s);
     DCatch(y,b,s) = Df(y,b,s) * NBar(y,b,s);
     BCatch(y,b,s) = Bf(y,b,s) * NBar(y,b,s);
     RCatch(y,b,s) = Rf(y,b,s) * NBar(y,b,s);
     PCatch(y,b,s) = Pf(y,b,s) * NBar(y,b,s);
     MCatch(y,b,s) =  M(y,s)   * NBar(y,b,s);
     //
     // Survivors. Plus group survivors are added to survivors
     // from the next younger age.
     if (b < LastAge)
        N(y+1,b+1,s) = N(y,b,s) * exp(-Z(y,b,s)); else
        N(y+1,b,s) = N(y+1,b,s) + N(y,b,s) * exp(-Z(y,b,s));
     }
  //
  // Done with core calculations.
  //
  // Unload numbers and removals into sex-specific arrays.
  for (y=FirstYear; y<=LastYear+1; y++)
  for (b=1; b<=LastAge; b++)
     {
     N_F(y,b) = N(y,b,1);
     N_M(y,b) = N(y,b,2);
     N_T(y,b) = N_F(y,b) + N_M(y,b);
     }
     //
  for (y=FirstYear; y<=LastYear; y++)
  for (b=1; b<=LastAge; b++)
     {
     NBar_F(y,b) = NBar(y,b,1);
     NBar_M(y,b) = NBar(y,b,2);
     NBar_T(y,b) = NBar_F(y,b) + NBar_M(y,b);
     //
     CCatch_F(y,b) = CCatch(y,b,1);
     CCatch_M(y,b) = CCatch(y,b,2);
     CCatch_T(y,b) = CCatch_F(y,b) + CCatch_M(y,b);
     //
     BCatch_F(y,b) = BCatch(y,b,1);
     BCatch_M(y,b) = BCatch(y,b,2);
     BCatch_T(y,b) = BCatch_F(y,b) + BCatch_M(y,b);
     //
     RCatch_F(y,b) = RCatch(y,b,1);
     RCatch_M(y,b) = RCatch(y,b,2);
     RCatch_T(y,b) = RCatch_F(y,b) + RCatch_M(y,b);
     //
     PCatch_F(y,b) = PCatch(y,b,1);
     PCatch_M(y,b) = PCatch(y,b,2);
     PCatch_T(y,b) = PCatch_F(y,b) + RCatch_M(y,b);
     //
     MCatch_F(y,b) = MCatch(y,b,1);
     MCatch_M(y,b) = MCatch(y,b,2);
     MCatch_T(y,b) = MCatch_F(y,b) + RCatch_M(y,b);
     }
  //
  // Compute model CCPUE, SurvCPUE, and SurvProp at true age,
  // called CCPUE_True, SurvCPUE_True, and SurvProp_True.
  //
  dvariable SurvCPUETot_True;
  for (y=FirstYear; y<=LastYear; y++)
     {
     SurvCPUETot_True = 0.0;
     for (b=1; b<=LastAge; b++)
        {
        CCPUE_True_F(y,b) = CQ(y,1) * CSel(y,b,1) * NBar(y,b,1);
        CCPUE_True_M(y,b) = CQ(y,2) * CSel(y,b,2) * NBar(y,b,2);
        CCPUE_True_T(y,b) = CCPUE_True_F(y,b) + CCPUE_True_M(y,b);
        //
        SurvCPUE_True_F(y,b) = SQ(y,1) * SSel(y,b,1) * NBar(y,b,1);
        SurvCPUE_True_M(y,b) = SQ(y,2) * SSel(y,b,2) * NBar(y,b,2);
        SurvCPUE_True_T(y,b) = SurvCPUE_True_F(y,b) + SurvCPUE_True_M(y,b);
        //
        // Compute predicted total survey CPUE in number on first pass.
        SurvCPUETot_True += SurvCPUE_True_T(y,b);
        } // End first loop on age.
     // Compute SurvProp_True on second pass.
     for (b=1; b<=LastAge; b++)
        {
        SurvProp_True_F(y,b) = SurvCPUE_True_F(y,b)/SurvCPUETot_True;
        SurvProp_True_M(y,b) = SurvCPUE_True_M(y,b)/SurvCPUETot_True;
        SurvProp_True_T(y,b) = SurvCPUE_True_T(y,b)/SurvCPUETot_True;
        } // End second loop on age.
     } // End loop on year.
  //
  //

FUNCTION Compute_SelSmooth_PSS
  //
  // Smooth the length-specific selectivities.
  // (Survey and commercial; NOT bycatch.)
  //
  SelSmooth_PSS = 0.0;
  //
  for (s=1; s<=2; s++)
  {
  //
  // Commercial.
  for (w=1; w<=CSelLEstN; w++)
     {
     for (v=1; v<=SLptN; v++) SelVec(v) = CSelLEst(w,v,s);
     SelSmooth_PSS += pen2diffs(SelVec, SelSmooth_SD);
     } // End loop on waypoint.
  // Survey.
  for (w=1; w<=SSelLEstN; w++)
     {
     for (v=1; v<=SLptN; v++) SelVec(v) = SSelLEst(w,v,s);
     SelSmooth_PSS += pen2diffs(SelVec, SelSmooth_SD);
     } // End loop on waypoint.
  // Compute penalty once only if males not split.
  if (SplitMalesFlag == 0) break;
  } // End loop on sex.
  //
  PSSTotal += SelSmooth_PSS;


FUNCTION Compute_Catch_RSS
  //
  // Compute RSS for catch at age/sex.
  //
  // Predict catch at observed age/sex.
  //
  SmearAges(CCatch_F, Catch_Pred_F, TrueToSurf, TrueToBurn, NA);
  SmearAges(CCatch_M, Catch_Pred_M, TrueToSurf, TrueToBurn, NA);
  SmearAges(CCatch_T, Catch_Pred_T, TrueToSurf, TrueToBurn, NA);
  //
  // Compute RSS. Inclusion in Total RSS depends on MultiFitFlag.
  //
  Catch_RSS = 0.0;
  Catch_n = 0.0;
  RSS_list = CalcRRSS(Catch_Pred_F, Catch_F, Catch_SE_F,
                      CatchTau_F, Robust,
                      CMinAge, CLastYearFit, NA);
  Catch_RSS_F = RSS_list(1);
  Catch_n_F = RSS_list(2);
  Catch_rms_F = sqrt(RSS_list(3));
  Catch_RSS += RSS_list(1);
  Catch_n += RSS_list(2);
  RSS_list = CalcRRSS(Catch_Pred_M, Catch_M, Catch_SE_M,
                      CatchTau_M, Robust,
                      CMinAge, CLastYearFit, NA);
  Catch_RSS_M = RSS_list(1);
  Catch_n_M = RSS_list(2);
  Catch_rms_M = sqrt(RSS_list(3));
  Catch_RSS += RSS_list(1);
  Catch_n += RSS_list(2);
  RSS_list = CalcRRSS(Catch_Pred_T, Catch_T, Catch_SE_T,
                      CatchTau_T, Robust,
                      CMinAge, CLastYearFit, NA);
  Catch_RSS_T = RSS_list(1);
  Catch_n_T = RSS_list(2);
  Catch_rms_T = sqrt(RSS_list(3));
  if (MultiFitFlag==1) Catch_RSS += RSS_list(1);
  if (MultiFitFlag==1) Catch_n += RSS_list(2);
  //
  RSSTotal += CatchLambda * Catch_RSS;
  if (CatchLambda > 0.0) RSSNobs += Catch_n;
  //

FUNCTION Compute_CCPUE_RSS
  //
  // Compute RSS for CCPUE at age/sex.
  //
  // Predict CCPUE at observed age/sex.
  //
  SmearAges(CCPUE_True_F, CCPUE_Pred_F, TrueToSurf, TrueToBurn, NA);
  SmearAges(CCPUE_True_M, CCPUE_Pred_M, TrueToSurf, TrueToBurn, NA);
  SmearAges(CCPUE_True_T, CCPUE_Pred_T, TrueToSurf, TrueToBurn, NA);
  //
  // Compute RSS. Inclusion of sex-specific RSS depends on SexFlag.
  //
  if (MultiFitFlag==1) {
  CCPUE_RSS = 0.0;
  CCPUE_n = 0.0;
  RSS_list = CalcRRSS(CCPUE_Pred_F, CCPUE_F, CCPUE_SE_F,
                      CCPUETau_F, Robust,
                      CMinAge, CLastYearFit, NA);
  CCPUE_RSS_F = RSS_list(1);
  CCPUE_n_F = RSS_list(2);
  CCPUE_rms_F = sqrt(RSS_list(3));
  CCPUE_RSS += RSS_list(1);
  CCPUE_n += RSS_list(2);
  RSS_list = CalcRRSS(CCPUE_Pred_M, CCPUE_M, CCPUE_SE_M,
                      CCPUETau_M, Robust,
                      CMinAge, CLastYearFit, NA);
  CCPUE_RSS_M = RSS_list(1);
  CCPUE_n_M = RSS_list(2);
  CCPUE_rms_M = sqrt(RSS_list(3));
  CCPUE_RSS += RSS_list(1);
  CCPUE_n += RSS_list(2);
  RSS_list = CalcRRSS(CCPUE_Pred_T, CCPUE_T, CCPUE_SE_T,
                      CCPUETau_T, Robust,
                      CMinAge, CLastYearFit, NA);
  CCPUE_RSS_T = RSS_list(1);
  CCPUE_n_T = RSS_list(2);
  CCPUE_rms_T = sqrt(RSS_list(3));
  CCPUE_RSS += RSS_list(1);
  CCPUE_n += RSS_list(2);
  //
  RSSTotal += CCPUELambda * CCPUE_RSS;
  if (CCPUELambda > 0.0) RSSNobs += CCPUE_n;
  //
  }
FUNCTION Compute_SurvProp_RSS
  //
  // Compute RSS for SurvProp at age/sex.
  //
  // Predict SurvProp at observed age/sex.
  //
  SmearAges(SurvProp_True_F, SurvProp_Pred_F, TrueToSurf, TrueToBurn, NA);
  SmearAges(SurvProp_True_M, SurvProp_Pred_M, TrueToSurf, TrueToBurn, NA);
  SmearAges(SurvProp_True_T, SurvProp_Pred_T, TrueToSurf, TrueToBurn, NA);
  //
  // The rows of SurvProp_Pred_T will not quite sum to one because a small
  // fraction of young fish will be smeared to ages less than e.g. 6. Normalize.
  dvariable SumProp;
  for (y=FirstYear; y<=LastYear; y++)
     {
     SumProp = 0.0;
     // Compute sum for this year.
     for (a=SMinAge; a<=SMaxAge; a++)
        {
        if (SurvProp_Pred_T(y,a) == NA) break;
        SumProp += SurvProp_Pred_T(y,a);
        }
     // Divide all proportions by this sum.
     for (a=SMinAge; a<=SMaxAge; a++)
        {
        if (SurvProp_Pred_T(y,a) == NA) break;
        SurvProp_Pred_F(y,a) /= SumProp;
        SurvProp_Pred_M(y,a) /= SumProp;
        SurvProp_Pred_T(y,a) /= SumProp;
        }
     } // End loop on year.
  //
  // Compute RSS. Inclusion of sex-specific RSS depends on SexFlag.
  //
  SurvProp_RSS = 0.0;
  SurvProp_n = 0.0;
  RSS_list = CalcRRSS(SurvProp_Pred_F, SurvProp_F, SurvProp_SE_F,
                      SurvPropTau_F, Robust,
                      SMinAge, SLastYearFit, NA);
  SurvProp_RSS_F = RSS_list(1);
  SurvProp_n_F = RSS_list(2);
  SurvProp_rms_F = sqrt(RSS_list(3));
  SurvProp_RSS += RSS_list(1);
  SurvProp_n += RSS_list(2);
  RSS_list = CalcRRSS(SurvProp_Pred_M, SurvProp_M, SurvProp_SE_M,
                      SurvPropTau_M, Robust,
                      SMinAge, SLastYearFit, NA);
  SurvProp_RSS_M = RSS_list(1);
  SurvProp_n_M = RSS_list(2);
  SurvProp_rms_M = sqrt(RSS_list(3));
  SurvProp_RSS += RSS_list(1);
  SurvProp_n += RSS_list(2);
  RSS_list = CalcRRSS(SurvProp_Pred_T, SurvProp_T, SurvProp_SE_T,
                      SurvPropTau_T, Robust,
                      SMinAge, SLastYearFit, NA);
  SurvProp_RSS_T = RSS_list(1);
  SurvProp_n_T = RSS_list(2);
  SurvProp_rms_T = sqrt(RSS_list(3));
  if (MultiFitFlag==1) SurvProp_RSS += RSS_list(1);
  if (MultiFitFlag==1) SurvProp_n += RSS_list(2);
  //
  RSSTotal += SurvPropLambda * SurvProp_RSS;
  if (SurvPropLambda > 0.0) RSSNobs += SurvProp_n;
  //

FUNCTION Compute_SurvCPUE_RSS
  //
  // Compute RSS for SurvCPUE at age/sex.
  //
  // Predict SurvCPUE at observed age/sex.
  //
  SmearAges(SurvCPUE_True_F, SurvCPUE_Pred_F, TrueToSurf, TrueToBurn, NA);
  SmearAges(SurvCPUE_True_M, SurvCPUE_Pred_M, TrueToSurf, TrueToBurn, NA);
  SmearAges(SurvCPUE_True_T, SurvCPUE_Pred_T, TrueToSurf, TrueToBurn, NA);
  //
  // Compute RSS. Inclusion of sex-specific RSS depends on SexFlag.
  //
  if (MultiFitFlag==1) {
  SurvCPUE_RSS = 0.0;
  SurvCPUE_n = 0.0;
  RSS_list = CalcRRSS(SurvCPUE_Pred_F, SurvCPUE_F, SurvCPUE_SE_F,
                      SurvCPUETau_F, Robust,
                      SMinAge, SLastYearFit, NA);
  SurvCPUE_RSS_F = RSS_list(1);
  SurvCPUE_n_F = RSS_list(2);
  SurvCPUE_rms_F = sqrt(RSS_list(3));
  SurvCPUE_RSS += RSS_list(1);
  SurvCPUE_n += RSS_list(2);
  RSS_list = CalcRRSS(SurvCPUE_Pred_M, SurvCPUE_M, SurvCPUE_SE_M,
                      SurvCPUETau_M, Robust,
                      SMinAge, SLastYearFit, NA);
  SurvCPUE_RSS_M = RSS_list(1);
  SurvCPUE_n_M = RSS_list(2);
  SurvCPUE_rms_M = sqrt(RSS_list(3));
  SurvCPUE_RSS += RSS_list(1);
  SurvCPUE_n += RSS_list(2);
  RSS_list = CalcRRSS(SurvCPUE_Pred_T, SurvCPUE_T, SurvCPUE_SE_T,
                      SurvCPUETau_T, Robust,
                      SMinAge, SLastYearFit, NA);
  SurvCPUE_RSS_T = RSS_list(1);
  SurvCPUE_n_T = RSS_list(2);
  SurvCPUE_rms_T = sqrt(RSS_list(3));
  SurvCPUE_RSS += RSS_list(1);
  SurvCPUE_n += RSS_list(2);
  //
  RSSTotal += SurvCPUELambda * SurvCPUE_RSS;
  if (SurvCPUELambda > 0.0) RSSNobs += SurvCPUE_n;
  }
  //

FUNCTION Compute_CPUE_Totals_RSS
  //
  // Compute RSS for commercial and survey total CPUE's in no. & weight.
  //
  // Compute the predictions.
  //
  for (y=FirstYear; y<=LastYear; y++)
     {
     if (y<=2001) AgePlus = AgePlusSurf; else AgePlus = AgePlusBurn;
     CCPUETot_Pred(y) = 0.0;
     CWPUETot_Pred(y) = 0.0;
     //
     SurvCPUETot_Pred(y) = 0.0;
     SurvWPUETot_Pred(y) = 0.0;
     for (a=SMinAge; a<=AgePlus; a++)
        {
        CCPUETot_Pred(y) += CCPUE_Pred_F(y,a) + CCPUE_Pred_M(y,a);
        CWPUETot_Pred(y) += CCPUE_Pred_F(y,a) * CatchW_F(y,a) +
                            CCPUE_Pred_M(y,a) * CatchW_M(y,a);
        SurvCPUETot_Pred(y) += SurvCPUE_Pred_F(y,a) + SurvCPUE_Pred_M(y,a);
        SurvWPUETot_Pred(y) += SurvCPUE_Pred_F(y,a) * SurvPropLegal_F(y,a) * SurvLegalW_F(y,a) +
                               SurvCPUE_Pred_M(y,a) * SurvPropLegal_M(y,a) * SurvLegalW_M(y,a);
        } // End loop on age.
     }  // End loop on year.
  //
  // Compute the RSS's.
  //
  RSS_list = CalcRRSS(SurvWPUETot_Pred, SurvWPUETot, SurvWPUETot_SE,
                      SurvWPUETotTau, Robust,
                      FirstYear, SLastYearFit, NA);
  SurvWPUETot_RSS = RSS_list(1);
  SurvWPUETot_n = RSS_list(2);
  SurvWPUETot_rms = sqrt(RSS_list(3));
  //
  RSS_list = CalcRRSS(CWPUETot_Pred, CWPUETot, CWPUETot_SE,
                      CWPUETotTau, Robust,
                      FirstYear, CLastYearFit, NA);
  CWPUETot_RSS = RSS_list(1);
  CWPUETot_n = RSS_list(2);
  CWPUETot_rms = sqrt(RSS_list(3));
  //
  if (MultiFitFlag==1) {
  RSS_list = CalcRRSS(CCPUETot_Pred, CCPUETot, CCPUETot_SE,
                      CCPUETotTau, Robust,
                      FirstYear, CLastYearFit, NA);
  CCPUETot_RSS = RSS_list(1);
  CCPUETot_n = RSS_list(2);
  CCPUETot_rms = sqrt(RSS_list(3));
  //
  RSS_list = CalcRRSS(SurvCPUETot_Pred, SurvCPUETot, SurvCPUETot_SE,
                      SurvCPUETotTau, Robust,
                      FirstYear, SLastYearFit, NA);
  SurvCPUETot_RSS = RSS_list(1);
  SurvCPUETot_n = RSS_list(2);
  SurvCPUETot_rms = sqrt(RSS_list(3));
  }
  //
  RSSTotal += CWPUETotLambda * CWPUETot_RSS +
               SurvWPUETotLambda * SurvWPUETot_RSS;
  if (MultiFitFlag==1) RSSTotal+= CCPUETotLambda * CCPUETot_RSS + SurvCPUETotLambda * SurvCPUETot_RSS;
  if (CWPUETotLambda > 0.0) RSSNobs += CWPUETot_n;
  if (SurvWPUETotLambda > 0.0) RSSNobs += SurvWPUETot_n;
  if (MultiFitFlag==1 & SurvCPUETotLambda > 0.0) RSSNobs += SurvCPUETot_n;
  if (MultiFitFlag==1 & CCPUETotLambda > 0.0) RSSNobs += CCPUETot_n;

  //

FUNCTION Compute_DiscardWt_RSS
  //
  // This is really simple because the average weight of sublegals
  // is 7.5 lb every year.
  //
  // Compute the predictions.
  for (y=FirstYear; y<= LastYear; y++)
     {
     DiscardWt_Pred(y) = 0.0;
     for (b=1; b<=LastAge; b++)
     for (s=1; s<=2; s++)
        DiscardWt_Pred(y) += 7.5 * DCatch(y,b,s);
     }
  //
  // Compute RSS.
  //
  RSS_list = CalcRRSS(DiscardWt_Pred, DiscardWt, XCatchWt_SE,
                     1.0, 0, // No variance scaling or robustification.
                     FirstYear, LastYear, NA);
  DiscardWt_RSS = RSS_list(1);
  DiscardWt_n = RSS_list(2);
  DiscardWt_rms = sqrt(RSS_list(3));
  RSSTotal += DiscardWt_RSS;
  RSSNobs += DiscardWt_n;
  //

FUNCTION Compute_Bycatch_RSS
  //
  // Compute RSS for bycatch in number at length.
  //
  for (y=FirstYear; y<=LastYear; y++)
     {
     // Zero out predicted bycatch nos. by length.
     for (v=1; v<=BLptN; v++) Bycatch_Pred(y,v) = 0.0;
     // Run though ages apportioning female bycatch.
     for (b=1; b<=LastAge; b++)
     for (v=1; v<=BLptN; v++)
        {
        LBot = BLptList(v);
        zBot = (LBot - BycatchLTrue_F(y,b)) / BycatchLSDTrue_F(y,b);
        if (v < BLptN) LTop = BLptList(v+1);
                   else  LTop = 250.0;
        zTop = (LTop - BycatchLTrue_F(y,b)) / BycatchLSDTrue_F(y,b);
        if (zTop < -3.0) continue;
        if (zBot >  3.0) break;
        LPorp = cumd_norm(zTop) - cumd_norm(zBot);
        Bycatch_Pred(y,v) += LPorp * BCatch_F(y,b);
        } // End age, lint loops for females.
     // Run though ages apportioning male bycatch.
     for (b=1; b<=LastAge; b++)
     for (v=1; v<=BLptN; v++)
        {
        LBot = BLptList(v);
        zBot = (LBot - BycatchLTrue_M(y,b)) / BycatchLSDTrue_M(y,b);
        if (v < BLptN) LTop = BLptList(v+1);
                   else  LTop = 250.0;
        zTop = (LTop - BycatchLTrue_M(y,b)) / BycatchLSDTrue_M(y,b);
        if (zTop < -3.0) continue;
        if (zBot >  3.0) break;
        LPorp = cumd_norm(zTop) - cumd_norm(zBot);
        Bycatch_Pred(y,v) += LPorp * BCatch_M(y,b);
        } // End age, lint loops for males.
     } // End loop on year.
  //
  // Compute RSS.
  RSS_list = CalcRRSS(Bycatch_Pred, Bycatch, Bycatch_SE,
                      BycatchTau, Robust, 1, LastYear, NA);
  Bycatch_RSS = RSS_list(1);
  Bycatch_n = RSS_list(2);
  Bycatch_rms = sqrt(RSS_list(3));
  if (FitBycatchNosFlag != 0) 
     {
     RSSTotal += BycatchLambda * Bycatch_RSS;
     RSSNobs += Bycatch_n;
     }
  //

FUNCTION Compute_BycatchWt_RSS
  //
  // Compute RSS for predicted total bycatch weight.
  //
  // Use weight at true age to compute bycatch prediction.
  //
  for (y=FirstYear; y<=LastYear; y++)
     {
     BycatchWt_Pred(y)=0.0;
     for (b=1; b<=LastAge; b++)
        BycatchWt_Pred(y) += BCatch_F(y,b) * BycatchWTrue_F(y,b) +
                             BCatch_M(y,b) * BycatchWTrue_M(y,b);
     }
  //
  // Compute RSS.
  //
  RSS_list = CalcRRSS(BycatchWt_Pred, BycatchWt, XCatchWt_SE,
                     1.0, 0, // No variance scaling or robustification.
                     FirstYear, LastYear, NA);
  BycatchWt_RSS = RSS_list(1);
  BycatchWt_n = RSS_list(2);
  BycatchWt_rms = sqrt(RSS_list(3));
  if (FitBycatchWtFlag != 0)
	{
	RSSTotal += BycatchLambda * BycatchWt_RSS;
	RSSNobs += BycatchWt_n;
	}
  //

FUNCTION Compute_SportCatchWt_RSS
  //
  // Smear ages and use survey weight at observed age for sport & persuse.
  //
  SmearAges(RCatch_F, SportCatch_Pred_F, TrueToSurf, TrueToBurn, NA);
  SmearAges(RCatch_M, SportCatch_Pred_M, TrueToSurf, TrueToBurn, NA);
  //
  // Compute predictions.
  for (y=FirstYear; y<=LastYear; y++)
     {
     SportCatchWt_Pred(y) = 0.0;
     if (y<=2001) AgePlus = AgePlusSurf; else AgePlus = AgePlusBurn;
     for (a=SMinAge; a<=AgePlus; a++)
        SportCatchWt_Pred(y) += SportCatch_Pred_F(y,a) * SurvW_F(y,a) +
                                SportCatch_Pred_M(y,a) * SurvW_M(y,a);
     } // End year loop.
  //
  // Compute RSS.
  //
  RSS_list = CalcRRSS(SportCatchWt_Pred, SportCatchWt, XCatchWt_SE,
                      1.0, 0, // No variance scaling or robustification.
                      FirstYear, LastYear, NA);
  SportCatchWt_RSS = RSS_list(1);
  SportCatchWt_n = RSS_list(2);
  SportCatchWt_rms = sqrt(RSS_list(3));
  RSSTotal += SportCatchWt_RSS;
  RSSNobs += SportCatchWt_n;
  //

FUNCTION Compute_PersUseWt_RSS
  //
  // Smear ages and use survey weight at observed age for sport & persuse.
  //
  SmearAges(PCatch_F, PersUse_Pred_F, TrueToSurf, TrueToBurn, NA);
  SmearAges(PCatch_M, PersUse_Pred_M, TrueToSurf, TrueToBurn, NA);
  //
  // Compute predictions.
  for (y=FirstYear; y<=LastYear; y++)
     {
     PersUseWt_Pred(y) = 0.0;
     if (y<=2001) AgePlus = AgePlusSurf; else AgePlus = AgePlusBurn;
     for (a=SMinAge; a<=AgePlus; a++)
        PersUseWt_Pred(y) += PersUse_Pred_F(y,a) * SurvW_F(y,a) +
                                PersUse_Pred_M(y,a) * SurvW_M(y,a);
     } // End year loop.
  //
  // Compute RSS.
  //
  RSS_list = CalcRRSS(PersUseWt_Pred, PersUseWt, XCatchWt_SE,
                     1.0, 0, // No variance scaling or robustification.
                     FirstYear, LastYear, NA);
  PersUseWt_RSS = RSS_list(1);
  PersUseWt_n = RSS_list(2);
  PersUseWt_rms = sqrt(RSS_list(3));
  RSSTotal += PersUseWt_RSS;
  RSSNobs += PersUseWt_n;
  //

FUNCTION Compute_CatchN_RSS
  //
  // Compute RSS for total annual catch in number.
  //
  for (y=FirstYear; y<=LastYear; y++)
     {
     if (y<=2001) AgePlus = AgePlusSurf; else AgePlus = AgePlusBurn;
     CatchN_Pred(y) = 0.0;
     for (a=CMinAge; a<=AgePlus; a++) CatchN_Pred(y) += Catch_Pred_T(y,a);
     }
  RSS_list = CalcRRSS(CatchN_Pred, CatchN, XCatchN_SE, 
                      1.0, 0, // No variance scaling or robustification.
                      FirstYear, LastYear, NA);
  CatchN_RSS = RSS_list(1);
  CatchN_n = RSS_list(2);
  CatchN_rms = sqrt(RSS_list(3));
  if (FitCatchNFlag != 0)
	{
	 RSSTotal += CatchN_RSS;
	 RSSNobs += CatchN_n;
	}
  //
  
FUNCTION Compute_outputs
  //
  // Compute values to be reported.
  //
  if (MultiFitFlag==0)
  Deviance = RSSTotal/2.0; else
	  Deviance = RSSTotal/4.0;

  //
  // Unload 3D selectivity arrays to a matrix for each sex.
  for (y=FirstYear; y<=LastYear; y++)
  for (v=1; v<=SLptN; v++)
     {
     CSelL_F(y,v) = CSelL(y,v,1);
     CSelL_M(y,v) = CSelL(y,v,2);
     SSelL_F(y,v) = SSelL(y,v,1);
     SSelL_M(y,v) = SSelL(y,v,2);
     }
  for (y=FirstYear; y<=LastYear; y++)
  for (v=1; v<=BLptN; v++)
     {
     BSelL_F(y,v) = BSelL(y,v,1);
     BSelL_M(y,v) = BSelL(y,v,2);
     }
  //
  for (y=FirstYear; y<=LastYear; y++)
  for (b=1; b<=LastAge; b++)
     {
     CSel_F(y,b) = CSel(y,b,1);
     CSel_M(y,b) = CSel(y,b,2);
     DSel_F(y,b) = DSel(y,b,1);
     DSel_M(y,b) = DSel(y,b,2);
     FSel_F(y,b) = FSel(y,b,1);
     FSel_M(y,b) = FSel(y,b,2);
     SSel_F(y,b) = SSel(y,b,1);
     SSel_M(y,b) = SSel(y,b,2);
     BSel_F(y,b) = BSel(y,b,1);
     BSel_M(y,b) = BSel(y,b,2);
     Z_F(y,b)    = Z(y,b,1);
     Z_M(y,b)    = Z(y,b,2);
     }
  //
  // Compute commercial and survey exploitable numbers at true age.
  int yp; // Parameter year yp=y except in LastYear+1 when it is LastYear.
  for (y=FirstYear; y<=LastYear+1; y++)
     {
     if (y<=LastYear) yp = y; else yp = LastYear;
     for (b=1; b<=LastAge; b++)
        {
        CEN_F(y,b) = CSel_F(yp,b) * N_F(y,b);
        CEN_M(y,b) = CSel_M(yp,b) * N_M(y,b);
        FEN_F(y,b) = FSel_F(yp,b) * N_F(y,b);
        FEN_M(y,b) = FSel_M(yp,b) * N_M(y,b);
        SEN_F(y,b) = SSel_F(yp,b) * N_F(y,b);
        SEN_M(y,b) = SSel_M(yp,b) * N_M(y,b);
        }
     }
  //
  // Compute the usual measures of estimated abundance using the smeared
  // population at age/sex NS and observed weight at age.
  //
  SmearAges(N_F, NS_F, TrueToSurf, TrueToBurn, NA);
  SmearAges(N_M, NS_M, TrueToSurf, TrueToBurn, NA);
  //
  SmearAges(CEN_F, CENS_F, TrueToSurf, TrueToBurn, NA);
  SmearAges(CEN_M, CENS_M, TrueToSurf, TrueToBurn, NA);
  //
  SmearAges(FEN_F, FENS_F, TrueToSurf, TrueToBurn, NA);
  SmearAges(FEN_M, FENS_M, TrueToSurf, TrueToBurn, NA);
  //
  SmearAges(SEN_F, SENS_F, TrueToSurf, TrueToBurn, NA);
  SmearAges(SEN_M, SENS_M, TrueToSurf, TrueToBurn, NA);
  //
  int yd; // Data year yd=year except in YearLastDat+1 it is YearLastDat.
  for (y=FirstYear; y<=LastYear+1; y++)
     {
     R6(y) = N_F(y,6) + N_M(y,6);
     R8(y) = N_F(y,8) + N_M(y,8);
     R8SD(y) = R8(y);
     //
     // Initialize accumulators.
     N8Plus_F(y) = 0.0;
     N8Plus_M(y) = 0.0;
     TBio8Plus_F(y) = 0.0;
     TBio8Plus_M(y) = 0.0;
     LBio(y) = 0.0;
     SBio(y) = 0.0;
     EBio(y) = 0.0;
     EBioFixed(y) = 0.0;
     EBioSurv(y) = 0.0;
     //
     // Run through observed ages.
     if (y<=2001) AgePlus = AgePlusSurf; else AgePlus = AgePlusBurn;
     if (y<=YearLastDat) yd = y; else yd = YearLastDat;
     for (a=SMinAge; a<=AgePlus; a++)
        {
        LBio(y) += NS_F(y,a) * SurvPropLegal_F(yd,a) * SurvLegalW_F(yd,a) +
                   NS_M(y,a) * SurvPropLegal_M(yd,a) * SurvLegalW_M(yd,a);
        SBio(y) += PropMat(a) * NS_F(y,a) * SurvW_F(yd,a);
        EBio(y) += CENS_F(y,a) * CatchW_F(yd,a) + CENS_M(y,a) * CatchW_M(yd,a);
        EBioFixed(y) += FENS_F(y,a) * CatchW_F(yd,a) + FENS_M(y,a) * CatchW_M(yd,a);
        EBioSurv(y) += SENS_F(y,a) * SurvW_F(yd,a) + SENS_M(y,a) * SurvW_M(yd,a);
        //
        if (a < 8) continue;  // Other measures are ages 8+.
        N8Plus_F(y) += NS_F(y,a);
        N8Plus_M(y) += NS_M(y,a);
        TBio8Plus_F(y) += NS_F(y,a) * SurvW_F(yd,a);
        TBio8Plus_M(y) += NS_M(y,a) * SurvW_M(yd,a);
        } // End loop on observed age.
     }  // End loop on year.
  LastEBioFixed = EBioFixed(LastYear+1);
  LastEBioFixedToo = EBioFixed(LastYear+1);
  LastSBio = SBio(LastYear+1);
  //

REPORT_SECTION
  //    
  // Report values to be reported.
  //
  Compute_outputs();
  ofstream rep("trouble.rep");
  //
  rep << "Area" << endl << Area << endl;
  rep << "SQ_bhat" << endl << SQ_bhat << endl;
  rep << "FirstYear" << endl << FirstYear << endl;
  rep << "LastYear" << endl << LastYear << endl;
  rep << "LastYearClass" << endl << LastYearClass << endl;
  rep << "MinusLogL" << endl << MinusLogL << endl;
  rep << "Deviance" << endl << Deviance << endl;
  rep << "RSSNobs" << endl << RSSNobs << endl;
  rep << "CatchWt" << endl << CatchWt/1.0E6 << endl;
  rep << "MPar" << endl << MPar << " " << MPar_M << endl;
  rep << "PropFem" << endl << PropFem << endl;
  rep << "R6" << endl << R6/1.0E6 << endl;
  rep << "R8" << endl << R8/1.0E6 << endl;
  rep << "EBio" << endl << EBio/1.0E6 << endl;
  rep << "EBioFixed" << endl << EBioFixed/1.0E6 << endl;
  rep << "EBioSurv" << endl << EBioSurv/1.0E6 << endl;
  rep << "SBio" << endl << SBio/1.0E6 << endl;
  rep << "LBio" << endl << LBio/1.0E6 << endl;
  rep << "N8Plus.F" << endl << N8Plus_F/1.0E6 << endl;
  rep << "N8Plus.M" << endl << N8Plus_M/1.0E6 << endl;
  rep << "TBio8Plus.F" << endl << TBio8Plus_F/1.0E6 << endl;
  rep << "TBio8Plus.M" << endl << TBio8Plus_M/1.0E6 << endl;
  rep << "M" << endl << M << endl;
  rep << "CF" << endl << CF << endl;
  rep << "Z.F" << endl << Z_F << endl;
  rep << "Z.M" << endl << Z_M << endl;
  rep << "CQMX" << endl << CQMX << endl;
  rep << "CQ" << endl << CQ << endl;
  rep << "CSelL.F" << endl << CSelL_F << endl;
  rep << "CSelL.M" << endl << CSelL_M << endl;
  rep << "CSel.F" << endl << CSel_F << endl;
  rep << "CSel.M" << endl << CSel_M << endl;
  rep << "SQMX" << endl << SQMX << endl;
  rep << "SQ" << endl << SQ << endl;
  rep << "SSelL.F" << endl << SSelL_F << endl;
  rep << "SSelL.M" << endl << SSelL_M << endl;
  rep << "SSel.F" << endl << SSel_F << endl;
  rep << "SSel.M" << endl << SSel_M << endl;
  rep << "DFPar" << endl << DFPar << endl;
  rep << "RFPar" << endl << RFPar << endl;
  rep << "PFPar" << endl << PFPar << endl;
  rep << "DSelL.F" << endl << DSel_F << endl;
  rep << "BF" << endl << BF << endl;
  rep << "BSelL.F" << endl << BSelL_F << endl;
  rep << "BSelL.M" << endl << BSelL_M << endl;
  rep << "BSel.F" << endl << BSel_F << endl;
  rep << "BSel.M" << endl << BSel_M << endl;
  rep << "FSel.F:" << endl << FSel_F << endl;
  rep << "FSel.M:" << endl << FSel_M << endl;
  rep << "R" << endl << R << endl;
  rep << "InitN" << endl << InitN << endl;
  //
  rep << "N.F" << endl << N_F << endl;
  rep << "N.M" << endl << N_M << endl;
  rep << "NS.F" << endl << NS_F << endl;
  rep << "NS.M" << endl << NS_M << endl;
  rep << "CENS.F" << endl << CENS_F << endl;
  rep << "CENS.M" << endl << CENS_M << endl;
  rep << "FENS.F" << endl << FENS_F << endl;
  rep << "FENS.M" << endl << FENS_M << endl;
  rep << "CatchW.F" << endl << CatchW_F << endl;
  rep << "CatchW.M" << endl << CatchW_M << endl;
  rep << "Catch.F" << endl << Catch_F << endl;
  rep << "Catch.Pred.F" << endl << Catch_Pred_F << endl;
  rep << "Catch.SE.F" << endl << Catch_SE_F << endl;
  rep << "Catch.M" << endl << Catch_M << endl;
  rep << "Catch.Pred.M" << endl << Catch_Pred_M << endl;
  rep << "Catch.SE.M" << endl << Catch_SE_M << endl;
  rep << "Catch.T" << endl << Catch_T << endl;
  rep << "Catch.Pred.T" << endl << Catch_Pred_T << endl;
  rep << "Catch.SE.T" << endl << Catch_SE_T << endl;
  rep << "CCPUE.F" << endl << CCPUE_F << endl;
  rep << "CCPUE.Pred.F" << endl << CCPUE_Pred_F << endl;
  rep << "CCPUE.SE.F" << endl << CCPUE_SE_F << endl;
  rep << "CCPUE.M" << endl << CCPUE_M << endl;
  rep << "CCPUE.Pred.M" << endl << CCPUE_Pred_M << endl;
  rep << "CCPUE.SE.M" << endl << CCPUE_SE_M << endl;
  rep << "CCPUE.T" << endl << CCPUE_T << endl;
  rep << "CCPUE.Pred.T" << endl << CCPUE_Pred_T << endl;
  rep << "CCPUE.SE.T" << endl << CCPUE_SE_T << endl;
  rep << "CCPUETot" << endl << CCPUETot << endl;
  rep << "CCPUETot.Pred" << endl << CCPUETot_Pred << endl;
  rep << "CCPUETot.SE" << endl << CCPUETot_SE << endl;
  rep << "CWPUETot" << endl << CWPUETot << endl;
  rep << "CWPUETot.Pred" << endl << CWPUETot_Pred << endl;
  rep << "CWPUETot.SE" << endl << CWPUETot_SE << endl;
  rep << "SurvProp.F" << endl << SurvProp_F << endl;
  rep << "SurvProp.Pred.F" << endl << SurvProp_Pred_F << endl;
  rep << "SurvProp.SE.F" << endl << SurvProp_SE_F << endl;
  rep << "SurvProp.M" << endl << SurvProp_M << endl;
  rep << "SurvProp.Pred.M" << endl << SurvProp_Pred_M << endl;
  rep << "SurvProp.SE.M" << endl << SurvProp_SE_M << endl;
  rep << "SurvProp.T" << endl << SurvProp_T << endl;
  rep << "SurvProp.Pred.T" << endl << SurvProp_Pred_T << endl;
  rep << "SurvProp.SE.T" << endl << SurvProp_SE_T << endl;
  rep << "SurvCPUE.True.F" << endl << SurvCPUE_True_F << endl;
  rep << "SurvCPUE.F" << endl << SurvCPUE_F << endl;
  rep << "SurvCPUE.Pred.F" << endl << SurvCPUE_Pred_F << endl;
  rep << "SurvCPUE.SE.F" << endl << SurvCPUE_SE_F << endl;
  rep << "SurvCPUE.True.M" << endl << SurvCPUE_True_M << endl;
  rep << "SurvCPUE.M" << endl << SurvCPUE_M << endl;
  rep << "SurvCPUE.Pred.M" << endl << SurvCPUE_Pred_M << endl;
  rep << "SurvCPUE.SE.M" << endl << SurvCPUE_SE_M << endl;
  rep << "SurvCPUE.True.T" << endl << SurvCPUE_True_T << endl;
  rep << "SurvCPUE.T" << endl << SurvCPUE_T << endl;
  rep << "SurvCPUE.Pred.T" << endl << SurvCPUE_Pred_T << endl;
  rep << "SurvCPUE.SE.T" << endl << SurvCPUE_SE_T << endl;
  rep << "SurvCPUETot" << endl << SurvCPUETot << endl;
  rep << "SurvCPUETot.Pred" << endl << SurvCPUETot_Pred << endl;
  rep << "SurvCPUETot.SE" << endl << SurvCPUETot_SE << endl;
  rep << "SurvWPUETot" << endl << SurvWPUETot << endl;
  rep << "SurvWPUETot.Pred" << endl << SurvWPUETot_Pred << endl;
  rep << "SurvWPUETot.SE" << endl << SurvWPUETot_SE << endl;
  rep << "DiscardWt" << endl << DiscardWt/1.0E6 << endl;
  rep << "DiscardWt.Pred" << endl << DiscardWt_Pred/1.0E6 << endl;
  rep << "Bycatch" << endl << Bycatch << endl;
  rep << "Bycatch.SE" << endl << Bycatch_SE << endl;
  rep << "Bycatch.Pred" << endl << Bycatch_Pred << endl;
  rep << "BycatchWt" << endl << BycatchWt/1.0E6 << endl;
  rep << "BycatchWt.Pred" << endl << BycatchWt_Pred/1.0E6 << endl;
  rep << "SportCatchWt" << endl << SportCatchWt/1.0E6 << endl;
  rep << "SportCatchWt.Pred" << endl << SportCatchWt_Pred/1.0E6 << endl;
  rep << "PersUseWt" << endl << PersUseWt/1.0E6 << endl;
  rep << "PersUseWt.Pred" << endl << PersUseWt_Pred/1.0E6 << endl;
  rep << "CatchN" << endl << CatchN/1.0E6 << endl;
  rep << "CatchN.Pred" << endl << CatchN_Pred/1.0E6 << endl;
  //
  ofstream ss("trouble.ss");
  //
  ss << "SSTotal " <<  SSTotal << endl;
  ss << endl;
  ss << "RSSTotal " << RSSTotal << endl;
  ss << "RSSNobs " << RSSNobs << endl;
  ss << "CatchLambda " << CatchLambda << endl;
  ss << "Catch_RSS_F & _M & _T " << Catch_RSS_F << " " << Catch_RSS_M << " " << Catch_RSS_T << endl;
  ss << "Catch_n_F & _M & _T " << Catch_n_F << " " << Catch_n_M << " " << Catch_n_T << endl;
  ss << "Catch_rms_F & _M & _T " << Catch_rms_F << " " << Catch_rms_M << " " << Catch_rms_T << endl;
  ss << "CCPUELambda " << CCPUELambda << endl;
  ss << "CCPUE_RSS_F & _M & _T " << CCPUE_RSS_F << " " << CCPUE_RSS_M << " " << CCPUE_RSS_T << endl;
  ss << "CCPUE_n_F & _M & _T " << CCPUE_n_F << " " << CCPUE_n_M << " " << CCPUE_n_T << endl;
  ss << "CCPUE_rms_F & _M & _T " << CCPUE_rms_F << " " << CCPUE_rms_M << " " << CCPUE_rms_T << endl;
  ss << "CCPUETotLambda " << CCPUETotLambda << endl;
  ss << "CCPUETot_RSS & _n & _rms " << CCPUETot_RSS << " " << CCPUETot_n << " " << CCPUETot_rms << endl;
  ss << "CWPUETotLambda " << CWPUETotLambda << endl;
  ss << "CWPUETot_RSS & _n & _rms " << CWPUETot_RSS << " " << CWPUETot_n << " " << CWPUETot_rms << endl;
  ss << "SurvCPUELambda " << SurvCPUELambda << endl;
  ss << "SurvCPUE_RSS_F & _M & _T " << SurvCPUE_RSS_F << " " << SurvCPUE_RSS_M << " " << SurvCPUE_RSS_T << endl;
  ss << "SurvCPUE_n_F & _M & _T " << SurvCPUE_n_F << " " << SurvCPUE_n_M << " " << SurvCPUE_n_T << endl;
  ss << "SurvCPUE_rms_F & _M & _T " << SurvCPUE_rms_F << " " << SurvCPUE_rms_M << " " << SurvCPUE_rms_T << endl;
  ss << "SurvCPUETotLambda " << SurvCPUETotLambda << endl;
  ss << "SurvCPUETot_RSS & _n & _rms " << SurvCPUETot_RSS << " " << SurvCPUETot_n << " " << SurvCPUETot_rms << endl;
  ss << "SurvWPUETotLambda " << SurvWPUETotLambda << endl;
  ss << "SurvWPUETot_RSS & _n & _rms " << SurvWPUETot_RSS << " " << SurvWPUETot_n << " " << SurvWPUETot_rms << endl;
  ss << "SurvPropLambda " << SurvPropLambda << endl;
  ss << "SurvProp_RSS_F & _M & _T " << SurvProp_RSS_F << " " << SurvProp_RSS_M << " " << SurvProp_RSS_T << endl;
  ss << "SurvProp_n_F & _M & _T " << SurvProp_n_F << " " << SurvProp_n_M << " " << SurvProp_n_T << endl;
  ss << "SurvProp_rms_F & _M & _T " << SurvProp_rms_F << " " << SurvProp_rms_M << " " << SurvProp_rms_T << endl;
  ss << "Bycatch_RSS " << BycatchLambda << endl;
  ss << "Bycatch_RSS  & _n & _rms " << Bycatch_RSS << " " << Bycatch_n << " " << Bycatch_rms << endl;
  ss << "BycatchWt_RSS " << BycatchLambda << endl;
  ss << "BycatchWt_RSS  & _n & _rms " << BycatchWt_RSS << " " << BycatchWt_n << " " << BycatchWt_rms << endl;
  ss << "DiscardWt_RSS  & _n & _rms " << DiscardWt_RSS << " " << DiscardWt_n << " " << DiscardWt_rms << endl;
  ss << "SportCatchWt_RSS  & _n & _rms " << SportCatchWt_RSS << " " << SportCatchWt_n << " " << SportCatchWt_rms << endl;
  ss << "PersUseWt_RSS  & _n & _rms " << PersUseWt_RSS << " " << PersUseWt_n << " " << PersUseWt_rms << endl;
  ss << "CatchN_RSS  & _n & _rms " << CatchN_RSS << " " << CatchN_n << " " << CatchN_rms << endl;
  ss << endl;
  ss << "PSSTotal " << PSSTotal << endl;
  ss << "SteadyCQ_PSS " << SteadyCQ_PSS << endl;
  ss << "SteadySQ_PSS " << SteadySQ_PSS << endl;
  ss << "TrendlessSQ_PSS " << TrendlessSQ_PSS << endl;
  ss << "SelSmooth_PSS " << SelSmooth_PSS << endl;
  ss << "UnevenSexRatio_PSS " << UnevenSexRatio_PSS << endl;
  ss << "WildYearClass_PSS " << WildYearClass_PSS << endl;
  ss << "WildR_PSS " << WildR_PSS << endl;
  ss << "LastF_PSS " << LastF_PSS << endl;
  //

RUNTIME_SECTION
  //
  maximum_function_evaluations 10000
  convergence_criteria 1.0E-4
  //
