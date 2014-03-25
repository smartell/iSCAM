/// @file iscam.tpl
/// @author Steve Martell, IPHC

///
/// \def REPORT(object)
/// \brief Prints name and value of \a object on ADMB report %ofstream file.
///

///
/// \def COUT(object)
/// \brief Screen dump using cout<<"object\n"<<object<<endl;
///

///
/// \def TINY
/// \brief A small number (1.e-08)`
///

// ----------------------------------------------------------------------------- //
//         integrated Statistical Catch Age Model (iSCAM)                        //
//                                                                               //
//                           VERSION 1.1                                         //
//               Tue Jul 19 22:23:58 PDT 2011                                    //
//                                                                               //
//                                                                               //
//           Created by Steven Martell on 2010-04-09                             //
//           Copyright (c) 2010. All rights reserved.                            //
//                                                                               //
// AUTHORS: SJDM Steven Martell                                                  //
//                                                                               //
// CONVENTIONS: Formatting conventions are based on the The                      //
//               Elements of C++ Style (Misfeldt et al. 2004)                    //
//                                                                               //
// NAMING CONVENTIONS:                                                           //
//             Macros       -> UPPERCASE                                         //
//             Constants    -> UpperCamelCase                                    //
//             Functions    -> lowerCamelCase                                    //
//             Variables    -> lowercase                                         //
//                                                                               //
// CHANGED add option for using empirical weight-at-age data                     //
// TODO:   ? add gtg options for length based fisheries                          //
// CHANGED add time varying natural mortality rate with splines                  //
// TODO:   ? add cubic spline interpolation for time varying M                   //
// CHANGED  Fix the type 6 selectivity implementation. not working.              //
// TODO:  fix cubic spline selectivity for only years when data avail            //
// CHANGED: fixed a bug in the simulation model log_ft_pars goes out             //
//        of bounds.                                                             //
// TODO: write a projection routine and verify equilibrium calcs                 //
// TODO: add DIC calculation for MCMC routines (in -mcveal phase)                //
// CHANGED: add SOK fishery a) egg fishing mort 2) bycatch for closed ponds      //
//                                                                               //
// TODO: correct the likelihood function as per Roberto's email. See Schnute     //
//       and Richards (1995) paper. Table 4 when 0<rho<1                         //
//                                                                               //
// TODO: Check the recruitment autocorrelation residuals & likelihood.           //
//                                                                               //
//                                                                               //
// ----------------------------------------------------------------------------- //
//-- CHANGE LOG:                                                               --//
//--  Nov 30, 2010 -modified survey biomass by the fraction of total           --//
//--                mortality that occurred during the time of the             --//
//--                survey. User specifies this fraction (0-1) in the          --//
//--                data file as the last column of the relative               --//
//--                abundance index.                                           --//
//--                                                                           --//
//--  Dec 6, 2010 -modified the code to allow for empiracle weight-            --//
//--               at-age data to be used.                                     --//
//--              -rescaled catch and relative abundance /1000, this           --//
//--               should be done in the data file and not here.               --//
//--                                                                           --//
//--  Dec 20, 2010-added prior to survey q's in control file                   --//
//--                                                                           --//
//--  Dec 24, 2010-added random walk for natural mortality.                    --//
//--                                                                           --//
//--  Jan 23, 2011-in Penticton Hospital with my mom in ICU, adopting          --//
//--               the naming conventions in The Elements of C++               --//
//--               style to keep my mind busy.                                 --//
//--                                                                           --//
//-- May 5, 2011- added logistic selectcitivty as a fucntion of                --//
//--              mean body weight.  3 parameter logistic.                     --//
//--              NOT WORKING YET                                              --//
//--                                                                           --//
//-- May 6, 2011- added pre-processor commands to determin PLATFORM            --//
//--              either "Windows" or "Linux"                                  --//
//--            - change April 10, 2013 to #if defined _WIN32 etc.             --//
//--                                                                           --//
//--                                                                           --//
//-- use -mcmult 1.5 for MCMC with log_m_nodes with SOG herrning               --//
//--                                                                           --//
//--                                                                           --//
//-- Dec 11, 2011- added halibut branch to local git repository aim is to      --//
//--               add gender dimension and stock dimension.                   --//
//--               This was created on the "twosex" branch in git merged       --//
//--                                                                           --//
//-- Dec 30, 2011- working on length-based selectivity for halibut.            --//
//--                                                                           --//
//-- Jan 5, 2012 - adding spawn on kelp fishery as catch_type ivector          --//
//--             - modified the following routines:                            --//
//--             - calcCatchAtAge                                              --//
//--             - calcTotalMortality                                          --//
//--                                                                           --//
//-- Oct 31,2012 - added penalty to time-varying changes in selex for          --//
//--             - isel_type 4 and 5 cases in the objective function.          --//
//--             - Requires and additional input in the control file.          --//
//--                                                                           --//
//--                                                                           --//
//--                                                                           --//
//-- TODO: add catch_type to equilibrium calculations for reference points     --//
//--                                                                           --//
//-- Feb 18, 2013 - Need to redesign the simulation selectivities.             --//
//--              - Should probably use a separate simulation control file.    --//
//--                                                                           --//
//-- April 16, - Created new IPHC branch for developing sex/area/group         --//
//--           - INDEXS:                                                       --//
//--             area     f                                                    --//
//--             group    g                                                    --//
//--             sex      h                                                    --//
//--             year     i                                                    --//
//--             age      j                                                    --//
//--             gear     k                                                    --//
//--                                                                           --//
//--  ToDo: add mirror selectivity option so one fishery can assume the same   --//
//--        Selectivity value as another fishery with informative data.        --//
//--                                                                           --//
//--                                                                           --//
//--                                                                           --//
//--                                                                           --//
// ----------------------------------------------------------------------------- //



DATA_SECTION
	// |---------------------------------------------------------------------------------|
	// | STRINGS FOR INPUT FILES                                                         |
	// |---------------------------------------------------------------------------------|
	/// | DataFile.dat           : data to condition the assessment model on     
	init_adstring DataFile;      ///< String for the input datafile name.
	/// | ControlFile.ctl        : controls for phases, selectivity options 
	init_adstring ControlFile;	 ///< String for the control file.
	/// | ProjectFileControl.pfc : used for stock projections under TAC
	init_adstring ProjectFileControl;  ///< String for the projection file.
	/// | BaseFileName           : file prefix used for all iSCAM model output
	!! BaseFileName = stripExtension(ControlFile);  ///< BaseName given by the control file
	/// | ReportFileName         : file name to copy report file to.
	!! ReportFileName = BaseFileName + adstring(".rep");
	!! cout<<BaseFileName<<endl;
	
	
	// |---------------------------------------------------------------------------------|
	// | READ IN PROJECTION FILE CONTROLS                                         
	// |---------------------------------------------------------------------------------|
	// | n_tac    : length of catch vector for decision table catch stream.
	// | tac      : vector of total catch values to be used in the decision table.
	// | pf_cntrl : vector of controls for the projections.
	// | Documentation for projection control file pf_cntrl
	// | 1) start year for m_bar calculation
	// | 2)   end year for m_bar calculation
	// | 3) start year for average fecundity/weight-at-age
	// | 4)   end year for average fecundity/weight-at-age
	// | 5) start year for recruitment period (not implemented yet)
	// | 6)   end year for recruitment period (not implemented yet)
	// |
	!! ad_comm::change_datafile_name(ProjectFileControl);
	/// | Number of catch options to explore in the decision table.
	init_int n_tac; ///< Number of catch options to explore in the decision table.
					///<
	!! COUT(n_tac);
	/// | Vector of catch options.
	init_vector tac(1,n_tac);
	init_int n_pfcntrl;
	init_vector pf_cntrl(1,n_pfcntrl);

	//init_vector mse_cntrl(1,1);

	init_int eof_pf;
	LOC_CALCS
		if(eof_pf!=-999)
		{
			cout<<"Error reading projection file."<<endl;
			cout<<"Last integer read is "<<eof_pf<<endl;
			cout<<"The file should end with -999.\n Aborting!"<<endl;
			exit(1);
		}
	END_CALCS
	
	
	// |---------------------------------------------------------------------------------|
	// | COMMAND LINE ARGUMENTS FOR SIMULATION & RETROSPECTIVE ANALYSIS
	// |---------------------------------------------------------------------------------|
	// | SimFlag    : if user specifies -sim, then turn SimFlag on.
	// | retro_yrs  : number of terminal years to remove.
	
	int SimFlag;  ///< Flag for simulation mode
	int mseFlag;  ///< Flag for management strategy evaluation mode
	int rseed;    ///< Random number seed for simulated data.
	int retro_yrs;///< Number of years to look back from terminal year.
	int NewFiles;
	int testMSY;
	LOC_CALCS
		SimFlag=0;
		rseed=999;
		int on,opt;
		//the following line checks for the "-SimFlag" command line option
		//if it exists the if statement retrieves the random number seed
		//that is required for the simulation model
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			SimFlag = 1;
			rseed   = atoi(ad_comm::argv[on+1]);
		}
		
		// Catarina implementing a new command for generating new data control and pfc file
		// for a new project.
		NewFiles = 0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-new",opt))>-1)
		{
			NewFiles = 1;
			NewFileName = ad_comm::argv[on+1];
		}


		// command line option for retrospective analysis. "-retro retro_yrs"
		retro_yrs=0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-retro",opt))>-1)
		{
			retro_yrs = atoi(ad_comm::argv[on+1]);
			cout<<"|____________________________________________________|\n";
			cout<<"| Implementing Retrospective analysis                |\n";
			cout<<"|____________________________________________________|\n";
			cout<<"| Number of retrospective years = "<<retro_yrs<<endl;
		}

		// Management strategy evaluation.
		mseFlag = 0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-mse",opt))>-1)
		{
			mseFlag = 1;
			rseed   = atoi(ad_comm::argv[on+1]);
			cout<<"|_________________________________________________|\n";
			cout<<"|Implementing Management Strategy Evaluation      |\n";
			cout<<"|_________________________________________________|\n";
		}

		// Test MSY
		testMSY = 0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-msy",opt))>-1)
		{
			cout<<"Testing MSY calculations with Spreadsheet MSF.xlsx"<<endl;
			testMSY = 1;
			
		}

	END_CALCS


	// |---------------------------------------------------------------------------------|
	// | MODEL DATA FROM DATA FILE
	// |---------------------------------------------------------------------------------|
	// |
	!! ad_comm::change_datafile_name(DataFile);

	// |---------------------------------------------------------------------------------|
	// | MODEL DIMENSIONS
	// |---------------------------------------------------------------------------------|
	// |
	// | area   f
	// | group  g
	// | sex    h
	// | year   i
	// | age    j
	// | gear   k  - number of gears with unique selectivity
	int f;
	int g; 
	int h;
	int i;
	int j;
	int k;

	init_int narea;			
	init_int ngroup;			
	init_int nsex;			
	init_int syr;			
	init_int nyr;				
	init_int sage;			
	init_int nage;			
	init_int ngear;				
	vector age(sage,nage);			

	// |---------------------------------------------------------------------------------|
	// | LINKS TO MANAGE ARRAY INDEXING
	// |---------------------------------------------------------------------------------|
	// | - n_ags: total number of areas * groups * sex
	// | - n_ag:  total number of areas * groups
	// | - n_gs:  total number of groups * sex
	// | - n_area:  vector of indexes for area for each sex & group combination.
	// | - n_group: vector of indexes for stock for each sex & area combination.
	// | - n_sex:   vector of indexes for sex foe each area & group combination.
	// | - pntr_ag: matrix of indices for area and group.
	// | - pntr_gs: matrix of indices for group and sex.
	// | - pntr_ags: d3_array of indices for area group sex.
	// |	
	int n_ags;
	!! n_ags = narea * ngroup * nsex;
	int n_ag;
	!! n_ag  = narea * ngroup;
	int n_gs;
	!! n_gs  = ngroup * nsex;
	ivector   n_area(1,n_ags);
	ivector  n_group(1,n_ags);
	ivector    n_sex(1,n_ags);
	imatrix  pntr_ag(1,narea,1,ngroup);
	imatrix  pntr_gs(1,ngroup,1,nsex);
	3darray pntr_ags(1,narea,1,ngroup,1,nsex);
	
	
	
	LOC_CALCS
		age.fill_seqadd(sage,1);
		int ig,ih,is;
		ig = 0;
		ih = 0;
		is = 0;
		for(f=1; f<=narea; f++)
		{
			for(g=1; g<=ngroup; g++)
			{
				ih ++;
				pntr_ag(f,g) = ih;
				for(h=1;h<=nsex;h++)
				{
					ig ++;
					n_area(ig)  = f;
					n_group(ig) = g;
					n_sex(ig)   = h;
					pntr_ags(f,g,h) = ig;
					if(f==1)
					{
						is ++;
						pntr_gs(g,h) = is;
					}
				}
			}
		}
		cout<<"| ----------------------- |"<<endl;
		cout<<"| MODEL DIMENSION         |"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<"| narea  \t"<<narea<<endl;
		cout<<"| ngroup \t"<<ngroup<<endl;
		cout<<"| nsex   \t"<<nsex<<endl;
		cout<<"| syr    \t"<<syr<<endl;
		cout<<"| nyr    \t"<<nyr<<endl;
		cout<<"| sage   \t"<<sage<<endl;
		cout<<"| nage   \t"<<nage<<endl;
		cout<<"| ngear  \t"<<ngear<<endl;
		cout<<"| n_area \t"<<n_area<<endl;
		cout<<"| n_group\t"<<n_group<<endl;
		cout<<"| n_sex  \t"<<n_sex<<endl;
		cout<<"| pntr_ag\n"<<pntr_ag<<endl;
		cout<<"| pntr_gs\n"<<pntr_gs<<endl;
		cout<<"| pntr_ags\n"<<pntr_ags(1)<<endl;
		cout<<"| ----------------------- |\n"<<endl;
		
		
		/* Check for dimension errors in projection control file. */
		if(pf_cntrl(1)<syr || pf_cntrl(3)<syr || pf_cntrl(5)<syr )
		{
			cout<<"ERROR: start year in projection file control is less than initial model year."<<endl;
			exit(1);
		}
		if(pf_cntrl(2)>nyr || pf_cntrl(4)>nyr || pf_cntrl(6)>nyr )
		{
			cout<<"ERROR: last year in projection file control is greater than last model year."<<endl;
			exit(1);
		}
	END_CALCS
	
	
	// |---------------------------------------------------------------------------------|
	// | Allocation for each gear in (ngear), use 0 for survey gears.
	// |---------------------------------------------------------------------------------|
	// | fsh_flag is used to determine which fleets should be in MSY-based referecen points
	// | If dAllocation >0 then set fish flag =1 else 0
	// | nfleet is the number of non-survey gear fleet with dAllocations > 0
	// |

	int nfleet;
	init_vector dAllocation(1,ngear);
	
	//init_ivector catch_sex_composition(1,ngear); 
	//init_ivector catch_type(1,ngear);

	ivector fsh_flag(1,ngear);
	LOC_CALCS
		dAllocation = dAllocation/sum(dAllocation);
		for(k=1;k<=ngear;k++)
		{
			if(dAllocation(k)>0)
				fsh_flag(k)=1;
			else
				fsh_flag(k)=0;
		}
		nfleet = sum(fsh_flag);
	END_CALCS
	
	ivector nFleetIndex(1,nfleet);
	LOC_CALCS
		j = 1;
		for(k=1; k<=ngear;k++)
		{
			if(fsh_flag(k)) nFleetIndex(j++) = k;
		}
		// cout<<"nFleetIndex index\t"<<nFleetIndex<<endl;
	END_CALCS
	
	
	// |---------------------------------------------------------------------------------|
	// | Growth and maturity parameters
	// |---------------------------------------------------------------------------------|
	// | n_ags -> number of areas * groups * sex
	// |

	init_vector  d_linf(1,n_ags);
	init_vector d_vonbk(1,n_ags);
	init_vector    d_to(1,n_ags);
	init_vector     d_a(1,n_ags);
	init_vector     d_b(1,n_ags);
	init_vector    d_ah(1,n_ags);
	init_vector    d_gh(1,n_ags);
	init_int 		n_MAT;
	int t1;
	int t2;
	LOC_CALCS
		if(n_MAT)
		{
			t1 = sage;
			t2 = nage;
		}
		else
		{
			t1 = 0;
			t2 = 0;
		}
	END_CALCS 
	!! COUT(n_MAT);
	!! COUT(t1);
	!! COUT(t2);
	init_vector 	d_maturityVector(t1,t2);
	!! COUT(d_maturityVector);

	matrix la(1,n_ags,sage,nage);		//length-at-age
	matrix wa(1,n_ags,sage,nage);		//weight-at-age
	matrix ma(1,n_ags,sage,nage);		//maturity-at-age
	LOC_CALCS
		cout<<setw(8)<<setprecision(4)<<endl;
	  	cout<<"| ----------------------- |"<<endl;
		cout<<"| GROWTH PARAMETERS       |"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<"| d_linf  \t"<<d_linf<<endl;
	  	cout<<"| d_vonbk \t"<<d_vonbk<<endl;
	  	cout<<"| d_to    \t"<<d_to<<endl;
	  	cout<<"| d_a     \t"<<d_a<<endl;
	  	cout<<"| d_b     \t"<<d_b<<endl;
	  	cout<<"| d_ah    \t"<<d_ah<<endl;
	  	cout<<"| d_gh    \t"<<d_gh<<endl;
	  	cout<<"| ----------------------- |\n"<<endl;

	  	// length & weight-at-age based on input growth pars
	  	ma.initialize();
	  	for(ig=1;ig<=n_ags;ig++)
	  	{
	  		la(ig) = d_linf(ig)*(1. - exp(-d_vonbk(ig)*(age-d_to(ig))));
	  		wa(ig) = d_a(ig) * pow(la(ig),d_b(ig));
	  		h = n_sex(ig);
	  		if(n_MAT==0)
	  		{
	  			ma(ig) = plogis(age,d_ah(ig),d_gh(ig));
	  		}
	  		else if( n_MAT>0 && h !=2 )
	  		{
	  			ma(ig) = d_maturityVector;
	  		}
	  	}
	END_CALCS
	
	// |---------------------------------------------------------------------------------|
	// | Historical removal
	// |---------------------------------------------------------------------------------|
	// | - Total catch in weight (type=1), numbers (type=2), or roe (type=3).
	// | - dCatchData matrix cols: (year gear area group sex type value).
	// | - If total catch is asexual (sex=0), pool predicted catch from nsex groups.
	// | - ft_count    -> Number of estimated fishing mortality rate parameters.
	// | - d3_Ct -> An array of observed catch in group(ig) year (row) by gear (col)
	// | - [?] - TODO: fix special case where nsex==2 and catch sex = 0 in catch array.
	init_int nCtNobs;
	!! COUT(nCtNobs)
	init_matrix dCatchData(1,nCtNobs,1,7);
	3darray d3_Ct(1,n_ags,syr,nyr,1,ngear);

	int ft_count;
 	

	LOC_CALCS
		ft_count = nCtNobs;
		cout<<"| ----------------------- |"<<endl;
		cout<<"| HEAD(dCatchData)        |"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<dCatchData.sub(1,3)<<endl;
		cout<<"| ----------------------- |\n"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<"| TAIL(dCatchData)        |"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<dCatchData.sub(nCtNobs-3,nCtNobs)<<endl;
		cout<<"| ----------------------- |\n"<<endl;
		d3_Ct.initialize();
		
		for(int ii=1;ii<=nCtNobs;ii++)
		{
			i = dCatchData(ii)(1);
			k = dCatchData(ii)(2);
			f = dCatchData(ii)(3);
			g = dCatchData(ii)(4);
			h = dCatchData(ii)(5);
			if( h==0 )
			{
				for(h=1;h<=nsex;h++)
				{
					ig = pntr_ags(f,g,h);
					d3_Ct(ig)(i)(k) = 1./nsex*dCatchData(ii)(7);
				}
			}
			else
			{
				ig = pntr_ags(f,g,h);
				d3_Ct(ig)(i)(k) = dCatchData(ii)(7);
			} 
			//if(verbose)  cout<<"Ok after reading catch data"<<ii<<" "<<ig <<" "<<n_ags<<endl;
		}
	END_CALCS
	



	// |---------------------------------------------------------------------------------|
	// | RELATIVE ABUNDANCE INDICIES (ragged array)
	// |---------------------------------------------------------------------------------|
	// | nItNobs     = number of independent surveys
	// | n_it_nobs   = number of survey observations
	// | n_survey_type = 1: survey is proportional to vulnerable numbers
	// | n_survey_type = 2: survey is proportional to vulnerable biomass
	// | n_survey_type = 3: survey is proportional to vulnerable spawning biomass
	// | d3_survey_data: (iyr index(it) gear area group sex wt timing)
	// | it_wt       = relative weights for each relative abundance normalized to have a
	// |               mean = 1 so rho = sig2/(sig^2+tau2) holds true in variance pars.
	// |

	init_int nItNobs;
	ivector        nSurveyIndex(1,nItNobs);
	init_ivector      n_it_nobs(1,nItNobs);
	init_ivector  n_survey_type(1,nItNobs);
	init_3darray d3_survey_data(1,nItNobs,1,n_it_nobs,1,8);
	matrix                it_wt(1,nItNobs,1,n_it_nobs);

// 	!! cout<<"Number of surveys "<<nItNobs<<endl;
	LOC_CALCS
		cout<<"| ----------------------- |"<<endl;
		cout<<"| TAIL(d3_survey_data)       |"<<endl;
		cout<<"| ----------------------- |"<<endl;
		cout<<d3_survey_data(nItNobs).sub(n_it_nobs(nItNobs)-3,n_it_nobs(nItNobs))<<endl;
		cout<<"| ----------------------- |\n"<<endl;
		for(k=1;k<=nItNobs;k++)
		{
			it_wt(k) = column(d3_survey_data(k),7) + 1.e-30;
			nSurveyIndex(k) = d3_survey_data(k)(1,3);
		}
		double tmp_mu = mean(it_wt);
		for(k=1;k<=nItNobs;k++)
		{
			it_wt(k) = it_wt(k)/tmp_mu;
		}
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | AGE COMPOSITION DATA (ragged object)
	// |---------------------------------------------------------------------------------|
	// | - nAgears    -> number of age-composition matrixes, one for each gear.
	// | - n_A_nobs   -> ivector for number of rows in age composition (A) matrix
	// | n_A_sage     -> imatrix for starting age in each row
	// | n_A_nage	  -> imatrix for plus group age in each row
	// | inp_nscaler  -> effective sample size for iterative re-weighting in multinomial.
	// | icol_A       -> number of columns for each row in A.
	// | A            -> array of data (year,gear,area,group,sex|Data...)
	// | d3_A_obs     -> array of catch-age data only.
	// |
	init_int nAgears
	init_ivector n_A_nobs(1,nAgears);
	init_ivector n_A_sage(1,nAgears);
	init_ivector n_A_nage(1,nAgears);
	init_vector  inp_nscaler(1,nAgears);
  // The 5 in the next command is to remove the first 5 columns
  // from the age comp 'data' because they are not the actual ages,
  // but the header data.
	init_3darray d3_A(1,nAgears,1,n_A_nobs,n_A_sage-5,n_A_nage);
	3darray d3_A_obs(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage);
	LOC_CALCS
		if( n_A_nobs(nAgears) > 0 )
			{
			cout<<"| ----------------------- |"<<endl;
			cout<<"| TAIL(A)       |"<<endl;
			cout<<"| ----------------------- |"<<endl;
			cout<<setw(4)<<d3_A(nAgears).sub(n_A_nobs(nAgears)-2,n_A_nobs(nAgears))<<endl;
			cout<<"| ----------------------- |\n"<<endl;
			for(k=1;k<=nAgears;k++)
			{
				dmatrix tmp = trans(trans(d3_A(k)).sub(n_A_sage(k),n_A_nage(k)));
				if(inp_nscaler(k) > 0)
				{
					for( i = 1; i <= n_A_nobs(k); i++ )
					{
						 tmp(i) = tmp(i)/sum(tmp(i)) * inp_nscaler(k);
					}
				}
				d3_A_obs(k) = tmp;
				//d3_A_obs(k) = trans(trans(d3_A(k)).sub(n_A_sage(k),n_A_nage(k)));
			}
		}
		else
		{
			cout<<"| ----------------------- |"<<endl;
			cout<<"| NO AGE DATA"<<endl;
			cout<<"| ----------------------- |"<<endl;
		}
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | EMPIRICAL WEIGHT_AT_AGE DATA
	// |---------------------------------------------------------------------------------|
	// | Mean weight-at-age data (kg) if nWtNobs > 0
	// | sage-5 = year
	// | sage-4 = gear
	// | sage-3 = area
	// | sage-2 = stock
	// | sage-1 = sex
	// | - construct and fill weight-at-age matrix for use in the model code  (d3_wt_avg)
	// | - construct and fill weight-at-age dev matrix for length-based selex (d3_wt_dev)
	// | - construct and fill fecundity-at-age matrix for ssb calculations.   (d3_wt_mat)
	// | [ ] - TODO fix h=0 option for weight-at-age data
	// | [ ] - TODO need to accomodate ragged arrays, or NA values, or partial d3_wt_avg.
	// |

	init_int nWtNobs;
	init_matrix inp_wt_avg(1,nWtNobs,sage-5,nage);

	matrix  dWt_bar(1,n_ags,sage,nage);
	3darray d3_wt_avg(1,n_ags,syr,nyr+1,sage,nage);
	3darray d3_wt_dev(1,n_ags,syr,nyr+1,sage,nage);
	3darray d3_wt_mat(1,n_ags,syr,nyr+1,sage,nage);
	3darray d3_len_age(1,n_ags,syr,nyr+1,sage,nage);
// 	vector fa_bar(sage,nage);				//average fecundity-at-age for all years.
// 	vector avg_fec(sage,nage);				//average fecundity-at-age
// 	vector avg_wt(sage,nage);				//average weight-at-age
	LOC_CALCS
		d3_wt_avg.initialize();
		d3_wt_dev.initialize();
		d3_wt_mat.initialize();
		d3_len_age.initialize();

		for(ig=1;ig<=n_ags;ig++)
		{
			for(int i=syr;i<=nyr;i++)
			{
				d3_wt_avg(ig)(i) = wa(ig);
				d3_wt_mat(ig)(i) = elem_prod(ma(ig),wa(ig));
				d3_len_age(ig)(i) = pow(wa(ig)/d_a(ig),1./d_b(ig));
			}
		}

		// the overwrite d3_wt_avg & d3_wt_mat with existing empirical data
		// SM Sept 6, 2013. Added option of using NA values (-99.0) for
		// missing weight-at-age data, or truncated age-data.
		int iyr;
		for(i=1;i<=nWtNobs;i++)
		{
			iyr = inp_wt_avg(i,sage-5);
			f   = inp_wt_avg(i,sage-3);
			g   = inp_wt_avg(i,sage-2);
			h   = inp_wt_avg(i,sage-1);

			// | SM Changed Sept 9, to accomodate NA's (-99) in empirical data.
			if( h )
			{
				ig                   = pntr_ags(f,g,h);
				dvector tmp          = inp_wt_avg(i)(sage,nage);
				ivector idx          = getIndex(age,tmp);
				for( int ii = 1; ii <= size_count(idx); ii++ )
				{
					d3_wt_avg(ig)(iyr)(idx(ii)) = inp_wt_avg(i)(idx(ii));
					d3_len_age(ig)(iyr)(idx(ii))= pow(d3_wt_avg(ig)(iyr)(idx(ii))
					                                  /d_a(ig),1./d_b(ig));
				}
				//d3_wt_avg(ig)(iyr)(idx) = inp_wt_avg(i)(idx);
				d3_wt_mat(ig)(iyr)      = elem_prod(ma(ig),d3_wt_avg(ig)(iyr));
				//cout<<"Yep \t"<<inp_wt_avg(i)(idx)<<endl;
				//cout<<"Yep \t"<<tmp(idx)<<endl;
				//cout<<"Yep \t"<<d3_wt_avg(ig)(iyr)(idx)<<endl;
			}
			else if( !h ) 
			{
					//cout<<h<<endl;

				for(int h=1;h<=nsex;h++)
				{
					ig                   = pntr_ags(f,g,h);
					dvector tmp          = inp_wt_avg(i)(sage,nage);
					ivector idx          = getIndex(age,tmp);
					// Problem, array indexed differ, must loop over idx;
					// d3_wt_avg(ig)(iyr)(idx) = inp_wt_avg(i)(idx);
					for( int ii = 1; ii <= size_count(idx); ii++)
					{
						d3_wt_avg(ig)(iyr)(idx(ii)) = inp_wt_avg(i)(idx(ii));
						d3_len_age(ig)(iyr)(idx(ii)) = pow(d3_wt_avg(ig)(iyr)(idx(ii))
						                               /d_a(ig),1./d_b(ig));
					}
					d3_wt_mat(ig)(iyr)      = elem_prod(ma(ig),d3_wt_avg(ig)(iyr));
				}
			}
		}
		

		// average weight-at-age in projection years
		for(ig=1;ig<=n_ags;ig++)
		{
			dWt_bar(ig)        = colsum(d3_wt_avg(ig).sub(pf_cntrl(3),pf_cntrl(4)));
			dWt_bar(ig)       /= pf_cntrl(4)-pf_cntrl(3)+1;
			d3_wt_avg(ig)(nyr+1) = dWt_bar(ig);
			d3_wt_mat(ig)(nyr+1) = elem_prod(dWt_bar(ig),ma(ig));
			d3_len_age(ig)(nyr+1) = pow(dWt_bar(ig)/d_a(ig),1./d_b(ig));
		}
		
		
		// deviations in mean weight-at-age
		for(ig=1;ig<=n_ags;ig++)
		{
			dmatrix mtmp = trans( d3_wt_avg(ig) );
			//COUT(mtmp);

			for(j=sage;j<=nage;j++)
			{
				//COUT(sum(first_difference(mtmp(j)(syr,nyr))));
				if( sum( first_difference(mtmp(j)(syr,nyr))) )
				{
					mtmp(j) = ( mtmp(j)-mean(mtmp(j)(syr,nyr)) ) 
							/ sqrt(var(mtmp(j)(syr,nyr)));
				}
				else
				{
					mtmp(j) = 0;
				}
			}
			d3_wt_dev(ig) = trans(mtmp);
		
			
			if( min(d3_wt_avg(ig))<=0.000 && min(d3_wt_avg(ig))!=NA )
			{
				cout<<"|-----------------------------------------------|"<<endl;
				cout<<"| ERROR IN INPUT DATA FILE FOR MEAN WEIGHT DATA |"<<endl;
				cout<<"|-----------------------------------------------|"<<endl;
				cout<<"| - Cannot have an observed mean weight-at-age  |"<<endl;
				cout<<"|   less than or equal to 0.  Please fix.       |"<<endl;
				cout<<"| - You are permitted to use '-99.0' for missing|"<<endl;
				cout<<"|   values in your weight-at-age data.          |"<<endl;
				cout<<"| - Aborting program!                           |"<<endl;
				cout<<"|-----------------------------------------------|"<<endl;
				ad_exit(1);
			}
		}
	END_CALCS
	
	
	// |---------------------------------------------------------------------------------|
	// | END OF DATA FILE
	// |---------------------------------------------------------------------------------|
	// |
	init_int eof;	
	LOC_CALCS
	  if(eof==999){
		cout<<"\n| -- END OF DATA SECTION -- |\n";
	  	cout<<"|         eof = "<<eof<<"         |"<<endl;
		cout<<"|___________________________|"<<endl;
	  }else{
		cout<<"\n *** ERROR READING DATA *** \n"<<endl; exit(1);
	  }
	END_CALCS

	
	
	// |---------------------------------------------------------------------------------|
	// | VARIABLES FOR MSY-BASED REFERENCE POINTS
	// |---------------------------------------------------------------------------------|
	// |
	vector fmsy(1,nfleet);			//Fishing mortality rate at Fmsy
	vector fall(1,nfleet);			//Fishing mortality based on dAllocation
	vector  msy(1,nfleet);			//Maximum sustainable yield
	number bmsy;					//Spawning biomass at MSY
 // number Umsy;					//Exploitation rate at MSY
	vector age_tau2(1,nAgears);	//MLE estimate of the variance for age comps
 // 	//catch-age for simulation model (could be declared locally 3d_array)
 // 	3darray d3C(1,ngear,syr,nyr,sage,nage);		
	
	
		
	
	// |---------------------------------------------------------------------------------|
	// | CONTROL FILE
	// |---------------------------------------------------------------------------------|
	// |
	!! ad_comm::change_datafile_name(ControlFile);
	

	// |---------------------------------------------------------------------------------|
	// | Leading Parameters
	// |---------------------------------------------------------------------------------|
	// | npar            -> number of leading parameters
	// | ipar_vector     -> integer vector based on the number of areas groups sexes
	// | -1) log_ro      - unfished sage recruitment
	// | -2) steepness   - steepness of the stock-recruitment relationship
	// | -3) log_m       - instantaneous natural mortality rate
	// | -4) log_avgrec  - average sage recruitment from syr+1 to nyr
	// | -5) log_recinit - average sage recruitment for initialization
	// | -6) rho         - proportion of total variance for observation errors
	// | -7) vartheta    - total precision (1/variance)
	init_int npar;
	init_matrix theta_control(1,npar,1,7);
	
	vector   theta_ival(1,npar);
	vector     theta_lb(1,npar);
	vector     theta_ub(1,npar);
	ivector   theta_phz(1,npar);
	ivector theta_prior(1,npar);
	ivector ipar_vector(1,npar);
	LOC_CALCS
		theta_ival  = column(theta_control,1);
		theta_lb    = column(theta_control,2);
		theta_ub    = column(theta_control,3);
		theta_phz   = ivector(column(theta_control,4));
		theta_prior = ivector(column(theta_control,5));
		ipar_vector(1,2) = ngroup;
		ipar_vector(6,7) = ngroup;
		ipar_vector(3)   = n_gs;
		ipar_vector(4,5) = n_ag;
	END_CALCS
	
	// |---------------------------------------------------------------------------------|
	// | CONTROLS PARAMETERS FOR AGE/SIZE COMPOSITION DATA FOR na_gears                  |
	// |---------------------------------------------------------------------------------|
	// |
	
	init_ivector nCompIndex(1,nAgears);
	init_ivector nCompLikelihood(1,nAgears);
	init_vector  dMinP(1,nAgears);
	init_vector  dEps(1,nAgears);
	init_ivector nPhz_age_tau2(1,nAgears);
	init_ivector nPhz_phi1(1,nAgears);
	init_ivector nPhz_phi2(1,nAgears);
	init_ivector nPhz_df(1,nAgears);
	init_int check;
	!! if(check != -12345) {COUT(check);cout<<"Error reading composition controls\n"<<endl; exit(1);}



	// |---------------------------------------------------------------------------------|
	// | CONTROLS FOR SELECTIVITY OPTIONS
	// |---------------------------------------------------------------------------------|
	// | - 12 different options for modelling selectivity which are summarized here:
	// | - isel_npar  -> ivector for # of parameters for each gear.
	// | - jsel_npar  -> ivector for the number of rows for time-varying selectivity.
	// | 
	// | SEL_TYPE  DESCRIPTION
	// |    1      age-based logistic function with 2 parameters.
	// |    2      age-based selectivity coefficients with nage-sage parameters.
	// |    3      cubic spline with age knots.
	// |    4      time-varying cubic spline with age knots.
	// |    5      time-varying bicubic spline with age and year knots.
	// |    6      logistic with fixed parameters.
	// |    7      logistic function of body weight with 2 parameters.
	// |    8      logistic 3 parameter function based on mean weight deviations.
	// |    11     length-based logistic function with 2 parametrs based on mean length.
	// |    12     length-based selectivity coefficients with cubic spline interpolation.
	// |
	// | selex_controls (1-10)
	// |  1  -> isel_type - switch for selectivity.
	// |  2  -> ahat      - age-at-50% vulnerbality for logistic function.
	// |  3  -> ghat      - std at 50% age of vulnerability for logistic function.
	// |  4  -> age_nodes - No. of age-nodes for bicubic spline.
	// |  5  -> yr_nodes  - No. of year-nodes for bicubic spline.
	// |  6  -> sel_phz   - phase for estimating selectivity parameters.
	// |  7  -> lambda_1  - penalty weight for 2nd difference in selectivity.
	// |  8  -> lambda_2  - penalty weight for dome-shaped selectivity.
	// |  9  -> lambda_3  - penalty weight for 2nd difference in time-varying selectivity.
	// |  10 -> Number of discrete selectivity blocks.
	// |

	init_matrix selex_controls(1,10,1,ngear);
	

	ivector    isel_npar(1,ngear);	
	ivector    jsel_npar(1,ngear);	
	ivector    isel_type(1,ngear);	
	ivector      sel_phz(1,ngear);	
	ivector n_sel_blocks(1,ngear);	

	vector      ahat(1,ngear);	
	vector      ghat(1,ngear);	
	vector age_nodes(1,ngear);	
	vector  yr_nodes(1,ngear);	
	vector  lambda_1(1,ngear);	
	vector  lambda_2(1,ngear);	
	vector  lambda_3(1,ngear);	
	
	LOC_CALCS
		ahat      = selex_controls(2);
		ghat      = selex_controls(3);
		age_nodes = selex_controls(4);
		yr_nodes  = selex_controls(5);
		lambda_1  = selex_controls(7);
		lambda_2  = selex_controls(8);
		lambda_3  = selex_controls(9);

		isel_type    = ivector(selex_controls(1));
		sel_phz      = ivector(selex_controls(6));
		n_sel_blocks = ivector(selex_controls(10));
	END_CALCS
	
	init_imatrix sel_blocks(1,ngear,1,n_sel_blocks);


	LOC_CALCS
		// | COUNT THE NUMBER OF ESTIMATED SELECTIVITY PARAMETERS TO ESTIMATE
		// | isel_npar number of columns for each gear.
		// | jsel_npar number of rows for each gear.
		isel_npar.initialize();
		for(i=1;i<=ngear;i++)
		{	
			jsel_npar(i)=1;
			switch(isel_type(i))
			{
				case 1:
					// logistic selectivity
					isel_npar(i) = 2;
					jsel_npar(i) = n_sel_blocks(i); 
					break;
					
				case 2:
					// age-specific coefficients
					isel_npar(i) = (nage-sage);
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				case 3:
				 	// cubic spline 
					isel_npar(i) = age_nodes(i);
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				case 4:	 
					// annual cubic splines
					isel_npar(i) = age_nodes(i);
					jsel_npar(i) = (nyr-syr-retro_yrs)+1;
					break;
					
				case 5:  
					// bicubic spline
					jsel_npar(i) = age_nodes(i);
					isel_npar(i) = yr_nodes(i);
					break;
				
				case 6:
					// fixed logistic (no parameters estimated)
					// ensure sel_phz is set to negative value.
					isel_npar(i) = 2;
					if(sel_phz(i)>0) sel_phz(i) = -1;
					break;
					
				case 7:
					// CHANGED: Now working, Vivian Haist fixed it.
					// logistic (3 parameters) with mean body 
					// weight deviations. 
					isel_npar(i) = 2;
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				case 8:
					// Alternative logistic selectivity with d3_wt_dev coefficients.
					isel_npar(i) = 3;
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				case 11:
					// Logistic length-based selectivity.
					isel_npar(i) = 2;
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				case 12:
					// Length-based selectivity coeffs with cubic spline interpolation
					isel_npar(i) = age_nodes(i);
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				default: break;
			}
		}
	END_CALCS
	





	// |---------------------------------------------------------------------------------|
	// | PRIOR FOR RELATIVE ABUNDANCE DATA
	// |---------------------------------------------------------------------------------|
	// | nits     -> number of relative abundance indices
	// | q_prior  -> type of prior to use, see legend
	// | mu_log_q -> mean q in log-space
	// | sd_log_q -> std of q prior in log-space
	// |
	// | q_prior type:
	// | 0 -> uninformative prior.
	// | 1 -> normal prior on q in log-space.
	// | 2 -> penalized random walk in q.



	init_int nits;					
	init_ivector q_prior(1,nits);
	init_vector mu_log_q(1,nits);
	init_vector sd_log_q(1,nits);
	
	// |---------------------------------------------------------------------------------|
	// | Miscellaneous controls                                                          |
	// |---------------------------------------------------------------------------------|
	// | 1 -> verbose
	// | 2 -> recruitment model (1=beverton-holt, 2=rickers)
	// | 3 -> std in catch first phase
	// | 4 -> std in catch in last phase
	// | 5 -> assumed unfished in first year (0=FALSE, 1=TRUE)
	// | 6 -> minimum proportion at age to consider in the dmvlogistic likelihood
	// | 7 -> mean fishing mortality rate to regularize the solution
	// | 8 -> standard deviation of mean F penalty in first phases
	// | 9 -> standard deviation of mean F penalty in last phase.
	// | 10-> phase for estimating deviations in natural mortality.
	// | 11-> std in natural mortality deviations.
	// | 12-> number of estimated nodes for deviations in natural mortality
	// | 13-> fraction of total mortality that takes place prior to spawning
	// | 14-> switch for age-composition likelihood (1=dmvlogistic,2=dmultinom)
	// | 15-> switch for generating selex based on IFD and cohort biomass
	init_vector d_iscamCntrl(1,15);
	int verbose;
	

	init_int eofc;
	LOC_CALCS
		verbose = d_iscamCntrl(1);
		if(verbose) COUT(d_iscamCntrl);
		
		if(eofc==999){
			cout<<"\n| -- END OF CONTROL SECTION -- |\n";
		  	cout<<"|          eofc = "<<eofc<<"          |"<<endl;
			cout<<"|______________________________|"<<endl;
		}else{
			cout<<"\n ***** ERROR CONTROL FILE ***** \n"<<endl; exit(1);
		}
	END_CALCS
	
	int nf;
	
	// |---------------------------------------------------------------------------------|
	// | VECTOR DIMENTIONS FOR NEGATIVE LOG LIKELIHOODS
	// |---------------------------------------------------------------------------------|
	// | ilvec[1,5,6,7] -> number of fishing gears (ngear)
	// | ilvec[2]       -> number of surveys       (nItNobs)
	// | ilvec[3]       -> number of age-compisition data sets (nAgears)
	// | ilvec[4]       -> container for recruitment deviations.
	ivector ilvec(1,7);
	!! ilvec    = ngear;
	!! ilvec(1) = 1;			
	!! ilvec(2) = nItNobs;			
	!! ilvec(3) = nAgears;		
	!! ilvec(4) = ngroup;
	

	// |---------------------------------------------------------------------------------|
	// | RETROSPECTIVE ADJUSTMENT TO nyrs
	// |---------------------------------------------------------------------------------|
	// | - Don't read any more input data from here on in. 
	// | - Modifying nyr to allow for retrospective analysis.
	// | - If retro_yrs > 0, then ensure that pf_cntrl arrays are not greater than nyr,
	// |   otherwise arrays for mbar will go out of bounds.
	// | 

	!! nyr = nyr - retro_yrs;
	
	LOC_CALCS
		if(retro_yrs)
		{
			if(pf_cntrl(2)>nyr) pf_cntrl(2) = nyr;
			if(pf_cntrl(4)>nyr) pf_cntrl(4) = nyr;
			if(pf_cntrl(6)>nyr) pf_cntrl(6) = nyr;
		}
	END_CALCS


	// |---------------------------------------------------------------------------------|
	// | MANAGEMENT STRATEGY EVALUATION INPUTS
	// |---------------------------------------------------------------------------------|
	// |
	

	LOC_CALCS
		ifstream ifile("Halibut2012.mse");
		if(ifile)
		{
			cout<<"Vader is happy"<<endl;
			readMseInputs();
			
			
			exit(1);
		}
	END_CALCS


	// END OF DATA_SECTION
	!! if(verbose) cout<<"||-- END OF DATA_SECTION --||"<<endl;

	

INITIALIZATION_SECTION
  theta theta_ival;
  phi1 0.01;
	
PARAMETER_SECTION
	// |---------------------------------------------------------------------------------|
	// | LEADING PARAMTERS
	// |---------------------------------------------------------------------------------|
	// | - Initialized in the INITIALIZATION_SECTION with theta_ival from control file.
	// | [ ] Change to init_bounded_vector_vector.
	// | theta[1] -> log_ro, or log_msy
	// | theta[2] -> steepness(h), or log_fmsy
	// | theta[3] -> log_m
	// | theta[4] -> log_avgrec
	// | theta[5] -> log_recinit
	// | theta[6] -> rho
	// | theta[7] -> vartheta
	// |
	init_bounded_vector_vector theta(1,npar,1,ipar_vector,theta_lb,theta_ub,theta_phz);
	
	// |---------------------------------------------------------------------------------|
	// | SELECTIVITY PARAMETERS
	// |---------------------------------------------------------------------------------|
	// | - This is a bounded matrix vector where the dimensions are defined by the 
	// | - selectivity options specified in the control file.
	// | - There are 1:ngear arrays, having jsel_npar rows and isel_npar columns.
	// | - If the user has not specified -ainp or -binp, the initial values are set
	// |   based on ahat and ghat in the control file for logistic selectivities.
	// | - Special case: if SimFlag=TRUE, then add some random noise to ahat.
	// | - NB  sel_par is in log space.
	// |
	init_bounded_matrix_vector sel_par(1,ngear,1,jsel_npar,1,isel_npar,-25.,25.,sel_phz);

	LOC_CALCS
		if ( !global_parfile )
		{
			for(int k=1; k<=ngear; k++)
			{
				if( isel_type(k)==1 || 
					isel_type(k)==6 || 
					(
					isel_type(k)>=7 && 
					isel_type(k) != 12 
					)
					)
				{
					for(int j = 1; j <= n_sel_blocks(k); j++ )
					{
						double uu = 0;
						if(SimFlag && j > 1)
						{
							uu = 0.05*randn(j+rseed);
						} 
						sel_par(k,j,1) = log(ahat(k)*exp(uu));
						sel_par(k,j,2) = log(ghat(k));
					}
				}
			}
		}
	END_CALCS
	

	// |---------------------------------------------------------------------------------|
	// | FISHING MORTALITY RATE PARAMETERS
	// |---------------------------------------------------------------------------------|
	// | - Estimate all fishing mortality rates in log-space.
	// | - If in simulation mode then initialize with F=0.1; Actual F is conditioned on 
	// |   the observed catch.
	// |
	
	init_bounded_vector log_ft_pars(1,ft_count,-30.,3.0,1);
	
	LOC_CALCS
		if(!SimFlag) log_ft_pars = log(0.10);
	END_CALCS
	
	

	// |---------------------------------------------------------------------------------|
	// | INITIAL AND ANNUAL RECRUITMENT 
	// |---------------------------------------------------------------------------------|
	// | - Estimate single mean initial recruitment and deviations for each initial 
	// |   cohort from sage+1 to nage. (Rinit + init_log_rec_devs)
	// | - Estimate mean overal recruitment and annual deviations from syr to nyr.
	// | - d_iscamCntrl(5) is a flag to initialize the model at unfished recruitment (ro),
	// |   if this is true, then do not estimate init_log_rec_devs
	// | [ ] - TODO add dev contstraint for rec_devs in calc_objective_function.

	!! int init_dev_phz = 2;
	!! if(d_iscamCntrl(5)) init_dev_phz = -1;
	init_bounded_matrix init_log_rec_devs(1,n_ag,sage+1,nage,-15.,15.,init_dev_phz);
	init_bounded_matrix log_rec_devs(1,n_ag,syr,nyr,-15.,15.,2);
	


	// |---------------------------------------------------------------------------------|
	// | DEVIATIONS FOR NATURAL MORTALITY BASED ON CUBIC SPLINE INTERPOLATION
	// |---------------------------------------------------------------------------------|
	// | - Estimating trends in natural mortality rates, where the user specified the 
	// |   number of knots (d_iscamCntrl(12)) and the std in M in the control file, and the phase
	// |   in which to estimate natural mortality devs (d_iscamCntrl(10)).  If the phase is neg.
	// |   then natural mortality rate deviates are not estimated and M is assumed const.
	// | - This model is implemented as a random walk, where M{t+1} = M{t} + dev.
	
	!! int m_dev_phz = -1;
	!!     m_dev_phz = d_iscamCntrl(10);
	!! int  n_m_devs = d_iscamCntrl(12);
	init_bounded_vector log_m_nodes(1,n_m_devs,-5.0,5.0,m_dev_phz);
	

	// |---------------------------------------------------------------------------------|
	// | CORRELATION COEFFICIENTS FOR AGE COMPOSITION DATA USED IN LOGISTIC NORMAL       |
	// |---------------------------------------------------------------------------------|
	// | log_age_tau2 is the variance of the composition errors.
	// | phi1 is the AR1 coefficient
	// | phi2 used in AR2 process.
	init_bounded_number_vector log_age_tau2(1,nAgears,-4.65,5.30,nPhz_age_tau2);
	init_bounded_number_vector phi1(1,nAgears,-1.0,1.0,nPhz_phi1);
	init_bounded_number_vector phi2(1,nAgears,0.0,1.0,nPhz_phi2);
	init_bounded_number_vector log_degrees_of_freedom(1,nAgears,0.70,10.0,nPhz_df);

	// |---------------------------------------------------------------------------------|
	// | AUTOCORRELATION IN RECRUITMENT DEVIATIONS                                       |
	// |---------------------------------------------------------------------------------|
	// | gamma_r: what fraction of the residual from year t-2 carries over to t-1.
	init_bounded_number gamma_r(0,1,-4);
	!!gamma_r = 0;

	// |---------------------------------------------------------------------------------|
	// | OBJECTIVE FUNCTION VALUE
	// |---------------------------------------------------------------------------------|
	// | - the value that ADMB will minimize, called objfun in iSCAM
	// |
	objective_function_value objfun;
	

    // |---------------------------------------------------------------------------------|
    // | POPULATION VARIABLES
    // |---------------------------------------------------------------------------------|
    // | - m_bar       -> Average natural mortality rate from syr to nyr.
    // | - rho         -> Proportion of total variance associated with obs error.
    // | - varphi      -> Total precision of CPUE and Recruitment deviations.
    // | - sig         -> STD of the observation errors in relative abundance data.
    // | - tau         -> STD of the process errors (recruitment deviations).
    // |
	number m_bar;	///< Average natural mortality rate.			
	number rho;					
	number varphi				
	number sig;					
	number tau; 				


	// |---------------------------------------------------------------------------------|
	// | POPULATION VECTORS
	// |---------------------------------------------------------------------------------|
    // | - ro          -> theoretical unfished age-sage recruits. 
    // | - bo          -> theoretical unfished spawning biomass (MSY-based ref point).
    // | - sbo         -> unfished spawning biomass at the time of spawning.
    // | - kappa       -> Goodyear recruitment compensation ratio K = 4h/(1-h); h=K/(4+K)
    // | - so          -> Initial slope (max R/S) of the stock-recruitment relationship.
    // | - beta        -> Density dependent term in the stock-recruitment relationship.
    // | - m           -> Instantaneous natural mortality rate by nsex
    // | - log_avgrec  -> Average sage recruitment(syr-nyr,area,group).
    // | - log_recinit -> Avg. initial recruitment for initial year cohorts(area,group).
	// | - log_m_devs  -> annual deviations in natural mortality.
	// | - q           -> conditional MLE estimates of q in It=q*Bt*exp(epsilon)
	// | - ct          -> predicted catch for each catch observation
	// | - eta         -> standardized log residual (log(obs_ct)-log(ct))/sigma_{ct}
	// |
	 
	vector        ro(1,ngroup);
	vector        bo(1,ngroup);
	vector       sbo(1,ngroup);
	vector     kappa(1,ngroup);
	vector steepness(1,ngroup);
	vector        so(1,ngroup);
	vector      beta(1,ngroup);
	vector           m(1,n_gs);	
	vector  log_avgrec(1,n_ag);			
	vector log_recinit(1,n_ag);			
	vector          q(1,nItNobs);
	vector         ct(1,nCtNobs);
	vector        eta(1,nCtNobs);	
	vector log_m_devs(syr+1,nyr);
	
	// |---------------------------------------------------------------------------------|
	// | MATRIX OBJECTS
	// |---------------------------------------------------------------------------------|
	// | - log_rt   -> age-sage recruitment for initial years and annual recruitment.
	// | - catch_df -> Catch data_frame (year,gear,area,group,sex,type,obs,pred,resid)
	// | - eta      -> log residuals between observed and predicted total catch.
	// | - nlvec    -> matrix for negative loglikelihoods.
	// | - epsilon  -> residuals for survey abundance index
	// | - it_hat   -> predicted survey index (no need to be differentiable)
	// | - qt       -> catchability coefficients (time-varying)
	// | - sbt      -> spawning stock biomass by group used in S-R relationship.
	// | - bt       -> average biomass by group used for stock projection
	// | - rt          -> predicted sage-recruits based on S-R relationship.
	// | - delta       -> residuals between estimated R and R from S-R curve (process err)
	// | 
	matrix  log_rt(1,n_ag,syr-nage+sage,nyr);
	matrix   nlvec(1,7,1,ilvec);	
	matrix epsilon(1,nItNobs,1,n_it_nobs);
	matrix  it_hat(1,nItNobs,1,n_it_nobs);
	matrix      qt(1,nItNobs,1,n_it_nobs);
	matrix     sbt(1,ngroup,syr,nyr+1);
	matrix      bt(1,ngroup,syr,nyr+1);
	matrix      rt(1,ngroup,syr+sage,nyr); 
	matrix   delta(1,ngroup,syr+sage,nyr);


	// |---------------------------------------------------------------------------------|
	// | THREE DIMENSIONAL ARRAYS
	// |---------------------------------------------------------------------------------|
	// | - ft       -> Mean fishing mortality rates for (area-sex, gear, year)
	// | F          -> Instantaneous fishing mortality rate for (group,year,age)
	// | M          -> Instantaneous natural mortality rate for (group,year,age)
	// | Z          -> Instantaneous total  mortalityr rate Z=M+F for (group,year,age)
	// | S          -> Annual survival rate exp(-Z) for (group,year,age)
	// | N          -> Numbers-at-age for (group,year+1,age)
	// | - A_hat    -> ragged matrix for predicted age-composition data.
	// | - A_nu		-> ragged matrix for age-composition residuals.
	// | 
	// |
	3darray  ft(1,n_ags,1,ngear,syr,nyr);
	3darray   F(1,n_ags,syr,nyr,sage,nage);
	3darray   M(1,n_ags,syr,nyr,sage,nage);
	3darray   Z(1,n_ags,syr,nyr,sage,nage);
	3darray   S(1,n_ags,syr,nyr,sage,nage);
	3darray   N(1,n_ags,syr,nyr+1,sage,nage);
	3darray  A_hat(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage);
	3darray   A_nu(1,nAgears,1,n_A_nobs,n_A_sage,n_A_nage);
	
	// //matrix jlog_sel(1,ngear,sage,nage);		//selectivity coefficients for each gear type.
	// //matrix log_sur_sel(syr,nyr,sage,nage);	//selectivity coefficients for survey.
	 
	// matrix Z(syr,nyr,sage,nage);
	// matrix S(syr,nyr,sage,nage);
	// matrix pit(1,nItNobs,1,n_it_nobs);			//predicted relative abundance index
	

	// 3darray Ahat(1,nAgears,1,n_A_nobs,n_A_sage-2,n_A_nage);		//predicted age proportions by gear & year
	// 3darray A_nu(1,nAgears,1,n_A_nobs,n_A_sage-2,n_A_nage);		//residuals for age proportions by gear & year
	
	// |---------------------------------------------------------------------------------|
	// | FOUR DIMENSIONAL ARRAYS
	// |---------------------------------------------------------------------------------|
	// | log_sel    -> Selectivity for (gear, group, year, age)
	// | Chat       -> Predicted catch-age array for (gear, group, year, age)
	// | 
	4darray log_sel(1,ngear,1,n_ags,syr,nyr,sage,nage);
	// 4darray    Chat(1,ngear,1,n_ags,syr,nyr,sage,nage);		
	

	// |---------------------------------------------------------------------------------|
	// | SDREPORT VARIABLES AND VECTORS
	// |---------------------------------------------------------------------------------|
	// | sd_depletion -> Predicted spawning biomass depletion level bt/Bo
	sdreport_vector sd_depletion(1,ngroup);	
	

PRELIMINARY_CALCS_SECTION
	// |---------------------------------------------------------------------------------|
	// | Run the model with input parameters to simulate real data.
	// |---------------------------------------------------------------------------------|
	// | - nf is a function evaluation counter.
 	// | - SimFlag comes from the -sim command line argument to simulate fake data.
 	// |
    nf=0;
  	if( testMSY )
  	{
  		testMSYxls();
  	}
	if( SimFlag ) 
	{
		initParameters();
		
		simulationModel(rseed);
	}
	
	if (NewFiles)
	{
		generate_new_files();	
	}
	
	if(verbose) cout<<"||-- END OF PRELIMINARY_CALCS_SECTION --||"<<endl;
	


RUNTIME_SECTION
    maximum_function_evaluations 100,  200,   500, 25000, 25000
    convergence_criteria        0.01, 0.01, 1.e-3, 1.e-4, 1.e-5


PROCEDURE_SECTION

	initParameters();
	
	calcSelectivities(isel_type);
	
	calcTotalMortality();
	
	calcNumbersAtAge();
	
	calcTotalCatch();
	
	calcAgeComposition();

	//calcLengthComposition();  // WIP: need to add variables for length comps
	
	calcSurveyObservations();
	
	calcStockRecruitment();
	
	calcObjectiveFunction();

	calcSdreportVariables();
	
	
	if(mc_phase())
	{
		mcmcPhase=1;
	}
	
	if(mceval_phase())
	{
		mcmcEvalPhase=1;
		mcmc_output();
	}
	
	if( verbose ) {cout<<"End of main function calls"<<endl;}


	/**
	Purpose:  This function calculates the sdreport variables.
	Author: Steven Martell
	
	Arguments:
		None
	
	NOTES:
		
	
	TODO list:
	  [?] - Calculate spawning biomass depletion for each group.
	*/
FUNCTION void calcSdreportVariables()
  {
	sd_depletion.initialize();

	for(g=1;g<=ngroup;g++)
	{
		sd_depletion(g) = sbt(g)(nyr)/sbo(g);
	}
	if( verbose ) { cout<<"**** Ok after calcSdreportVariables ****"<<endl;}
  }


  	/**
  	Purpose: This function extracts the specific parameter values from the theta vector
  	       to initialize the leading parameters in the model.
  	Author: Steven Martell
  	
  	Arguments:
  		None
  	
  	NOTES:
  		- You must call this routine before running the simulation model to generate 
  		  fake data, otherwise you'll have goofy initial values for your leading parameters.
  		- Variance partitioning:
  	  Estimating total variance as = 1/precision
  	  and partition variance by rho = sig^2/(sig^2+tau^2).
  	  
  	  E.g. if sig = 0.2 and tau =1.12 then
  	  rho = 0.2^2/(0.2^2+1.12^2) = 0.03090235
  	  the total variance is kappa^2 = sig^2 + tau^2 = 1.2944
  	
  	TODO list:
  	[ ] - Alternative parameterization using MSY and FMSY as leading parameters (Martell).
  	[*] - avg recruitment limited to area, may consider ragged object for area & stock.
  	
  	*/
FUNCTION void initParameters()
  {

	
  	
  	int ih;
  	
	ro        = mfexp(theta(1));
	steepness = theta(2);
	m         = mfexp(theta(3));
	rho       = theta(6,1);
	varphi    = sqrt(1.0/theta(7,1));
	sig       = sqrt(rho) * varphi;
	tau       = sqrt(1.0-rho) * varphi;

	for(ih=1;ih<=n_ag;ih++)
	{
		log_avgrec(ih)  = theta(4,ih);
		log_recinit(ih) = theta(5,ih);
	}
	

	
	switch(int(d_iscamCntrl(2)))
	{
		case 1:
			//Beverton-Holt model
			kappa = elem_div(4.*steepness,(1.-steepness));
			break;
		case 2:
			//Ricker model
			kappa = pow((5.*steepness),1.25);
		break;
	}
	
	if(verbose)cout<<"**** Ok after initParameters ****"<<endl;
	
  }
	
FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs)
  {
	RETURN_ARRAYS_INCREMENT();
	int nodes=size_count(spline_coffs);
	dvector ia(1,nodes);
	dvector fa(sage,nage);
	ia.fill_seqadd(0,1./(nodes-1));
	fa.fill_seqadd(0,1./(nage-sage));
	vcubic_spline_function ffa(ia,spline_coffs);
	RETURN_ARRAYS_DECREMENT();
	
	//some testing here
	/*dvar_vector spline_nodes(1,nodes);
		spline_nodes.fill_seqadd(-0.5,1./(nodes-1));
		cout<<spline_nodes<<endl;
		vcubic_spline_function test_ffa(ia,spline_nodes);
		cout<<test_ffa(fa)<<endl;
		exit(1);*/
	return(ffa(fa));
  }

FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs, const dvector& la)
  {
	/*interplolation for length-based selectivity coefficeients*/
	RETURN_ARRAYS_INCREMENT();
	int nodes=size_count(spline_coffs);
	dvector ia(1,nodes);
	ia.fill_seqadd(0,1./(nodes-1));
	dvector fa = (la-min(la))/(max(la)-min(la));
	vcubic_spline_function ffa(ia,spline_coffs);
	RETURN_ARRAYS_DECREMENT();
	return(ffa(fa));
  }

  /**
   * @brief cubic spline interpolation
   * @details Uses cubic spline interpolatoin for data type variables based on a 
   * vector of spline coefficients, or nodes, and independent points.  
   * The nodes are rescaled to 0-1.  This function does not extrapolate beyond the 
   * independent points.
   * 
   * @param spline_coffs a data vector of spline coefficients (nodes)
   * @param la a vector of independent points for use in interpolation.
   * 
   * @return A data vector containing the interpolated points.
   */
  
FUNCTION dvector cubic_spline(const dvector& spline_coffs, const dvector& la)
  {
	/*interplolation for length-based selectivity coefficeients*/
	//RETURN_ARRAYS_INCREMENT();
	int nodes=size_count(spline_coffs);
	dvector ia(1,nodes);
	ia.fill_seqadd(0,1./(nodes-1));
	dvector fa = (la-min(la))/(max(la)-min(la));
	vcubic_spline_function ffa(ia,spline_coffs);
	//RETURN_ARRAYS_DECREMENT();
	return(value(ffa(fa)));
	//return(1.0*la);
  }

// FUNCTION dvar_matrix cubic_spline_matrix(const dvar_matrix& spline_coffs)
//   {
// 	RETURN_ARRAYS_INCREMENT();
// 	int nodes= spline_coffs.colmax()-spline_coffs.colmin()+1;
// 	int rmin = spline_coffs.rowmin();
// 	int rmax = spline_coffs.rowmax();
	
// 	dvector ia(1,nodes);
// 	dvector fa(sage,nage);
// 	ia.fill_seqadd(0,1./(nodes-1));
// 	//fa.fill_seqadd(sage,1);
// 	fa.fill_seqadd(0,1./(nage-sage));
// 	vcubic_spline_function_array fna(rmin,rmax,ia,spline_coffs);
// 	RETURN_ARRAYS_DECREMENT();
// 	return(fna(fa));
	
//   }


  	/**
  	Purpose: This function loops over each of ngears and calculates the corresponding
  	         selectivity coefficients for that gear in each year.  It uses a switch 
  	         statement based on isel_type to determine which selectivty function to use
  	         for each particular gear that is specified in the control file.  See NOTES
  	         below for more information on selectivity models.

  	Author: Steven Martell
  	
  	Arguments:
  		isel_type -> an ivector with integers that determine what selectivity model to use.
  	
  	NOTES:
  		- The following is a list of the current selectivity models that are implemented:
		1)  Logistic selectivity with 2 parameters.
		2)  Age-specific selectivity coefficients with (nage-sage) parameters.
		    and the last two age-classes are assumed to have the same selectivity.
		3)  A reduced age-specific parameter set based on a bicubic spline.
		4)  Time varying cubic spline.
		5)  Time varying bicubic spline (2d version).
		6)  Fixed logistic.
		7)  Logistic selectivity based on relative changes in mean weight at age
		8)  Time varying selectivity based on logistic with deviations in 
		    weights at age (3 estimated parameters).
		11) Logistic selectivity with 2 parameters based on mean length.
		12) Length-based selectivity using cubic spline interpolation.
  		
  		- The bicubic_spline function is located in stats.cxx library.
  	
  	TODO list:
  	[*] add an option for length-based selectivity.  Use inverse of
		allometric relationship w_a = a*l_a^b; to get mean length-at-age from
		empirical weight-at-age data, then calculate selectivity based on 
		mean length. IMPLEMENTED IN CASE 11

	[*] change index for gear loop from j to k, and be consistent with year (i) and
	    age (j), and sex (h) indexing.

  	*/
FUNCTION void calcSelectivities(const ivector& isel_type)
  {
	

	int ig,i,j,k,byr,bpar,kgear;
	double tiny=1.e-10;
	dvariable p1,p2,p3;
	dvar_vector age_dev=age;
	dvar_matrix t1;
	dvar_matrix   tmp(syr,nyr-1,sage,nage);
	dvar_matrix  tmp2(syr,nyr,sage,nage);
	dvar_matrix ttmp2(sage,nage,syr,nyr);
	
	// Selex cSelex(age);
	// logistic_selectivity cLogisticSelex(age);
	log_sel.initialize();

	for(kgear=1; kgear<=ngear; kgear++)
	{
		// The following is used to mirror another gear-type
		// based on the absolute value of sel_phz.
		k  = kgear;
		if(sel_phz(k) < 0)
		{
			k = abs(sel_phz(kgear));
		}

		for( ig = 1; ig <= n_ags; ig++ )
		{
			tmp.initialize(); tmp2.initialize();
			dvector iy(1,yr_nodes(k));
			dvector ia(1,age_nodes(k));
			byr  = 1;
			bpar = 0; 
			switch(isel_type(k))
			{
				case 1: //logistic selectivity (2 parameters)
					for(i=syr; i<=nyr; i++)
					{
						if( i == sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < n_sel_blocks(k) ) byr++;
						}

						// cout<<"Testing selex class"<<endl;
						// log_sel(k)(ig)(i) = log( cSelex.logistic(sel_par(k)(bpar)) );
						// log_sel(k)(ig)(i) = log( cLogisticSelex(sel_par(k)(bpar)) );
						p1 = mfexp(sel_par(k,bpar,1));
						p2 = mfexp(sel_par(k,bpar,2));
						log_sel(kgear)(ig)(i) = log( plogis<dvar_vector>(age,p1,p2)+tiny );
					}
					break;
				
				case 6:	// fixed logistic selectivity
					p1 = mfexp(sel_par(k,1,1));
					p2 = mfexp(sel_par(k,1,2));
					for(i=syr; i<=nyr; i++)
					{
						log_sel(kgear)(ig)(i) = log( plogis<dvar_vector>(age,p1,p2) );
						// log_sel(k)(ig)(i) = log( cLogisticSelex(sel_par(k)(1)) );
					}
					break;
					
				case 2:	// age-specific selectivity coefficients
					for(i=syr; i<=nyr; i++)
					{
						if( i == sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < n_sel_blocks(k) ) byr++;
						}
						for(j=sage;j<=nage-1;j++)
						{
							log_sel(k)(ig)(i)(j)   = sel_par(k)(bpar)(j-sage+1);
						}
						log_sel(kgear)(ig)(i,nage) = log_sel(k)(ig)(i,nage-1);
					}
					break;
					
				case 3:	// cubic spline 
					for(i=syr; i<nyr; i++)
					{
						if( i==sel_blocks(k,byr) )
						{
							bpar ++;	
							log_sel(k)(ig)(i)=cubic_spline( sel_par(k)(bpar) );
							if( byr < n_sel_blocks(k) ) byr++;
						}
						log_sel(kgear)(ig)(i+1) = log_sel(k)(ig)(i);
					}
					break;
					
				case 4:	// time-varying cubic spline every year				
					for(i=syr; i<=nyr; i++)
					{
						log_sel(kgear)(ig)(i) = cubic_spline(sel_par(k)(i-syr+1));
					}
					break;
					
				case 5:	// time-varying bicubic spline
					ia.fill_seqadd( 0,1./(age_nodes(k)-1) );
					iy.fill_seqadd( 0,1./( yr_nodes(k)-1) );	
					bicubic_spline( iy,ia,sel_par(k),tmp2 );
					log_sel(kgear)(ig) = tmp2; 
					break;
					
				case 7:
					// time-varying selectivity based on deviations in weight-at-age
					// CHANGED This is not working and should not be used. (May 5, 2011)
					// SkDM:  I was not able to get this to run very well.
					// AUG 5, CHANGED so it no longer has the random walk component.
					p1 = mfexp(sel_par(k,1,1));
					p2 = mfexp(sel_par(k,1,2));
					
					for(i = syr; i<=nyr; i++)
					{
						dvar_vector tmpwt=log(d3_wt_avg(ig)(i)*1000)/mean(log(d3_wt_avg(ig)*1000.));
						log_sel(kgear)(ig)(i) = log( plogis(tmpwt,p1,p2)+tiny );
					}	 
					break;
					
				case 8:
					//Alternative time-varying selectivity based on weight 
					//deviations (d3_wt_dev) d3_wt_dev is a matrix(syr,nyr+1,sage,nage)
					//p3 is the coefficient that describes variation in log_sel.
					p1 = mfexp(sel_par(k,1,1));
					p2 = mfexp(sel_par(k,1,2));
					p3 = sel_par(k,1,3);
					
					for(i=syr; i<=nyr; i++)
					{
						tmp2(i) = p3*d3_wt_dev(ig)(i);
						log_sel(kgear)(ig)(i) = log( plogis<dvar_vector>(age,p1,p2)+tiny ) + tmp2(i);
					}
					break;
					
				case 11: // logistic selectivity based on mean length-at-age
					for(i=syr; i<=nyr; i++)
					{
						if( i == sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < n_sel_blocks(k) ) byr++;
						}
						p1 = mfexp(sel_par(k,bpar,1));
						p2 = mfexp(sel_par(k,bpar,2));

						dvector len = pow(d3_wt_avg(ig)(i)/d_a(ig),1./d_b(ig));

						log_sel(kgear)(ig)(i) = log( plogis<dvar_vector>(len,p1,p2) );
						//log_sel(kgear)(ig)(i) = log( plogis(len,p1,p2) );
					}	
					break;
					
				case 12: // cubic spline length-based coefficients.
					for(i=syr; i<=nyr; i++)
					{
						if( i == sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < n_sel_blocks(k) ) byr++;
						}
					
						dvector len = pow(d3_wt_avg(ig)(i)/d_a(ig),1./d_b(ig));
						log_sel(kgear)(ig)(i)=cubic_spline( sel_par(k)(bpar), len );
					}
					break;
					
					
				default:
					log_sel(kgear)(ig)=0;
					break;
					
			}  // switch

			//subtract mean to ensure mean(exp(log_sel))==1
			for(i=syr;i<=nyr;i++)
			{
				log_sel(kgear)(ig)(i) -= log( mean(mfexp(log_sel(kgear)(ig)(i))) );
				// log_sel(k)(ig)(i) -= log( max(mfexp(log_sel(k)(ig)(i))) );
			}
		}
	}  //end of gear k
	
	if(verbose)cout<<"**** Ok after calcSelectivities ****"<<endl;
	
  }	

  	
	
  	/**
  	Purpose: This function calculates fishing mortality, total mortality and annual
  	         surivival rates S=exp(-Z) for each age and year based on fishing mortality
  	         and selectivity coefficients.  Z also is updated with time-varying 
  	         natural mortality rates if specificed by user.
  	Author: Steven Martell
  	
  	Arguments:
  		None
  	
  	NOTES:
  		- Jan 5, 2012 Added catch_type to allow for catch in numbers, weight or spawn.
          In the case of spawn on kelp (roe fisheries), the Fishing mortality does not
          occur on the adult component.  
        - Added if(catch_type(k)!=3) //exclude roe fisheries
  		- F(group,year,age)
  		- Exclude type = 3, roe fisheries harvesting eggs only, not adults.
		- if dCatchData$sex is male & female combined, then allocate to both sexes.

  	TODO list:
  	[*] Dec 24, 2010.  Adding time-varying natural mortality.
  	[*] May 20, 2011.  Add cubic spline to the time-varying natural mortality.
	[ ] Calculate average M for reference point calculations based on pfc file.
	[ ] 
  	*/
FUNCTION calcTotalMortality
  {

	int ig,ii,i,k,l;
	dvariable ftmp;
	F.initialize(); 
	ft.initialize();
	
	
	// |---------------------------------------------------------------------------------|
	// | FISHING MORTALITY
	// |---------------------------------------------------------------------------------|
	// |

	for(ig=1;ig<=nCtNobs;ig++)
	{
		i  = dCatchData(ig)(1);	 //year
		k  = dCatchData(ig)(2);  //gear
		f  = dCatchData(ig)(3);  //area
		g  = dCatchData(ig)(4);  //group
		h  = dCatchData(ig)(5);  //sex
		l  = dCatchData(ig)(6);  //type

		if( i > nyr ) continue;
		if( h )
		{
			ii = pntr_ags(f,g,h);    
			ftmp = mfexp(log_ft_pars(ig));
			ft(ii)(k,i) = ftmp;
			if( l != 3 )
			{
				F(ii)(i) += ftmp*mfexp(log_sel(k)(ii)(i));
			}
		}
		else if( !h ) // h=0 case for asexual catch
		{
			for(h=1;h<=nsex;h++)
			{
				ii = pntr_ags(f,g,h);    
				ftmp = mfexp(log_ft_pars(ig));
				ft(ii)(k,i) = ftmp;
				if( l != 3 )
				{
					F(ii)(i) += ftmp*mfexp(log_sel(k)(ii)(i));
				}		
			}
		}
	}

	// |---------------------------------------------------------------------------------|
	// | NATURAL MORTALITY
	// |---------------------------------------------------------------------------------|
	// | - uses cubic spline to interpolate time-varying natural mortality
	M.initialize();
	log_m_devs.initialize();
	for(ig=1;ig<=n_ags;ig++)
	{
		g = n_group(ig);
		h = n_sex(ig);
		M(ig) = m( pntr_gs(g,h) );
		if( active( log_m_nodes) )
		{
			int nodes = size_count(log_m_nodes);
			dvector im(1,nodes);
			dvector fm(syr+1,nyr);
			im.fill_seqadd(0,1./(nodes-1));
			fm.fill_seqadd(0,1./(nyr-syr));
			vcubic_spline_function m_spline(im,log_m_nodes);
			log_m_devs = m_spline( fm );
		}

		for(i=syr+1; i<=nyr; i++)
		{
			M(ig)(i) = M(ig)(i-1) * mfexp(log_m_devs(i));
		}
		// TODO fix for reference point calculations
		// m_bar = mean( M_tot.sub(pf_cntrl(1),pf_cntrl(2)) );
	}

	// |---------------------------------------------------------------------------------|
	// | TOTAL MORTALITY
	// |---------------------------------------------------------------------------------|
	// |
	for(ig=1;ig<=n_ags;ig++)
	{
		Z(ig) = M(ig) + F(ig);
		S(ig) = mfexp(-Z(ig));
	}

	if(verbose) cout<<"**** OK after calcTotalMortality ****"<<endl;	
  }
	
	
  	/**
  	Purpose: This function initializes the numbers-at-age matrix in syr
  	         based on log_rinit and log_init_rec_devs, the annual recruitment
  	         based on log_rbar and log_rec_devs, and updates the number-at-age
  	         over time based on the survival rate calculated in calcTotalMortality.

  	Author: Steven Martell
  	
  	Arguments:
  		None
  	
  	NOTES:
		- Aug 9, 2012.  Made a change here to initialize the numbers
		  at age in syr using the natural mortality rate at age in syr. 
		  Prior to this the average (m_bar) rate was used, since this 
		  has now changed with new projection control files.  Should only
		  affect models that were using time varying natural mortality.
  		- d_iscamCntrl(5) is a flag to start at unfished conditions, so set N(syr,sage) = ro
  	
  	TODO list:
  	[ ] - Restrict log_avgrec and rec_devs to area and group dimensions (remove sex).
  	[ ] - Initialize from unfished conditions (d_iscamCntrl 5 flag is true then rt(syr) = ro)
  	*/
FUNCTION calcNumbersAtAge
  {
	int ig,ih;

	N.initialize();
	bt.initialize();
	for(ig=1;ig<=n_ags;ig++)
	{
		f  = n_area(ig);
		g  = n_group(ig);
		ih = pntr_ag(f,g);

		dvar_vector lx(sage,nage);
		dvar_vector tr(sage,nage);
		lx(sage) = 1.0;
		for(j=sage;j< nage;j++)
		{
			lx(j+1) = lx(j) * exp( -M(ig)(syr)(j) );
		}
		lx(nage) /= (1.-exp(-M(ig)(syr,nage)));
		
		if( d_iscamCntrl(5) ) // initialize at unfished conditions.
		{
			tr =  log( ro(g) ) + log(lx);
		}
		else if ( !d_iscamCntrl(5) )
		{
			tr(sage)        = ( log_avgrec(ih)+log_rec_devs(ih)(syr));
			tr(sage+1,nage) = (log_recinit(ih)+init_log_rec_devs(ih));
			tr(sage+1,nage) = tr(sage+1,nage)+log(lx(sage+1,nage));
		}
		N(ig)(syr)(sage,nage) = 1./nsex * mfexp(tr);
		log_rt(ih)(syr-nage+sage,syr) = tr.shift(syr-nage+sage);


		for(i=syr;i<=nyr;i++)
		{
			if( i>syr )
			{
				log_rt(ih)(i) = (log_avgrec(ih)+log_rec_devs(ih)(i));
				N(ig)(i,sage) = 1./nsex * mfexp( log_rt(ih)(i) );				
			}

			N(ig)(i+1)(sage+1,nage) =++elem_prod(N(ig)(i)(sage,nage-1)
			                                     ,S(ig)(i)(sage,nage-1));
			N(ig)(i+1,nage)        +=  N(ig)(i,nage)*S(ig)(i,nage);

			// average biomass for group in year i
			bt(g)(i) += N(ig)(i) * d3_wt_avg(ig)(i);
		}
		N(ig)(nyr+1,sage) = 1./nsex * mfexp( log_avgrec(ih));
		bt(g)(nyr+1) += N(ig)(nyr+1) * d3_wt_avg(ig)(nyr+1);

	}
	if(verbose)cout<<"**** Ok after calcNumbersAtAge ****"<<endl;	
  }





  	/**
  	Purpose:  This function calculates the predicted age-composition samples (A) for 
  	          both directed commercial fisheries and survey age-composition data. For 
  	          all years of data specified in the A matrix, calculated the predicted 
  	          proportions-at-age in the sampled catch-at-age.  If no catch-age data exist
  	          for a particular year i, for gear k (i.e. no directed fishery or from a 
  	          survey sample process that does not have an appreciable F), the calculate 
  	          the predicted proportion based on log(N) + log_sel(group,gear,year)
  	Author: Steven Martell
  	
  	Arguments:
  		None
  	
  	NOTES:
  		- Adapted from iSCAM 1.5.  
  		- No longer using ragged arrays for gear, the ragged matrix is indexed by:
  		  year gear area, group, sex | age columns
  		- For the residuals, note that each gear is weigthed by the conditional MLE
  		  of the variance.
  	
  	TODO list:
  	[x] - Merge redundant code from calcCatchAtAge
  	[*] - Add case where Chat data do not exsist.
	[x] - Calculate residuals A_nu; gets done automatically in dmvlogistic
	[?] - add plus group if n_A_nage < nage;  Aug 7, 2013

  	*/
  	
FUNCTION calcAgeComposition
  {
  	int ii,ig,kk;
  	dvar_vector va(sage,nage);
  	dvar_vector fa(sage,nage);
  	dvar_vector sa(sage,nage);
  	dvar_vector za(sage,nage);
  	dvar_vector ca(sage,nage);
  	dvar_vector na(sage,nage);
  	A_hat.initialize();

  	 for(kk=1;kk<=nAgears;kk++)
  	 {
  	 	for(ii=1;ii<=n_A_nobs(kk);ii++)
  	 	{
	  		i = d3_A(kk)(ii)(n_A_sage(kk)-5);
	  		k = d3_A(kk)(ii)(n_A_sage(kk)-4);
	  		f = d3_A(kk)(ii)(n_A_sage(kk)-3);
	  		g = d3_A(kk)(ii)(n_A_sage(kk)-2);
	  		h = d3_A(kk)(ii)(n_A_sage(kk)-1);
	  		// | trap for retrospecitve analysis.
	  		if(i > nyr) continue;

	  		if( h )  // age comps are sexed (h > 0)
	  		{
				ig = pntr_ags(f,g,h);
				va = mfexp(log_sel(k)(ig)(i));
				za = Z(ig)(i);
				sa = S(ig)(i);
				na = N(ig)(i);
				if( ft(ig)(k)(i)==0 )
				{
					ca = elem_prod(elem_prod(na,va),0.5*sa);
				}
				else
				{
					fa = ft(ig)(k)(i) * va;
					ca = elem_prod(elem_prod(elem_div(fa,za),1.-sa),na);					
				}
				A_hat(kk)(ii) = ca(n_A_sage(kk),n_A_nage(kk));

				// | +group if n_A_nage(kk) < nage
				if( n_A_nage(kk) < nage )
				{
					A_hat(kk)(ii)(n_A_nage(kk)) += sum( ca(n_A_nage(kk)+1,nage) );
				}
	  		}
	  		else if( !h )  // age-comps are unsexed
	  		{
	  			for(h=1;h<=nsex;h++)
	  			{
					ig = pntr_ags(f,g,h);
					va = mfexp(log_sel(k)(ig)(i));
					za = Z(ig)(i);
					sa = S(ig)(i);
					na = N(ig)(i);
					if( ft(ig)(k)(i)==0 )
					{
						ca = elem_prod(na,0.5*sa);
					}
					else
					{
						fa = ft(ig)(k)(i) * va;
						ca = elem_prod(elem_prod(elem_div(fa,za),1.-sa),na);					
					}
					A_hat(kk)(ii) += ca(n_A_sage(kk),n_A_nage(kk));

					// | +group if n_A_nage(kk) < nage
					if( n_A_nage(kk) < nage )
					{
						A_hat(kk)(ii)(n_A_nage(kk)) += sum( ca(n_A_nage(kk)+1,nage) );
					}
		  		}
	  		}
	  		A_hat(kk)(ii) /= sum( A_hat(kk)(ii) );
  	 	}
  	}
  	
	if(verbose)cout<<"**** Ok after calcAgeComposition ****"<<endl;

  }	

FUNCTION calcLengthComposition
  {
  	/*
	Purpose:  This function calculates the predicted length-composition samples (B) for 
  	          both directed commercial fisheries and survey age-composition data. For 
  	          all years of data specified in the B matrix, calculated the predicted 
  	          proportions-at-length in the sampled catch-at-length.  If no catch-length 
  	          data exist for a particular year i, for gear k (i.e. no directed fishery or 
  	          from a survey sample process that does not have an appreciable F), then 
  	          calculate the predicted proportion based on:
  	          							N(i) * sel(group,gear,year) * ALK(i)
  	          where ALK is the Age-Length Transition Key.
  	Author: Steven Martell

  	Arguments:  none

  	NOTES: 	To calculate the length composition, the catch-at-age vector is multiplied by
  			an Age_length transition key.  In theory, A_hat from the calcAgeComposition 
  			routine can be multiplied by the ALK.

  	TODO List:
  	[ ] - write this function. and Let Robyn Forrest know about it.
  	*/

  	int ii,ig,kk;
  	dvar_vector va(sage,nage);
  	dvar_vector fa(sage,nage);
  	dvar_vector sa(sage,nage);
  	dvar_vector za(sage,nage);
  	dvar_vector ca(sage,nage);
  	dvar_vector na(sage,nage);
  	A_hat.initialize();

  	 for(kk=1;kk<=nAgears;kk++)
  	 {
  	 	for(ii=1;ii<=n_A_nobs(kk);ii++)
  	 	{
	  		i = d3_A(kk)(ii)(n_A_sage(kk)-5);
	  		k = d3_A(kk)(ii)(n_A_sage(kk)-4);
	  		f = d3_A(kk)(ii)(n_A_sage(kk)-3);
	  		g = d3_A(kk)(ii)(n_A_sage(kk)-2);
	  		h = d3_A(kk)(ii)(n_A_sage(kk)-1);
	  		
	  		// | trap for retrospecitve analysis.
	  		if(i > nyr) continue;

	  		if( h )
	  		{
				ig = pntr_ags(f,g,h);
				va = mfexp(log_sel(k)(ig)(i));
				za = Z(ig)(i);
				sa = S(ig)(i);
				na = N(ig)(i);
				if( ft(ig)(k)(i)==0 )
				{
					ca = elem_prod(na,0.5*sa);
				}
				else
				{
					fa = ft(ig)(k)(i) * va;
					ca = elem_prod(elem_prod(elem_div(fa,za),1.-sa),na);					
				}
				A_hat(kk)(ii) = ca(n_A_sage(kk),n_A_nage(kk));

				// | +group if n_A_nage(kk) < nage
				if( n_A_nage(kk) < nage )
				{
					A_hat(kk)(ii)(n_A_nage(kk)) += sum( ca(n_A_nage(kk)+1,nage) );
				}
	  		}
	  		else if( !h )
	  		{
	  			for(h=1;h<=nsex;h++)
	  			{
					ig = pntr_ags(f,g,h);
					va = mfexp(log_sel(k)(ig)(i));
					za = Z(ig)(i);
					sa = S(ig)(i);
					na = N(ig)(i);
					if( ft(ig)(k)(i)==0 )
					{
						ca = elem_prod(na,0.5*sa);
					}
					else
					{
						fa = ft(ig)(k)(i) * va;
						ca = elem_prod(elem_prod(elem_div(fa,za),1.-sa),na);					
					}
					A_hat(kk)(ii) += ca(n_A_sage(kk),n_A_nage(kk));

					// | +group if n_A_nage(kk) < nage
					if( n_A_nage(kk) < nage )
					{
						A_hat(kk)(ii)(n_A_nage(kk)) += sum( ca(n_A_nage(kk)+1,nage) );
					}
		  		}
	  		}
	  		A_hat(kk)(ii) /= sum( A_hat(kk)(ii) );
  	 	}
  	}
  	
	if(verbose)cout<<"**** Ok after calcLengthComposition ****"<<endl;

  }

FUNCTION calcTotalCatch
  {
  	/*
  	Purpose:  This function calculates the total catch.  
  	Dependencies: Must call calcCatchAtAge function first.
  	Author: Steven Martell
  	
  	Arguments:
  		None
  	
  	NOTES:
  		
  	
  	TODO list:
  	[ ] get rid of the obs_ct, ct, eta array structures, inefficient, better to use
  	    a matrix, then cbind the predicted catch and residuals for report. (ie. an R
  	    data.frame structure and use melt to ggplot for efficient plots.)
  	*/
  	int ii,l,ig;
  	double d_ct;

  	ct.initialize();
  	eta.initialize();
  	
  	dvar_vector     fa(sage,nage);
  	dvar_vector     ca(sage,nage);
  	dvar_vector     sa(sage,nage);
  	dvar_vector     za(sage,nage);
  	
  	

  	for(ii=1;ii<=nCtNobs;ii++)
	{
		i    = dCatchData(ii,1);
		k    = dCatchData(ii,2);
		f    = dCatchData(ii,3);
		g    = dCatchData(ii,4);
		h    = dCatchData(ii,5);
		l    = dCatchData(ii,6);
		d_ct = dCatchData(ii,7);
  		
  		// | trap for retro year
  		if( i>nyr ) continue;


		switch(l)
		{
			case 1:  // catch in weight
				if( h )
				{
					ig     = pntr_ags(f,g,h);
					fa     = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
					za     = Z(ig)(i);
					sa     = S(ig)(i);
					ca     = elem_prod(elem_prod(elem_div(fa,za),1.-sa),N(ig)(i));
					ct(ii) = ca * d3_wt_avg(ig)(i);
				}
				else if( !h )
				{
					for(h=1;h<=nsex;h++)
					{
						ig     = pntr_ags(f,g,h);
						fa     = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
						za     = Z(ig)(i);
						sa     = S(ig)(i);
						ca     = elem_prod(elem_prod(elem_div(fa,za),1.-sa),N(ig)(i));
						ct(ii)+= ca * d3_wt_avg(ig)(i);		
					}
				}
			break;

			case 2:  // catch in numbers
				if( h )
				{
					ig     = pntr_ags(f,g,h);
					fa     = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
					za     = Z(ig)(i);
					sa     = S(ig)(i);
					ca     = elem_prod(elem_prod(elem_div(fa,za),1.-sa),N(ig)(i));
					ct(ii) = sum( ca );
				}
				else if( !h )
				{
					for(h=1;h<=nsex;h++)
					{
						ig     = pntr_ags(f,g,h);
						fa     = ft(ig)(k)(i) * mfexp(log_sel(k)(ig)(i));
						za     = Z(ig)(i);
						sa     = S(ig)(i);
						ca     = elem_prod(elem_prod(elem_div(fa,za),1.-sa),N(ig)(i));
						ct(ii)+= sum( ca );
					}
				}
			break;

			case 3:  // roe fisheries, special case
				if( h )
				{
					ig            = pntr_ags(f,g,h);
					dvariable ssb = N(ig)(i) * d3_wt_mat(ig)(i);
					ct(ii)        = (1.-exp(-ft(ig)(k)(i))) * ssb;
				}
				else if( !h )
				{
					for(h=1;h<=nsex;h++)
					{
						ig            = pntr_ags(f,g,h);
						dvariable ssb = N(ig)(i) * d3_wt_mat(ig)(i);
						ct(ii)       += (1.-exp(-ft(ig)(k)(i))) * ssb;
					}
				}
			break;
		}	// end of switch

		// | catch residual
		eta(ii) = log(d_ct+TINY) - log(ct(ii)+TINY);
	}
	if(verbose)cout<<"**** Ok after calcTotalCatch ****"<<endl;
  }


	
  	/**
  	Purpose:  This function computes the mle for survey q, calculates the survey 
  	          residuals (epsilon).
  	Author: Steven Martell
  	
  	Arguments:
  		None
  	
  	NOTES:
  		- Oct 31, 2010, added retrospective counter.
  		- Nov 22, 2010, adding multiple surveys. 
  		Still need to check with retrospective option
  		- Nov 30, 2010, adjust the suvery biomass by the fraction of Z that has occurred 
  		when the survey was conducted. For herring spawning biomass this would be 
  		after the fishery has taken place.
  		- Dec 6, 2010, modified predicted survey biomass to accomodate empirical
  		weight-at-age data (d3_wt_avg).
  		- May 11, 2011.  Vivian Haist pointed out an error in survey biomass comparison.
  		The spawning biomass was not properly calculated in this routine. I.e. its 
  		different than the spawning biomass in the stock-recruitment routine. (Based on 
  		fecundity which changes with time when given empirical weight-at-age data.)
  		- Jan 6, 2012.  CHANGED corrected spawn survey observations to include a roe 
  		fishery that would remove potential spawn that would not be surveyed.
  		- d3_survey_data: (iyr index(it) gear area group sex wt timing)
  		- for MLE of survey q, using weighted mean of zt to calculate q.

  	TODO list:
  	    [?] - add capability to accomodate priors for survey q's.
  	    [ ] - verify q_prior=2 option for random walk in q.
  	    [ ] - For sel_type==3, may need to reduce abundance by F on spawning biomass (herring)
 
  	TODO LIST:
	  [ ] - add capability to accompodate priors for survey catchabiliyt coefficients.

  */
FUNCTION calcSurveyObservations
  {
	
	int ii,kk,ig,nz;
	double di;
	dvariable ftmp;
	dvar_vector Na(sage,nage);
	dvar_vector va(sage,nage);
	dvar_vector sa(sage,nage);
	epsilon.initialize();
	it_hat.initialize();

	for(kk=1;kk<=nItNobs;kk++)
	{
		// | Vulnerable number-at-age to survey.
		dvar_matrix V(1,n_it_nobs(kk),sage,nage);
		V.initialize();
		nz = 0;
		for(ii=1;ii<=n_it_nobs(kk);ii++)
		{
			i    = d3_survey_data(kk)(ii)(1);
			k    = d3_survey_data(kk)(ii)(3);
			f    = d3_survey_data(kk)(ii)(4);
			g    = d3_survey_data(kk)(ii)(5);
			h    = d3_survey_data(kk)(ii)(6);
			di   = d3_survey_data(kk)(ii)(8);

			// | trap for retrospective nyr change
			if( i > nyr ) continue;

			nz ++;

			// h ==0?h=1:NULL;
			Na.initialize();
			for(h=1;h<=nsex;h++)
			{
				ig  = pntr_ags(f,g,h);
				va  = mfexp( log_sel(k)(ig)(i) );
				sa  = mfexp( -Z(ig)(i)*di );
				Na  = elem_prod(N(ig)(i),sa);
				switch(n_survey_type(kk))
				{
					case 1:
						V(ii) += elem_prod(Na,va);
					break; 
					case 2:
						V(ii) += elem_prod(elem_prod(Na,va),d3_wt_avg(ig)(i));
					break;
					case 3:
						V(ii) += elem_prod(Na,d3_wt_mat(ig)(i));
					break;
				}
			}
		
		} // end of ii loop
		dvector     it = trans(d3_survey_data(kk))(2)(1,nz);
		dvector     wt = trans(d3_survey_data(kk))(7)(1,nz);
		            wt = wt/sum(wt);
		dvar_vector t1 = rowsum(V);
		dvar_vector zt = log(it) - log(t1(1,nz));
		dvariable zbar = sum(elem_prod(zt,wt));
				 q(kk) = mfexp(zbar);

		// | survey residuals
		epsilon(kk).sub(1,nz) = zt - zbar;
		 it_hat(kk).sub(1,nz) = q(kk) * t1(1,nz);

		// | SPECIAL CASE: penalized random walk in q.
		if( q_prior(kk)==2 )
		{
			epsilon(kk).initialize();
			dvar_vector fd_zt     = first_difference(zt);
			dvariable  zw_bar     = sum(elem_prod(fd_zt,wt(1,nz-1)));
			epsilon(kk).sub(1,nz-1) = fd_zt - zw_bar;
			qt(kk)(1) = exp(zt(1));
			for(ii=2;ii<=nz;ii++)
			{
				qt(kk)(ii) = qt(kk)(ii-1) * exp(fd_zt(ii-1));
			}
			it_hat(kk).sub(1,nz) = elem_prod(qt(kk)(1,nz),t1(1,nz));
		}
	}
	if(verbose)cout<<"**** Ok after calcSurveyObservations ****"<<endl;
	
  }

  
  	/**
	Purpose:  
		This function is used to derive the underlying stock-recruitment 
		relationship that is ultimately used in determining MSY-based reference 
		points.  The objective of this function is to determine the appropriate 
		Ro, Bo and steepness values of either the Beverton-Holt or Ricker  Stock-
		Recruitment Model:

		Beverton-Holt Model
		\f$ Rt=k*Ro*St/(Bo+(k-1)*St)*exp(delta-0.5*tau*tau) \f$

		Ricker Model
		\f$ Rt=so*St*exp(-beta*St)*exp(delta-0.5*tau*tau) \f$
			
		The definition of a stock is based on group only. At this point, spawning biomass
		from all areas for a given group is the assumed stock, and the resulting
		recruitment is compared with the estimated recruits|group for all areas.

  	Author: Steven Martell
  	
  	Arguments:
  		None
  	
  	NOTES:
		Psuedocode:
		-1) Get average natural mortality rate at age.
		-2) Calculate survivorship to time of spawning.
		-3) Calculate unfished spawning biomass per recruit.
		-4) Compute spawning biomass vector & substract roe fishery
		-5) Project spawning biomass to nyr+1 under natural mortality.
		-6) Calculate stock recruitment parameters (so, beta);
		-7) Calculate predicted recruitment
		-8) Compute residuals from estimated recruitments.
  			
  	
  	TODO list:
	  [?] - Change step 3 to be a weighted average of spawning biomass per recruit by area.
	  [?] - Increase dimensionality of ro, sbo, so, beta, and steepness to ngroup.
	  [?] - Add autocorrelation in recruitment residuals with parameter \f$ \gamma_r \f$.

  	*/
FUNCTION void calcStockRecruitment()
  {

  	int ig,ih;
  	rt.initialize();
  	sbt.initialize();
  	delta.initialize();
	
	dvariable phib;//,so,beta;
	dvector         fa(sage,nage);
	dvar_vector   stmp(sage,nage);
	dvar_vector     ma(sage,nage);
	dvar_vector tmp_rt(syr+sage,nyr);
	dvar_vector     lx(sage,nage); 
	dvar_vector     lw(sage,nage); 
	
	
	for(g=1;g<=ngroup;g++)
	{
		lx.initialize();
		lw.initialize();
		lx(sage) = 1.0;
		lw(sage) = 1.0;
		phib = 0;
		for(f=1;f<=narea;f++)
		{
			for(h=1;h<=nsex;h++)
			{
				ig = pntr_ags(f,g,h);

				// | Step 1. average natural mortality rate at age.
				// | Step 2. calculate survivorship
				for(j=sage;j<=nage;j++)
				{
					ma(j) = mean(trans(M(ig))(j));
					fa(j) = mean( trans(d3_wt_mat(ig))(j) );
					if(j > sage)
					{
						lx(j) = lx(j-1) * mfexp(-ma(j-1));
					}
					lw(j) = lx(j) * mfexp(-ma(j)*d_iscamCntrl(13));
				}
				lx(nage) /= 1.0 - mfexp(-ma(nage));
				lw(nage) /= 1.0 - mfexp(-ma(nage));
				
				// | Step 3. calculate average spawing biomass per recruit.
				phib += 1./(narea*nsex) * lw*fa;

				// | Step 4. compute spawning biomass at time of spawning.
				for(i=syr;i<=nyr;i++)
				{
					stmp      = mfexp(-Z(ig)(i)*d_iscamCntrl(13));
					sbt(g,i) += elem_prod(N(ig)(i),d3_wt_mat(ig)(i)) * stmp;
				}

				// | Step 5. spawning biomass projection under natural mortality only.
				stmp          = mfexp(-M(ig)(nyr)*d_iscamCntrl(13));
				sbt(g,nyr+1) += elem_prod(N(ig)(nyr),d3_wt_mat(ig)(i)) * stmp;
			}

			// | Estimated recruits
			ih     = pntr_ag(f,g);
			rt(g) += mfexp(log_rt(ih)(syr+sage,nyr));
		}

		// | Step 6. calculate stock recruitment parameters (so, beta, sbo);
		so(g)  = kappa(g)/phib;
		sbo(g) = ro(g) * phib;

		// | Step 7. calculate predicted recruitment.
		dvar_vector tmp_st = sbt(g)(syr,nyr-sage).shift(syr+sage);
		switch(int(d_iscamCntrl(2)))
		{
			case 1:  // | Beverton Holt model
				beta(g)   = (kappa(g)-1.)/sbo(g);
				tmp_rt    = elem_div(so(g)*tmp_st,1.+beta(g)*tmp_st);
			break;

			case 2:  // | Ricker model
				beta(g)   = log(kappa(g))/sbo(g);
				tmp_rt    = elem_prod(so(g)*tmp_st,exp(-beta(g)*tmp_st));
			break;
		}
		
		// | Step 8. // residuals in stock-recruitment curve with gamma_r = 0
		delta(g) = log(rt(g))-log(tmp_rt)+0.5*tau*tau;

		// Autocorrelation in recruitment residuals.
		// if gamma_r > 0 then 
		if( active(gamma_r) )
		{
			int byr = syr+sage+1;
			delta(g)(byr,nyr) 	= log(rt(g)(byr,nyr)) 
									- (1.0-gamma_r)*log(tmp_rt(byr,nyr)) 
									- gamma_r*log(++rt(g)(byr-1,nyr-1))
									+ 0.5*tau*tau;			
		}


	}
	
	if(verbose)cout<<"**** Ok after calcStockRecruitment ****"<<endl;
	
  }
	
FUNCTION calcObjectiveFunction
  {
  	/*
  	Purpose:  This function computes the objective function that ADMB will minimize.
  	Author: Steven Martell
  	
  	Arguments:
  		None
  	
  	NOTES:
		There are several components to the objective function
		Likelihoods (nlvec):
			-1) likelihood of the catch data
			-2) likelihood of the survey abundance index
			-3) likelihood of age composition data 
			-4) likelihood for stock-recruitment relationship
			-5) penalized likelihood for fishery selectivities
			-6) penalized likelihood for fishery selectivities
  		
  	
  	TODO list:
	[*]	- Dec 20, 2010.  SJDM added prior to survey qs.
		  q_prior is an ivector with current options of 0 & 1 & 2.
		  0 is a uniform density (ignored) and 1 is a normal
		  prior density applied to log(q), and 2 is a random walk in q.
  	[ ] - Allow for annual sig_c values in catch data likelihood.
  	[ ] - Increase dimensionality of sig and tau to ngroup.
  	[ ] - Correct likelihood for cases when rho > 0 (Schnute & Richards, 1995)
  	*/



// 	int i,j,k;
// 	double o=1.e-10;

	nlvec.initialize();
	
	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR CATCH DATA
	// |---------------------------------------------------------------------------------|
	// | - This likelihood changes between phases n-1 and n:
	// | - Phase (n-1): standard deviation in the catch based on user input d_iscamCntrl(3)
	// | - Phase (n)  : standard deviation in the catch based on user input d_iscamCntrl(4)
	// | 

	double sig_c =d_iscamCntrl(3);
	if(last_phase())
	{
		sig_c=d_iscamCntrl(4);
	}
	if( active(log_ft_pars) )
	{
		nlvec(1) = dnorm(eta,0.0,sig_c);
	}


	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR RELATIVE ABUNDANCE INDICES
	// |---------------------------------------------------------------------------------|
	// | - sig_it     -> vector of standard deviations based on relative wt for survey.
	// |
	for(k=1;k<=nItNobs;k++)
	{
		dvar_vector sig_it = sig/it_wt(k);
		nlvec(2,k)=dnorm(epsilon(k),sig_it);
	}
	
	
	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD FOR AGE-COMPOSITION DATA
	// |---------------------------------------------------------------------------------|
	// | - Two options based on d_iscamCntrl(14):
	// | - 	1 -> multivariate logistic using conditional MLE of the variance for weight.
	// | -  2 -> multnomial, assumes input sample size as n in n log(p)
	// | -  Both likelihoods pool pmin (d_iscamCntrl(16)) into adjacent yearclass.
	// | -  PSEUDOCODE:
	// | -    => first determine appropriate dimensions for each of nAgears arrays (naa)
	// | -    => second extract sub arrays into obs (O) and predicted (P)
	// | -    => Compute either dmvlogistic, or dmultinom negative loglikehood.
	// | 
	// | TODO:
	// | [ ] - change A_nu to data-type variable, does not need to be differentiable.
	// |
	A_nu.initialize();
	for(k=1;k<=nAgears;k++)
	{	
		if( n_A_nobs(k)>0 )
		{
			int naa=0;
			int iyr;
			//retrospective counter
			for(i=1;i<=n_A_nobs(k);i++)
			{
				iyr = d3_A(k)(i)(n_A_sage(k)-5);	//index for year
				if(iyr<=nyr) naa++;
			}
			
			dmatrix     O = trans(trans(d3_A_obs(k)).sub(n_A_sage(k),n_A_nage(k))).sub(1,naa);
			dvar_matrix P = trans(trans(A_hat(k)).sub(n_A_sage(k),n_A_nage(k))).sub(1,naa);
			dvar_matrix nu(O.rowmin(),O.rowmax(),O.colmin(),O.colmax()); 
			nu.initialize();
			
			// | Choose form of the likelihood based on d_iscamCntrl(14) switch
			//switch(int(d_iscamCntrl(14)))
			logistic_normal cLN_Age( O,P,dMinP(k),dEps(k) );
			logistic_student_t cLST_Age( O,P,dMinP(k),dEps(k) );
			switch( int(nCompLikelihood(k)) )
			{
				case 1:
					nlvec(3,k) = dmvlogistic(O,P,nu,age_tau2(k),dMinP(k));
				break;
				case 2:
					nlvec(3,k) = dmultinom(O,P,nu,age_tau2(k),dMinP(k));
				break;
				case 3:
					if( !active(log_age_tau2(k)) )                 // LN1 Model
					{
						nlvec(3,k)  = cLN_Age();	
					}
					else
					{
						nlvec(3,k) = cLN_Age( exp(log_age_tau2(k)) );
					}

					// Residual
					if(last_phase())
					{
						nu          = cLN_Age.get_standardized_residuals();
						age_tau2(k) = cLN_Age.get_sigma2();
					}
				break;

				case 4:
					//logistic_normal cLN_Age( O,P,dMinP(k),dEps(k) );
					if( active(phi1(k)) && !active(phi2(k)) )  // LN2 Model
					{
						nlvec(3,k)   = cLN_Age(exp(log_age_tau2(k)),phi1(k));	
					}
					if( active(phi1(k)) && active(phi2(k)) )   // LN3 Model
					{
						nlvec(3,k)   = cLN_Age(exp(log_age_tau2(k)),phi1(k),phi2(k));	
					}

					// Residual
					if(last_phase())
					{
						nu          = cLN_Age.get_standardized_residuals();
						age_tau2(k) = cLN_Age.get_sigma2();
					}

				break;

				case 5: // Logistic-normal with student-t
					if( !active(log_degrees_of_freedom(k)) )
					{
						nlvec(3,k) = cLST_Age();
					}
					else
					{
						nlvec(3,k) = cLST_Age(exp(log_degrees_of_freedom(k)));
					}

					// Residual
					if(last_phase())
					{
						nu          = cLST_Age.get_standardized_residuals();
						age_tau2(k) = cLST_Age.get_sigma2();
					}
				break;
			}
			
			// | Extract residuals.
			for(i=1;i<=naa;i++)
			{
				A_nu(k)(i)(n_A_sage(k),n_A_nage(k))=nu(i);
			}
		}
	}
	
	// |---------------------------------------------------------------------------------|
	// | STOCK-RECRUITMENT LIKELIHOOD COMPONENT
	// |---------------------------------------------------------------------------------|
	// | - tau is the process error standard deviation.
	if( active(theta(1)) || active(theta(2)) )
	{
		for(g=1;g<=ngroup;g++)
		{
			nlvec(4,g) = dnorm(delta(g),tau);
		}
	}

	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD COMPONENTS FOR SELECTIVITY PARAMETERS
	// |---------------------------------------------------------------------------------|
	// | - lambda_1  -> penalty weight for smoothness
	// | - lambda_2  -> penalty weight for dome-shape
	// | - lambda_3  -> penalty weight for inter-annual variation.
	dvar_vector lvec(1,7); 
	lvec.initialize();
	int ig;
	for(k=1;k<=ngear;k++)
	{
		if(active(sel_par(k)))
		{
			//if not using logistic selectivity then
			//CHANGED from || to &&  May 18, 2011 Vivian
			if( isel_type(k)!=1 && 
				isel_type(k)!=7 && 
				isel_type(k)!=8 &&
				isel_type(k)!=11 )  
			{
				for(ig=1;ig<=n_ags;ig++)
				{
				for(i=syr;i<=nyr;i++)
				{
					//curvature in selectivity parameters
					dvar_vector df2 = first_difference(first_difference(log_sel(k)(ig)(i)));
					nlvec(5,k)     += lambda_1(k)/(nage-sage+1)*df2*df2;

					//penalty for dome-shapeness
					for(j=sage;j<=nage-1;j++)
						if(log_sel(k,ig,i,j)>log_sel(k,ig,i,j+1))
							nlvec(6,k)+=lambda_2(k)
										*square( log_sel(k,ig,i,j)-log_sel(k,ig,i,j+1) );
				}
				}
			}
			
			/*
			Oct 31, 2012 Halloween! Added 2nd difference penalty on time 
			for isel_type==(4)

			Mar 13, 2013, added 2nd difference penalty on isel_type==5 
			*/
			if( isel_type(k)==4 || isel_type(k)==5 || n_sel_blocks(k) > 1 )
			{
				for(ig=1;ig<=n_ags;ig++)
				{
					
				dvar_matrix trans_log_sel = trans( log_sel(k)(ig) );
				for(j=sage;j<=nage;j++)
				{
					dvar_vector df2 = first_difference(first_difference(trans_log_sel(j)));
					nlvec(7,k)     +=  lambda_3(k)/(nage-sage+1)*norm2(df2);
				}
				}
			}
			
		}
	}
	
	
	// |---------------------------------------------------------------------------------|
	// | CONSTRAINTS FOR SELECTIVITY DEVIATION VECTORS
	// |---------------------------------------------------------------------------------|
	// | [?] - TODO for isel_type==2 ensure mean 0 as well.
	// |
	for(k=1;k<=ngear;k++)
	{
		if( active(sel_par(k)) &&
			isel_type(k)!=1    &&
			isel_type(k)!=7    &&
			isel_type(k)!=8    &&
			isel_type(k)!=11 )
		{
			dvariable s = 0;
			if(isel_type(k)==5)  //bicubic spline version ensure column mean = 0
			{
				dvar_matrix tmp = trans(sel_par(k));
				for(j=1;j<=tmp.rowmax();j++)
				{
					s=mean(tmp(j));
					lvec(1)+=10000.0*s*s;
				}
			}
			if( isel_type(k)==2 ||
			    isel_type(k)==3 ||
			 	isel_type(k)==4 || 
				isel_type(k)==12 )
			{
				dvar_matrix tmp = sel_par(k);
				for(j=1;j<=tmp.rowmax();j++)
				{
					s=mean(tmp(j));
					lvec(1)+=10000.0*s*s;
				}
			}
		}
	}
	
	
	// |---------------------------------------------------------------------------------|
	// | PRIORS FOR LEADING PARAMETERS p(theta)
	// |---------------------------------------------------------------------------------|
	// | - theta_prior is a switch to determine which statistical distribution to use.
	// |
	dvariable ptmp; 
	dvar_vector priors(1,npar);
	priors.initialize();
	for(i=1;i<=npar;i++)
	{
		ptmp = 0;
		for(j=1;j<=ipar_vector(i);j++)
		{
			if( active(theta(i)) )
			{
				switch(theta_prior(i))
				{
				case 1:		//normal
					ptmp += dnorm(theta(i,j),theta_control(i,6),theta_control(i,7));
					break;
					
				case 2:		//lognormal CHANGED RF found an error in dlnorm prior. rev 116
					ptmp += dlnorm(theta(i,j),theta_control(i,6),theta_control(i,7));
					break;
					
				case 3:		//beta distribution (0-1 scale)
					double lb,ub;
					lb=theta_lb(i);
					ub=theta_ub(i);
					ptmp += dbeta((theta(i,j)-lb)/(ub-lb),theta_control(i,6),theta_control(i,7));
					break;
					
				case 4:		//gamma distribution
					ptmp += dgamma(theta(i,j),theta_control(i,6),theta_control(i,7));
					break;
					
				default:	//uniform density
					ptmp += log(1./(theta_control(i,3)-theta_control(i,2)));
					break;
				}
			}
		}
		priors(i) = ptmp;	
	}
	
	// |---------------------------------------------------------------------------------|
	// | PRIOR FOR SURVEY Q
	// |---------------------------------------------------------------------------------|
	// |
	dvar_vector qvec(1,nits);
	qvec.initialize();
	for(k=1;k<=nits;k++)
	{
		if(q_prior(k) == 1 )
		{
			qvec(k) = dnorm( log(q(k)), mu_log_q(k), sd_log_q(k) );
		}
	}
	
	
// 	//** Legacy **  By accident took Rick Methot's bag from Nantes.
// 	//301 787 0241  Richard Methot cell phone.
// 	//ibis charles du gaulle at
// 	//01 49 19 19 20
	

	// |---------------------------------------------------------------------------------|
	// | LIKELIHOOD PENALTIES TO REGULARIZE SOLUTION
	// |---------------------------------------------------------------------------------|
	// | - pvec(1)  -> penalty on mean fishing mortality rate.
	// | - pvec(2)  -> penalty on first difference in natural mortality rate deviations.
	// | - pvec(4)  -> penalty on recruitment deviations.
	// | - pvec(5)  -> penalty on initial recruitment vector.
	// | - pvec(6)  -> constraint to ensure sum(log_rec_dev) = 0
	// | - pvec(7)  -> constraint to ensure sum(init_log_rec_dev) = 0
	// |
	dvar_vector pvec(1,7);
	pvec.initialize();
	
	dvariable log_fbar = mean(log_ft_pars);
	if(last_phase())
	{
		
		pvec(1) = dnorm(log_fbar,log(d_iscamCntrl(7)),d_iscamCntrl(9));
		
		// | Penalty for log_rec_devs (large variance here)
		for(g=1;g<=n_ag;g++)
		{
			pvec(4) += dnorm(log_rec_devs(g),2.0);
			pvec(5) += dnorm(init_log_rec_devs(g),2.0);
			dvariable s = 0;
			s = mean(log_rec_devs(g));
			pvec(6) += 1.e5 * s*s;
			s = mean(init_log_rec_devs(g));
			pvec(7) += 1.e5 * s*s;
		}
	}
	else
	{

		pvec(1) = dnorm(log_fbar,log(d_iscamCntrl(7)),d_iscamCntrl(8));
		
		//Penalty for log_rec_devs (CV ~ 0.0707) in early phases
		for(g=1;g<=n_ag;g++)
		{
			pvec(4) += 100.*norm2(log_rec_devs(g));
			pvec(5) += 100.*norm2(init_log_rec_devs(g));
			dvariable s = 0;
			s = mean(log_rec_devs(g));
			pvec(6) += 1.e5 * s*s;
			s = mean(init_log_rec_devs(g));
			pvec(7) += 1.e5 * s*s;
		}
	}
	
	if(active(log_m_nodes))
	{
		double std_mdev = d_iscamCntrl(11);
		dvar_vector fd_mdevs=first_difference(log_m_devs);
		pvec(2)  = dnorm(fd_mdevs,std_mdev);
		pvec(2) += 0.5*norm2(log_m_nodes);
	}
	
	
	if(verbose)
	{
		COUT(nlvec);
		COUT(lvec);
		COUT(priors);
		COUT(pvec);
		COUT(qvec);
	}
	// COUT(nlvec);
	objfun  = sum(nlvec);
	objfun += sum(lvec);
	objfun += sum(priors);
	objfun += sum(pvec);
	objfun += sum(qvec);
	nf++;
	if(verbose)cout<<"**** Ok after calcObjectiveFunction ****"<<endl;
	
  }






// FUNCTION void equilibrium(const double& fe, const dvector& ak, const double& ro, const double& kap, const double& m, const dvector& age, const dvector& wa, const dvector& fa, const dmatrix& va,double& re,double& ye,double& be,double& ve,double& dye_df,double& d2ye_df2)//,double& phiq,double& dphiq_df, double& dre_df)
//   {
// 	/*
// 	Equilibrium age-structured model used to determin Fmsy and MSY based reference points.
// 	Author: Steven Martell
	
// 	Comments: 
// 	This code uses a numerical approach to determine a vector of fe_multipliers
// 	to ensure that the dAllocation is met for each gear type.
	
// 	args:
// 	fe	-steady state fishing mortality
// 	ak	-dAllocation of total ye to gear k.
// 	ro	-unfished sage recruits
// 	kap	-recruitment compensation ration
// 	m	-instantaneous natural mortality rate
// 	age	-vector of ages
// 	wa	-mean weight at age
// 	fa	-mean fecundity at age
// 	va	-mean vulnerablity at age for fe gear.

	
// 	Modified args:
// 	re	-steady state recruitment
// 	ye	-steady state yield
// 	be	-steady state spawning biomass
// 	phiq		-per recruit yield
// 	dre_df		-partial of recruitment wrt fe
// 	dphiq_df	-partial of per recruit yield wrt fe
	
// 	LUCIE'S RULE: the derivative of a sum is the sum of its derivatives.
// 	Lucie says: there is some nasty calculus in here daddy, are you sure
// 	you've got it right?
	
// 	I've got it pretty close. 
	
// 	DEPRECATE THIS FUNCTION.  NOW DONE IN THE MSY CLASS
	
// 	*/
// 	int i,j,k;
// 	int nage    = max(age);
// 	int sage    = min(age);
// 	double  dre_df;
// 	double  phif;
// 	dvector lx(sage,nage);
// 	dvector lz(sage,nage);
// 	dvector lambda(1,ngear);        //F-multiplier
// 	dvector phix(1,ngear);
// 	dvector phiq(1,ngear);
// 	dvector dphiq_df(1,ngear);
// 	dvector dyek_df(1,ngear);
// 	dvector d2yek_df2(1,ngear);
// 	dvector yek(1,ngear);
// 	dmatrix qa(1,ngear,sage,nage);
// 	dmatrix xa(1,ngear,sage,nage);  //vulnerable numbers per recruit
	
// 	lx          = pow(exp(-m),age-double(sage));
// 	lx(nage)   /=(1.-exp(-m));
// 	double phie = lx*fa;		// eggs per recruit
// 	double so   = kap/phie;
// 	double beta = (kap-1.)/(ro*phie);
// 	lambda      = ak/mean(ak);	// multiplier for fe for each gear
	
	
// 	/* Must iteratively solve for f-multilier */
// 	for(int iter=1;iter<=30;iter++)
// 	{
// 		/* Survivorship under fished conditions */
// 		lz(sage)    = 1.0;
// 		lambda     /= mean(lambda);
// 		dvector fk  = fe*lambda;
// 		dvector ra  = lambda*va;
// 		dvector za  = m + fe*ra;
// 		dvector sa  = mfexp(-za);
// 		dvector oa  = 1.0 - sa;
		
		
// 		for(k=1;k<=ngear;k++)
// 		{
// 			qa(k) = elem_prod(elem_div(lambda(k)*va(k),za),oa);
// 			xa(k) = elem_prod(elem_div(va(k),za),oa);
// 		}
		
// 		double dlz_df = 0, dphif_df = 0;
// 		dphiq_df.initialize();
// 		dre_df   = 0;
// 		for(j=sage;j<=nage;j++)
// 		{
// 			if(j>sage) lz(j)  = lz(j-1) * sa(j-1);
// 			if(j>sage) dlz_df = dlz_df  * sa(j-1) - lz(j-1)*ra(j-1)*sa(j-1);
			
// 			if(j==nage)
// 			{
// 				lz(j)  = lz(j) / oa(j);
				
// 				double t4 = (-ra(j-1)+ra(j))*sa(j)+ra(j-1);
// 				dlz_df = dlz_df/oa(j) - (lz(j-1)*sa(j-1)*t4) / square(oa(j));
// 			}
			
// 			dphif_df   = dphif_df+fa(j)*dlz_df;
// 			for(k=1;k<=ngear;k++)
// 			{
// 				double t1   = lambda(k) * wa(j) *va(k,j) * ra(j) * lz(j);
// 				double t3   = -1. + (1.+za(j)) * sa(j);
// 				double t9   = square(za(j));
// 				dphiq_df(k)+= wa(j)*qa(k,j)*dlz_df + t1 * t3 / t9; 
// 			}
// 		} 
		
// 		phif   = elem_prod(lz,exp(-za*d_iscamCntrl(13)))*fa;
// 		re     = ro*(kap-phie/phif)/(kap-1.);
// 		dre_df = ro/(kap-1.0)*phie/square(phif)*dphif_df;
		
// 		/* Equilibrium yield */
// 		for(k=1;k<=ngear;k++)
// 		{
// 			phix(k)      = sum(elem_prod(elem_prod(lz,wa),xa(k)));
// 			phiq(k)      = sum(elem_prod(elem_prod(lz,wa),qa(k)));
// 			yek(k)       = fe*re*phiq(k);
// 			dyek_df(k)   = re*phiq(k) + fe*phiq(k)*dre_df + fe*re*dphiq_df(k);
// 			d2yek_df2(k) = phiq(k)*dre_df + re*dphiq_df(k);
// 		}
		
// 		/* Iterative soln for lambda */
// 		dvector pk = yek/sum(yek);
// 		dvector t1 = elem_div(ak,pk+1.e-30);
// 		lambda     = elem_prod(lambda,t1);
// 		if(abs(sum(ak-pk))<1.e-6) break;
// 	} // end of iter
// 	ve       = re*sum(elem_prod(ak,phix));
// 	be       = re*phif;
// 	ye       = sum(yek);
// 	dye_df   = sum(dyek_df);
// 	d2ye_df2 = sum(d2yek_df2);

// 	// cout<<"EQUILIBRIUM CODE "<<setprecision(4)<<setw(2)<<fe<<setw(3)<<" "
// 	// <<ye<<setw(5)<<" "<<dye_df<<"  "<<dyek_df(1,3)<<endl;
//   }



	
// FUNCTION void equilibrium(const double& fe,const double& ro, const double& kap, const double& m, const dvector& age, const dvector& wa, const dvector& fa, const dvector& va,double& re,double& ye,double& be,double& phiq,double& dphiq_df, double& dre_df)
//   {
// 	/*
// 	This is the equilibrium age-structured model that is 
// 	conditioned on fe (the steady state fishing mortality rate).
	
// 	In the case of multiple fisheries, fe is to be considered as the
// 	total fishing mortality rate and each fleet is given a specified
// 	dAllocation based on its selectivity curve.  The dAllocation to 
// 	each fleet must be specified a priori.
	
// 	args:
// 	fe	-steady state fishing mortality
// 	ro	-unfished sage recruits
// 	kap	-recruitment compensation ration
// 	m	-instantaneous natural mortality rate
// 	age	-vector of ages
// 	wa	-mean weight at age
// 	fa	-mean fecundity at age
// 	va	-mean vulnerablity at age for fe gear.
// 	ak	-dAllocation of total ye to gear k.
	
// 	Modified args:
// 	re	-steady state recruitment
// 	ye	-steady state yield
// 	be	-steady state spawning biomass
// 	phiq		-per recruit yield
// 	dre_df		-partial of recruitment wrt fe
// 	dphiq_df	-partial of per recruit yield wrt fe
	
// 	FIXME add Ricker model to reference points calculations.
// 	FIXME partial derivatives for dphif_df need to be fixed when d_iscamCntrl(13)>0.
// 	*/
// 	int i;
	
// 	int nage=max(age);
// 	int sage=min(age);
// 	dvector lx=pow(exp(-m),age-double(sage));
// 	lx(nage)/=(1.-exp(-m));
// 	dvector lz=lx;
// 	dvector za=m+fe*va;
// 	dvector sa=1.-exp(-za);
// 	dvector qa=elem_prod(elem_div(va,za),sa);
	
// 	double phie = lx*fa;		//eggs per recruit
// 	double so = kap/phie;
// 	double beta = (kap-1.)/(ro*phie);
	
	
// 	double dlz_df = 0, dphif_df = 0;
// 	dphiq_df=0; dre_df=0;
// 	for(i=sage; i<=nage; i++)
// 	{
// 		if(i>sage) lz[i]=lz[i-1]*exp(-za[i-1]);
// 		if(i>sage) dlz_df=dlz_df*exp(-za[i-1]) - lz[i-1]*va[i-1]*exp(-za[i-1]);
// 		if(i==nage){ //6/11/2007 added plus group.
// 					lz[i]/=(1.-mfexp(-za[i]));
					
// 					dlz_df=dlz_df/(1.-mfexp(-za[i]))
// 							-lz[i-1]*mfexp(-za[i-1])*va[i]*mfexp(-za[i])
// 					/((1.-mfexp(-za[i]))*(1.-mfexp(-za[i])));
// 				}
// 		dphif_df=dphif_df+fa(i)*dlz_df;
// 		dphiq_df=dphiq_df+wa(i)*qa(i)*dlz_df+(lz(i)*wa(i)*va(i)*va(i))/za(i)*(exp(-za[i])-sa(i)/za(i));
// 	}
// 	//CHANGED need to account for fraction of mortality that occurs
// 	//before the spawning season in the recruitment calculation.
// 	//cout<<"lz\t"<<elem_prod(lz,exp(-za*d_iscamCntrl(13)))<<endl;
// 	//exit(1);
// 	//double phif = lz*fa;
// 	double phif = elem_prod(lz,exp(-za*d_iscamCntrl(13)))*fa;
// 	phiq=sum(elem_prod(elem_prod(lz,wa),qa));
// 	re=ro*(kap-phie/phif)/(kap-1.);
// 	//cout<<fe<<" spr ="<<phif/phie<<endl;
// 	if(re<=0) re=0;
// 	dre_df=(ro/(kap-1.))*phie/square(phif)*dphif_df;
// 	ye=fe*re*phiq;
// 	be=re*phif;	//spawning biomass
	
// 	//cout<<"Equilibrium\n"<<ro<<"\n"<<re<<"\n"<<ye<<endl;
	
//   }
	
FUNCTION void calcReferencePoints()
  {
  	/*
  	Purpose:  This function calculates the MSY-based reference points, and also loops
  	          over values of F and F-multipliers and calculates the equilibrium yield
  	          for each fishing gear.
  	Author: Steven Martell
  	
  	Arguments:
  		None
  	
  	NOTES:
  		- This function is based on the msyReferencePoint class object written by 
  		  Steve Martell on the Island of Maui while on Sabbatical leave from UBC.
  		- The msyReferencePoint class uses Newton-Raphson method to iteratively solve
  		  for the Fmsy values that maximize catch for each fleet. You can compare 
  		  MSY-reference points with the numerical values calculated at the end of this
  		  function.
  		- Code check: appears to find the correct value of MSY
		  in terms of maximizing ye.  Check to ensure rec-devs
		  need a bias correction term to get this right.
	
		- Modification for multiple fleets:
    	  	Need to pass a weighted average vector of selectivities
    	  	to the equilibrium routine, where the weights for each
    	  	selectivity is based on the dAllocation to each fleet.
		
	 		Perhaps as a default, assign an equal dAllocation to each
	 		fleet.  Eventually,user must specify dAllocation in 
	 		control file.  DONE
		
	 	- Use selectivity in the terminal year to calculate reference
	 	  points.  See todo, this is something that should be specified by the user.
	
		- June 8, 2012.  SJDM.  Made the following changes to this routine.
			1) changed reference points calculations to use the average
			   weight-at-age and fecundity-at-age.
			2) change equilibrium calculations to use the catch dAllocation
			   for multiple gear types. Not the average vulnerablity... this was wrong.
	
		- July 29, 2012.  SJDM Issue1.  New routine for calculating reference points
		  for multiple fleets. In this case, finds a vector of Fmsy's that simultaneously 
		  maximizes the total catch for each of the fleets respectively.  See
		  iSCAMequil_soln.R for an example.
	
		- August 1, 2012.  SJDM, In response to Issue1. A major overhaul of this routine.
		  Now using the new Msy class to calculate reference points. This greatly simplifies
		  the code in this routine and makes other routines (equilibrium) redundant.  Also
		  the new Msy class does a much better job in the case of multiple fleets.
	
		- Aug 11, 2012.
		  For the Pacific herring branch omit the get_fmsy calculations and use only the 
		  Bo calculuation for the reference points.  As there are no MSY based reference
		  points required for the descision table. 

		- May 8, 2013.
		  Starting to modify this code to allow for multiple areas, sex and stock-specific
		  reference points.  Key to this modification is the definition of a stock. For
		  the purposes of reference points, a stock is assumed to be distributed over all
		  areas, can be unisex or two sex, and can be fished by all fleets given the stock
		  exists in an area where the fishing gear operates.

		- Aug, 3-6, 2013.
		  Major effort went into revising this routine as well as the Msy class to 
		  calculate MSY-based reference points for multiple fleets and multiple sexes.  
		  I had found a significant bug in the dye calculation where I need to use the 
		  proper linear algebra to calculate the vector of dye values. This appears to 
		  be working properly and I've commented out the lines of code where I numerically
		  checked the derivatives of the catch equation.  This is a major acomplishment.

		- Mar, 2013.
		  A major new development here with the use of msy.hpp and a template class for 
		  calculating MSY-based reference points.  The user can now calculate reference
		  points for each gear based on fixed allocation, and optimum allocations based
		  on relative differences in selectivities among the gears landing fish. Uses the
		  name space "rfp".

   	PSEUDOCODE: 
   		(1) : Construct array of selectivities (potentially sex based log_sel)
   		(2) : Construct arrays of d3_wt_avg and d3_wt_mat for reference years.
	  	(3) : Come up with a reasonable guess for fmsy for each gear in nfleet.
	  	(4) : Instantiate an Msy class object and get_fmsy.
	  	(5) : Use Msy object to get reference points.

		
	slx::Selex<dvar_vector> * ptr;  //Pointer to Selex base class
  ptr = new slx::LogisticCurve<dvar_vector,dvariable>(mu,sd);
  log_sel = ptr->logSelectivity(age);
  delete ptr;

  	TODO list:
  	[ ] - allow user to specify which selectivity years are used in reference point
  	      calculations. This should probably be done in the projection File Control.
  	*/
	int kk,ig;
	
	

	// | (1) : Matrix of selectivities for directed fisheries.
	// |     : log_sel(gear)(n_ags)(year)(age)
	// |     : ensure dAllocation sums to 1.
	dvector d_ak(1,nfleet);
	d3_array  d_V(1,n_ags,1,nfleet,sage,nage);
	dvar3_array  dvar_V(1,n_ags,1,nfleet,sage,nage);
	for(k=1;k<=nfleet;k++)
	{
		kk      = nFleetIndex(k);
		d_ak(k) = dAllocation(kk);
		for(ig=1;ig<=n_ags;ig++)
		{
			d_V(ig)(k) = value( exp(log_sel(kk)(ig)(nyr)) );
			dvar_V(ig)(k) =( exp(log_sel(kk)(ig)(nyr)) );

		}
	}
	d_ak /= sum(d_ak);

	// | (2) : Average weight and mature spawning biomass for reference years
	// |     : dWt_bar(1,n_ags,sage,nage)
	dmatrix fa_bar(1,n_ags,sage,nage);
	dmatrix  M_bar(1,n_ags,sage,nage);
	for(ig=1;ig<=n_ags;ig++)
	{
		fa_bar(ig) = elem_prod(dWt_bar(ig),ma(ig));
		M_bar(ig)  = colsum(value(M(ig).sub(pf_cntrl(3),pf_cntrl(4))));
		M_bar(ig) /= pf_cntrl(4)-pf_cntrl(3)+1;	
	}
	
	// | (3) : Initial guess for fmsy for each fleet
	// |     : set fmsy = 2/3 of M divided by the number of fleets
	fmsy.initialize();
	fall.initialize();
	msy.initialize();
	bmsy = 0;

	dvar_vector dftry(1,nfleet);

	dvector ftry(1,nfleet);
	ftry  = 0.6/nfleet * mean(M_bar);
	fmsy  = ftry;
	dftry = ftry;
	


	// | (4) : Instantiate msy class for each stock
	for(g=1;g<=ngroup;g++)
	{
		//double d_ro = value(ro(g));
		//double  d_h = value(steepness(g));
		double d_rho = d_iscamCntrl(13);

		dvector d_mbar = M_bar(g);
		dvector   d_wa = dWt_bar(g);
		dvector   d_fa = fa_bar(g);

		//Pointer to the base class
		//rfp::referencePoints<dvariable,dvar_vector,dvar_matrix> * pMSY; 
		//pMSY = new rfp::msy<dvariable,dvar_vector,dvar_matrix,dvar3_array>
		//(ro(g),steepness(g),d_rho,M_bar,dWt_bar,fa_bar,dvar_V);
		//dvar_vector dfmsy = pMSY->getFmsy(dftry);
		//delete pMSY;
		

		// cout<<"Initial Fe "<<dftry<<endl;
		rfp::msy<dvariable,dvar_vector,dvar_matrix,dvar3_array> 
		c_MSY(ro(g),steepness(g),d_rho,M_bar,dWt_bar,fa_bar,dvar_V);

		dvar_vector dfmsy = c_MSY.getFmsy(dftry);
		bo  = c_MSY.getBo();
		dvariable dbmsy = c_MSY.getBmsy();
		dvar_vector dmsy = c_MSY.getMsy();
		bmsy = value(dbmsy);
		msy  = value(dmsy);

		c_MSY.print();

		
	}


	if(verbose)cout<<"**** Ok after calcReferencePoints ****"<<endl;
  }

  /**
   * This is a simple test routine for comparing the MSY class output to the 
   * MSF.xlsx spreadsheet that was used to develop the multiple fleet msy 
   * method.  Its a permanent feature of the iscam code for testing.
   */	
FUNCTION void testMSYxls()

	double ro = 1.0;
	double steepness = 0.75;
	double d_rho = 0.0;
	dmatrix m_bar(1,1,1,20);
	m_bar = 0.30;
	dmatrix dWt_bar(1,1,1,20);
	dWt_bar(1).fill("{0.005956243,0.035832542,0.091848839,0.166984708,0.252580458,0.341247502,0.427643719,0.508367311,0.581557922,0.646462315,0.703062177,0.751788795,0.793319224,0.828438536,0.857951642,0.882630331,0.903184345,0.920248191,0.934377843,0.946053327}");
	dmatrix fa_bar(1,1,1,20);
	fa_bar(1).fill("{0,0,0,0,0.252580458,0.341247502,0.427643719,0.508367311,0.581557922,0.646462315,0.703062177,0.751788795,0.793319224,0.828438536,0.857951642,0.882630331,0.903184345,0.920248191,0.934377843,0.946053327}");
	d3_array dvar_V(1,1,1,2,1,20);
	dvar_V(1)(1).fill("{0.001271016,0.034445196,0.5,0.965554804,0.998728984,0.999954602,0.99999838,0.999999942,0.999999998,1,1,1,1,1,1,1,1,1,1,1}");
	dvar_V(1)(2).fill("{0.000552779,0.006692851,0.07585818,0.5,0.92414182,0.993307149,0.999447221,0.999954602,0.999996273,0.999999694,0.999999975,0.999999998,1,1,1,1,1,1,1,1}");
	//dvar_V(1)(2).fill("{0.000189406,0.000789866,0.003287661,0.013576917,0.054313266,0.19332137,0.5,0.80667863,0.945686734,0.986423083,0.996712339,0.999210134,0.999810594,0.999954602,0.99998912,0.999997393,0.999999375,0.99999985,0.999999964,0.999999991}");

	dvector dftry(1,2);
	dftry = 0.1 ;
	cout<<"Initial Fe "<<dftry<<endl;
	rfp::msy<double,dvector,dmatrix,d3_array> 
	c_MSY(ro,steepness,d_rho,m_bar,dWt_bar,fa_bar,dvar_V);
	dvector dfmsy = c_MSY.getFmsy(dftry);
	cout<<"Fmsy = "<<dfmsy<<endl;


	dvector ak(1,2);
	ak = 0.3;
	ak(2) = 1-ak(1);
	rfp::msy<double,dvector,dmatrix,d3_array>
	c_MSYk(ro,steepness,d_rho,m_bar,dWt_bar,fa_bar,dvar_V);
	dvector dkmsy = c_MSYk.getFmsy(dftry,ak);
	cout<<"Fmsy_k ="<<dkmsy<<endl;

	c_MSYk.print();

	dvector akmsy = c_MSYk.getFmsy(dftry);
	c_MSYk.print();

	exit(1);





  	/**
  	Purpose:  This routine gets called from the PRELIMINARY_CALCS_SECTION if the 
  	          user has specified the -sim command line option.  The random seed
  	          is specifed at the command line.

  	Author: Steven Martell
  	
  	Arguments:
  		seed -> a random seed for generating a unique, repeatable, sequence of 
  		        random numbers to be used as observation and process errors.
  	
  	NOTES:
		- This routine will over-write the observations in memory
		  with simulated data, where the true parameter values are
		  the initial values.  Change the standard deviations of the 
		  random number vectors epsilon (observation error) or 
 		  recruitment devs wt (process error).
 		- At the end of DATA_SECTION nyrs is modified by retro_yrs if -retro.
 		- Add back the retro_yrs to ensure same random number sequence for 
		  observation errors.
 
 	PSUEDOCODE:
 		1)  calculate selectivities to be used in the simulations.
 		2)  calculate mortality rates (M), F is conditioned on catch.
 		3)  generate random numbers for observation & process errors.
 		4)  calculate survivorship and stock-recruitment pars based on average M & fec.
 		5)  initialize state variables.
 		6)  population dynamics with F conditioned on catch.
 		7)  compute catch-at-age samples.
 		8)  compute total catch.
 		9)  compute the relative abundance indices.
 		10) write simulated data to file.

  	TODO list:
	[?] - March 9, 2013.  Fix simulation model to generate consistent data when 
		  doing retrospective analyses on simulated datasets.
	[ ] - TODO: Stock-Recruitment model.
	[ ] - TODO: switch statement for catch-type to get Fishing mortality rate.
  	[ ] 
  	*/
FUNCTION void simulationModel(const long& seed)
  {
	cout<<global_parfile<<endl;
	bool pinfile = 0;
	cout<<"___________________________________________________\n"<<endl;
	cout<<"  **Implementing Simulation--Estimation trial**    "<<endl;
	cout<<"___________________________________________________"<<endl;
	//if(norm(log_rec_devs)!=0)
	if( global_parfile )
	{
		cout<<"\tUsing pin file for simulation"<<endl;
		pinfile = 1;
	}
	cout<<"\tRandom Seed No.:\t"<< rseed<<endl;
	cout<<"\tNumber of retrospective years: "<<retro_yrs<<endl;
	cout<<"___________________________________________________\n"<<endl;
	
	

	int ii,ig,ih;
	// |---------------------------------------------------------------------------------|
	// | 1) SELECTIVITY
	// |---------------------------------------------------------------------------------|
	// |
    calcSelectivities(isel_type);
    

    // |---------------------------------------------------------------------------------|
    // | 2) MORTALITY
    // |---------------------------------------------------------------------------------|
    // | - NOTE only natural mortality is computed at this time.
    // | [ ] - add simulated random-walk in natural mortality rate here.
    // |
    calcTotalMortality();
    F.initialize();
    Z.initialize();
    S.initialize();



    // |---------------------------------------------------------------------------------|
    // | 3) GENERATE RANDOM NUMBERS
    // |---------------------------------------------------------------------------------|
    // | - epsilon -> Observation errors
    // | - rec_dev -> Process errors
    // | - init_rec_dev
    // | [ ] - add other required random numbers if necessary.
    // |
	random_number_generator rng(seed);
	dmatrix      epsilon(1,nItNobs,1,n_it_nobs);
	dmatrix      rec_dev(1,n_ag,syr,nyr+retro_yrs);
	dmatrix init_rec_dev(1,n_ag,sage+1,nage);
	dvector      eta(1,nCtNobs);
    

	epsilon.fill_randn(rng);
	rec_dev.fill_randn(rng);
	init_rec_dev.fill_randn(rng);
	eta.fill_randn(rng);

    // | Scale survey observation errors
    double std;
    for(k=1;k<=nItNobs;k++)
    {
    	for(i=1;i<=n_it_nobs(k);i++)
    	{
    		std = 1.0e3;
    		if( it_wt(k,i)>0 )
    		{
    			std = value(sig/it_wt(k,i));
    		}
    		epsilon(k,i) = epsilon(k,i)*std - 0.5*std*std;
    	}
    }

    // | Scale process errors
    for(ih=1;ih<=n_ag;ih++)
    {
		std              = value(tau);
		rec_dev(ih)      = rec_dev(ih) * std - 0.5*std*std;
		init_rec_dev(ih) = init_rec_dev(ih)*std - 0.5*std*std;
    }

    // | Scale total catch errors
    std = d_iscamCntrl(4);
    for(ii=1;ii<=nCtNobs;ii++)
    {
    	eta(ii) = eta(ii)* std  - 0.5*std*std;
    }

    // |---------------------------------------------------------------------------------|
    // | 4) SURVIVORSHIP & STOCK-RECRUITMENT PARAMETERS BASED ON AVERAGE M & FECUNDITY
    // |---------------------------------------------------------------------------------|
    // | -> Loop over each group/stock and compute survivorship, phib, so and beta.
    // | - fa is the average mature weight-at-age
    // |
    double phib;
    dvector ma(sage,nage);
    dvector fa(sage,nage);
    dvector lx(sage,nage);
    dvector lw(sage,nage);

    for(g=1;g<=ngroup;g++)
    {
    	lx.initialize();
		lw.initialize();
		lx(sage) = 1.0;
		lw(sage) = 1.0;
		phib = 0;
		for(f=1;f<=narea;f++)
		{
			for(h=1;h<=nsex;h++)
			{
				ig = pntr_ags(f,g,h);
				for(j=sage;j<=nage;j++)
				{
					ma(j) = mean( trans(value(M(ig)))(j)  );
					fa(j) = mean( trans(d3_wt_mat(ig))(j) );
					if(j>sage)
					{
						lx(j) = lx(j-1) * exp(-ma(j-1));
					}
					lw(j) = lx(j) * exp( -ma(j)*d_iscamCntrl(13) );
				}
				lx(nage) /= 1.0 - exp(-ma(nage));
				lw(nage) /= 1.0 - exp(-ma(nage));

				phib += 1./(narea*nsex) * lw * fa;
			}
		}
		so(g)  = kappa(g)/phib;
		sbo(g) = ro(g) * phib;
		switch(int(d_iscamCntrl(2)))
		{
			case 1:
				beta(g) = (kappa(g)-1.0)/sbo(g);
			break;
			case 2:
				beta(g) = log(kappa(g))/sbo(g);
			break;
		}
    }


   // |---------------------------------------------------------------------------------|
   // | 5) INITIALIZE STATE VARIABLES
   // |---------------------------------------------------------------------------------|
   // |
	N.initialize();
	for(ig=1;ig<=n_ags;ig++)
	{
		dvector tr(sage,nage);
		f  = n_area(ig);
		g  = n_group(ig);
		ih = pntr_ag(f,g);
		
		lx.initialize();
		lx(sage) = 1.0;
		for(j=sage;j<nage;j++)
		{
			lx(j+1) = lx(j) * exp(-value(M(ig)(syr)(j)));
		}
		lx(nage) /= 1.0 - exp(-value(M(ig)(syr)(nage)));
		if( d_iscamCntrl(5) )
		{
			tr = log( value(ro(g)) ) + log(lx);
		}
		else if( !d_iscamCntrl(5) )
		{
			tr(sage)        = value(log_avgrec(ih)+rec_dev(ih)(syr));;
			tr(sage+1,nage) = value(log_recinit(ih)+init_rec_dev(ih));
			tr(sage+1,nage) = tr(sage+1,nage)+log(lx(sage+1,nage));
		}
		N(ig)(syr)(sage,nage) = 1./nsex * exp(tr);
		log_rt(ih)(syr-nage+sage,syr) = tr.shift(syr-nage+sage);

		for(i=syr+1;i<=nyr;i++)
		{
			log_rt(ih)(i) = (log_avgrec(ih)+rec_dev(ih)(i));
			N(ig)(i,sage) = 1./nsex * exp( log_rt(ih)(i) );
		}
		N(ig)(nyr+1,sage) = 1./nsex * exp(log_avgrec(ih));
	}

	// |---------------------------------------------------------------------------------|
	// | 6) POPULATION DYNAMICS WITH F CONDITIONED ON OBSERVED CATCH
	// |---------------------------------------------------------------------------------|
	// | - va  -> matrix of fisheries selectivity coefficients.
	// | - [ ] TODO: switch statement for catch-type to get Fishing mortality rate.
	// | - [ ] TODO: Stock-Recruitment model (must loop over area sex for each group).
	// | - bug! ft is a global variable used in calcCatchAtAge and calcTotalCatch. FIXED
	// | - [ ] TODO: fishing mortality rate with sex ratio of catch is unknown and assume
	// |             the same F value for both male and females.
	dmatrix va(1,ngear,sage,nage);
	dmatrix zt(syr,nyr,sage,nage);
	dmatrix st(syr,nyr,sage,nage);
	BaranovCatchEquation cBaranov;
	for(ig=1;ig<=n_ags;ig++)
	{
		dmatrix tmp_ft(syr,nyr,1,ngear);

		for(i=syr;i<=nyr;i++)
		{
			dvector ba = elem_prod(value(N(ig)(i)),d3_wt_avg(ig)(i));
			dvector ct = d3_Ct(ig)(i);

			// | Selectivity modifications if necessary
			for(k=1;k<=ngear;k++)
			{
				va(k) = exp(value(log_sel(k)(ig)(i)));
				if( d_iscamCntrl(15) == 1 && dAllocation(k) > 0 )
				{
					va(k)             = ifdSelex(va(k),ba,0.25);
					log_sel(k)(ig)(i) = log(va(k));
					// dlog_sel(k)(i) = log(va(k));
				}
			}

			// | [ ] TODO switch statement for catch_type to determine F.
			tmp_ft(i) = cBaranov.getFishingMortality(ct, value(M(ig)(i)), va, value(N(ig)(i)),d3_wt_avg(ig)(i));
			
			zt(i) = value(M(ig)(i));
			for(k=1;k<=ngear;k++)
			{
				ft(ig)(k)(i) = tmp_ft(i,k);
				F(ig)(i) += tmp_ft(i,k) * va(k);
				zt(i) += tmp_ft(i,k) * va(k);
			}
			st(i) = exp(-zt(i));
			Z(ig)(i) = M(ig)(i) + F(ig)(i);
			S(ig)(i) = exp(-Z(ig)(i));

			// | [ ] TODO: Stock-Recruitment model
			// | [ ] TODO: Add autocorrelation.

			if( !pinfile ) COUT("Add stock recruitment Model")
			/* 
			sbt(ig)(i) = (elem_prod(N(ig)(i),exp(-zt(i)*d_iscamCntrl(13)))*d3_wt_mat(ig)(i));
			if(i>=syr+sage-1 && !pinfile)
			{
				double rt,et;
				et = value(sbt(ig)(i-sage+1));
				if(d_iscamCntrl(2)==1)
				{
					rt = value(so(g)*et/(1.+beta(g)*et));
				}
				else
				{
					rt = value(so(g)*et*exp(-beta(g)*et));
				}
				N(ig)(i)(sage) = rt * exp(rec_dev(ih)(i)-0.5*tau*tau);
			}
			*/

			// | Update state variables
			N(ig)(i+1)(sage+1,nage) =++ elem_prod(N(ig)(i)(sage,nage-1),st(i)(sage,nage-1));
			N(ig)(i+1)(nage) += N(ig)(i)(nage)*st(i)(nage);
		}
	}


	// |---------------------------------------------------------------------------------|
	// | 7) CATCH-AT-AGE
	// |---------------------------------------------------------------------------------|
	// | - A is the matrix of observed catch-age data.
	// | - A_hat is the predicted matrix from which to draw samples.
	// |
	int kk,aa,AA;
	double age_tau = value(sig);
	
	calcAgeComposition();
	for(kk=1;kk<=nAgears;kk++)
	{
		aa = n_A_sage(kk);
		AA = n_A_nage(kk);
		dvector pa(aa,AA);
		for(ii=1;ii<=n_A_nobs(kk);ii++)
		{
			pa = value(A_hat(kk)(ii));
			d3_A(kk)(ii)(aa,AA)=rmvlogistic(pa,age_tau,i+seed);
		}
	}
	
	// |---------------------------------------------------------------------------------|
	// | 8) TOTAL CATCH
	// |---------------------------------------------------------------------------------|
	// | - dCatchData is the matrix of observations
	// | - need to over-write the d3_Ct with the new errors.
	
	calcTotalCatch();
	d3_Ct.initialize();
	for(ii=1;ii<=nCtNobs;ii++)
	{
		dCatchData(ii,7) = value(ct(ii)) * exp(eta(ii));
		i = dCatchData(ii)(1);
		k = dCatchData(ii)(2);
		f = dCatchData(ii)(3);
		g = dCatchData(ii)(4);
		h = dCatchData(ii)(5);
		if( h==0 )
		{
			for(h=1;h<=nsex;h++)
			{
				ig = pntr_ags(f,g,h);
				d3_Ct(ig)(i)(k) = 1./nsex*dCatchData(ii)(7);
			}
		}
		else
		{
			ig = pntr_ags(f,g,h);
			d3_Ct(ig)(i)(k) = dCatchData(ii)(7);
		} 
	}
	// cout<<d3_Ct(1)<<endl;


	// |---------------------------------------------------------------------------------|
	// | 9) RELATIVE ABUNDANCE INDICES
	// |---------------------------------------------------------------------------------|
	// | - d3_survey_data is the matrix of input data.
	// |

	calcSurveyObservations();
	for(kk=1;kk<=nItNobs;kk++)
	{
		for(ii=1;ii<=n_it_nobs(kk);ii++)
		{
			d3_survey_data(kk)(ii)(2) *= exp(epsilon(kk)(ii));	
		}
	}
	
	// |---------------------------------------------------------------------------------|
	// | 10) WRITE SIMULATED DATA TO FILE
	// |---------------------------------------------------------------------------------|
	// |
	writeSimulatedDataFile();

	
// 	calcReferencePoints();
// 	//cout<<"	OK after reference points\n"<<fmsy<<endl;
// 	//exit(1);
// 	//	REPORT(fmsy);
// 	//	REPORT(msy);
// 	//	REPORT(bmsy);
	
	
// 	cout<<"___________________________________________________"<<endl;
// 	ofstream ofs("iscam.sim");
// 	ofs<<"fmsy\n"<<fmsy<<endl;
// 	ofs<<"msy\n"<<msy<<endl;
// 	ofs<<"bmsy\n"<<bmsy<<endl;
// 	ofs<<"bo\n"<<bo<<endl;
// 	ofs<<"va\n"<<va<<endl;
// 	ofs<<"sbt\n"<<sbt<<endl;//<<rowsum(elem_prod(N,fec))<<endl;
// 	ofs<<"log_rec_devs\n"<<log_rec_devs<<endl;
// 	ofs<<"rt\n"<<rt<<endl;
// 	ofs<<"ct\n"<<obs_ct<<endl;
// 	ofs<<"ft\n"<<trans(ft)<<endl;
// 	//ofs<<"ut\n"<<elem_div(colsum(obs_ct),N.sub(syr,nyr)*wa)<<endl;
// 	ofs<<"iyr\n"<<iyr<<endl;
// 	ofs<<"it\n"<<it<<endl;
// 	ofs<<"N\n"<<N<<endl;
// 	ofs<<"A\n"<<A<<endl;
// 	ofs<<"dlog_sel\n"<<dlog_sel<<endl;
// 	cout<<"  -- Simuation results written to iscam.sim --\n";
// 	cout<<"___________________________________________________"<<endl;
	
// 	//cout<<N<<endl;
// 	//exit(1);
  }






















  	/**
  	Purpose:  This function writes a simulated data file based on the simulation
  		  model output when the user specifies the -sim option.  This is only
  	      necessary if the user wishes to perform a retrospecrtive analysis on
  	      simulated data. 

  	Author: Steven Martell
  	
  	Arguments:
  	\param	seed -> the random number seed that is concatenated into the file name.
  	
  	NOTES:
  		
  	
  	TODO list:
  	[ ] 
  	*/
FUNCTION writeSimulatedDataFile
  {
  	adstring sim_datafile_name = "Simulated_Data_"+str(rseed)+".dat";
  	ofstream dfs(sim_datafile_name);
  	dfs<<"#Model dimensions"<<endl;
  	dfs<< narea 		<<endl;
  	dfs<< ngroup		<<endl;
  	dfs<< nsex			<<endl;
  	dfs<< syr   		<<endl;
  	dfs<< nyr   		<<endl;
  	dfs<< sage  		<<endl;
  	dfs<< nage  		<<endl;
  	dfs<< ngear 		<<endl;
 
  	dfs<<"#Allocation"	<<endl;
  	dfs<< dAllocation 	<<endl;
  	

  	dfs<<"#Age-schedule and population parameters"<<endl;
  	dfs<< d_linf  		<<endl;
  	dfs<< d_vonbk  		<<endl;
  	dfs<< d_to  			<<endl;
  	dfs<< d_a  			<<endl;
  	dfs<< d_b  			<<endl;
  	dfs<< d_ah  			<<endl;
  	dfs<< d_gh  			<<endl;

  	dfs<<"#Observed catch data"<<endl;
  	dfs<< nCtNobs 		<<endl;
  	dfs<< dCatchData    <<endl;

  	dfs<<"#Abundance indices"	<<endl;
  	dfs<< nItNobs 					<<endl;
  	dfs<< n_it_nobs 				<<endl;
  	dfs<< n_survey_type 			<<endl;
  	dfs<< d3_survey_data 			<<endl;

  	dfs<<"#Age composition"		<<endl;
  	dfs<< nAgears				<<endl;
  	dfs<< n_A_nobs				<<endl;
  	dfs<< n_A_sage				<<endl;
  	dfs<< n_A_nage				<<endl;
  	dfs<< d3_A						<<endl;

  	dfs<<"#Empirical weight-at-age data"	<<endl;
  	dfs<< nWtNobs				<<endl;
	dfs<< inp_wt_avg			<<endl;

	dfs<<"#EOF"	<<endl;
	dfs<< 999	<<endl;
	
	// | END OF WRITING SIMULATED DATAFILE.
  }



  	/**
  	Purpose:  This function returns a modified selectivity vector (va) based on
  			  the assumption that age-based selectivity will operate on the principle
  	          of ideal free distribution.
  	Author: Steven Martell
  	
  	Arguments:
  	\param	va -> age-specific vulnerability
  	\param	ba -> age-specific biomass (relative abundance is fine)
  	
  	NOTES:
  		
  	
  	TODO list:
  	[ ] 
  	*/
FUNCTION dvector ifdSelex(const dvector& va, const dvector& ba, const double& mpow)
  {

  	dvector pa(sage,nage);

  	pa = (elem_prod(va,pow(ba,mpow)));
  	pa = pa/sum(pa);
  	pa = exp( log(pa) - log(mean(pa)) );
  	return (pa);
  }

REPORT_SECTION
  
  	cout<<"You got here"<<endl;
	if(verbose)cout<<"Start of Report Section..."<<endl;
	report<<DataFile<<endl;
	report<<ControlFile<<endl;
	report<<ProjectFileControl<<endl;
	REPORT(objfun);
	REPORT(nlvec);
	REPORT(ro);
	dvector rbar=value(exp(log_avgrec));
	REPORT(rbar);
	dvector rinit=value(exp(log_recinit));
	REPORT(rinit);
	REPORT(sbo);
	REPORT(kappa);
	dvector steepness=value(theta(2));
	REPORT(steepness);
	REPORT(m);
	// double tau = value(sqrt(1.-rho)*varphi);
	// double sig = value(sqrt(rho)*varphi);
	REPORT(tau);
	REPORT(sig);
	REPORT(age_tau2);
	
	// |---------------------------------------------------------------------------------|
	// | MODEL DIMENSIONS & AGE-SCHEDULE INFORMATION ON GROWTH AND MATURITY
	// |---------------------------------------------------------------------------------|
	// |
	REPORT(narea);
	REPORT(ngroup);
	REPORT(nsex);
	REPORT(syr);
	REPORT(nyr);
	REPORT(sage);
	REPORT(nage);
	REPORT(ngear);

	ivector yr(syr,nyr);
	ivector yrs(syr,nyr+1);
	yr.fill_seqadd(syr,1); 
	yrs.fill_seqadd(syr,1); 
	REPORT(yr);
	REPORT(yrs);
	// REPORT(iyr);  //DEPRECATE, old survey years
	REPORT(age);
	REPORT(la);
	REPORT(wa);
	REPORT(ma);

	// |---------------------------------------------------------------------------------|
	// | OBSERVED AND PREDICTED DATA AND RESIDUALS
	// |---------------------------------------------------------------------------------|
	// | - Catch data
	// | - Survey data
	// | - Age composition data
	// | - Empirical weight-at-age data
	REPORT(dCatchData);
	REPORT(ct);
	REPORT(eta);

	REPORT(q);
	REPORT(qt);
	REPORT(d3_survey_data);
	REPORT(it_hat);
	REPORT(epsilon);

	if(n_A_nobs(nAgears) > 0)
	{
		REPORT(n_A_sage);
		REPORT(n_A_nage);
		REPORT(d3_A);
		REPORT(A_hat);
		REPORT(A_nu);


		/// The following is a total hack job to get the effective sample size
		/// for the multinomial distributions.

		// TODO Fix the retrospective bug here near line 4507 (if iyr<=nyr)
		report<<"Neff"<<endl;
		dvector nscaler(1,nAgears);
		nscaler.initialize();
		for(k = 1; k<=nAgears; k++)
		{
			if( int(nCompLikelihood(k)) )
			{
				int naa=0;
				int iyr;
				//retrospective counter
				for(i=1;i<=n_A_nobs(k);i++)
				{
					iyr = d3_A(k)(i)(n_A_sage(k)-5);	//index for year
					if(iyr<=nyr) naa++; else continue;
				}
				
				dmatrix     O = trans(trans(d3_A_obs(k)).sub(n_A_sage(k),n_A_nage(k))).sub(1,naa);
				dvar_matrix P = trans(trans(A_hat(k)).sub(n_A_sage(k),n_A_nage(k))).sub(1,naa);
				
				for(j = 1; j<= naa; j++)
				{
					double effectiveN = neff(O(j)/sum(O(j)),P(j));
					report<<sum(O(j))<<"\t"<<effectiveN<<endl;
					nscaler(k) += effectiveN;
				}	
				
				nscaler(k) /= naa;
			}
		}
		REPORT(nscaler);
	}

	// d3_wt_avg(1,n_ags,syr,nyr+1,sage,nage);
	adstring tt = "\t";
	REPORT(inp_wt_avg);
	REPORT(dWt_bar);
	// REPORT(d3_wt_avg);
	REPORT(d3_wt_mat);
	REPORT(d3_wt_dev);

	report<<"d3_wt_avg"<<endl;
	for(int ig=1;ig<=n_ags;ig++)
	{
		f = n_area(ig);
		g = n_group(ig);
		h = n_sex(ig);
		
		for(i=syr;i<=nyr;i++)
		{
			//year area stock sex |age columns (sage, nage) of weight at age data |
			report<<i<<tt;
			report<<f<<tt;
			report<<g<<tt;
			report<<h<<tt;
			report<<d3_wt_avg(ig)(i)<<endl;
		}
	
	}


	// |---------------------------------------------------------------------------------|
	// | SELECTIVITIES (4darray)
	// |---------------------------------------------------------------------------------|
	// |
	report<<"log_sel"<<endl;
	for(k=1;k<=ngear;k++)
	{
		for(int ig=1;ig<=n_ags;ig++)
		{
			for(i=syr;i<=nyr;i++)
			{
				report<<k<<"\t"<<ig<<"\t"<<i<<"\t"<<log_sel(k)(ig)(i)<<endl;	
			}
		}
	}
	// |---------------------------------------------------------------------------------|
	// | MORTALITY
	// |---------------------------------------------------------------------------------|
	// |
	// REPORT(ft);
	report<<"ft"<<endl;
	for(int ig = 1; ig <= n_ags; ig++ )
	{
		report<<ft(ig)<<endl;
	}
	REPORT(M);
	REPORT(F);
	REPORT(Z);

	// |---------------------------------------------------------------------------------|
	// | STOCK-RECRUITMENT
	// |---------------------------------------------------------------------------------|
	// |
	int rectype=int(d_iscamCntrl(2));
	REPORT(rectype);
	REPORT(so);
	REPORT(beta);
	REPORT(sbt);
	REPORT(bt);
	REPORT(rt);
	REPORT(delta);
	
	dmatrix rep_rt = value( exp(trans(trans(log_rt).sub(syr,nyr))) );
	for(int ig = 1; ig <= n_ag; ig++ )
	{
		rep_rt(ig)(syr) = value( exp( log_rt(ig)(syr-nage+sage) ) );
	}
	REPORT(rep_rt);

	// |---------------------------------------------------------------------------------|
	// | ABUNDANCE IN NUMBERS 
	// |---------------------------------------------------------------------------------|
	// |
	REPORT(N);

	// |---------------------------------------------------------------------------------|
	// | MSY-BASED REFERENCE POINTS
	// |---------------------------------------------------------------------------------|
	// |
	if( last_phase() )
	{
		calcReferencePoints();
		cout<<"Finished calcReferencePoints"<<endl;
		//exit(1);
		REPORT(bo);
		REPORT(fmsy);
		REPORT(msy);
		REPORT(bmsy);
		// REPORT(Umsy);
	}

	cout<<"You got to the end of the report section"<<endl;
	// |---------------------------------------------------------------------------------|
	// | OUTPUT FOR OPERATING MODEL
	// |---------------------------------------------------------------------------------|
	// | Move to final section?
	if( last_phase() )
	{
		ofstream ofs("iSCAM.res");
		ofs<<bo<<endl;
		ofs<<fmsy<<endl;
		ofs<<msy<<endl;
		ofs<<bmsy<<endl;
		for( g = 1; g <= ngroup; g++ )
		{
			ofs<<sbt(g)(nyr+1)<<"\t";
		}
		ofs<<endl;

		// projected biomass
		// The total biomass for each stock
		for( g = 1; g <= ngroup; g++ )
		{
			ofs<<bt(g)(nyr+1)<<"\t";
		}
		ofs<<endl;
	}


// 	/*
// 	Stock status info
// 	Bstatus = sbt/bmsy;
// 	Fstatus = ft/fmsy; If fmsy > 0 
// 	*/
// 	if(bmsy>0)
// 	{
// 		dvector Bstatus=value(sbt/bmsy);
// 		REPORT(Bstatus);
// 	}
	
// 	dmatrix Fstatus(1,ngear,syr,nyr);
// 	Fstatus.initialize();
// 	for(k = 1; k <= nfleet; k++)
// 	{
// 		if(fmsy(k) >0 )
// 		{
// 			j    = nFleetIndex(k);
// 			Fstatus(j) = value(ft(j)/fmsy(k));
// 		}
// 	}
// 	REPORT(Fstatus);
	
// 	//Parameter controls
// 	dmatrix ctrl=theta_control;
// 	REPORT(ctrl);
	
	
// 	if(last_phase()) decision_table();
	
	
// 	dvector rt3(1,3);
// 	if(last_phase())
// 	{
// 		dvector rt3 = age3_recruitment(value(column(N,3)),d3_wt_avg(nyr+1,3),value(M_tot(nyr,3)));
// 		REPORT(rt3);
// 	}
	
// 	//dvector future_bt = value(elem_prod(elem_prod(N(nyr+1),exp(-M_tot(nyr))),d3_wt_avg(nyr+1)));
// 	dvector future_bt = value(elem_prod(N(nyr+1)*exp(-m_bar),d3_wt_avg(nyr+1)));
// 	REPORT(future_bt);
// 	double future_bt4 = sum(future_bt(4,nage));
// 	REPORT(future_bt4);
	

	if(verbose)cout<<"END of Report Section..."<<endl;
	
	/*IN the following, I'm renaming the report file
	in the case where retrospective analysis is occurring*/
	//#if defined __APPLE__ || defined __linux
	//if( retro_yrs && last_phase() )
	//{
	//	//adstring rep="iscam.ret"+str(retro_yrs);
	//	//rename("iscam.rep",rep);
	//	adstring copyrep = "cp iscam.rep iscam.ret"+str(retro_yrs);
	//	system(copyrep);
	//}
	//// cout<<"Ya hoo"<<endl;
	////exit(1);
	//#endif
//
//	//#if defined _WIN32 || defined _WIN64
//	//if( retro_yrs && last_phase() )
//	//{
//	//	//adstring rep="iscam.ret"+str(retro_yrs);
//	//	//rename("iscam.rep",rep);
//	//	adstring copyrep = "copy iscam.rep iscam.ret"+str(retro_yrs);
//	//	system(copyrep);
//	//}
//	//// cout<<"Windows"<<endl;
//	//// exit(1);
	//#endif
	
   // REPORT_SECTION END
	
// FUNCTION decision_table
//   {
// 	/*
// 	This function takes a vector of projected catches and computes the following
// 	Reference points: Bmsy, Bo, Fmsy, Umsy.
	
// 	Biomass Metrics for the decision table:
// 	1) P(SB_{t+1} < SB_{t})
// 	2) P(SB_{t+1} < 0.25 B_{0})
// 	3) P(SB_{t+1} < 0.75 B_{0})
// 	4) P(SB_{t+1} < 0.40 B_{MSY})
// 	5) P(SB_{t+1} < 0.80 B_{MSY})
	
// 	Harvest Metrics for the decision table:
// 	1) P(U_{t+1} > Target harvest rate)
// 	2) P(U_{t+1} > 1/2 Fmsy)
// 	3) P(U_{t+1} > 2/3 Fmsy)
// 	4) P(tac/3+  > 20%)
	
// 	Key to the harvest metric is the definition of Umsy and dAllocation to fleets.
	
// 	Pseudocode:
// 		1) Calculate reference points (Fmsy, Bmsy)
// 		2) Loop over vector of proposed catches
// 		3) Evaluate biomass metrics for each posterior sample
// 		4) Evaluate harvest metrics for each posterior sample
	
// 	*/
// 	int i;
// 	// 1) Calculate reference pionts.
// 	//calcReferencePoints();  //redundant b/c its called in mcmc_output?
	
// 	// 2) Loop over vector of proposed catches
// 	//    This vector should is now read in from the projection file control (pfc).
// 	for(i=1;i<=n_tac;i++)
// 	{
// 		cout<<i<<" "<<n_tac<<endl;
// 		projection_model(tac(i));
// 	}
// 	// cout<<"Ok to here"<<endl;
//   }
	
FUNCTION generate_new_files
   {

   	ofstream rd("RUN.dat");
		rd<<NewFileName + ".dat"<<endl;
		rd<<NewFileName + ".ctl"<<endl;
		rd<<NewFileName + ".pfc"<<endl;



	#if defined __APPLE__ || defined __linux

	adstring bscmddat = "cp ../lib/iscam.dat" + NewFileName +".dat";
		system(bscmddat);

	adstring bscmdctl = "cp ../lib/ iscam.ctl" + NewFileName +".ctl";
		system(bscmdctl);

	adstring bscmdpfc = "cp ../lib/ iscam.PFC" + NewFileName +".pfc";
		system(bscmdpfc);	

	#endif

	#if defined _WIN32 || defined _WIN64

	adstring bscmddat = "copy ../lib/iscam.dat" + NewFileName +".dat";
		system(bscmddat);

	adstring bscmdctl = "copy ../lib/ iscam.ctl" + NewFileName +".ctl";
		system(bscmdctl);

	adstring bscmdpfc = "copy ../lib/ iscam.PFC" + NewFileName +".pfc";
		system(bscmdpfc);	

	#endif

  }

FUNCTION mcmc_output
//   {
// 	if(nf==1){
// 		ofstream ofs("iscam.mcmc");
// 		ofs<<"       log.ro";
// 		ofs<<"        log.h";
// 		ofs<<"        log.m";
// 		ofs<<"     log.rbar";
// 		ofs<<"    log.rinit";
// 		ofs<<"          rho";
// 		ofs<<"     vartheta";
// 		ofs<<"           bo";
// 		ofs<<"         bmsy";
// 		for(int k=1;k<=nfleet;k++) ofs<<"         msy"<<k;
// 		for(int k=1;k<=nfleet;k++) ofs<<"        fmsy"<<k;
// 		ofs<<"          SSB";
// 		ofs<<"        Age-4";
// 		ofs<<"         Poor";
// 		ofs<<"      Average";
// 		ofs<<"         Good";
// 		for(int i=1;i<=nItNobs;i++) ofs<<"         lnq"<<i;
// 		ofs<<"            f";
// 		ofs<<endl;
		
// 		ofstream of1("sbt.mcmc");
// 		ofstream of2("rt.mcmc");
		
// 	}
	
// 	// leading parameters & reference points
// 	calcReferencePoints();
// 	// decision table output
// 	//dvector future_bt = value(elem_prod(elem_prod(N(nyr+1),exp(-M_tot(nyr))),d3_wt_avg(nyr+1)));
// 	dvector future_bt = value(elem_prod(N(nyr+1)*exp(-m_bar),d3_wt_avg(nyr+1)));
// 	double future_bt4 = sum(future_bt(4,nage));
// 	dvector rt3 = age3_recruitment(value(column(N,3)),d3_wt_avg(nyr+1,3),value(M_tot(nyr,3)));	
	
	
// 	//if( bmsy > 0 && min(fmsy) >= 0 )
// 	{
// 		ofstream ofs("iscam.mcmc",ios::app);
// 		ofs<<setw(12)<<theta;
// 		ofs<<setw(13)<< bo;
// 		ofs<<setw(13)<< bmsy;
// 		ofs<<setw(12)<< msy;
// 		ofs<<setw(12)<< fmsy;
// 		ofs<<setw(13)<< sbt(nyr);
// 		ofs<<setw(13)<< future_bt4;
// 		ofs<<setw(12)<< future_bt4+rt3;
// 		ofs<<setw(12)<< log(q);
// 		ofs<<setw(13)<< objfun;
// 		ofs<<endl;
		
// 		/* June 12, 2012.  SJDM Call decision table. */
// 		decision_table();  
// 	}
// 	// output spawning stock biomass
// 	ofstream of1("sbt.mcmc",ios::app);
// 	of1<<setw(10)<<sbt(syr,nyr)<<endl;
	
// 	// output age-1 recruits
// 	ofstream of2("rt.mcmc",ios::app);
// 	of2<<setw(10)<<rt<<endl;
	
	
//   }

// FUNCTION dvector age3_recruitment(const dvector& rt, const double& wt,const double& M)
//   {
// 	/*
// 	This routine returns the poor average and good age-3 recruits
// 	that is used in constructing the decision table for the pacific
// 	herring fisheries.
	
// 	-1) sort the rt vector from small to large
// 	-2) find the 33rd and 66th percentiles
// 	-3) take the averge of the 0-33, 34-66, 67-100
// 	-4) multiply by the average weight
// 	-5) return the age-3 recruitment biomass
// 	*/
	
// 	dvector s_rt = sort(rt);
// 	dvector rbar(1,3);
	
// 	double idx = floor((nyr-syr+1.)/3.);
// 	int ix1 = syr+int(idx);
// 	int ix2 = syr+int(2.*idx);
// 	rbar(1) = mean(s_rt(syr,ix1));
// 	rbar(2) = mean(s_rt(ix1+1,ix2));
// 	rbar(3) = mean(s_rt(ix2+1,nyr));
// 	rbar = rbar*wt*exp(-M);
// 	//cout<<rbar<<endl;
// 	return(rbar);
//   }

// FUNCTION void projection_model(const double& tac);
//   {
// 	/*
// 	This routine conducts population projections based on 
// 	the estimated values of theta.  Note that all variables
// 	in this routine are data type variables.
	
// 	Arguments:
// 	tac is the total allowable catch that must be allocated 
// 	to each gear type based on dAllocation(k)
	
// 	theta(1) = log_ro
// 	theta(2) = h
// 	theta(3) = log_m
// 	theta(4) = log_avgrec
// 	theta(5) = log_recinit
// 	theta(6) = rho
// 	theta(7) = vartheta
	
// 	** NOTES **
// 	* Projections are based on average natural mortality and fecundity.
// 	* Selectivity is based on selectivity in terminal year.
// 	* Average weight-at-age is based on mean weight in the last 5 years.
	
// 	Aug 20, 2012 Found a bug, See issue 2 on github.
// 	*/
// 	static int runNo=0;
// 	runNo ++;
// 	int i,j,k;
// 	int pyr = nyr+1;	//projection year.
	
// 	// --derive stock recruitment parameters
// 	// --survivorship of spawning biomass
// 	dvector lx(sage,nage);
// 	double   tau = value(sqrt(1.-rho)*varphi); 
// 	double   m_M = value(m_bar); 
// 	double m_rho = d_iscamCntrl(13);
// 	lx(sage)     = 1;
// 	for(i=sage; i<=nage; i++)
// 	{
// 		lx(i) = exp( -m_M*(i-sage) -m_rho*m_M );
// 		if(i==nage) 
// 			lx(i) /= 1.0 - exp( -m_M );
// 	}	
// 	double phib = lx*avg_fec;
// 	double so   = value(kappa)/phib;
// 	double bo   = value(ro)*phib;
	
// 	double beta;
// 	switch(int(d_iscamCntrl(2)))
// 	{
// 		case 1:  // Beverton-Holt
// 			beta = (value(kappa)-1.)/bo;
// 		break;
// 		case 2:  // Ricker
// 			beta = log(value(kappa))/bo;
// 		break;
// 	}
	
// 	/* Fill arrays with historical values */
// 	dvector p_sbt(syr,pyr);
// 	dvector  p_ct(1,ngear);
// 	dmatrix  p_ft(nyr+1,pyr,1,ngear);
// 	dmatrix   p_N(syr,pyr+1,sage,nage);
// 	dmatrix   p_Z(syr,pyr,sage,nage);
// 	p_N.initialize();
// 	p_sbt.initialize();
// 	p_Z.initialize();
// 	p_N.sub(syr,nyr)   = value( N.sub(syr,nyr) );
// 	p_sbt(syr,nyr)     = value( sbt(syr,nyr)   );
// 	p_Z.sub(syr,nyr)   = value( Z.sub(syr,nyr) );
	
	
// 	/* Selectivity and dAllocation to gears */
// 	dmatrix va_bar(1,ngear,sage,nage);
// 	for(k=1;k<=ngear;k++)
// 	{
// 		p_ct(k)   = dAllocation(k)*tac;
// 		va_bar(k) = exp(value(log_sel(k)(nyr)));
// 	}
		
// 	/* Simulate population into the future under constant tac policy. */
// 	for(i = nyr; i<=pyr; i++)
// 	{
		
// 		if(i > nyr)
// 		{
// 			// get_ft is defined in the Baranov.cxx file
// 			p_ft(i) = getFishingMortality(p_ct, value(m_bar), va_bar, p_N(i),avg_wt);
// 			// p_ft(i)(1,nfleet) = fmsy(1,nfleet);
			
// 			// calculate total mortality in future years
// 			p_Z(i) = value(m_bar);
// 			for(k=1;k<=ngear;k++)
// 			{
// 				p_Z(i)+=p_ft(i,k)*va_bar(k);
// 			}
// 		}
		
		
// 		// spawning biomass
// 		p_sbt(i) = elem_prod(p_N(i),exp(-p_Z(i)*d_iscamCntrl(13))) * avg_fec;
		
		
// 		// sage recruits with random deviate xx
// 		// note the random number seed is repeated for each tac level.
// 		double  xx = randn(nf+i)*tau;
// 		if(i>=syr+sage-1)
// 		{
// 			double rt;
// 			double et = p_sbt(i-sage+1);		// lagged spawning biomass
			
// 			if(d_iscamCntrl(2)==1)						// Beverton-Holt model
// 			{
// 				rt=(so*et/(1.+beta*et));
// 			}
// 			if(d_iscamCntrl(2)==2)						// Ricker model
// 			{
// 				rt=(so*et*exp(-beta*et));
// 			}
			
// 			p_N(i+1,sage)=rt*exp(xx-0.5*tau*tau); 
// 		}
		
// 		/* Update numbers at age in future years */
// 		p_N(i+1)(sage+1,nage) =++ elem_prod(p_N(i)(sage,nage-1),exp(-p_Z(i)(sage,nage-1)));
// 		p_N(i+1,nage)        +=   p_N(i,nage)*exp(-p_Z(i,nage));
		
// 		//Predicted catch for checking calculations
// 		//for(k=1;k<=nfleet;k++)
// 		//{
// 		//	dvector ba = elem_prod(p_N(i),avg_wt);
// 		//	cout<<k<<" tac = "<<tac<<"\t ct = ";
// 		//	cout<<sum(elem_div(elem_prod(elem_prod(ba,p_ft(i,k)*va_bar(k)),1.-exp(-p_Z(i))),p_Z(i)));
// 		//	cout<<" fmsy = "<<fmsy<<" ft = "<<p_ft(i,k)<<endl;
// 		//}
// 	}	
// 	//cout<<"fmsy\n"<<fmsy<<endl;
// 	//cout<<"Spawning biomass\n"<<p_sbt<<endl;
// 	//exit(1);
// 	/* 
// 	  Write output to *.proj file for constructing decision tables. 
	
// 	  Biomass Metrics for the decision table:
// 	  1) P(SB_{t+1} < SB_{t})
// 	  2) P(SB_{t+1} < 0.25 B_{0})
// 	  3) P(SB_{t+1} < 0.75 B_{0})
// 	  4) P(SB_{t+1} < 0.40 B_{MSY})
// 	  5) P(SB_{t+1} < 0.80 B_{MSY})
	  
// 	  Harvest Metrics for the decision table:
// 	  1) P(U_{t+1} > Umsy)
// 	  2) P(U_{t+1} > 1/2 Umsy)
// 	  3) P(U_{t+1} > 2/3 Umsy)
// 	  4) P(tac/2+  > 20%)
	
// 	  Metric for spawning depletion and spawning trends:
// 	  1) P(5-year decline)
	  
// 	  Defn: Stock status is based on spawning biomass
// 	  Defn: Removal rate is based on removals/spawning biomass
// 	  Defn: Harvest rate    : U_{t+1}=TAC/SBio
// 	  Defn: MSY harvest rate: Umsy=MSY/SBmsy
	  
	
// 	*/
// 	if(nf==1 && runNo==1)
// 	{
// 		ofstream ofs(BaseFileName + ".proj");
// 		ofs<<" tac";
// 		ofs<<"      P(SB1)"; 
// 		ofs<<"      P(SB2)";
// 		ofs<<"      P(SB3)";
// 		ofs<<"      P(SB4)";
// 		ofs<<"      P(SB5)";
// 		ofs<<"       P(U1)";
// 		ofs<<"       P(U2)";
// 		ofs<<"       P(U3)";
// 		ofs<<"       P(U4)";
// 		ofs<<"       P(D5)";
// 		ofs<<endl;
// 		cout<<"Bo when nf==1 \t"<<bo<<endl;
// 	}
	
// 	double  ut  = tac / p_sbt(pyr);
// 	double u20  = tac / ( (p_N(pyr)(3,nage)*exp(-value(M_tot(nyr,3))))* avg_wt(3,nage) );
	
// 	/* Average rate of change in spawning biomass in last 5 years */
// 	double dSb5 = mean(log(p_sbt(pyr-5,pyr)) - log(p_sbt(pyr-6,pyr-1).shift(pyr-5)));
	
// 	ofstream ofs(BaseFileName + ".proj",ios::app);
// 	ofs<< setprecision(4)               <<setw(4) 
// 	   << tac                           <<setw(12)
// 	   << p_sbt(pyr-1)/p_sbt(pyr)       <<setw(12)
// 	   << 0.25*bo/p_sbt(pyr)            <<setw(12)
// 	   << 0.75*bo/p_sbt(pyr)            <<setw(12)
// 	   << 0.40*bmsy/p_sbt(pyr)          <<setw(12)
// 	   << 0.80*bmsy/p_sbt(pyr)          <<setw(12)
// 	   << ut/Umsy                       <<setw(12)
// 	   << ut/(0.5*Umsy)                 <<setw(12)
// 	   << ut/(2./3.*Umsy)               <<setw(12)
// 	   << u20/0.2                       <<setw(12)
// 	   << dSb5+1                        <<setw(12)
// 	   << endl;
// 	// cout<<"Finished projection model"<<endl;
//   }

///
/// \brief Run the managment strategy Evaluation routine
///


FUNCTION void runMSE()
//    cout<<"Top of runMSE"<<endl;//

//	s_iSCAMdata      s_mseData;
//	s_iSCAMvariables s_mseVars;//

//	/* SCENARIO DATA */
//	s_mseData.nStock      = ngroup;
//	s_mseData.nArea       = narea;
//	s_mseData.nSex        = nsex;
//	s_mseData.nSyr        = syr;
//	s_mseData.nNyr        = nyr;
//	s_mseData.nPyr        = nyr+10;
//	s_mseData.nSage       = sage;
//	s_mseData.nNage       = nage;
//	s_mseData.nGear       = ngear;
//	s_mseData.nFleet      = nfleet;
//	s_mseData.nFleetIndex = nFleetIndex;
//	s_mseData.dAllocation = dAllocation;
//	s_mseData.d_linf      = d_linf;
//	s_mseData.d_vonbk     = d_vonbk;
//	s_mseData.d_to        = d_to;
//	s_mseData.d_a         = d_a;
//	s_mseData.d_b         = d_b;
//	s_mseData.d_ah        = d_ah;
//	s_mseData.d_gh        = d_gh;
//	s_mseData.nCtNobs     = nCtNobs;
//	s_mseData.dCatchData  = dCatchData;
//	s_mseData.nIt         = nItNobs;
//	s_mseData.nItNobs     = n_it_nobs;
//	s_mseData.n_survey_type = n_survey_type;
//	s_mseData.dSurveyData = &d3_survey_data;
//	
//	s_mseData.nAgears     = nAgears;
//	s_mseData.nAnobs      = n_A_nobs;
//	s_mseData.nAsage      = n_A_sage;
//	s_mseData.nAnage      = n_A_nage;
//	s_mseData.dA          = &d3_A;
//	s_mseData.nWtNobs     = nWtNobs;
//	s_mseData.d3_wt_avg     = &d3_wt_avg;
//	s_mseData.d3_wt_mat     = &d3_wt_mat;//

//	s_mseData.d_iscamCntrl        = d_iscamCntrl;//

//    //

//    /* SCENARIO PARAMETERS */
//	d3_array d3_selpar(1,ngear,1,jsel_npar,1,isel_npar);
//	d3_array d3_M(1,n_ags,syr,nyr,sage,nage);
//	d3_array d3_S(1,n_ags,syr,nyr,sage,nage);
//	d3_array d3_ft(1,n_ags,1,ngear,syr,nyr);
//	d4_array d4_log_sel(1,ngear,1,n_ags,syr,nyr,sage,nage);
//	for(k=1;k<=ngear;k++)
//	{
//		d3_selpar(k) = value(sel_par(k));
//		d4_log_sel(k) = value(log_sel(k));
//	}
//	for( g = 1; g <= n_ags; g++ )
//	{
//		d3_M(g) = value(M(g));
//		d3_S(g) = value(S(g));
//		d3_ft(g) = value(ft(g));
//	}//

//	s_mseVars.d_log_ro        = value(theta(1));
//	s_mseVars.d_steepness     = value(theta(2));
//	s_mseVars.d_log_m         = value(theta(3));
//	s_mseVars.d_log_rbar      = value(theta(4));
//	s_mseVars.d_log_rinit     = value(theta(5));
//	s_mseVars.d_rho           = value(theta(6));
//	s_mseVars.d_varphi        = value(theta(7));
//	s_mseVars.dLog_M_devs     = value(log_m_devs);
//	s_mseVars.dLog_Rbar_devs  = value(log_rec_devs);
//	s_mseVars.dLog_Rinit_devs = value(init_log_rec_devs);
//	s_mseVars.nSel_type       = isel_type;
//	s_mseVars.nSel_block      = sel_blocks;
//	s_mseVars.dSelPars        = &d3_selpar;
//	s_mseVars.d4_log_sel      = &d4_log_sel;
//	s_mseVars.d3_Ft           = &d3_ft;
//	s_mseVars.d3_Mt			  = &d3_M;
//	s_mseVars.d3_St           = &d3_S;//
//

//	// Instantiate the Operating Model Class
//    //OperatingModel om(s_mseData,s_mseVars);
//    cout<<"This is red leader, Im going in"<<endl;
//    OperatingModel om(s_mseData,s_mseVars,argc,argv);
//    cout<<"Use the force Luke"<<endl;
//    // Methods for the class om.
//    om.runScenario(rseed);//
//
//
//
//


TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
    
	#undef COUT
	#define COUT(object) cout << #object "\n" << object <<endl;

	#undef TINY
	#define TINY 1.e-08

	#undef NA
	#define NA -99.0

	#include <admodel.h>
	#include <time.h>
	#include <string.h>
	#include "lib/msy.h"
	#include "lib/msy.hpp"
	#include "lib/baranov.h"
    #include "lib/LogisticNormal.h"
	#include "Selex.h"
	#include "lib/msy.cpp"
	#include "lib/baranov.cpp"
	#include "lib/LogisticNormal.cpp"
	#include "lib/LogisticStudentT.cpp"
	#include "OpMod.h"

	ivector getIndex(const dvector& a, const dvector& b)
	{
		int i,j,n;
		n = 0;
		j = 1;
		for( i = a.indexmin(); i <= a.indexmax(); i++ )
		{
			// cout<<"Luke I am your father"<<endl;
			// cout<<i<<endl;
			// exit(1);
			 if(b(i) != NA) n++;
		}
		ivector tmp(1,n);
		for( i = a.indexmin(); i <= a.indexmax(); i++ )
		{
			if(b(i) != NA )
			{
				tmp(j++) = a(i);
			}
		}
		
		return(tmp);
	}

	void readMseInputs()
	  {
	  	cout<<"yep this worked"<<endl;
	  }

	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	bool mcmcPhase = 0;
	bool mcmcEvalPhase = 0;
	
	adstring BaseFileName;
	adstring ReportFileName;
	adstring NewFileName;

	adstring stripExtension(adstring fileName)
	{
		/*
		This function strips the file extension
		from the fileName argument and returns
		the file name without the extension.
		*/
		const int length = fileName.size();
		for (int i=length; i>=0; --i)
		{
			if (fileName(i)=='.')
			{
				return fileName(1,i-1);
			}
		}
		return fileName;
	}
	
	

	// class selex_vector
	// {
	// 	private:
	// 		int m_length;
			

	// 	public:

	// 	selex_vector()
	// 	{
	// 		m_length = 0;
			
	// 	}

	// 	~selex_vector()
	// 	{
			
	// 	}

		
	// };


	
	// #ifdef __GNUDOS__
	//   #include <gccmanip.h>
	// #endif
	// Variables to store results from DIC calculations.
	double dicNoPar = 0;
	double dicValue = 0;




// 	void function_minimizer::mcmc_eval(void)
// 	{
// 		// |---------------------------------------------------------------------------|
// 		// | Added DIC calculation.  Martell, Jan 29, 2013                             |
// 		// |---------------------------------------------------------------------------|
// 		// | DIC = pd + dbar
// 		// | pd  = dbar - dtheta  (Effective number of parameters)
// 		// | dbar   = expectation of the likelihood function (average f)
// 		// | dtheta = expectation of the parameter sample (average y) 

// 	  gradient_structure::set_NO_DERIVATIVES();
// 	  initial_params::current_phase=initial_params::max_number_phases;
// 	  uistream * pifs_psave = NULL;

// 	#if defined(USE_LAPLACE)
// 	#endif

// 	#if defined(USE_LAPLACE)
// 	    initial_params::set_active_random_effects();
// 	    int nvar1=initial_params::nvarcalc(); 
// 	#else
// 	  int nvar1=initial_params::nvarcalc(); // get the number of active parameters
// 	#endif
// 	  int nvar;
	  
// 	  pifs_psave= new
// 	    uistream((char*)(ad_comm::adprogram_name + adstring(".psv")));
// 	  if (!pifs_psave || !(*pifs_psave))
// 	  {
// 	    cerr << "Error opening file "
// 	            << (char*)(ad_comm::adprogram_name + adstring(".psv"))
// 	       << endl;
// 	    if (pifs_psave)
// 	    {
// 	      delete pifs_psave;
// 	      pifs_psave=NULL;
// 	      return;
// 	    }
// 	  }
// 	  else
// 	  {     
// 	    (*pifs_psave) >> nvar;
// 	    if (nvar!=nvar1)
// 	    {
// 	      cout << "Incorrect value for nvar in file "
// 	           << "should be " << nvar1 << " but read " << nvar << endl;
// 	      if (pifs_psave)
// 	      {
// 	        delete pifs_psave;
// 	        pifs_psave=NULL;
// 	      }
// 	      return;
// 	    }
// 	  }
	  
// 	  int nsamp = 0;
// 	  double sumll = 0;
// 	  independent_variables y(1,nvar);
// 	  independent_variables sumy(1,nvar);

// 	  do
// 	  {
// 	    if (pifs_psave->eof())
// 	    {
// 	      break;
// 	    }
// 	    else
// 	    {
// 	      (*pifs_psave) >> y;
// 	      sumy = sumy + y;
// 	      if (pifs_psave->eof())
// 	      {
// 	      	double dbar = sumll/nsamp;
// 	      	int ii=1;
// 	      	y = sumy/nsamp;
// 	      	initial_params::restore_all_values(y,ii);
// 	        initial_params::xinit(y);   
// 	        double dtheta = 2.0 * get_monte_carlo_value(nvar,y);
// 	        double pd     = dbar - dtheta;
// 	        double dic    = pd + dbar;
// 	        dicValue      = dic;
// 	        dicNoPar      = pd;

// 	        cout<<"Number of posterior samples    = "<<nsamp    <<endl;
// 	        cout<<"Expectation of log-likelihood  = "<<dbar     <<endl;
// 	        cout<<"Expectation of theta           = "<<dtheta   <<endl;
// 	        cout<<"Number of estimated parameters = "<<nvar1    <<endl;
// 		    cout<<"Effective number of parameters = "<<dicNoPar <<endl;
// 		    cout<<"DIC                            = "<<dicValue <<endl;
// 	        break;
// 	      }
// 	      int ii=1;
// 	      initial_params::restore_all_values(y,ii);
// 	      initial_params::xinit(y);   
// 	      double ll = 2.0 * get_monte_carlo_value(nvar,y);
// 	      sumll    += ll;
// 	      nsamp++;
// 	      // cout<<sumy(1,3)/nsamp<<" "<<get_monte_carlo_value(nvar,y)<<endl;
// 	    }
// 	  }
// 	  while(1);
// 	  if (pifs_psave)
// 	  {
// 	    delete pifs_psave;
// 	    pifs_psave=NULL;
// 	  }
// 	  return;
// 	}
	
	
FINAL_SECTION
	cout<<"Here I am "<<endl;
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"--Number of function evaluations: "<<nf<<endl;
	cout<<"--Results are saved with the base name:\n"<<"\t"<<BaseFileName<<endl;
	cout<<"*******************************************"<<endl;

		

	if(mseFlag) runMSE();
	cout<<"End of class testing"<<endl;


	//Make copies of the report file using the ReportFileName
	//to ensure the results are saved to the same directory 
	//that the data file is in. This should probably go in the 
	//FINAL_SECTION
	
	//CHANGED only copy over the mcmc files if in mceval_phase()
	
	#if defined __APPLE__ || defined __linux
	if(last_phase() && !retro_yrs)
	{
		adstring bscmd = "cp iscam.rep " +ReportFileName;
		system(bscmd);
		
		bscmd = "cp iscam.par " + BaseFileName + ".par";
		system(bscmd); 
		
		bscmd = "cp iscam.std " + BaseFileName + ".std";
		system(bscmd);
		
		bscmd = "cp iscam.cor " + BaseFileName + ".cor";
		system(bscmd);
		
		if( SimFlag )
		{
			bscmd = "cp iscam.sim " + BaseFileName + ".sim";
			system(bscmd);
		}

		if( mcmcPhase )
		{
			bscmd = "cp iscam.psv " + BaseFileName + ".psv";
			system(bscmd);
			
			cout<<"Copied binary posterior sample values"<<endl;
		}
		
		if( mcmcEvalPhase )
		{		
			bscmd = "cp iscam.mcmc " + BaseFileName + ".mcmc";
			system(bscmd);
		
			bscmd = "cp sbt.mcmc " + BaseFileName + ".mcst";
			system(bscmd);
		
			bscmd = "cp rt.mcmc " + BaseFileName + ".mcrt";
			system(bscmd);
			
			ofstream mcofs(ReportFileName,ios::app);
			mcofs<<"ENpar\n"<<dicNoPar<<endl;
			mcofs<<"DIC\n"<<dicValue<<endl;
			mcofs.close();
			cout<<"Copied MCMC Files"<<endl;
		}
	}

	if( last_phase() && retro_yrs )
	{
		//copy report file with .ret# extension for retrospective analysis
		adstring bscmd = "cp iscam.rep " + BaseFileName + ".ret" + str(retro_yrs);
		system(bscmd);
	}
	#endif

	#if defined _WIN32 || defined _WIN64
	if(last_phase() && !retro_yrs)
	{
		adstring bscmd = "copy iscam.rep " +ReportFileName;
		system(bscmd);
		
		bscmd = "copy iscam.par " + BaseFileName + ".par";
		system(bscmd); 
		
		bscmd = "copy iscam.std " + BaseFileName + ".std";
		system(bscmd);
		
		bscmd = "copy iscam.cor " + BaseFileName + ".cor";
		system(bscmd);
		
		if( mcmcPhase )
		{
			bscmd = "copy iscam.psv " + BaseFileName + ".psv";
			system(bscmd);
			
			cout<<"Copied binary posterior sample values"<<endl;
		}
		
		if( mcmcEvalPhase )
		{		
			bscmd = "copy iscam.mcmc " + BaseFileName + ".mcmc";
			system(bscmd);
		
			bscmd = "copy sbt.mcmc " + BaseFileName + ".mcst";
			system(bscmd);
		
			bscmd = "copy rt.mcmc " + BaseFileName + ".mcrt";
			system(bscmd);
		
			cout<<"Copied MCMC Files"<<endl;
		}
	}

	if( last_phase() && retro_yrs )
	{
		//copy report file with .ret# extension for retrospective analysis
		adstring bscmd = "copy iscam.rep " + BaseFileName + ".ret" + str(retro_yrs);
		system(bscmd);
	}
	#endif


