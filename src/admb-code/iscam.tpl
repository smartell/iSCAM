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
// TODO:    add gtg options for length based fisheries                          //
// CHANGED add time varying natural mortality rate with splines                  //
// TODO:    add cubic spline interpolation for time varying M                   //
// CHANGED  Fix the type 6 selectivity implementation. not working.              //
// TODO:  fix cubic spline selectivity for only years when data avail            //
// CHANGED: fixed a bug in the simulation model log_ft_pars goes out             //
//        of bounds.                                                             //
// TODO: write a projection routine and verify equilibrium calcs                 //
// TODO: add DIC calculation for MCMC routines (in -mcveal phase)                //
// CHANGED: add SOK fishery a) egg fishing mort 2) bycatch for closed ponds      //
//                                                                               //
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
//--             - calcFisheryObservations                                     --//
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
//--                                                                           --//
//--                                                                           --//
// ----------------------------------------------------------------------------- //



DATA_SECTION
	// ------------------------------------------------------------------------- //
	// In the DATA_SECTION 3 separate files are read in:                         //
	// 1) ProjectFileControl.pfc (used for stock projections under TAC)          //
	// 2) DataFile.dat           (data to condition the assessment model on)     //
	// 3) ControlFile.ctl        (controls for phases, selectivity options )     //
	//                                                                           //
	// NOTES:                                                                    //
	//                                                                           //
	// ------------------------------------------------------------------------- //

	!! cout<<"iSCAM has detected that you are on a "<<PLATFORM<<" box"<<endl;
	
	// |---------------------------------------------------------------------------------|
	// | STRINGS FOR INPUT FILES                                                         |
	// |---------------------------------------------------------------------------------|
	// |
	init_adstring DataFile;
	init_adstring ControlFile;
	init_adstring ProjectFileControl;

	
	// ------------------------------------------------------------------------- //
	// READ IN PROJECTION FILE CONTROLS                                          //
	// ------------------------------------------------------------------------- //	
	!! ad_comm::change_datafile_name(ProjectFileControl);
	init_int n_tac;
	init_vector tac(1,n_tac);
	
	// Documentation for projection control file pf_cntrl
	// 1) start year for m_bar calculation
	// 2)   end year for m_bar calculation
	// 3) start year for average fecundity/weight-at-age
	// 4)   end year for average fecundity/weight-at-age
	// 5) start year for recruitment period (not implemented yet)
	// 6)   end year for recruitment period (not implemented yet)
	init_int n_pfcntrl;
	init_vector pf_cntrl(1,n_pfcntrl);
	
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
	
	
	!! BaseFileName=stripExtension(ControlFile);
	!! cout<<BaseFileName<<endl;
	!! ReportFileName = BaseFileName + adstring(".rep");
	
	
	
	int SimFlag;
	int rseed;
	int retro_yrs;
	LOC_CALCS
		SimFlag=0;
		rseed=999;
		int on,opt;
		//the following line checks for the "-SimFlag" command line option
		//if it exists the if statement retreives the random number seed
		//that is required for the simulation model
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			SimFlag=1;
			rseed=atoi(ad_comm::argv[on+1]);
			//if(SimFlag)exit(1);
		}
		
		// command line option for retrospective analysis. "-retro retro_yrs"
		retro_yrs=0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-retro",opt))>-1)
		{
			retro_yrs=atoi(ad_comm::argv[on+1]);
			cout<<"|____________________________________________________|\n";
			cout<<"| Implementing Retrospective analysis                |\n";
			cout<<"|____________________________________________________|\n";
			cout<<"| Number of retrospective years = "<<retro_yrs<<endl;
		}
	END_CALCS

	// ************************************************************************* //
	// ** READ IN MODEL DATA FROM  DataFile                                   ** //
	// ************************************************************************* //
	!! ad_comm::change_datafile_name(DataFile);

	
	// ------------------------------------------------------------------------- //
	// MODEL DIMENSIONS                                                          //
	// ------------------------------------------------------------------------- //
	init_int syr;
	init_int nyr;	
	init_int sage;
	init_int nage;
	vector age(sage,nage);
	!! age.fill_seqadd(sage,1);
	
	//number of gear types with unique selectivities
	init_int ngear;	
	
	LOC_CALCS
		cout<<"** __MODEL DIMENSION__ **"<<endl;
		cout<<"  syr\t"<<syr<<endl;
		cout<<"  nyr\t"<<nyr<<endl;
		cout<<"  sage\t"<<sage<<endl;
		cout<<"  nage\t"<<nage<<endl;
		cout<<"  ngear\t"<<ngear<<endl;
		cout<<"** ___________________ **"<<endl;
		
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
	
	
	
	// ------------------------------------------------------------------------- //
	// Allocation for each gear in (ngear), use 0 for survey gears.              //
	// ------------------------------------------------------------------------- //
	init_vector allocation(1,ngear);
	init_ivector catch_type(1,ngear);
	int nfleet;
	ivector fsh_flag(1,ngear);
	LOC_CALCS
		//If allocation >0 then set fish flag =1 else 0
		int k;
		allocation = allocation/sum(allocation);
		for(k=1;k<=ngear;k++)
		{
			if(allocation(k)>0)
				fsh_flag(k)=1;
			else
				fsh_flag(k)=0;
		}
		nfleet = sum(fsh_flag);
	END_CALCS
	
	ivector ifleet(1,nfleet);
	LOC_CALCS	
		int j=1;
		for(k=1; k<=ngear;k++)
		{
			if(fsh_flag(k)) ifleet(j++) = k;
		}
		cout<<"ifleet index\t"<<ifleet<<endl;
	END_CALCS
	
	//The following code has been deprecated
	//ivector ft_phz(1,ngear);
	//LOC_CALCS
	//	int k;
	//	ft_phz=1;
	//	for(k=1;k<=ngear;k++)
	//	{
	//		if(!fsh_flag(k))
	//			ft_phz(k)=-1;
	//	}
	//END_CALCS
	
	init_number fixed_m;		//FIXME: depricate this from data files
	init_number linf;
	init_number vonbk;
	init_number to;
	init_number a;
	init_number b;
	init_number ah;
	init_number gh;
	
	vector la(sage,nage);		//length-at-age
	vector wa(sage,nage);		//weight-at-age
	LOC_CALCS
	  cout<<"_________________________"<<endl;
	  COUT(linf);
	  COUT(vonbk);
	  COUT(to);
	  COUT(a);
	  COUT(b);
	  COUT(ah);
	  COUT(gh);
	  cout<<"_________________________"<<endl;
	  la=linf*(1.-exp(-vonbk*(age-to)));
	  wa=a*pow(la,b);
	  cout<<setprecision(2);		//2 decimal places for output
	  cout<<"la\n"<<la<<endl;
	  cout<<"wa\n"<<wa<<endl;
	  cout<<setprecision(5);
	END_CALCS
	
	//Time series data
	init_matrix catch_data(syr,nyr,1,ngear+1);
	matrix obs_ct(1,ngear,syr,nyr);
	int ft_count;
	int i;
	LOC_CALCS
		ft_count=0;
		for(k=1;k<=ngear;k++)
			obs_ct(k)=column(catch_data,k+1);
		
		for(k=1;k<=ngear;k++)
		{	
			for(i=syr;i<=nyr;i++)
				if( obs_ct(k,i)>0 ) ft_count++;
		}
		COUT(catch_data);
		cout<<"ft_count\n"<<ft_count<<endl;
		cout<<"last row of catch \n"<<catch_data(nyr)<<endl;
		cout<<"Ok after catch extraction"<<endl;
	END_CALCS
	
	init_int nit;
	!! cout<<"Number of surveys "<<nit<<endl;
	init_ivector nit_nobs(1,nit);
	//#survey type 
	//## 1 = survey is proportional to vulnerable numbers
	//## 2 = survey is proportional to vulnerable biomass
	//## 3 = survey is proportional to spawning biomass (e.g., herring spawn survey)
	init_ivector survey_type(1,nit);
	init_3darray survey_data(1,nit,1,nit_nobs,1,5);
	//init_matrix survey_data(1,nit,1,4);
	imatrix iyr(1,nit,1,nit_nobs);
	imatrix igr(1,nit,1,nit_nobs);
	matrix it(1,nit,1,nit_nobs);
	matrix it_wt(1,nit,1,nit_nobs);		//relative weight
	matrix it_timing(1,nit,1,nit_nobs);	//timing of the survey (0-1)
	LOC_CALCS
		for(i=1;i<=nit;i++)
		{
			iyr(i)=ivector(column(survey_data(i),1));
			igr(i)=ivector(column(survey_data(i),3));
			it(i)=column(survey_data(i),2);
			it_wt(i)=column(survey_data(i),4)+1.e-10;//add a small constant to allow 0 weight
			it_timing(i)=column(survey_data(i),5);
		}
		cout<<"Last row of the relative abundance data\n"<<survey_data(nit)(nit_nobs(nit))<<endl;
		cout<<"OK after relative abundance index"<<endl;
		/*Normalize survey weights so estimated rho parameter is consistent with sd(epsilon)/sd(rec_devs)*/
		double mean_it_wt;
		mean_it_wt = mean(it_wt);
		for(i=1;i<=nit;i++)
		{
			it_wt(i) = it_wt(i)/mean_it_wt;
		}
	END_CALCS
	
	
	//Age comps for all gears.
	init_int na_gears;	//total number of aging observations
	init_ivector na_nobs(1,na_gears);
	init_ivector a_sage(1,na_gears);	//youngest age in the ageing matrix
	init_ivector a_nage(1,na_gears);	//oldest age in the ageing matrix

	init_3darray A(1,na_gears,1,na_nobs,a_sage-2,a_nage);
	
	//Mean weight-at-age data (units are kg) (if exists)
	init_int n_wt_nobs;
	init_matrix tmp_wt_obs(1,n_wt_nobs,sage-1,nage);
	
	matrix wt_obs(syr,nyr+1,sage,nage);		//weight-at-age
	matrix wt_dev(syr,nyr+1,sage,nage);		//standardized deviations in weight-at-age
	matrix fec(syr,nyr+1,sage,nage);		//fecundity-at-age
	vector fa_bar(sage,nage);				//average fecundity-at-age for all years.
	vector avg_fec(sage,nage);				//average fecundity-at-age
	vector avg_wt(sage,nage);				//average weight-at-age
	LOC_CALCS
		int iyr;
		avg_fec.initialize();
		for(i=syr;i<=nyr+1;i++)
		{
			wt_obs(i) = wa;			
			fec(i)    = elem_prod(plogis(age,ah,gh),wt_obs(i));
		}
		//if empiracle weight-at-age data exist, the overwrite wt_obs & fec.
		for(i=1;i<=n_wt_nobs;i++)
		{
			iyr         = tmp_wt_obs(i,sage-1);  //index for year
			wt_obs(iyr) = tmp_wt_obs(i)(sage,nage);
			fec(iyr)    = elem_prod(plogis(age,ah,gh),wt_obs(iyr));
		}
		//CHANGED SM Deprecated Aug 9, 2012 for new projection file control.
		//int nfec = fec.rowmax()-fec.rowmin()+1;
		//avg_fec=colsum(fec)/nfec;
		
		// SM Aug 9, 2012 calculate averge weight-at-age and
		// fecundity-at-age based on years specificed in the
		// projection control file (pf_cntrl(3-4)). Also set
		// the average weight-at-age in nyr+1 to the same.
		// Deprecate the 5-year average that Jake suggested.
		avg_wt   = colsum(wt_obs.sub(pf_cntrl(3),pf_cntrl(4)));
		avg_wt  /= pf_cntrl(4)-pf_cntrl(3)+1;
		wt_obs(nyr+1) = avg_wt;
		
		avg_fec  = colsum(fec.sub(pf_cntrl(3),pf_cntrl(4)));
		avg_fec /= pf_cntrl(4)-pf_cntrl(3)+1;
		
		fa_bar   = colsum(fec.sub(syr,nyr))/(nyr-syr+1);
		cout<<avg_wt<<endl;
		cout<<"avg_fec\n"<<avg_fec<<endl;
		cout<<"fa_bar\n"<<fa_bar<<endl;
		//DEPRECATED SM AUG 9, 2012
		//from Jake Schweigert: use mean-weight-at-age data
		//from the last 5 years for the projected mean wt.
		//dvector tmp=colsum(wt_obs.sub(nyr-5,nyr))/6.;
		//wt_obs(nyr+1) = tmp;
		//
		///*June 8, 2012, average wt at age for all years*/
		//tmp = colsum(wt_obs.sub(syr,nyr))/(nyr-syr+1.);
		//avg_wt = tmp;
		
		cout<<"n_wt_nobs\t"<<n_wt_nobs<<endl;
		cout<<"Ok after empiracle weight-at-age data"<<endl;
		
		//July 14, 2011 Error handler to ensure non-zero wt_obs
		//in the data. Otherwise causes and error with weight-based
		//selectivities.
		if(min(wt_obs)==0)
		{
			cout<<"Cannont have a observed 0 mean weight at age\n";
			cout<<"in the data file.  Please fix.\n Aborting program!"<<endl;
			exit(2);
		}
		
		
		
		//May 5, 2011 SJDM: Calculating standardized deviates
		//in mean weights-at-age from empiracle data to use
		//with selectivity as a function of mean weight at age.
		//Idea borrowed from Vivian Haist in HCAM model.
		/*
			CHANGED: selectivity function based on average weight.
			Based on the logistic function;
			1) calculate a matrix of standardized deviates
			wt_dev = (wt_obs-mean(wt_obs))/sd(wt_obs);
			
			Not implemented, using the version provided by Vivian.
			
			Aug 5, 2011, noticed that Vivian's method was implemented
			in an awkward way where GN selectivities were a random walk
			sel(t+1)=sel(t) + tmp(2);
			
			Trying a compromize where estimating an additional parameter
			that attemps to explain residual variation in age-comps
			via changes in standardized mean weights at age. This is 
			implemented as:
			
			log_sel = log(plogis(age,s1,s2)) + s3*delta
			where delta is a matrix of standardized deviations
			(mu=0, std=1) of weights at age.  delta is calculated
			based on the wt_dev matrix above.
		*/
		wt_dev.initialize();
		dmatrix mtmp = trans(wt_obs);
		for(i=sage;i<=nage;i++)
		{
			dvector wa_dev = (mtmp(i)-mean(mtmp(i)))/sqrt(var(mtmp(i)));
			mtmp(i) = wa_dev;
		}
		wt_dev = trans(mtmp);	//each column has mean=0 sd=1
		
	END_CALCS
	
	//End of data file
	init_int eof;	
	LOC_CALCS
	  cout<<"eof = "<<eof<<endl;
	  if(eof==999){
		cout<<"\n -- END OF DATA SECTION -- \n"<<endl;
	  }else{
		cout<<"\n *** ERROR READING DATA *** \n"<<endl; exit(1);
	  }
	END_CALCS

	
	
	
	vector fmsy(1,nfleet);			//Fishing mortality rate at Fmsy
	vector fall(1,nfleet);			//Fishing mortality based on allocation
	vector  msy(1,nfleet);			//Maximum sustainable yield
	number bmsy;					//Spawning biomass at MSY
	number Umsy;					//Exploitation rate at MSY
	vector age_tau2(1,na_gears);	//MLE estimate of the variance for age comps
	//catch-age for simulation model (could be declared locally 3d_array)
	3darray d3C(1,ngear,syr,nyr,sage,nage);		
	
	
	
	
	
	
	
	
	
	
	// ***************************************************
	// ** Read parameter controls from ControlFile      **
	// ***************************************************
	!! ad_comm::change_datafile_name(ControlFile);
	
	init_int npar;
	init_matrix theta_control(1,npar,1,7);
	!! COUT(npar);
	!! COUT(theta_control);
	
	vector theta_ival(1,npar);
	vector theta_lb(1,npar);
	vector theta_ub(1,npar);
	ivector theta_phz(1,npar);
	ivector theta_prior(1,npar);
	LOC_CALCS
		theta_ival = column(theta_control,1);
		theta_lb = column(theta_control,2);
		theta_ub = column(theta_control,3);
		theta_phz = ivector(column(theta_control,4));
		theta_prior=ivector(column(theta_control,5));
	END_CALCS
	
	
	// ***************************************************
	// ** Read selectivity parameter options
	// ***************************************************
	// type 1 = logistic (2pars)
	// type 2 = selcoffs (A-1 pars)
	// type 3 = cubic spline (age_nodes)
	// type 4 = time varying cubic spline (age_nodes for each year)
	// type 5 = bicubic spline with age_nodes adn yr_nodes
	// type 6 = fixed logistic by turning sel_phz to (-ve)
	// type 7 = logistic (3pars) as a function of body weight.
	init_matrix selex_controls(1,10,1,ngear);
	!! COUT(selex_controls);

	ivector isel_npar(1,ngear);				//ivector for # of parameters for each gear.
	ivector jsel_npar(1,ngear);				//ivector for the number of rows for time-varying selectivity.
	ivector isel_type(1,ngear);  	   	 	//Switch for selectivity
	vector ahat(1,ngear);			   	 	//age-at-50% vulnerbality for logistic function
	vector ghat(1,ngear);			   	 	//std at 50% age of vulnerability for logistic funciton
	vector age_nodes(1,ngear);		   	 	//No. of age-nodes for bicubic spline.
	vector yr_nodes(1,ngear);		   	 	//No. of year-nodes for bicubic spline.
	ivector sel_phz(1,ngear);		   	 	//Phase for estimating selectivity parameters.
	vector sel_2nd_diff_wt(1,ngear);	   	//Penalty weight for 2nd difference in selectivity.
	vector sel_dome_wt(1,ngear);		   	//Penalty weight for dome-shaped selectivity.
	vector sel_2nd_diff_wt_time(1,ngear); 	//Penalty weight for 2nd difference in time-varying selectivity.
	ivector n_sel_blocks(1,ngear);			//Number of discrete selectivity blocks.
	
	LOC_CALCS
		isel_type = ivector(selex_controls(1));
		ahat      = selex_controls(2);
		ghat      = selex_controls(3);
		age_nodes = selex_controls(4);
		yr_nodes  = selex_controls(5);

		sel_phz   = ivector(selex_controls(6));
		sel_2nd_diff_wt = selex_controls(7);
		sel_dome_wt     = selex_controls(8);
		sel_2nd_diff_wt_time = selex_controls(9);
		n_sel_blocks    = ivector(selex_controls(10));
	END_CALCS
	!! COUT(isel_type);
	!! COUT(sel_dome_wt);
	// Selectivity blocks for each gear.
	init_imatrix sel_blocks(1,ngear,1,n_sel_blocks);
	!! COUT(sel_blocks);

	LOC_CALCS
		//cout up the number of selectivity parameters
		//depending on the value of isel_type
		//testing();
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
					// Alternative logistic selectivity with wt_dev coefficients.
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
					isel_npar(i) = (nage-sage);
					jsel_npar(i) = n_sel_blocks(i);
					break;
					
				default: break;
			}
		}
		//cout<<"Number of estimated selectivity parameters\n"<<isel_npar<<endl;
	END_CALCS
	
	//Controls for prior on survey q.
	init_int nits;					//FIXME (redundant with nit, could be deprecated)
	init_ivector q_prior(1,nits);
	init_vector q_mu(1,nits);
	init_vector q_sd(1,nits);
	!! cout<<"nits\n"<<nits<<endl;
	!! cout<<"q Prior\n"<<q_mu<<endl<<q_sd<<endl;
	
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
	// | 
	init_vector cntrl(1,14);
	int verbose;
	

	// |---------------------------------------------------------------------------------|
	// | SIMULATION CONTROLS
	// |---------------------------------------------------------------------------------|
	// | Defn' for sim_ctrl
	// | 1) Flag for modelling selectivity using ideal free distribution.
	// | 2) Index for simulated sel_type
	// | 3) Number of selectivity blocks for each gear type.

	init_matrix sim_ctrl(1,3,1,ngear);
	!! COUT(sim_ctrl);

	ivector sim_isel_type(1,ngear);
	ivector nsim_sel_blocks(1,ngear);
	
	LOC_CALCS
		sim_isel_type   = ivector(sim_ctrl(2));
		nsim_sel_blocks = ivector(sim_ctrl(3));
	END_CALCS
	init_imatrix sim_sel_blocks(1,ngear,1,nsim_sel_blocks);
	// TODO modify calcSelectivities to include sim_sel_blocks.

	init_int eofc;
	LOC_CALCS
		verbose = cntrl(1);
		cout<<"cntrl\n"<<cntrl<<endl;
		cout<<"eofc\t"<<eofc<<endl;
		if(eofc==999){
			cout<<"\n -- END OF CONTROL FILE -- \n"<<endl;
		}else{
			cout<<"\n ***** ERROR CONTROL FILE ***** \n"<<endl; exit(1);
		}
	END_CALCS
	
	
	int nf;
	
	// |
	ivector ilvec(1,7);
	!!ilvec    = ngear;			//number of fisheries
	!!ilvec(2) = nit;			//number of surveys
	!!ilvec(3) = na_gears;		//number of age-comps
	!!ilvec(4) = 1;
	
	// SM Oct 31, 2010.  Implementing retrospective analysis.
	//Dont read in any more data below the retrospective reset of nyr
	!! nyr = nyr - retro_yrs;
	
	// SM Sept 2, 2012. If in retrospective analysis, make sure arrays from pfc 
	// are ajusted downwards, otherwise arrays for m_bar will go out of bounds.
	LOC_CALCS
		if(retro_yrs)
		{
			if(pf_cntrl(2)>nyr) pf_cntrl(2) = nyr;
			if(pf_cntrl(4)>nyr) pf_cntrl(4) = nyr;
			if(pf_cntrl(6)>nyr) pf_cntrl(6) = nyr;
		}
		if(verbose) COUT(pf_cntrl);
	END_CALCS

	
	// END OF DATA_SECTION
	!! if(verbose) cout<<"||-- END OF DATA_SECTION --||"<<endl;
	
INITIALIZATION_SECTION
 theta theta_ival;
	
PARAMETER_SECTION
	//Leading parameters
	//theta[1]		log_ro, or log_msy
	//theta[2]		steepness(h), or log_fmsy
	//theta[3]		log_m
	//theta[4]		log_avgrec
	//theta[5]		log_recinit
	//theta[6]		rho
	//theta[7]		vartheta
	
	
	init_bounded_number_vector theta(1,npar,theta_lb,theta_ub,theta_phz);
	//!! for(int i=1;i<=npar;i++) theta(i)=theta_ival(i);
	
	//Selectivity parameters (A very complicated ragged array)
	init_bounded_matrix_vector sel_par(1,ngear,1,jsel_npar,1,isel_npar,-15.,15.,sel_phz);
	LOC_CALCS
		//initial values for logistic selectivity parameters
		//set phase to -1 for fixed selectivity.
		if (!global_parfile)
		{
			for(int k=1;k<=ngear;k++)
			{
				if( isel_type(k)==1 || 
					isel_type(k)==6 || 
					isel_type(k)>=7 ||
					sim_isel_type(k)==1 ||
					sim_isel_type(k)>=7 )
				{
					for(int j = 1; j <= n_sel_blocks(k); j++ )
					{
						double uu = 0;
						if(sim_isel_type(k)==1 && j > 1)
						{
							uu = 0.05*randn(j+rseed);
						} 

						// SPECIAL CASE for simulation.
						// Special case if sim with sel_type=1 and estimate with sel_type11
						if(sim_isel_type(k)==1 && isel_type(k)>=11 && j==1)
						{
							// convert length to age
							ahat(k) = -log(-(ahat(k)-linf)/linf)/vonbk;
						}

							
						sel_par(k,j,1) = log(ahat(k)*exp(uu));
						sel_par(k,j,2) = log(ghat(k));
					}
				}
			}
		}
	END_CALCS
	
	//Fishing mortality rate parameters
	init_bounded_vector log_ft_pars(1,ft_count,-30.,3.0,1);
	
	LOC_CALCS
		if(!SimFlag) log_ft_pars = log(0.1);
	END_CALCS
	
	
	//CHANGED Vivian Haist: possible bias in Bo due to dev_vector for log_rec_devs
	//	-Try and compare estimates using init_bounded_vector version.
	//	-Now estimating Rinit, Rbar and Ro.
	
	// May 16, Just looked at Ianelli's code, he used bounded_vector, and 
	// separates the problem into init_log_rec_devs and log_rec_devs
	
	
	//Annual recruitment deviations
	//!! int ii;
	//!! ii=syr-nage+sage;
	//!! if(cntrl(5)) ii = syr;  //if initializing with ro
	//init_bounded_vector log_rec_devs(ii,nyr,-15.,15.,2);
	
	!! int init_dev_phz = 2;
	!! if(cntrl(5)) init_dev_phz = -1;
	init_bounded_vector init_log_rec_devs(sage+1,nage,-15.,15.,init_dev_phz);
	init_bounded_vector log_rec_devs(syr,nyr,-15.,15.,2);
	
	//Deviations for natural mortality
	!! int m_dev_phz = -1;
	!! m_dev_phz = cntrl(10);
	!! int n_m_devs = cntrl(12);
	init_bounded_vector log_m_nodes(1,n_m_devs,-5.0,5.0,m_dev_phz);
	//init_bounded_vector log_m_devs(syr+1,nyr,-5.0,5.0,m_dev_phz);
	
	objective_function_value f;
    
	number ro;					//unfished age-1 recruits
	number bo;					//unfished spawning stock biomass (reference point)
	number sbo;					//unfished spawning biomass at time of spawning from SR curve
	number kappa;				//Goodyear compensation ratio
	number m;					//initial natural mortality rate
	number m_bar;				//average natural mortality rate
	number log_avgrec;			//log of average recruitment.
	number log_recinit;			//log of initial recruitment in syr.
	//number log_avg_f;			//log of average fishing mortality DEPRICATED
	number rho;					//proportion of the observation error
	number varphi				//total precision in the CPUE & Rec anomalies.
	number so;
	number beta;
	
	vector log_rt(syr-nage+sage,nyr);
	vector vax(sage,nage);		//survey selectivity coefficients
	vector q(1,nit);			//survey catchability coefficients
	
	vector sbt(syr,nyr+1);		//spawning stock biomass
	vector rt(syr+sage,nyr); 	//predicted sage recruits from S-R curve
	
	vector delta(syr+sage,nyr);	//residuals for stock recruitment
	vector avg_log_sel(1,ngear);//conditional penalty for objective function
	vector log_m_devs(syr+1,nyr);// log deviations in natural mortality
	
	matrix nlvec(1,7,1,ilvec);	//matrix for negative loglikelihoods
	
	//matrix jlog_sel(1,ngear,sage,nage);		//selectivity coefficients for each gear type.
	//matrix log_sur_sel(syr,nyr,sage,nage);	//selectivity coefficients for survey.
	 
	matrix N(syr,nyr+1,sage,nage);			//Numbers at age
	matrix F(syr,nyr,sage,nage);			//Age-specific fishing mortality
	matrix M_tot(syr,nyr,sage,nage);		//Age-specific natural mortality
	matrix ft(1,ngear,syr,nyr);				//Gear specific fishing mortality rates
	matrix log_ft(1,ngear,syr,nyr);			//Gear specific log fishing mortlity rates
	matrix Z(syr,nyr,sage,nage);
	matrix S(syr,nyr,sage,nage);
	matrix ct(1,ngear,syr,nyr);				//predicted catch biomass
	matrix eta(1,ngear,syr,nyr);			//residuals for catch
	matrix epsilon(1,nit,1,nit_nobs);		//residuals for survey abundance index
	matrix pit(1,nit,1,nit_nobs);			//predicted relative abundance index
	matrix qt(1,nit,1,nit_nobs);			//catchability coefficients (time-varying)
	

	3darray Ahat(1,na_gears,1,na_nobs,a_sage-2,a_nage);		//predicted age proportions by gear & year
	3darray A_nu(1,na_gears,1,na_nobs,a_sage-2,a_nage);		//residuals for age proportions by gear & year
	
	3darray log_sel(1,ngear,syr,nyr,sage,nage);		//selectivity coefficients for each gear type.
	3darray Chat(1,ngear,syr,nyr,sage,nage);		//predicted catch-at-age
	
	sdreport_number sd_depletion;	
	
PRELIMINARY_CALCS_SECTION
	// |---------------------------------------------------------------------------------|
	// | Run the model with input parameters to simulate real data.
	// |---------------------------------------------------------------------------------|
	// | 1) get Simulation controls from *.sim
 	// | 
  nf=0;
  cout<<"OK TO HERE"<<endl;
  if(SimFlag) 
  {
    initParameters();
    simulationModel(rseed);
  }
  if(verbose) cout<<"||-- END OF PRELIMINARY_CALCS_SECTION --||"<<endl;

RUNTIME_SECTION
    maximum_function_evaluations 100,200,500,25000,25000
    convergence_criteria 0.01,0.01,1.e-5,1.e-5


PROCEDURE_SECTION
	initParameters();
	
	calcSelectivities(isel_type);
	
	calcTotalMortality();
	
	calcNumbersAtAge();
	
	calcFisheryObservations();
	
	calcAgeProportions();
	
	calcSurveyObservations();
	
	calcStockRecruitment();
	
	calc_objective_function();

	sd_depletion=sbt(nyr)/sbo;
	
	if(mc_phase())
	{
		mcmcPhase=1;
	}
	
	if(mceval_phase())
	{
		mcmcEvalPhase=1;
		mcmc_output();
	}
	
	//The following causes a linker error
	//duplicate symbol in libdf1b2o.a
	//dvariable a=3.0;
	//cout<<"testing gammln(dvariable)"<<gammln(a)<<endl;
	

FUNCTION initParameters
  {
	/*
	This function is used to extract the specific parameter values
	from the init_bounded_number_vector to the specific variables
	used in the code.
	
	Note that you must call this routine before runnning the 
	simulation model to generate fake data.
	*/
	
	ro          = mfexp(theta(1));
	dvariable h = theta(2);
	m           = mfexp(theta(3));
	log_avgrec  = theta(4);
	log_recinit = theta(5);
	rho         = theta(6);
	varphi      = theta(7);
	
	
	switch(int(cntrl(2)))
	{
		case 1:
			//Beverton-Holt model
			kappa = (4.*h/(1.-h));
			break;
		case 2:
			//Ricker model
			kappa = pow((5.*h),1.25);
		break;
	}
	
	
	//TODO Alternative parameterization using MSY and FMSY as leading parameters
	
	

	
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

FUNCTION dvar_matrix cubic_spline_matrix(const dvar_matrix& spline_coffs)
  {
	RETURN_ARRAYS_INCREMENT();
	int nodes= spline_coffs.colmax()-spline_coffs.colmin()+1;
	int rmin = spline_coffs.rowmin();
	int rmax = spline_coffs.rowmax();
	
	dvector ia(1,nodes);
	dvector fa(sage,nage);
	ia.fill_seqadd(0,1./(nodes-1));
	//fa.fill_seqadd(sage,1);
	fa.fill_seqadd(0,1./(nage-sage));
	vcubic_spline_function_array fna(rmin,rmax,ia,spline_coffs);
	RETURN_ARRAYS_DECREMENT();
	return(fna(fa));
	
  }


FUNCTION void calcSelectivities(const ivector& isel_type)
  {
	/*
		This function loops over each ngear and calculates the corresponding
		selectivity for that gear type. It first uses a switch statement 
		to calculate selectivity curves based on isel_type where:
		1) logistic selectivity with 2 parameters
		2) age-specific selectivity coefficients with (nage-sage) parameters
		   and the last two age-classes are assumed to have the same selectivity.
		3) a reduced age-specific parameter set based on a bicubic spline.
		4) Time varying cubic spline.
		5) Time varying bicubic spline (2d version)
		6) Fixed logistic
		7) logistic selectivity based on relative changes in mean weight at age
		8) Time varying selectivity based on logistic with deviations in 
		   weights at age (3 estimated parameters).
		11) logistic selectivity with 2 parameters based on mean length
		12) length-based selectivity using cubic spline interpolation
		
		Following the initialization of the selectivity curves, time-varying 
		considerations are implemented.
		
		CHANGED: Add penality (10.*square(avg_log_sel)) to objective function 
		in cases where estimating sel_coffs to mimic an init_bounded_dev vector.
		
		CHANGED: Problem with case 7: turns out to be a random walk, so there
		is changes in selectivity even with constant growth.  Remove the
		random walk outside the switch statement.
		
		TODO: add an option for length-based selectivity.  Use inverse of
		allometric relationship w_a = a*l_a^b; to get mean length-at-age from
		empirical weight-at-age data, then calculate selectivity based on 
		mean length.
	*/

	int i,j,k,byr,bpar;
	double tiny=1.e-10;
	dvariable p1,p2,p3;
	dvar_vector age_dev=age;
	dvar_matrix t1;
	dvar_matrix tmp(syr,nyr-1,sage,nage);
	dvar_matrix tmp2(syr,nyr,sage,nage);
	dvar_matrix ttmp2(sage,nage,syr,nyr);
	//jlog_sel.initialize();
	log_sel.initialize();
	avg_log_sel.initialize();
	
	for(j=1;j<=ngear;j++)
	{
		tmp.initialize(); tmp2.initialize();
		dvector iy(1,yr_nodes(j));
		dvector ia(1,age_nodes(j));
		byr = 1;
		bpar = 0;
		switch(isel_type(j))
		{
			case 1:
				// logistic selectivity for case 1 or 6
				for(i=syr; i<=nyr; i++)
				{
					if( i == sel_blocks(j,byr) )
					{
						bpar ++;
						if( byr < n_sel_blocks(j) ) byr++;
					}

					p1 = mfexp(sel_par(j,bpar,1));
					p2 = mfexp(sel_par(j,bpar,2));
					log_sel(j)(i) = log( plogis(age,p1,p2)+tiny );
				}
				break;
			
			case 6:
				// logistic selectivity for case 1 or 6
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				for(i=syr; i<=nyr; i++)
				{
					log_sel(j)(i) = log( plogis(age,p1,p2) );
				}
				break;
				
			case 2:		
				// age-specific selectivity coefficients
				for(i=syr; i<=nyr; i++)
				{
					if( i == sel_blocks(j,byr) )
					{
						bpar ++;
						if( byr < n_sel_blocks(j) ) byr++;
					}
					for(k=sage;k<=nage-1;k++)
					{
						log_sel(j)(i)(k)   = sel_par(j)(bpar)(k-sage+1);
					}
					log_sel(j)(i,nage) = log_sel(j)(i,nage-1);
				}
				break;
				
			case 3:		
				// cubic spline
				// log_sel(j)(syr)=cubic_spline( sel_par(j)(1) );
				for(i=syr; i<nyr; i++)
				{
					if( i==sel_blocks(j,byr) )
					{
						bpar ++;	
						log_sel(j)(i)=cubic_spline( sel_par(j)(bpar) );
						if( byr < n_sel_blocks(j) ) byr++;
					}
					log_sel(j)(i+1) = log_sel(j)(i);
				}
				break;
				
			case 4:		
				// time-varying cubic spline every year
				for(i=syr; i<=nyr; i++)
				{
					log_sel(j)(i) = cubic_spline(sel_par(j)(i-syr+1));
				}
				
				//jlog_sel(j) = cubic_spline(sel_par(j)(1));
				//t1 = cubic_spline_matrix(sel_par(j).sub(2,jsel_npar(j)));
				//for(i = t1.indexmin(); i <= t1.indexmax(); i++)
				//{
				//	tmp( syr+(i-t1.indexmin()) ) = t1(i);
				//}
				break;
				
			case 5:		
				// time-varying bicubic spline
				ia.fill_seqadd( 0,1./(age_nodes(j)-1) );
				iy.fill_seqadd( 0,1./(yr_nodes(j)-1) );
				// bicubic_spline function is located in stats.cxx library
				bicubic_spline( iy,ia,sel_par(j),tmp2 );
				log_sel(j) = tmp2; 
				break;
				
			case 7:
				// time-varying selectivity based on deviations in weight-at-age
				// CHANGED This is not working and should not be used. (May 5, 2011)
				// SJDM:  I was not able to get this to run very well.
				// AUG 5, CHANGED so it no longer has the random walk component.
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				
				for(i = syr; i<=nyr; i++)
				{
					dvar_vector tmpwt=log(wt_obs(i)*1000)/mean(log(wt_obs*1000.));
					log_sel(j)(i) = log( plogis(tmpwt,p1,p2)+tiny );
				}	 
				break;
				
			case 8:
				//Alternative time-varying selectivity based on weight 
				//deviations (wt_dev) wt_dev is a matrix(syr,nyr+1,sage,nage)
				//p3 is the coefficient that describes variation in log_sel.
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				p3 = sel_par(j,1,3);
				
				for(i=syr; i<=nyr; i++)
				{
					tmp2(i) = p3*wt_dev(i);
					log_sel(j)(i) = log( plogis(age,p1,p2)+tiny ) + tmp2(i);
				}
				break;
				
			case 11:
				//logistic selectivity based on mean length-at-age				
				for(i=syr; i<=nyr; i++)
				{
					if( i == sel_blocks(j,byr) )
					{
						bpar ++;
						if( byr < n_sel_blocks(j) ) byr++;
					}
					p1 = mfexp(sel_par(j,bpar,1));
					p2 = mfexp(sel_par(j,bpar,2));

					dvector len = pow(wt_obs(i)/a,1./b);

					log_sel(j)(i) = log( plogis(len,p1,p2) );
					// COUT(exp(log_sel(j)(i)- log(mean(exp(log_sel(j)(i))))) );
				}
				break;
				
			case 12:
				//length-specific selectivity coefficients
				//based on cubic spline interpolation
				for(i=syr; i<=nyr; i++)
				{
					if( i == sel_blocks(j,byr) )
					{
						bpar ++;
						if( byr < n_sel_blocks(j) ) byr++;
					}
					dvector len = pow(wt_obs(i)/a,1./b);
					log_sel(j)(i)=cubic_spline( sel_par(j)(bpar), len );
				}
				break;
				
				
			default:
				log_sel(j)=0;
				break;
				
		}  // switch
		

		//log_sel(j)(syr) = jlog_sel(j)+tmp2(syr);
		//for(i=syr;i<nyr;i++)
		//{
		//	log_sel(j)(i+1)(sage,nage) = log_sel(j)(i)(sage,nage)+tmp(i)+tmp2(i+1);			
		//}
		
		//subtract mean to ensure mean(exp(log_sel))==1
		//substract max to ensure exp(log_sel) ranges from 0-1
		for(i=syr;i<=nyr;i++)
		{
			log_sel(j)(i) -= log( mean(mfexp(log_sel(j)(i))) );
			// log_sel(j)(i) -= log(max(mfexp(log_sel(j)(i))));
		}
			
		//cout<<"log_sel \t"<<j<<"\n"<<log_sel(j)<<"\n \n"<<endl;
		//testing bicubic spline  (SM Checked OCT 25,2010.  Works on the example below.)
		/*ia.fill_seqadd(0,1./(age_nodes(j)-1));
				iy.fill_seqadd(0,1./(yr_nodes(j)-1));
				dvar_matrix tn(1,age_nodes(j),1,yr_nodes(j));
				tn.colfill_seqadd(1,-.5,0.1);
				tn.colfill_seqadd(2,-.4,0.1);
				bicubic_spline(iy,ia,tn,tmp2);
				cout<<ia<<endl;
				cout<<iy<<endl;
				cout<<tn<<endl;
				cout<<tmp2<<endl;
				exit(1);*/
	
	}
	
	if(verbose)cout<<"**** Ok after calcSelectivities ****"<<endl;
	
  }	
	
FUNCTION calcTotalMortality
  {
	/*
	This routine calculates fishing mortality, total mortality
	and annaul survival rates (exp(-Z)) for each age in each
	year.
	
	There is a complication in that if there is a survey gear
	then the fishing mortality rate is set to an extremely small
	value: exp(-70.)~3.975e-31, which is effectively 0.
	
	
	The above issue is no longer an issue b/c now estimating Fs
	for years where there is non-zero catch.
	
	SJDM.  Dec 24, 2010.  Adding time-varying natural mortality
	
	CHANGED May 20, 2011  Add cubic spline to the time-varying natural mortality
	
	Jan 5, 2012 Adding catch_type to allow for catch in numbers, weight or spawn.
	In the case of spawn on kelp (roe fisheries), the Fishing mortality does not
	occur on the adult component.  Added if(catch_type(k)!=3) //exclude roe fisheries
	*/
	int i,k,ki;
	dvariable ftmp;
	F.initialize();
	ft.initialize();
	log_ft.initialize();
	
	//Fishing mortality
	ki=1;
	for(k=1;k<=ngear;k++)
	{
		
		for(i=syr;i<=nyr;i++)
		{	
			ftmp=0;
			if( obs_ct(k,i)>0  )
				ftmp = mfexp(log_ft_pars(ki++));
			
			ft(k,i)=ftmp;
			
			if(catch_type(k)!=3){	//exclude roe fisheries
				F(i)+=ftmp*mfexp(log_sel(k)(i));
				//cout<<obs_ct(k,i)<<endl;
			}
		}
	}
	
	//Natural mortality (year and age specific)
	//M_tot(syr,nyr,sage,nage);
	M_tot = m;


	// Cubic spline to interpolate log_m_devs (log_m_nodes)
	log_m_devs = 0.;
	if(active(log_m_nodes))
	{
		int nodes = size_count(log_m_nodes);
		dvector im(1,nodes);
		dvector fm(syr+1,nyr);
		im.fill_seqadd(0,1./(nodes-1));
		fm.fill_seqadd(0,1./(nyr-syr));
		vcubic_spline_function m_spline(im,log_m_nodes);
		//m_spline(fm);
		log_m_devs = m_spline(fm);
	}
	
	// Random walk in natural mortality.
	for(i=syr;i<=nyr;i++)
	{
		// if(active(log_m_devs)&&i>syr)
		if(active(log_m_nodes)&&i>syr)
		{
			M_tot(i)=M_tot(i-1)*exp(log_m_devs(i));
		}
	}
	m_bar = mean( M_tot.sub(pf_cntrl(1),pf_cntrl(2)) );
	
	Z=M_tot+F;
	S=mfexp(-Z);
	if(verbose) cout<<"**** OK after calcTotalMortality ****"<<endl;
	
  }
	
	
FUNCTION calcNumbersAtAge
  {
	/*
		Aug 9, 2012.  Made a change here to initialize the numbers
		at age in syr using the natural mortality rate at age in syr. 
		Prior to this the average (m_bar) rate was used, since this 
		has now changed with new projection control files.  Should only
		affect models that were using time varying natural mortality.
	*/
	
	int i,j;
	N.initialize();
	dvariable avg_M = mean(M_tot(syr));
	dvar_vector lx(sage,nage);
	lx(sage)=1;
	for( j=sage+1; j<=nage;j++) 
	{
		lx(j) = lx(j-1)*exp(-avg_M);
	}
	lx(nage) /= (1.-exp(-avg_M));
	
	if(cntrl(5))    //If initializing in at unfished conditions
	{	
		log_rt(syr) = log(ro);
		N(syr)      = ro * lx;
		//SM Deprecated Aug 9, 2012
		//for(j=sage;j<=nage;j++)
		//{
		//	N(syr,j)=ro*exp(-m_bar*(j-1.));    
		//}
	}
	else            //If starting at unfished conditions
	{
		log_rt(syr)         = log_avgrec+log_rec_devs(syr);
		N(syr,sage)         = mfexp(log_rt(syr));
		dvar_vector tmpr    = mfexp(log_recinit + init_log_rec_devs(sage+1,nage));
		N(syr)(sage+1,nage) = elem_prod(tmpr,lx(sage+1,nage));
		//SM Deprecated Aug 9, 2012
		//for(j=sage+1;j<=nage;j++)
		//{
		//	N(syr,j)=mfexp(log_recinit+init_log_rec_devs(j))*exp(-m_bar*(j-sage));
		//}
	}
	// SM Depreceated Aug 9, 2012
	//N(syr,nage)/=(1.-exp(-m_bar));
	
	//initial number of sage recruits from year syr+1, nyr;
	for(i=syr+1;i<=nyr;i++){
		log_rt(i)=log_avgrec+log_rec_devs(i);
		N(i,sage)=mfexp(log_rt(i));
	}
	N(nyr+1,sage)=mfexp(log_avgrec);
	
	/* Dynamic state variables */
	for(i=syr;i<=nyr;i++)
	{
		N(i+1)(sage+1,nage)=++elem_prod(N(i)(sage,nage-1),S(i)(sage,nage-1))+1.e-10;
		N(i+1,nage)+=N(i,nage)*S(i,nage);
	}
	if(verbose)cout<<"**** Ok after calcNumbersAtAge ****"<<endl;
	
  }

FUNCTION calcAgeProportions
  {
	/*This function loops over each gear and year
	and calculates the predicted proportions at age
	sampled based on the selectivity of that gear and
	the numbers-at-age in the population.*/
	
	int i,k,iyr,ig;
	for(k=1;k<=na_gears;k++)
	{
		for(i=1;i<=na_nobs(k);i++)
		{
			iyr=A(k,i,a_sage(k)-2);	//index for year
			ig=A(k,i,a_sage(k)-1);	//index for gear
			if(iyr>nyr)break;		//trap for retrospective analysis
			
			A_nu(k,i,a_sage(k)-2)=iyr;
			A_nu(k,i,a_sage(k)-1)=ig;
			Ahat(k,i,a_sage(k)-2)=iyr;
			Ahat(k,i,a_sage(k)-1)=ig;
			Ahat(k)(i)(a_sage(k),a_nage(k))=Chat(k)(iyr)(a_sage(k),a_nage(k))
										/sum(Chat(k)(iyr)(a_sage(k),a_nage(k)));
		}
	}
	if(verbose)cout<<"**** Ok after calcAgeProportions ****"<<endl;

  }	

FUNCTION calcFisheryObservations
  {
	/*
	Dec 6, 2010.  Modified ct calculations to include
				  empirical weight at age data (wt_obs);
				
	Jan 16, 2011. 	modified this code to get age-comps from surveys, rather than 
					computing the age-comps in calc_fisheries_observations
	
	Jan 6, 2012. 	modified code to allow for catch observations in numbers,
					biomass, and harvest of roe.  
	
	Jun 22, 2012.	added eta variable for catch residuals.
	
	*/
	
	/*
		FIXED Reconcile the difference between the predicted catch 
		here and in the simulation model.
	*/
	int i,k;
	ct.initialize();
	eta.initialize();
	for(i=syr;i<=nyr;i++)
	{
		for(k=1;k<=ngear;k++)
		{
			
			dvar_vector log_va=log_sel(k)(i);
			
			
			//SJDM Jan 16, 2011 Modification as noted above.
			//SJDM Jan 06, 2012 Modification as noted above.
			if(obs_ct(k,i)>0)
			{
				dvar_vector fa=ft(k,i)*mfexp(log_va);
				Chat(k,i)=elem_prod(elem_prod(elem_div(fa,Z(i)),1.-S(i)),N(i));
				switch(catch_type(k))
				{
					case 1:	//catch in weight
						ct(k,i) = Chat(k,i)*wt_obs(i);
					break;
					case 2:	//catch in numbers
						ct(k,i) = sum(Chat(k,i));
					break;
					case 3:	//catch in roe that does not contribute to SSB
						dvariable ssb = elem_prod(N(i),exp(-Z(i)*cntrl(13)))*fec(i);
						ct(k,i) = ( 1.-mfexp(-ft(k,i)) )*ssb;
					break;
				}
				eta(k,i) = log(obs_ct(k,i))-log(ct(k,i));
			}
			else
			{/*If there is no commercial fishery the set Chat equal to 
			   the expected proportions at age.*/
				dvar_vector fa = mfexp(log_va);
				Chat(k,i)=elem_prod(elem_prod(elem_div(fa,Z(i)),1.-S(i)),N(i));
				
				/*
				dvar_matrix alk(sage,nage,1,nbins) = ALK(la, cv*la, xbin)
				
				Lhat(i) = Chat *ALK; 
				*/
				
			}
			
			
			/*
			Changed this on Jan 16, 2011 as it was preventing
			convergence for simuation with zero error due to the tiny
			constant added to F.
			
			need to add a tiny constant to deal with catch-age data 
			for non-extractive survey age comps
			dvar_vector fa=ft(k,i)*mfexp(log_va)+1.e-30;
			
			//Catch-at-age by commercial gear
			if(fsh_flag(k))
				Chat(k,i)=elem_prod(elem_prod(elem_div(fa,Z(i)),1.-S(i)),N(i));
			
			//Catch weight by gear
			//ct(k,i)=Chat(k,i)*wa;  
			ct(k,i)=Chat(k,i)*wt_obs(i);  //SM Dec 6, 2010*/
		}
	}
	if(verbose)cout<<"**** Ok after calcFisheryObservations ****"<<endl;

  }	
	
FUNCTION calcSurveyObservations
  {
	/*This code needs to be modified to accomodate
	multiple surveys or block changes in survey q.
	
	Oct 31, 2010, added retrospective counter.
	
	Nov 22, 2010, adding multiple surveys. Still need to check with retrospective option
	
	Nov 30, 2010, adjust the suvery biomass by the fraction of Z that has occurred 
	when the survey was conducted. For herring spawning biomass this would be after the 
	fishery has taken place.
	
	Dec 6, 2010, modified predicted survey biomass to accomodate empirical weight-at-age 
	data (wt_obs).
	
	May 11, 2011.  CHANGED Vivian Haist pointed out an error in survey biomass comparison.
	The spawning biomass was not properly calculated in this routine. I.e. its different 
	than the spawning biomass in the stock-recruitment routine. (Based on fecundity which
	changes with time when given empirical weight-at-age data.)
	
	Jan 6, 2012.  CHANGED corrected spawn survey observations to include a roe 
	fishery that would remove potential spawn that would not be surveyed.
	
	*/
	/*
		CHANGED add capability to accomodate priors for survey q's.
		DONE
	*/
	
	int i,j,ii,k,kk;
	
	
	//survey abudance index residuals
	epsilon.initialize();
	pit.initialize();
	
	for(i=1;i<=nit;i++)
	{	
		int nx=0;		//counter for retrospective analysis
		dvar_matrix V(1,nit_nobs(i),sage,nage);  //vulnerable units for survey comparison
		V.initialize();
		//dvar_matrix VB(1,nit_nobs(i),sage,nage); //vulnerable biomass
		for(j=1;j<=nit_nobs(i);j++)
		{
			ii=iyr(i,j);
			k=igr(i,j);
			if(ii>nyr) break;	//trap for retrospective analysis.
			dvar_vector log_va=log_sel(k)(ii);
			//Adjust survey biomass by the fraction of the mortality that 
			//occurred during the time of the survey.
			dvar_vector Np = elem_prod(N(ii),exp( -Z(ii)*it_timing(i,j) ));
			
			//get fishing mortality rate on spawn.
			dvariable ftmp = 0;
			for(kk=1;kk<=ngear;kk++)
				if(catch_type(kk)==3)
					ftmp += ft(kk,ii);
					
			switch(survey_type(i))
			{
				case 1:
					V(j)=elem_prod(Np,mfexp(log_va));
				break;
				case 2:
					V(j)=elem_prod(elem_prod(Np,mfexp(log_va)),wt_obs(ii));
				break;
				case 3:
					//SJDM Jan 6, 2012 Modified for roe fishery
					V(j)=elem_prod(Np,fec(ii))*exp(-ftmp);
				break;
			}
			
			//VB(j)=elem_prod(V(j),wt_obs(ii));		//SM Dec 6, 2010
			
			//If the survey is a spawn index, then need to account
			//for changes in fecundity.
			
			if(iyr(i,j)<=nyr) nx++;
		}
		dvar_vector t1 = rowsum(V);//V*wa;
		//cout<<"V\n"<<V<<endl;
		
		//See Ludwig & Walters 1994
		//Note this is incorrect if each survey has different weights.
		dvar_vector zt=log(it(i).sub(1,nx))-log(t1(1,nx));
		//cout<<"zt\n"<<t1(1,nx)<<endl;
		epsilon(i).sub(1,nx) = zt-mean(zt);
		q(i) = exp(mean(zt));
		pit(i).sub(1,nx)=t1(1,nx)*q(i);	//predicted index
		
		
		//TODO, this might not be working correctly, simulation test it.
		// The following allows for a penalized random walk in survey q.
		if(q_prior(i)==2)
		{
			//random walk in q
			epsilon(i)=0;
			dvar_vector fd_zt=first_difference(zt);
			epsilon(i).sub(1,nx-1) = fd_zt-mean(fd_zt);
			//dvar_vector qt(1,nx);
			qt(i,1) = exp(zt(1));
			for(j=2;j<=nx;j++)
				qt(i,j) = qt(i,j-1)*exp(fd_zt(j-1));
				
			pit(i).sub(1,nx)=elem_prod(t1(1,nx),qt(i)(1,nx));
			//cout<<sum(epsilon(i))<<endl;
			//exit(1);
		}
	
	}
	
	
	
	
	if(verbose)cout<<"**** Ok after calcSurveyObservations ****"<<endl;
	
  }
	
FUNCTION calcStockRecruitment
  {
	/*
	The following code is used to derive unfished
	spawning stock biomass bo and the stock-
	recruitment parameters for the:
	Beverton-Holt Model
		-Rt=k*Ro*St/(Bo+(k-1)*St)*exp(delta-0.5*tau*tau)
	Ricker Model
		-Rt=so*St*exp(-beta*St)*exp(delta-0.5*tau*tau)
		
	Dec 6, 2010.  Modified code to allow empirical weight-at-age data
				This required a fecundity-at-age matrix.  Need to 
				project spawning biomass into the future.
	CHANGED bo should be a function of the average natural mortality
	TODO update phib calculation if age-specific M's are used.
	
	May 6, 2010.  Changed phib calculation based on average M 
	in cases where M varies over time. Otherwise bo is biased based
	on the initial value of M.
	
	CHANGED Need to adjust spawning biomass to post fishery numbers.
	CHANGED Need to adjust spawners per recruit (phib) to average fecundity.
	
	Jan 6, 2012.  Need to adjust stock-recruitment curve for reductions 
	in fecundity associated with removal of roe from a spawn on kelp fishery.
	
	Aug 9, 2012. Revised routine so that the slope of the stock recruitment
	relationship is based on all of the average fecundity over the whole series.
	Bo is no longer calculated in this routine. This is in response to recent 
	issue with the herring assessment, where reference points (Bo included) are 
	based on recent trends in mean weight-at-age/fecundity-at-age.
	
	Psuedocode:
		-1) Get average natural mortality rate at age.
		-2) Calculate survivorship to time of spawning.
		-3) Calculate unfished spawning biomass per recruit.
		-4) Compute spawning biomass vector & substract roe fishery
		-5) Project spawning biomass to nyr+1 under natural mortality.
		-6) Calculate stock recruitment parameters (so, beta);
		-7) Calculate predicted recruitment
		-8) Compute residuals from estimated recruitments.
	
	*/ 
	int i,j,k;
	
	dvariable   phib,so,beta;
	dvariable   tau = (1.-rho)/varphi;
	dvar_vector     ma(sage,nage);
	dvar_vector tmp_rt(syr+sage,nyr);
	dvar_vector     lx(sage,nage); lx(sage) = 1.0;
	dvar_vector     lw(sage,nage); lw(sage) = 1.0;

	
	// -steps (1),(2),(3)
	dvar_matrix t_M_tot = trans(M_tot);
	for(j=sage; j<=nage;j++)
	{
		ma(j) = mean(t_M_tot(j));
		if(j>sage)
		{
			lx(j) = lx(j-1)*mfexp(-ma(j-1));
		} 
		lw(j) = lx(j) * mfexp(-ma(j)*cntrl(13));
	}
	lx(nage) /= 1.0 - mfexp(-ma(nage));
	lw(nage) /= 1.0 - mfexp(-ma(nage));
	phib      = lw * fa_bar;
	
	
	// step (4)
	for(i=syr;i<=nyr;i++)
	{
		sbt(i) = elem_prod(N(i),exp(-Z(i)*cntrl(13)))*fec(i);
		//Adjustment to spawning biomass for roe fisheries
		for(k=1;k<=ngear;k++)
		{
			if(catch_type(k)==3)
			{
				sbt(i) *= mfexp(-ft(k,i));
			}
		}
			
	}
	sbt(nyr+1) = elem_prod(N(nyr+1),exp(-M_tot(nyr))) * fec(nyr+1);	
	dvar_vector tmp_st=sbt(syr,nyr-sage).shift(syr+sage);
	
	sbo = ro*phib;
	// steps (6),(7)
	so = kappa/phib;			//max recruits per spawner
	switch(int(cntrl(2)))
	{
		case 1:
			//Beverton-Holt model
			beta   = (kappa-1.)/sbo;
			tmp_rt = elem_div(so*tmp_st,1.+beta*tmp_st);
			break;
		case 2:
			//Ricker model
			beta   = log(kappa)/sbo;
			tmp_rt = elem_prod(so*tmp_st,exp(-beta*tmp_st));
		break;
	}
	
	// step (8) residuals in stock-recruitment curve
	rt    = mfexp(log_rt(syr+sage,nyr));
	delta = log(rt)-log(tmp_rt)+0.5*tau*tau;
	
	if(verbose)cout<<"**** Ok after calcStockRecruitment ****"<<endl;
	
  }
	
FUNCTION calc_objective_function
  {
	//Dec 20, 2010.  SJDM added prior to survey qs.
	/*q_prior is an ivector with current options of 0 & 1.
	0 is a uniform density (ignored) and 1 is a normal
	prior density applied to log(q).*/

	/*
	There are several components to the objective function
	Likelihoods:
		-1) likelihood of the catch data
		-2) likelihood of the survey abundance index
		-3) likelihood of age composition data 
		-4) likelihood for stock-recruitment relationship
		-5) penalized likelihood for fishery selectivities
		-6) penalized likelihood for fishery selectivities
	*/
	int i,j,k;
	double o=1.e-10;

	dvar_vector lvec(1,7); 
	lvec.initialize();
	nlvec.initialize();
	
	//1) likelihood of the catch data (retro)
	double sig_c =cntrl(3);
	if(last_phase())sig_c=cntrl(4);
	if(active(log_ft_pars))
	for(k=1;k<=ngear;k++){
		for(i=syr;i<=nyr;i++)
		{
			if(obs_ct(k,i)!=0)
				nlvec(1,k)+=dnorm(eta(k,i),0.0,sig_c);
				//nlvec(1,k)+=dnorm(log(ct(k,i)),log(obs_ct(k,i)),sig_c);
		}
		//if(active(log_ft_pars))
		//	nlvec(1,k)=dnorm(log(obs_ct(k).sub(syr,nyr)+o)-log(ct(k).sub(syr,nyr)+o),sig_c);
	}
	
	
	//2) likelihood of the survey abundance index (retro)
	for(k=1;k<=nit;k++)
	{
		dvar_vector sig = (rho/varphi)/it_wt(k);
		nlvec(2,k)=dnorm(epsilon(k),sig);
	}
	
	
	
	//3) likelihood for age-composition data
	for(k=1;k<=na_gears;k++)
	{	
		if(na_nobs(k)>0){
			int naa=0;
			int iyr;
			//retrospective counter
			for(i=1;i<=na_nobs(k);i++)
			{
				iyr=A(k,i,a_sage(k)-2);	//index for year
				if(iyr<=nyr)naa++;
			}
			
			dmatrix O=trans(trans(A(k)).sub(a_sage(k),a_nage(k))).sub(1,naa);
			dvar_matrix P=trans(trans(Ahat(k)).sub(a_sage(k),a_nage(k))).sub(1,naa);
			//dvar_matrix nu=trans(trans(Ahat(k)).sub(a_sage(k),a_nage(k))).sub(1,naa); //residuals
			dvar_matrix nu(O.rowmin(),O.rowmax(),O.colmin(),O.colmax()); //residuals
			nu.initialize();
			
			//CHANGED add a switch statement here to choose form of the likelihood
			switch(int(cntrl(14)))
			{
				case 1:
					nlvec(3,k) = dmvlogistic(O,P,nu,age_tau2(k),cntrl(6));
				break;
				case 2:
					nlvec(3,k) = dmultinom(O,P,nu,age_tau2(k),cntrl(6));
				break;
			}
			
			for(i=1;i<=naa/*na_nobs(k)*/;i++)
			{
				iyr=A(k,i,a_sage(k)-2);	//index for year
				A_nu(k)(i)(a_sage(k),a_nage(k))=nu(i);
			}
		}
	}
	
	
	
	//4) likelihood for stock-recruitment relationship
	dvariable tau = (1.-rho)/varphi;
	if(active(theta(1)))
		nlvec(4,1)=dnorm(delta,tau);
	
	
	
	//5-6) likelihood for selectivity paramters
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

				for(i=syr;i<=nyr;i++)
				{
					//curvature in selectivity parameters
					dvar_vector df2 = first_difference(first_difference(log_sel(k)(i)));
					nlvec(5,k)     += sel_2nd_diff_wt(k)/(nage-sage+1)*df2*df2;

					//penalty for dome-shapeness
					for(j=sage;j<=nage-1;j++)
						if(log_sel(k,i,j)>log_sel(k,i,j+1))
							nlvec(6,k)+=sel_dome_wt(k)
										*square(log_sel(k,i,j)-log_sel(k,i,j+1));
				}
			}
			
			/*Oct 31, 2012 Halloween! Added 2nd difference penalty on time for isel_type==(4)*/
			if( isel_type(k)==4 )
			{
				dvar_matrix trans_log_sel = trans( log_sel(k) );
				for(j=sage;j<=nage;j++)
				{
					dvar_vector df2 = first_difference(first_difference(trans_log_sel(j)));
					nlvec(7,k)     += sel_2nd_diff_wt_time(k)/(nage-sage+1)*norm2(df2);
				}
			}
			
		}
	}
	
	
	// CONSTRAINT FOR SELECTIVITY DEV VECTORS
	// Ensure vector of sel_par sums to 0. (i.e., a dev_vector)
	// TODO for isel_type==2 ensure mean 0 as well (ie. a dev_vector)
	
	for(k=1;k<=ngear;k++)
	{
		if( active(sel_par(k)) &&
			isel_type(k)!=1 &&
			isel_type(k)!=7 &&
			isel_type(k)!=8 &&
			isel_type(k)!=11 )
		{
			dvariable s=0;
			if(isel_type(k)==5)  //bicubic spline version ensure column mean = 0
			{
				dvar_matrix tmp = trans(sel_par(k));
				for(j=1;j<=tmp.rowmax();j++)
				{
					s=mean(tmp(j));
					lvec(1)+=10000.0*s*s;
				}
			}
			if( isel_type(k)==4 ||
			 	isel_type(k)==3 || 
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
	
	
	/*
	PRIORS for estimated model parameters from the control file.
	*/
	dvar_vector priors(1,npar);
	dvariable ptmp; priors.initialize();
	for(i=1;i<=npar;i++)
	{
		if(active(theta(i)))
		{
			switch(theta_prior(i))
			{
			case 1:		//normal
				ptmp=dnorm(theta(i),theta_control(i,6),theta_control(i,7));
				break;
				
			case 2:		//lognormal CHANGED RF found an error in dlnorm prior. rev 116
				ptmp=dlnorm(theta(i),theta_control(i,6),theta_control(i,7));
				break;
				
			case 3:		//beta distribution (0-1 scale)
				double lb,ub;
				lb=theta_lb(i);
				ub=theta_ub(i);
				ptmp=dbeta((theta(i)-lb)/(ub-lb),theta_control(i,6),theta_control(i,7));
				break;
				
			case 4:		//gamma distribution
				ptmp=dgamma(theta(i),theta_control(i,6),theta_control(i,7));
				break;
				
			default:	//uniform density
				ptmp=1./(theta_control(i,3)-theta_control(i,2));
				break;
			}
			priors(i)=ptmp;	
		}
	}
	
	//Priors for suvey q based on control file.
	dvar_vector qvec(1,nits);

	qvec.initialize();
	for(i=1;i<=nits;i++)
	{
		if(q_prior(i)==1) 
		{
			qvec(i)=dnorm(log(q(i)),q_mu(i),q_sd(i));
		}
	}
	
	
	//** Legacy **  By accident took Rick Methot's bag from Nantes.
	//301 787 0241  Richard Methot cell phone.
	//ibis charles du gaulle at
	//01 49 19 19 20
	
	/*
	The following are penalties that are invoked in early
	phases to ensure reasonable solutions are obtained,
	and to reduce the sensitivity of initial starting
	conditions.  Penalties include:
		-1) keep average fishing mortality rate near 
			0.2 and in the last phase relax this constraint.
		-3) normal prior for log rec devs with std=50.
	*/
	
	dvar_vector pvec(1,7);
	pvec.initialize();
	
	//Penalties to regularize the solution for fishing mortality rates
	dvariable log_fbar = mean(log_ft_pars);
	if(last_phase())
	{
		pvec(1) = dnorm(log_fbar,log(cntrl(7)),cntrl(9));
		//Penalty for log_rec_devs (large variance here)
		pvec(4) = dnorm(log_rec_devs,2.0);
		pvec(5) = dnorm(init_log_rec_devs,2.0);
	}
	else
	{
		pvec(1) = dnorm(log_fbar,log(cntrl(7)),cntrl(8));
		//Penalty for log_rec_devs (CV ~ 0.0707) in early phases
		pvec(4)=100.*norm2(log_rec_devs);
		pvec(5)=100.*norm2(init_log_rec_devs);
	}
	
	//Priors for deviations in natural mortality rates
	//if(active(log_m_devs))
	if(active(log_m_nodes))
	{
		double std_mdev = cntrl(11);
		dvar_vector fd_mdevs=first_difference(log_m_devs);
		pvec(2) = dnorm(fd_mdevs,std_mdev);
		pvec(2) += 0.5*norm2(log_m_nodes);
	}
	
	
	if(verbose)
	{
		COUT(nlvec);
		cout<<"lvec\t"<<lvec<<endl;
		cout<<"priors\t"<<priors<<endl;
		cout<<"penalties\t"<<pvec<<endl;
	}
	f=sum(nlvec)+sum(lvec)+sum(priors)+sum(pvec)+sum(qvec);
	//cout<<f<<endl;
	nf++;
	if(verbose)cout<<"**** Ok after calc_objective_function ****"<<endl;
	
  }

FUNCTION void equilibrium(const double& fe, const dvector& ak, const double& ro, const double& kap, const double& m, const dvector& age, const dvector& wa, const dvector& fa, const dmatrix& va,double& re,double& ye,double& be,double& ve,double& dye_df,double& d2ye_df2)//,double& phiq,double& dphiq_df, double& dre_df)
  {
	/*
	Equilibrium age-structured model used to determin Fmsy and MSY based reference points.
	Author: Steven Martell
	
	Comments: 
	This code uses a numerical approach to determine a vector of fe_multipliers
	to ensure that the allocation is met for each gear type.
	
	args:
	fe	-steady state fishing mortality
	ak	-allocation of total ye to gear k.
	ro	-unfished sage recruits
	kap	-recruitment compensation ration
	m	-instantaneous natural mortality rate
	age	-vector of ages
	wa	-mean weight at age
	fa	-mean fecundity at age
	va	-mean vulnerablity at age for fe gear.

	
	Modified args:
	re	-steady state recruitment
	ye	-steady state yield
	be	-steady state spawning biomass
	phiq		-per recruit yield
	dre_df		-partial of recruitment wrt fe
	dphiq_df	-partial of per recruit yield wrt fe
	
	LUCIE'S RULE: the derivative of a sum is the sum of its derivatives.
	Lucie says: there is some nasty calculus in here daddy, are you sure
	you've got it right?
	
	I've got it pretty close. 
	
	DEPRECATE THIS FUNCTION.  NOW DONE IN THE MSY CLASS
	
	*/
	int i,j,k;
	int nage    = max(age);
	int sage    = min(age);
	double  dre_df;
	double  phif;
	dvector lx(sage,nage);
	dvector lz(sage,nage);
	dvector lambda(1,ngear);        //F-multiplier
	dvector phix(1,ngear);
	dvector phiq(1,ngear);
	dvector dphiq_df(1,ngear);
	dvector dyek_df(1,ngear);
	dvector d2yek_df2(1,ngear);
	dvector yek(1,ngear);
	dmatrix qa(1,ngear,sage,nage);
	dmatrix xa(1,ngear,sage,nage);  //vulnerable numbers per recruit
	
	lx          = pow(exp(-m),age-double(sage));
	lx(nage)   /=(1.-exp(-m));
	double phie = lx*fa;		// eggs per recruit
	double so   = kap/phie;
	double beta = (kap-1.)/(ro*phie);
	lambda      = ak/mean(ak);	// multiplier for fe for each gear
	
	
	/* Must iteratively solve for f-multilier */
	for(int iter=1;iter<=30;iter++)
	{
		/* Survivorship under fished conditions */
		lz(sage)    = 1.0;
		lambda     /= mean(lambda);
		dvector fk  = fe*lambda;
		dvector ra  = lambda*va;
		dvector za  = m + fe*ra;
		dvector sa  = mfexp(-za);
		dvector oa  = 1.0 - sa;
		
		
		for(k=1;k<=ngear;k++)
		{
			qa(k) = elem_prod(elem_div(lambda(k)*va(k),za),oa);
			xa(k) = elem_prod(elem_div(va(k),za),oa);
		}
		
		double dlz_df = 0, dphif_df = 0;
		dphiq_df.initialize();
		dre_df   = 0;
		for(j=sage;j<=nage;j++)
		{
			if(j>sage) lz(j)  = lz(j-1) * sa(j-1);
			if(j>sage) dlz_df = dlz_df  * sa(j-1) - lz(j-1)*ra(j-1)*sa(j-1);
			
			if(j==nage)
			{
				lz(j)  = lz(j) / oa(j);
				
				double t4 = (-ra(j-1)+ra(j))*sa(j)+ra(j-1);
				dlz_df = dlz_df/oa(j) - (lz(j-1)*sa(j-1)*t4) / square(oa(j));
			}
			
			dphif_df   = dphif_df+fa(j)*dlz_df;
			for(k=1;k<=ngear;k++)
			{
				double t1   = lambda(k) * wa(j) *va(k,j) * ra(j) * lz(j);
				double t3   = -1. + (1.+za(j)) * sa(j);
				double t9   = square(za(j));
				dphiq_df(k)+= wa(j)*qa(k,j)*dlz_df + t1 * t3 / t9; 
			}
		} 
		
		phif   = elem_prod(lz,exp(-za*cntrl(13)))*fa;
		re     = ro*(kap-phie/phif)/(kap-1.);
		dre_df = ro/(kap-1.0)*phie/square(phif)*dphif_df;
		
		/* Equilibrium yield */
		for(k=1;k<=ngear;k++)
		{
			phix(k)      = sum(elem_prod(elem_prod(lz,wa),xa(k)));
			phiq(k)      = sum(elem_prod(elem_prod(lz,wa),qa(k)));
			yek(k)       = fe*re*phiq(k);
			dyek_df(k)   = re*phiq(k) + fe*phiq(k)*dre_df + fe*re*dphiq_df(k);
			d2yek_df2(k) = phiq(k)*dre_df + re*dphiq_df(k);
		}
		
		/* Iterative soln for lambda */
		dvector pk = yek/sum(yek);
		dvector t1 = elem_div(ak,pk+1.e-30);
		lambda     = elem_prod(lambda,t1);
		if(abs(sum(ak-pk))<1.e-6) break;
	} // end of iter
	ve       = re*sum(elem_prod(ak,phix));
	be       = re*phif;
	ye       = sum(yek);
	dye_df   = sum(dyek_df);
	d2ye_df2 = sum(d2yek_df2);

	// cout<<"EQUILIBRIUM CODE "<<setprecision(4)<<setw(2)<<fe<<setw(3)<<" "
	// <<ye<<setw(5)<<" "<<dye_df<<"  "<<dyek_df(1,3)<<endl;
  }



	
FUNCTION void equilibrium(const double& fe,const double& ro, const double& kap, const double& m, const dvector& age, const dvector& wa, const dvector& fa, const dvector& va,double& re,double& ye,double& be,double& phiq,double& dphiq_df, double& dre_df)
  {
	/*
	This is the equilibrium age-structured model that is 
	conditioned on fe (the steady state fishing mortality rate).
	
	In the case of multiple fisheries, fe is to be considered as the
	total fishing mortality rate and each fleet is given a specified
	allocation based on its selectivity curve.  The allocation to 
	each fleet must be specified a priori.
	
	args:
	fe	-steady state fishing mortality
	ro	-unfished sage recruits
	kap	-recruitment compensation ration
	m	-instantaneous natural mortality rate
	age	-vector of ages
	wa	-mean weight at age
	fa	-mean fecundity at age
	va	-mean vulnerablity at age for fe gear.
	ak	-allocation of total ye to gear k.
	
	Modified args:
	re	-steady state recruitment
	ye	-steady state yield
	be	-steady state spawning biomass
	phiq		-per recruit yield
	dre_df		-partial of recruitment wrt fe
	dphiq_df	-partial of per recruit yield wrt fe
	
	FIXME add Ricker model to reference points calculations.
	FIXME partial derivatives for dphif_df need to be fixed when cntrl(13)>0.
	*/
	int i;
	
	int nage=max(age);
	int sage=min(age);
	dvector lx=pow(exp(-m),age-double(sage));
	lx(nage)/=(1.-exp(-m));
	dvector lz=lx;
	dvector za=m+fe*va;
	dvector sa=1.-exp(-za);
	dvector qa=elem_prod(elem_div(va,za),sa);
	
	double phie = lx*fa;		//eggs per recruit
	double so = kap/phie;
	double beta = (kap-1.)/(ro*phie);
	
	
	double dlz_df = 0, dphif_df = 0;
	dphiq_df=0; dre_df=0;
	for(i=sage; i<=nage; i++)
	{
		if(i>sage) lz[i]=lz[i-1]*exp(-za[i-1]);
		if(i>sage) dlz_df=dlz_df*exp(-za[i-1]) - lz[i-1]*va[i-1]*exp(-za[i-1]);
		if(i==nage){ //6/11/2007 added plus group.
					lz[i]/=(1.-mfexp(-za[i]));
					
					dlz_df=dlz_df/(1.-mfexp(-za[i]))
							-lz[i-1]*mfexp(-za[i-1])*va[i]*mfexp(-za[i])
					/((1.-mfexp(-za[i]))*(1.-mfexp(-za[i])));
				}
		dphif_df=dphif_df+fa(i)*dlz_df;
		dphiq_df=dphiq_df+wa(i)*qa(i)*dlz_df+(lz(i)*wa(i)*va(i)*va(i))/za(i)*(exp(-za[i])-sa(i)/za(i));
	}
	//CHANGED need to account for fraction of mortality that occurs
	//before the spawning season in the recruitment calculation.
	//cout<<"lz\t"<<elem_prod(lz,exp(-za*cntrl(13)))<<endl;
	//exit(1);
	//double phif = lz*fa;
	double phif = elem_prod(lz,exp(-za*cntrl(13)))*fa;
	phiq=sum(elem_prod(elem_prod(lz,wa),qa));
	re=ro*(kap-phie/phif)/(kap-1.);
	//cout<<fe<<" spr ="<<phif/phie<<endl;
	if(re<=0) re=0;
	dre_df=(ro/(kap-1.))*phie/square(phif)*dphif_df;
	ye=fe*re*phiq;
	be=re*phif;	//spawning biomass
	
	//cout<<"Equilibrium\n"<<ro<<"\n"<<re<<"\n"<<ye<<endl;
	
  }
	
FUNCTION void calc_reference_points()
  {
	/**
	\file iscam.tpl
	\author Steven Martell
	Uses Newton_Raphson method to determine Fmsy and MSY 
	based reference points.  
	
	Code check: appears to find the correct value of MSY
	in terms of maximizing ye.  Check to ensure rec-devs
	need a bias correction term to get this right.
	
	Modification for multiple fleets:
		Need to pass a weighted average vector of selectivities
		to the equilibrium routine, where the weights for each
		selectivity is based on the allocation to each fleet.
		
		Perhaps as a default, assign an equal allocation to each
		fleet.  Eventually,user must specify allocation in 
		control file.  DONE
		
		Use selectivity in the terminal year to calculate reference
		points.
	
	June 8, 2012.  SJDM.  Made the following changes to this routine.
		1) changed reference points calculations to use the average
		   weight-at-age and fecundity-at-age.
		2) change equilibrium calculations to use the catch allocation
		   for multiple gear types. Not the average vulnerablity... this was wrong.
	
	July 29, 2012.  SJDM Issue1.  New routine for calculating reference points
	for multiple fleets. In this case, finds a vector of Fmsy's that simultaneously 
	maximizes the total catch for each of the fleets respectively.  See
	iSCAMequil_soln.R for an example.
	
	August 1, 2012.  SJDM, In response to Issue1. A major overhaul of this routine.
	Now using the new Msy class to calculate reference points. This greatly simplifies
	the code in this routine and makes other routines (equilibrium) redundant.  Also
	the new Msy class does a much better job in the case of multiple fleets.
	
	The algorithm is as follows:
		(2) Construct a matrix of selectivities for the directed fleets.
		(3) Come up with a reasonable guess for fmsy for each gear.
		(4) Instantiate an Msy class object and get_fmsy.
		(5) Use Msy object to get reference points.
		
		
	Aug 11, 2012.
	For the Pacific herring branch omit the get_fmsy calculations and use only the 
	Bo calculuation for the reference points.  As there are no MSY based reference
	points required for the descision table. 
	*/
	int i,j,k;
	
	/* (1) Determine which fleets are directed fishing fleets. */
	/* This is done in the data section with nfeet and ifleet. */
	
	
	/* (2) Matrix of selectivities for directed fleets */
	dmatrix d_V(1,nfleet,sage,nage);
	dvector d_ak(1,nfleet);
	for(k = 1; k<= nfleet; k++)
	{
		j    = ifleet(k);
		d_V(k) = value(exp(log_sel(j)(nyr)));
		d_ak(k)= allocation(j);
	}
	d_ak /= sum(d_ak);
	
	/* 
	(3) Come up with a reasonable estimate of Fmsy 
	In practice seems to  work best if start with 
	extreme low values.
	*/
	dvector ftry(1,nfleet);
	ftry = 0.6*value(m_bar);    // initial guess for Fmsy
	fmsy = ftry;
	fall = ftry;
	
	
	/* (4) Instantiate an Msy class object and get_fmsy */
	double  d_ro  = value(ro);
	double  d_h   = value(theta(2));
	double  d_m   = value(m_bar);
	double  d_rho = cntrl(13);
	dvector d_wa  = (avg_wt);
	dvector d_fa  = (avg_fec);
	Msy cMSY(d_ro,d_h,d_m,d_rho,d_wa,d_fa,d_V);
	
	bo   = cMSY.getBo();
	cMSY.get_fmsy(fall,d_ak);
	cMSY.get_fmsy(fmsy);
	if(cMSY.getFail())
	cout<<"FAILED TO CONVERGE"<<endl;
	
	
	msy  = cMSY.getYe();
	bmsy = cMSY.getBe();
	Umsy = sum(cMSY.getYe())/cMSY.getBi();
	
	
	
	cout<<"|------------------------------------------|" <<endl;
	cout<<"| Bo   = "<<setw(10)<<bo                      <<endl;
	cout<<"| Bmsy = "<<setw(10)<<bmsy                    <<endl;
	cout<<"| Fmsy ="<<setw(10)<<fmsy                     <<endl;
	cout<<"| Fall ="<<setw(10)<<fall                     <<endl;
	cout<<"| MSY  ="<<setw(10)<<msy                      <<endl;
	cout<<"| dYe  = "<<setw(10)<<sum(cMSY.getdYe())      <<endl;
	cout<<"|------------------------------------------|" <<endl;
	
	
	
	//fall = ftry;
	//
	////fmsy = fall;
	//cMSY.get_fmsy(fmsy);
	//bmsy = cMSY.getBmsy();
	//msy  = cMSY.getMsy();
	//bo   = cMSY.getBo();  //Spawning biomass just prior to spawning.
	//
	////if(nf==1) ftry = fmsy;
	//
	//cout<<"------------------------"<<endl;
	//cout<<"Ftry      \t"<<ftry<<endl;
	//cout<<"Fmsy      \t"<<fmsy<<endl;
	//cout<<"MSY       \t"<<msy<<endl;
	//cout<<"dYe       \t"<<cMSY.getdYe()<<endl;
	//cout<<"Bo        \t"<<bo<<endl;
	//cout<<"Bmsy      \t"<<bmsy<<endl;
	//cout<<"Bi        \t"<<cMSY.getBi()<<endl;
	//cout<<"SPR at MSY\t"<<cMSY.getSprMsy()<<endl;
	//cout<<"phiB      \t"<<cMSY.getPhie()<<endl;
	//cout<<"------------------------"<<endl;
	//
	///* (4) Now do it with allocation */
	////cout<<"\nAllocation"<<allocation(ifleet)<<endl;
	//fall = ftry;
	//cMSY.get_fmsy(fall,d_ak);
	////bmsy = cMSY.getBmsy();
	////msy  = cMSY.getMsy();
	////bo   = cMSY.getBo();  //Spawning biomass just prior to spawning.
	//
	///* 
	//I've defined Umsy as the sum of catches divided 
	//by spawning biomass at the start of the year.
	//*/
	//Umsy = sum(cMSY.getYe())/cMSY.getBi();
	//
	//
	//cout<<"------------------------"<<endl;
	//cout<<"Fall      \t"<<fall<<endl;
	//cout<<"Yield     \t"<<cMSY.getYe()<<endl;
	//cout<<"Be        \t"<<cMSY.getBe()<<endl;
	//cout<<"Spr       \t"<<cMSY.getSpr()<<endl;
	//cout<<"Umsy      \t"<<Umsy<<endl;
	//cout<<"------------------------"<<endl;
	
	//The following code should be deprecated, along with the two equilibrium functions
	//as this reference point material is now hanlded by the Msy class.
	
	//dmatrix va(1,ngear,sage,nage);
	//dvector va_bar(sage,nage);
	//va_bar.initialize();
	//
	//
	///*CHANGED Allow for user to specify allocation among gear types.*/
	///*FIXME:  this allocation should be on the catch on the vulnerabilities*/
	//for(j=1;j<=ngear;j++)
	//{
	//	va_bar+=allocation(j)*value(exp(log_sel(j)(nyr)));
	//	va(j) = value(exp(log_sel(j)(nyr)));
	//}
	//dmatrix V(1,ngear-2,sage,nage);
	//dvector fk(1,ngear-2);
	//fk = 0.6*value(m_bar);
	//for(j=1;j<=ngear-2;j++)
	//{
	//	V(j) = value(exp(log_sel(j)(nyr)));
	//}
	//
	//double h = value(theta(2));
	//cout<<"Declaring class"<<endl;
	//Msy cMSY(value(ro),h,value(m_bar),avg_wt,avg_fec,V);
	//cout<<"About to call get_fmsy"<<endl;
	//fk = cMSY.get_fmsy(fk);
	
	/*CHANGED: SJDM June 8, 2012 fixed average weight-at-age for reference points
	           and average fecundity-at-age.
	*/
	
	//#if defined(USE_NEW_EQUILIBRIUM)
	//	/* Newton-Raphson method to determine MSY-based reference points. */
	//	for(i=1;i<=15;i++)
	//	{
	//		equilibrium(fe,allocation,value(ro),value(kappa),value(m_bar),age,avg_wt,
	//				avg_fec,va,re,ye,be,ve,dye_df,d2ye_df2);
	//	
	//		fe = fe - dye_df/d2ye_df2;
	//		if(square(dye_df)<1e-12)break;
	//	}
	//	fmsy=fe;
	//	equilibrium(fe,allocation,value(ro),value(kappa),value(m_bar),age,avg_wt,
	//			avg_fec,va,re,ye,be,ve,dye_df,d2ye_df2);
	//#endif
	//
	//#if !defined(USE_NEW_EQUILIBRIUM)
	//	for(i=1;i<=20;i++)
	//	{
	//		//equilibrium(fe,value(ro),value(kappa),value(m),age,wa,fa,value(exp(log_sel(1)(nyr))),re,ye,be,phiq,dphiq_df,dre_df);
	//		//equilibrium(fe,value(ro),value(kappa),value(m_bar),age,wt_obs(nyr),
	//		//			fec(nyr),va_bar,re,ye,be,phiq,dphiq_df,dre_df);
	//		equilibrium(fe,value(ro),value(kappa),value(m_bar),age,avg_wt,
	//					avg_fec,va_bar,re,ye,be,phiq,dphiq_df,dre_df);
	//	
	//		dye_df = re*phiq+fe*phiq*dre_df+fe*re*dphiq_df;
	//		ddye_df = phiq*dre_df + re*dphiq_df;
	//		fe = fe - dye_df/ddye_df;
	//		if(verbose) cout<<"fe\t"<<fe<<"\t"<<dye_df<<"\t"<<ye<<endl;
	//		if(sfabs(dye_df)<1.e-5)break;
	//	}
	//	fmsy=fe;
	//	equilibrium(fmsy,value(ro),value(kappa),value(m_bar),age,avg_wt,
	//				avg_fec,va_bar,re,ye,be,phiq,dphiq_df,dre_df);
	//#endif
	
	//msy=ye;
	//bmsy=be;
	//Umsy=msy/Vmsy;
	
	/*TODO print this to the REPORT file for plotting.*/
	/*SM Loop over discrete value of fe and ensure above code is 
	finding the correct value of msy.*/
	
	//if(!mceval_phase())
	//{
	//	ofstream report_file("iscam.eql");
	//
	//	if(report_file.is_open())
	//	{
	//		report_file<<"index\t fe \t ye \t be \t ve \t re \t spr\n";
	//	
	//		fe = 0; i=0;
	//		while(i < 1500)
	//		{
	//			#if !defined(USE_NEW_EQUILIBRIUM)
	//			equilibrium(fe,value(ro),value(kappa),value(m_bar),age,wt_obs(nyr),
	//						fec(nyr),va_bar,re,ye,be,phiq,dphiq_df,dre_df);
	//			#endif
	//		
	//			#if defined(USE_NEW_EQUILIBRIUM)
	//			equilibrium(fe,allocation,value(ro),value(kappa),value(m_bar),age,avg_wt,
	//					avg_fec,va,re,ye,be,ve,dye_df,d2ye_df2);
	//			#endif
	//			if(re<=0)break;
	//		
	//			double spr = value(-ro/((kappa-1)*re-ro*kappa));
	//			report_file<<i++<<"\t"<<fe<<"\t"<<ye<<"\t"<<be<<"\t"<<ve<<"\t";
	//			report_file<<re<<"\t"<<spr<<endl;
	//		
	//			fe += 0.01;
	//		}
	//	}
	//}//exit(1);
	
	
	//cout<<"Ro     "<<d_ro<<endl;
	//cout<<"h      "<<d_h<<endl;
	//cout<<"m      "<<d_m<<endl;
	//cout<<"wa     "<<d_wa<<endl;
	//cout<<"fa     "<<d_fa<<endl;
	//cout<<"V      "<<d_V<<endl;
	/*
		TODO Need to rethink this, should call equibrium with calc_equilirbium(fe,allocation)
		Then loop over values of F and return fe=lambda*F to satisfy allocation scheme.
	*/
	if(!mceval_phase())
	{
		Msy cRFP(d_ro,d_h,d_m,d_rho,d_wa,d_fa,d_V);
		double fmult;
		dvector fe(1,nfleet);
		dvector fadj(1,nfleet);
		
		ofstream report_file("iscam.eql");
		if(report_file.is_open())
		{
			report_file<<"      index";
			for(k=1;k<=nfleet;k++) report_file<<"        fe"<<k;
			for(k=1;k<=nfleet;k++) report_file<<"        ye"<<k;
			for(k=1;k<=nfleet;k++) report_file<<"       dye"<<k;
			report_file<<"         be";
			report_file<<"         re";
			report_file<<"        spr";
			report_file<<endl;
			
			fmult = 0; i=1;
			while(i<300)
			{
				fe = fmult*fmsy;
				cRFP.calc_equilibrium(fe);
				report_file<<setw(11)<<i++;
				report_file<<setw(10)<<fe;
				report_file<<setw(10)<<cRFP.getYe();
				report_file<<setw(10)<<cRFP.getdYe();
				report_file<<setw(11)<<cRFP.getBe();
				report_file<<setw(11)<<cRFP.getRe();
				report_file<<setw(11)<<cRFP.getSpr();
				report_file<<endl;

				fmult += 0.01;
			}
		}
	}
	//exit(1);
	if(verbose)cout<<"**** Ok after calc_reference_points ****"<<endl;
  }
	
FUNCTION void simulationModel(const long& seed)
  {
	/*
	Call this routine to simulate data for simulation testing.
	The random number seed can be used to repeat the same 
	sequence of random number for simulation testing.
	
	Implemented using the "-SimFlag 99" command line option where
	99 is the random number seed.
	
	-SimFlag 99 is a special case used for the manuscript for case 1.
	-SimFlag 000 is a special case with 0 error (exact data)
	
	-This routine will over-write the observations in memory
	with simulated data, where the true parameter values are
	the initial values.  Change the standard deviations of the 
	random number vectors epsilon (observation error) or 
	recruitment devs wt (process error).

	Jan 30, 2013  SM
	- added some simulation controls to the control file for simulating
	  fake data.  These controls include:
	  - sim_control(1) -> flag modeling selectivity using Ideal Free Distribution.
	  - sim_control(2) -> 
	*/
	
	bool pinfile = 0;
	cout<<"___________________________________________________\n"<<endl;
	cout<<"  **Implementing Simulation--Estimation trial**    "<<endl;
	cout<<"___________________________________________________"<<endl;
	if(norm(log_rec_devs)!=0)
	{
		cout<<"\tUsing pin file for simulation"<<endl;
		pinfile = 1;
	}
	cout<<"\tRandom Seed No.:\t"<< rseed<<endl;
	cout<<"___________________________________________________\n"<<endl;
	
	
	//Indexes:
	int i,j,k,ii,ki;

	// Initialize selectivity based on model parameters.
    calcSelectivities(sim_isel_type);
    
    // Initialize natural mortality rates.  Should add Random-walk in M parameters here.
    calcTotalMortality();

	
	/*----------------------------------*/
	/*	-- Generate random numbers --	*/
	/*----------------------------------*/
	random_number_generator rng(seed);
	dvector wt(syr-nage-1,nyr);			//recruitment anomalies if !pinfile
	dmatrix epsilon(1,nit,1,nit_nobs);  //observation errors in survey
	
	double sig = value(rho/varphi);
	double tau = value((1.-rho)/varphi);
	COUT(varphi);
	if(seed==000)
	{
		cout<<"No Error\n";
		sig=0;
		tau=0;
	}
	wt.fill_randn(rng); wt *= tau;
	epsilon.fill_randn(rng); 
	//now loop over surveys and scale the observation errors
	for(k=1;k<=nit;k++)
	{
		for(j=1;j<=nit_nobs(k);j++)
		{
			if(it_wt(k,j)!=0)
			epsilon(k,j) *= sig/it_wt(k,j);
		}
	}
	COUT(wt);
	COUT(epsilon);
	cout<<"	OK after random numbers\n";
	/*----------------------------------*/
	
	
	
	
	
	/*----------------------------------*/
    /*		--Initialize model--		*/
	/*CHANGED now calculating phie based on m_bar and avg_fec*/
	/*----------------------------------*/
	dvector lx=pow(exp(-value(m_bar)),age-min(age));
	lx(nage)/=(1.-exp(-value(m_bar)));
	double phie=(lx*exp(-value(m_bar)*cntrl(13)))*avg_fec;//fec(syr);
	so=kappa/phie;
	
	
	if(cntrl(2)==1) beta=(kappa-1.)/(ro*phie);
	if(cntrl(2)==2) beta=log(kappa)/(ro*phie);
	
	//Initial numbers-at-age with recruitment devs
	N.initialize();
	if(cntrl(5))    //If initializing in at unfished conditions
	{	
		log_rt(syr) = log(ro);
		N(syr)      = ro * lx;
	}
	else            //If starting at fished conditions
	{
		log_rt(syr)         = log_avgrec+log_rec_devs(syr);
		N(syr,sage)         = mfexp(log_rt(syr));
		dvar_vector tmpr    = mfexp(log_recinit + init_log_rec_devs(sage+1,nage));
		N(syr)(sage+1,nage) = elem_prod(tmpr,lx(sage+1,nage));
	}
	
	for(i=syr+1;i<=nyr;i++){
		log_rt(i)=log_avgrec+log_rec_devs(i);
		N(i,sage)=mfexp(log_rt(i));
	}
	N(nyr+1,sage)=mfexp(log_avgrec);
	cout<<"	Ok after initialize model\n";
	/*----------------------------------*/
	
	
	
	/*----------------------------------*/
    /*		--    Selectivity   --		*/
	/*----------------------------------*/
	/*
		-Based on values in the control file.
		-Add autocorrelated random numbers
		for time varying or something to that
		effect.
		
		-If seed==99 then set up a special case
		for the cubic spline manunscript using
		the eplogistic function where g goes from
		strongly domed to asymptotic, e.g.,
		g = 0.2 * (nyr-i)/(nyr-syr);

		-If simulation blocks differ from estimation blocks,
		then need to loop over years and change selectivity
		to correspond to the simulation blocks.
		
	*/

	
	dmatrix va(1,ngear,sage,nage);			//fishery selectivity
	d3_array dlog_sel(1,ngear,syr,nyr,sage,nage);

	// Check to see if Selectivity blocks differ from estimation blocks.
	int byr;
	double p1,p2,xx;
	
	for(k=1;k<=ngear;k++)
	{
		byr = 1;
		if( nsim_sel_blocks(k) > 1 )
		{
			p1 = mfexp(value(sel_par(k,1,1)));
			p2 = mfexp(value(sel_par(k,1,2)));

			for(i=syr; i<=nyr; i++)
			{
				if( i == sim_sel_blocks(k,byr) )
				{
					if( byr <= nsim_sel_blocks(k) )
					{
						byr++;

						// Get new selectivity parameters
						xx  = randn(rng);
						p1 *= exp(0.1 * xx);
						xx  = randn(rng);
						p2 *= exp( 0.05 * xx);
					}
				}				
				log_sel(k)(i) = log( plogis(age,p1,p2) );
			}
		}
	}
	
	dlog_sel=value(log_sel);

	
	/*
		for(i=syr;i<=nyr;i++)
		{
			//sel(k)(i)=plogis(age,ahat(k),ghat(k));
			log_sel(k)(i)=log(plogis(age,ahat(k),ghat(k)));
		log_sel(k)(i) -= log(mean(exp(log_sel(k)(i))));
		}
	//log_sel(j)(i) -= log(mean(mfexp(log_sel(j)(i))));
	*/
	cout<<"	Ok after selectivity\n";

	/*----------------------------------*/
	
	
	/*----------------------------------*/
    /*	--  Population dynamics  --		*/
	/*----------------------------------*/
	
	dmatrix zt(syr,nyr,sage,nage);			//total mortality
	zt.initialize();
	dmatrix ft(syr,nyr,1,ngear);
	ft.initialize();
	dvector sbt(syr,nyr+1);
	sbt.initialize();
	
	
	for(i=syr;i<=nyr;i++)
	{   
		
		//total biomass at age
		//dvector bt = elem_prod(value(N(i)),wa);
		dvector bt = elem_prod(value(N(i)),wt_obs(i));

		/*calculate instantaneous fishing mortalities
		based on Baranov's catch equation and the 
		observed catch from each fleet.*/
		dvector oct = trans(obs_ct)(i);
		
		/*
		Feb 1, 2013.  SM
		Added simulation options for selectivity in:
		sim_ctrl(1) -> uses Ideal Free Distribution to modify dlog_sel

		*/
		for(k=1;k<=ngear;k++)
		{
			va(k)=exp(dlog_sel(k)(i));
			if( sim_ctrl(1)(k) )
			{
				va(k) = ifdSelex(va(k),bt);
				dlog_sel(k)(i) = log(va(k));
			}
		}
		
		//get_ft is defined in the Baranov.cxx file
		//CHANGED these ft are based on biomass at age, should be numbers at age
		//ft(i) = get_ft(oct,value(m),va,bt);
		//ft(i) = get_ft(oct,value(m),va,value(N(i)),wt_obs(i));
		ft(i) = getFishingMortality(oct, value(m), va, value(N(i)),wt_obs(i));
		//cout<<"ft\t"<<ft(i)<<endl;
		//cout<<trans(obs_ct)(i)<<"\t"<<oct<<endl;
		
		// overwrite observed catch incase it was modified by get_ft
		// SJDM Deprecated Sept 10, 2012
		//for(k=1;k<=ngear;k++)
		//	obs_ct(k,i)=oct(k);
		
		//total age-specific mortality
		//dvector zt(sage,nage);
		zt(i)=value(m);
		for(k=1;k<=ngear;k++){
			zt(i)+= ft(i,k)*exp(dlog_sel(k)(i));
		}
		
		
		//CHANGED definition of spawning biomass based on ctrl(13)
		sbt(i) = value(elem_prod(N(i),exp(-zt(i)*cntrl(13)))*fec(i));
		
		//Update numbers at age
		// SM Changed to allow for pinfile to over-ride S_R relationship.
		if(i>=syr+sage-1 && !pinfile)
		{
			double rt;
			//double et=value(N(i-sage+1))*fec(i-sage+1);
			double et=sbt(i-sage+1);
			if(cntrl(2)==1)rt=value(so*et/(1.+beta*et));
			if(cntrl(2)==2)rt=value(so*et*exp(-beta*et));
			N(i+1,sage)=rt*exp(wt(i)-0.5*tau*tau);
			
			/*CHANGED The recruitment calculation above is incosistent
			  with the assessment model.  Below recruitment is based on
			  rt=exp(log_avgrec + wt + rt_dev), where the rt_dev calculation
			is based on the BH or Ricker model.*/
			//double rt_dev = log(rt)-value(log_avgrec);
			//N(i+1,sage)=exp(log_avgrec+wt(i));
			
		}


		N(i+1)(sage+1,nage)=++elem_prod(N(i)(sage,nage-1),exp(-zt(i)(sage,nage-1)));
		N(i+1,nage)+=N(i,nage)*exp(-zt(i,nage));
		
		
		//Catch & Catch-at-age
		for(k=1;k<=ngear;k++)
		{
			if(ft(i,k)>0)
			{
				dvector sel = exp(dlog_sel(k)(i));
				d3C(k)(i)=elem_prod(elem_div(ft(i,k)*sel,zt(i)),elem_prod(1.-exp(-zt(i)),value(N(i))));
				obs_ct(k,i)=d3C(k)(i)*wt_obs(i);
			}
			else	//if this is a survey
			{
				dvector sel = exp(dlog_sel(k)(i));
				d3C(k)(i)=elem_prod(elem_div(sel,zt(i)),elem_prod(1.-exp(-zt(i)),value(N(i))));
			}
		}
		
	}
	// COUT(N);
	
	//initial values of log_ft_pars set to true values
	ki=1;
	for(k=1;k<=ngear;k++)
		for(i=syr;i<=nyr;i++)
			if(obs_ct(k,i)>0){
				log_ft_pars(ki++)=log(ft(i,k));
			}
	
	// Error handler to inform user population went extinct.
	if(min(sbt(syr,nyr))<=1.e-5)
	{
		cout<<"---------------------------------------\n";
		cout<<"Simulated population went extinct, try\n";
		cout<<"increasing steepness, Ro and Rbar\n";
		cout<<sbt<<endl;
		cout<<"Minimum spawning biomass="<<min(sbt(syr,nyr))<<endl;
		cout<<"---------------------------------------\n";
		exit(1);
	}
	
	//Average recruitment calculation
	
	
	cout<<"	log(mean(column(N,sage))) = "<<mean(log(column(N,sage)))<<endl;
	cout<<"	log_avgrec = "<<log_avgrec<<endl;
	cout<<"	Ok after population dynamics\n";
	/*----------------------------------*/
	
	/*----------------------------------*/
    /*	--  Observation models  --		*/
	/*----------------------------------*/
	
	//Simulated Age-compositions
	int ig;
	double age_tau = sqrt(value(rho/varphi));
	for(k=1;k<=na_gears;k++)
	{
		for(i=1;i<=na_nobs(k);i++)
		{
			ii=A(k,i,a_sage(k)-2);	//index for year
			ig=A(k,i,a_sage(k)-1);	//index for gear
			dvector pa = d3C(ig)(ii);	//
			pa/=sum(pa);
			
			
			dvector t1=pa(a_sage(k),a_nage(k));
			t1/=sum(t1);
			A(k)(i)(a_sage(k),a_nage(k))=rmvlogistic(t1,age_tau,i+seed);
			if(seed==000)
			{
				A(k)(i)(a_sage(k),a_nage(k))=t1;
			}
			//cout<<iyr<<"\t"<<k<<endl;
		}
	}

	//cout<<Ahat<<endl;
	
	//Relative abundance indices
	//CHANGED fixed this to reflect survey timing etc & survey_type
	for(k=1;k<=nit;k++)
	{   
		for(i=1;i<=nit_nobs(k);i++)
		{
			ii=iyr(k,i);
			ig=igr(k,i);
			dvector sel = exp(dlog_sel(ig)(ii));
			dvector Np = value(elem_prod(N(ii),exp(-zt(ii)*it_timing(k,i))));
			switch(survey_type(k))
			{
				case 1: //survey based on numbers
					Np = elem_prod(Np,sel);
				break;
				case 2: //survey based on biomass
					Np = elem_prod(elem_prod(Np,sel),wt_obs(ii));
				break;
				case 3: //survey based on spawning biomass
					Np = elem_prod(Np,fec(ii));
				break;
			}
			it(k,i) = sum(Np) * exp(epsilon(k,i));
		}
	}
	
	

	cout<<"	OK after observation models\n";
	/*----------------------------------*/
	
	//CHANGED Fixed bug in reference points calc call from simulation model,
	//had to calculate m_bar before running this routine.
	
	calc_reference_points();
	//cout<<"	OK after reference points\n"<<fmsy<<endl;
	//exit(1);
	//	REPORT(fmsy);
	//	REPORT(msy);
	//	REPORT(bmsy);
	
	
	cout<<"___________________________________________________"<<endl;
	ofstream ofs("iscam.sim");
	ofs<<"fmsy\n"<<fmsy<<endl;
	ofs<<"msy\n"<<msy<<endl;
	ofs<<"bmsy\n"<<bmsy<<endl;
	ofs<<"va\n"<<va<<endl;
	ofs<<"sbt\n"<<sbt<<endl;//<<rowsum(elem_prod(N,fec))<<endl;
	ofs<<"rt\n"<<rt<<endl;
	ofs<<"ct\n"<<obs_ct<<endl;
	ofs<<"ft\n"<<trans(ft)<<endl;
	ofs<<"ut\n"<<elem_div(colsum(obs_ct),N.sub(syr,nyr)*wa)<<endl;
	ofs<<"iyr\n"<<iyr<<endl;
	ofs<<"it\n"<<it<<endl;
	ofs<<"N\n"<<N<<endl;
	ofs<<"A\n"<<A<<endl;
	ofs<<"dlog_sel\n"<<dlog_sel<<endl;
	cout<<"  -- Simuation results written to iscam.sim --\n";
	cout<<"___________________________________________________"<<endl;
	
	//cout<<N<<endl;
	//exit(1);
  }
FUNCTION dvector ifdSelex(const dvector& va, const dvector& ba)
  {
  	/*
  	This function returns a modified selectivity vector (va) based on
  	the assumption that age-based selectivity will operate on the principle
  	of ideal free distribution.  

  	va -> is a vector representing the selectivity of the gear
  	ba -> is a vector of relative biomasses-at-age.  
  	*/
  	dvector pa(sage,nage);

  	pa = (elem_prod(va,pow(ba,0.25)));
  	// pa = (pa - mean(pa))/sqrt(var(pa));
  	// pa = va * exp(0.1*pa);
  	pa = pa/sum(pa);
  	pa = exp( log(pa) - log(mean(pa)) );
  	return (pa);
  }

REPORT_SECTION
  {
	if(verbose)cout<<"Start of Report Section..."<<endl;
	int i,j,k;
	report<<DataFile<<endl;
	report<<ControlFile<<endl;
	report<<ProjectFileControl<<endl;
	REPORT(f);
	REPORT(nlvec);
	REPORT(ro);
	double rbar=value(exp(log_avgrec));
	REPORT(rbar);
	double rinit=value(exp(log_recinit));
	REPORT(rinit);
	REPORT(sbo);
	REPORT(kappa);
	double steepness=value(theta(2));
	REPORT(steepness);
	REPORT(m);
	double tau = value((1.-rho)/varphi);
	double sig = value(rho/varphi);
	REPORT(tau);
	REPORT(sig);
	REPORT(age_tau2);
	
	ivector yr(syr,nyr);
	ivector yrs(syr,nyr+1);
	yr.fill_seqadd(syr,1); 
	yrs.fill_seqadd(syr,1); 
	REPORT(ngear);
	REPORT(yr);
	REPORT(yrs);
	REPORT(iyr);
	REPORT(age);
	REPORT(la);
	REPORT(wa);
	REPORT(fec);
	//Selectivity
	report<<"log_sel"<<endl;
	for(k=1;k<=ngear;k++)
		for(i=syr;i<=nyr;i++)
			report<<k<<"\t"<<log_sel(k)(i)<<endl;
	//REPORT(log_sel);
	REPORT(vax);
	REPORT(obs_ct);
	REPORT(ct);
	REPORT(eta);
	REPORT(ft);
	/*FIXED small problem here with array bounds if using -retro option*/
	report<<"ut\n"<<elem_div(colsum(obs_ct)(syr,nyr),N.sub(syr,nyr)*wa)<<endl;
	report<<"bt\n"<<rowsum(elem_prod(N,wt_obs))<<endl;
	report<<"sbt\n"<<sbt<<endl;

	int rectype=int(cntrl(2));
	REPORT(rectype);
	REPORT(rt);
	dvector ln_rt=value(log_rt(syr,nyr));
	REPORT(ln_rt);
	REPORT(delta);
	REPORT(q);
	REPORT(qt);
	REPORT(it);
	REPORT(pit);
	REPORT(it_wt);
	REPORT(epsilon);
	REPORT(F);
	REPORT(M_tot);
	
	REPORT(a_sage);
	REPORT(a_nage);
	REPORT(A); 
	REPORT(Ahat);
	REPORT(A_nu);
	REPORT(N);
	REPORT(wt_obs);
	
	if(last_phase())
	{
		calc_reference_points();
		REPORT(bo);
		REPORT(fmsy);
		REPORT(msy);
		REPORT(bmsy);
		REPORT(Umsy);
	}
	/*
	Stock status info
	Bstatus = sbt/bmsy;
	Fstatus = ft/fmsy; If fmsy > 0 
	*/
	if(bmsy>0)
	{
		dvector Bstatus=value(sbt/bmsy);
		REPORT(Bstatus);
	}
	
	dmatrix Fstatus(1,ngear,syr,nyr);
	Fstatus.initialize();
	for(k = 1; k <= nfleet; k++)
	{
		if(fmsy(k) >0 )
		{
			j    = ifleet(k);
			Fstatus(j) = value(ft(j)/fmsy(k));
		}
	}
	REPORT(Fstatus);
	
	//Parameter controls
	dmatrix ctrl=theta_control;
	REPORT(ctrl);
	
	
	if(last_phase()) decision_table();
	
	
	dvector rt3(1,3);
	if(last_phase())
	{
		dvector rt3 = age3_recruitment(value(column(N,3)),wt_obs(nyr+1,3),value(M_tot(nyr,3)));
		REPORT(rt3);
	}
	
	//dvector future_bt = value(elem_prod(elem_prod(N(nyr+1),exp(-M_tot(nyr))),wt_obs(nyr+1)));
	dvector future_bt = value(elem_prod(N(nyr+1)*exp(-m_bar),wt_obs(nyr+1)));
	REPORT(future_bt);
	double future_bt4 = sum(future_bt(4,nage));
	REPORT(future_bt4);
	

	
	if(verbose)cout<<"END of Report Section..."<<endl;
	
	/*IN the following, I'm renaming the report file
	in the case where retrospective analysis is occurring*/
	if(retro_yrs && last_phase() && PLATFORM =="Linux")
	{
		//adstring rep="iscam.ret"+str(retro_yrs);
		//rename("iscam.rep",rep);
		adstring copyrep = "cp iscam.rep iscam.ret"+str(retro_yrs);
		system(copyrep);
	}

	if(retro_yrs && last_phase() && PLATFORM =="Windows")
	{
		//adstring rep="iscam.ret"+str(retro_yrs);
		//rename("iscam.rep",rep);
		adstring copyrep = "copy iscam.rep iscam.ret"+str(retro_yrs);
		system(copyrep);
	}
	
	
  }
	
FUNCTION decision_table
  {
	/*
	This function takes a vector of projected catches and computes the following
	Reference points: Bmsy, Bo, Fmsy, Umsy.
	
	Biomass Metrics for the decision table:
	1) P(SB_{t+1} < SB_{t})
	2) P(SB_{t+1} < 0.25 B_{0})
	3) P(SB_{t+1} < 0.75 B_{0})
	4) P(SB_{t+1} < 0.40 B_{MSY})
	5) P(SB_{t+1} < 0.80 B_{MSY})
	
	Harvest Metrics for the decision table:
	1) P(U_{t+1} > Target harvest rate)
	2) P(U_{t+1} > 1/2 Fmsy)
	3) P(U_{t+1} > 2/3 Fmsy)
	4) P(tac/3+  > 20%)
	
	Key to the harvest metric is the definition of Umsy and allocation to fleets.
	
	Pseudocode:
		1) Calculate reference points (Fmsy, Bmsy)
		2) Loop over vector of proposed catches
		3) Evaluate biomass metrics for each posterior sample
		4) Evaluate harvest metrics for each posterior sample
	
	*/
	int i;
	// 1) Calculate reference pionts.
	//calc_reference_points();  //redundant b/c its called in mcmc_output?
	
	// 2) Loop over vector of proposed catches
	//    This vector should is now read in from the projection file control (pfc).
	for(i=1;i<=n_tac;i++)
	{
		projection_model(tac(i));
	}
  }
	
FUNCTION mcmc_output
  {
	if(nf==1){
		ofstream ofs("iscam.mcmc");
		ofs<<"       log.ro";
		ofs<<"        log.h";
		ofs<<"        log.m";
		ofs<<"     log.rbar";
		ofs<<"    log.rinit";
		ofs<<"          rho";
		ofs<<"     vartheta";
		ofs<<"           bo";
		ofs<<"         bmsy";
		for(int k=1;k<=nfleet;k++) ofs<<"         msy"<<k;
		for(int k=1;k<=nfleet;k++) ofs<<"        fmsy"<<k;
		ofs<<"          SSB";
		ofs<<"        Age-4";
		ofs<<"         Poor";
		ofs<<"      Average";
		ofs<<"         Good";
		for(int i=1;i<=nit;i++) ofs<<"         lnq"<<i;
		ofs<<"            f";
		ofs<<endl;
		
		ofstream of1("sbt.mcmc");
		ofstream of2("rt.mcmc");
		
	}
	
	// leading parameters & reference points
	calc_reference_points();
	// decision table output
	//dvector future_bt = value(elem_prod(elem_prod(N(nyr+1),exp(-M_tot(nyr))),wt_obs(nyr+1)));
	dvector future_bt = value(elem_prod(N(nyr+1)*exp(-m_bar),wt_obs(nyr+1)));
	double future_bt4 = sum(future_bt(4,nage));
	dvector rt3 = age3_recruitment(value(column(N,3)),wt_obs(nyr+1,3),value(M_tot(nyr,3)));	
	
	
	//if( bmsy > 0 && min(fmsy) >= 0 )
	{
		ofstream ofs("iscam.mcmc",ios::app);
		ofs<<setw(12)<<theta;
		ofs<<setw(13)<< bo;
		ofs<<setw(13)<< bmsy;
		ofs<<setw(12)<< msy;
		ofs<<setw(12)<< fmsy;
		ofs<<setw(13)<< sbt(nyr);
		ofs<<setw(13)<< future_bt4;
		ofs<<setw(12)<< future_bt4+rt3;
		ofs<<setw(12)<< log(q);
		ofs<<setw(13)<< f;
		ofs<<endl;
		
		/* June 12, 2012.  SJDM Call decision table. */
		decision_table();  
	}
	// output spawning stock biomass
	ofstream of1("sbt.mcmc",ios::app);
	of1<<setw(10)<<sbt(syr,nyr)<<endl;
	
	// output age-1 recruits
	ofstream of2("rt.mcmc",ios::app);
	of2<<setw(10)<<rt<<endl;
	
	
  }

FUNCTION dvector age3_recruitment(const dvector& rt, const double& wt,const double& M)
  {
	/*
	This routine returns the poor average and good age-3 recruits
	that is used in constructing the decision table for the pacific
	herring fisheries.
	
	-1) sort the rt vector from small to large
	-2) find the 33rd and 66th percentiles
	-3) take the averge of the 0-33, 34-66, 67-100
	-4) multiply by the average weight
	-5) return the age-3 recruitment biomass
	*/
	
	dvector s_rt = sort(rt);
	dvector rbar(1,3);
	
	double idx = floor((nyr-syr+1.)/3.);
	int ix1 = syr+int(idx);
	int ix2 = syr+int(2.*idx);
	rbar(1) = mean(s_rt(syr,ix1));
	rbar(2) = mean(s_rt(ix1+1,ix2));
	rbar(3) = mean(s_rt(ix2+1,nyr));
	rbar = rbar*wt*exp(-M);
	//cout<<rbar<<endl;
	return(rbar);
  }

FUNCTION void projection_model(const double& tac);
  {
	/*
	This routine conducts population projections based on 
	the estimated values of theta.  Note that all variables
	in this routine are data type variables.
	
	Arguments:
	tac is the total allowable catch that must be allocated 
	to each gear type based on allocation(k)
	
	theta(1) = log_ro
	theta(2) = h
	theta(3) = log_m
	theta(4) = log_avgrec
	theta(5) = log_recinit
	theta(6) = rho
	theta(7) = vartheta
	
	** NOTES **
	* Projections are based on average natural mortality and fecundity.
	* Selectivity is based on selectivity in terminal year.
	* Average weight-at-age is based on mean weight in the last 5 years.
	
	Aug 20, 2012 Found a bug, See issue 2 on github.
	*/
	static int runNo=0;
	runNo ++;
	int i,j,k;
	int pyr = nyr+1;	//projection year.
	
	// --derive stock recruitment parameters
	// --survivorship of spawning biomass
	dvector lx(sage,nage);
	double   tau = value((1.-rho)/varphi); 
	double   m_M = value(m_bar); 
	double m_rho = cntrl(13);
	lx(sage)     = 1;
	for(i=sage; i<=nage; i++)
	{
		lx(i) = exp( -m_M*(i-sage) -m_rho*m_M );
		if(i==nage) 
			lx(i) /= 1.0 - exp( -m_M );
	}	
	double phib = lx*avg_fec;
	double so   = value(kappa)/phib;
	double bo   = value(ro)*phib;
	
	double beta;
	switch(int(cntrl(2)))
	{
		case 1:  // Beverton-Holt
			beta = (value(kappa)-1.)/bo;
		break;
		case 2:  // Ricker
			beta = log(value(kappa))/bo;
		break;
	}
	
	/* Fill arrays with historical values */
	dvector p_sbt(syr,pyr);
	dvector  p_ct(1,ngear);
	dmatrix  p_ft(nyr+1,pyr,1,ngear);
	dmatrix   p_N(syr,pyr+1,sage,nage);
	dmatrix   p_Z(syr,pyr,sage,nage);
	p_N.initialize();
	p_sbt.initialize();
	p_Z.initialize();
	p_N.sub(syr,nyr)   = value( N.sub(syr,nyr) );
	p_sbt(syr,nyr)     = value( sbt(syr,nyr)   );
	p_Z.sub(syr,nyr)   = value( Z.sub(syr,nyr) );
	
	
	/* Selectivity and allocation to gears */
	dmatrix va_bar(1,ngear,sage,nage);
	for(k=1;k<=ngear;k++)
	{
		p_ct(k)   = allocation(k)*tac;
		va_bar(k) = exp(value(log_sel(k)(nyr)));
	}
		
	/* Simulate population into the future under constant tac policy. */
	for(i = nyr; i<=pyr; i++)
	{
		
		if(i > nyr)
		{
			// get_ft is defined in the Baranov.cxx file
			p_ft(i) = getFishingMortality(p_ct, value(m_bar), va_bar, p_N(i),avg_wt);
			// p_ft(i)(1,nfleet) = fmsy(1,nfleet);
			
			// calculate total mortality in future years
			p_Z(i) = value(m_bar);
			for(k=1;k<=ngear;k++)
			{
				p_Z(i)+=p_ft(i,k)*va_bar(k);
			}
		}
		
		
		// spawning biomass
		p_sbt(i) = elem_prod(p_N(i),exp(-p_Z(i)*cntrl(13))) * avg_fec;
		
		
		// sage recruits with random deviate xx
		// note the random number seed is repeated for each tac level.
		double  xx = randn(nf+i)*tau;
		if(i>=syr+sage-1)
		{
			double rt;
			double et = p_sbt(i-sage+1);		// lagged spawning biomass
			
			if(cntrl(2)==1)						// Beverton-Holt model
			{
				rt=(so*et/(1.+beta*et));
			}
			if(cntrl(2)==2)						// Ricker model
			{
				rt=(so*et*exp(-beta*et));
			}
			
			p_N(i+1,sage)=rt*exp(xx-0.5*tau*tau); 
		}
		
		/* Update numbers at age in future years */
		p_N(i+1)(sage+1,nage) =++ elem_prod(p_N(i)(sage,nage-1),exp(-p_Z(i)(sage,nage-1)));
		p_N(i+1,nage)        +=   p_N(i,nage)*exp(-p_Z(i,nage));
		
		//Predicted catch for checking calculations
		//for(k=1;k<=nfleet;k++)
		//{
		//	dvector ba = elem_prod(p_N(i),avg_wt);
		//	cout<<k<<" tac = "<<tac<<"\t ct = ";
		//	cout<<sum(elem_div(elem_prod(elem_prod(ba,p_ft(i,k)*va_bar(k)),1.-exp(-p_Z(i))),p_Z(i)));
		//	cout<<" fmsy = "<<fmsy<<" ft = "<<p_ft(i,k)<<endl;
		//}
	}	
	//cout<<"fmsy\n"<<fmsy<<endl;
	//cout<<"Spawning biomass\n"<<p_sbt<<endl;
	//exit(1);
	/* 
	  Write output to *.proj file for constructing decision tables. 
	
	  Biomass Metrics for the decision table:
	  1) P(SB_{t+1} < SB_{t})
	  2) P(SB_{t+1} < 0.25 B_{0})
	  3) P(SB_{t+1} < 0.75 B_{0})
	  4) P(SB_{t+1} < 0.40 B_{MSY})
	  5) P(SB_{t+1} < 0.80 B_{MSY})
	  
	  Harvest Metrics for the decision table:
	  1) P(U_{t+1} > Umsy)
	  2) P(U_{t+1} > 1/2 Umsy)
	  3) P(U_{t+1} > 2/3 Umsy)
	  4) P(tac/2+  > 20%)
	
	  Metric for spawning depletion and spawning trends:
	  1) P(5-year decline)
	  
	  Defn: Stock status is based on spawning biomass
	  Defn: Removal rate is based on removals/spawning biomass
	  Defn: Harvest rate    : U_{t+1}=TAC/SBio
	  Defn: MSY harvest rate: Umsy=MSY/SBmsy
	  
	
	*/
	if(nf==1 && runNo==1)
	{
		ofstream ofs(BaseFileName + ".proj");
		ofs<<" tac";
		ofs<<"      P(SB1)"; 
		ofs<<"      P(SB2)";
		ofs<<"      P(SB3)";
		ofs<<"      P(SB4)";
		ofs<<"      P(SB5)";
		ofs<<"       P(U1)";
		ofs<<"       P(U2)";
		ofs<<"       P(U3)";
		ofs<<"       P(U4)";
		ofs<<"       P(D5)";
		ofs<<endl;
		cout<<"Bo when nf==1 \t"<<bo<<endl;
	}
	
	double  ut  = tac / p_sbt(pyr);
	double u20  = tac / ( (p_N(pyr)(3,nage)*exp(-value(M_tot(nyr,3))))* avg_wt(3,nage) );
	
	/* Average rate of change in spawning biomass in last 5 years */
	double dSb5 = mean(log(p_sbt(pyr-5,pyr)) - log(p_sbt(pyr-6,pyr-1).shift(pyr-5)));
	
	ofstream ofs(BaseFileName + ".proj",ios::app);
	ofs<< setprecision(4)               <<setw(4) 
	   << tac                           <<setw(12)
	   << p_sbt(pyr-1)/p_sbt(pyr)       <<setw(12)
	   << 0.25*bo/p_sbt(pyr)            <<setw(12)
	   << 0.75*bo/p_sbt(pyr)            <<setw(12)
	   << 0.40*bmsy/p_sbt(pyr)          <<setw(12)
	   << 0.80*bmsy/p_sbt(pyr)          <<setw(12)
	   << ut/Umsy                       <<setw(12)
	   << ut/(0.5*Umsy)                 <<setw(12)
	   << ut/(2./3.*Umsy)               <<setw(12)
	   << u20/0.2                       <<setw(12)
	   << dSb5+1                        <<setw(12)
	   << endl;
	
  }

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

	#if defined(_WIN32) && !defined(__linux__)
		const char* PLATFORM = "Windows";
	#else
		const char* PLATFORM = "Linux";
	#endif

	#include <admodel.h>
	#include <time.h>
	#include <string.h>
	// #include <statsLib.h>
	#include "msy.cpp"
	//#include "stats.cxx"
	#include "baranov.cxx"
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	bool mcmcPhase = 0;
	bool mcmcEvalPhase = 0;
	
	adstring BaseFileName;
	adstring ReportFileName;
	
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
	
	
	// #ifdef __GNUDOS__
	//   #include <gccmanip.h>
	// #endif
	// Variables to store results from DIC calculations.
	double dicNoPar;
	double dicValue;

    void testing(void)
    {
		cout<<"Hello World"<<endl;
		exit(1);
    }


	void function_minimizer::mcmc_eval(void)
	{
		// |---------------------------------------------------------------------------|
		// | Added DIC calculation.  Martell, Jan 29, 2013                             |
		// |---------------------------------------------------------------------------|
		// | DIC = pd + dbar
		// | pd  = dbar - dtheta  (Eeffective number of parameters)
		// | dbar   = expectation of the likelihood function (average f)
		// | dtheta = expectation of the parameter sample (average y) 

	  gradient_structure::set_NO_DERIVATIVES();
	  initial_params::current_phase=initial_params::max_number_phases;
	  uistream * pifs_psave = NULL;

	#if defined(USE_LAPLACE)
	#endif

	#if defined(USE_LAPLACE)
	    initial_params::set_active_random_effects();
	    int nvar1=initial_params::nvarcalc(); 
	#else
	  int nvar1=initial_params::nvarcalc(); // get the number of active parameters
	#endif
	  int nvar;
	  
	  pifs_psave= new
	    uistream((char*)(ad_comm::adprogram_name + adstring(".psv")));
	  if (!pifs_psave || !(*pifs_psave))
	  {
	    cerr << "Error opening file "
	            << (char*)(ad_comm::adprogram_name + adstring(".psv"))
	       << endl;
	    if (pifs_psave)
	    {
	      delete pifs_psave;
	      pifs_psave=NULL;
	      return;
	    }
	  }
	  else
	  {     
	    (*pifs_psave) >> nvar;
	    if (nvar!=nvar1)
	    {
	      cout << "Incorrect value for nvar in file "
	           << "should be " << nvar1 << " but read " << nvar << endl;
	      if (pifs_psave)
	      {
	        delete pifs_psave;
	        pifs_psave=NULL;
	      }
	      return;
	    }
	  }
	  
	  int nsamp = 0;
	  double sumll = 0;
	  independent_variables y(1,nvar);
	  independent_variables sumy(1,nvar);

	  do
	  {
	    if (pifs_psave->eof())
	    {
	      break;
	    }
	    else
	    {
	      (*pifs_psave) >> y;
	      sumy = sumy + y;
	      if (pifs_psave->eof())
	      {
	      	double dbar = sumll/nsamp;
	      	int ii=1;
	      	y = sumy/nsamp;
	      	initial_params::restore_all_values(y,ii);
	        initial_params::xinit(y);   
	        double dtheta = 2.0 * get_monte_carlo_value(nvar,y);
	        double pd     = dbar - dtheta;
	        double dic    = pd + dbar;
	        dicValue      = dic;
	        dicNoPar      = pd;

	        cout<<"Number of posterior samples    = "<<nsamp  <<endl;
	        cout<<"Expectation of log-likelihood  = "<<dbar   <<endl;
	        cout<<"Expectation of theta           = "<<dtheta <<endl;
	        cout<<"Number of estimated parameters = "<<nvar1  <<endl;
		    cout<<"Effective number of parameters = "<<pd     <<endl;
		    cout<<"DIC = "<<dic<<endl;
	        break;
	      }
	      int ii=1;
	      initial_params::restore_all_values(y,ii);
	      initial_params::xinit(y);   
	      double ll = 2.0 * get_monte_carlo_value(nvar,y);
	      sumll    += ll;
	      nsamp++;
	      //cout<<get_monte_carlo_value(nvar,y)<<endl;
	    }
	  }
	  while(1);
	  if (pifs_psave)
	  {
	    delete pifs_psave;
	    pifs_psave=NULL;
	  }
	  return;
	}
	
	
FINAL_SECTION
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

	//Make copies of the report file using the ReportFileName
	//to ensure the results are saved to the same directory 
	//that the data file is in. This should probably go in the 
	//FINAL_SECTION
	
	//CHANGED only copy over the mcmc files if in mceval_phase()
	
	if(last_phase() && PLATFORM =="Linux" && !retro_yrs)
	{
		adstring bscmd = "cp iscam.rep " +ReportFileName;
		system(bscmd);
		
		bscmd = "cp iscam.par " + BaseFileName + ".par";
		system(bscmd); 
		
		bscmd = "cp iscam.std " + BaseFileName + ".std";
		system(bscmd);
		
		bscmd = "cp iscam.cor " + BaseFileName + ".cor";
		system(bscmd);
		
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

	if( last_phase() && PLATFORM =="Linux" && retro_yrs )
	{
		//copy report file with .ret# extension for retrospective analysis
		adstring bscmd = "cp iscam.rep " + BaseFileName + ".ret" + str(retro_yrs);
		system(bscmd);
	}

	if(last_phase() && PLATFORM =="Windows" && !retro_yrs)
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

	if( last_phase() && PLATFORM =="Windows" && retro_yrs )
	{
		//copy report file with .ret# extension for retrospective analysis
		adstring bscmd = "copy iscam.rep " + BaseFileName + ".ret" + str(retro_yrs);
		system(bscmd);
	}



