// ------------------------------------------------------------------------- //
//         integrated Statistical Catch Age Model (iSCAM)                    //
//                                                                           //
//                           VERSION 1.1                                     //
//               Tue Jul 19 22:23:58 PDT 2011                                //
//                                                                           //
//                                                                           //
//           Created by Steven Martell on 2010-04-09                         //
//           Copyright (c) 2010. All rights reserved.                        //
//                                                                           //
// AUTHORS: SJDM Steven Martell                                              //
//                                                                           //
// CONVENTIONS: Formatting conventions are based on the The                  //
//               Elements of C++ Style (Misfeldt et al. 2004)                //
//                                                                           //
// NAMING CONVENTIONS:                                                       //
//             Macros       -> UPPERCASE                                     //
//             Constants    -> UpperCamelCase                                //
//             Functions    -> lowerCamelCase                                //
//             Variables    -> lowercase                                     //
//                                                                           //
// INDEXES:                                                                  //
//             h -> sex                                                      //
//             i -> year                                                     //
//             j -> age                                                      //
//             k -> gear                                                     //
//                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
// CHANGED add option for using empirical weight-at-age data                 //
// TODO:    add gtg options for length based fisheries                      //
// CHANGED add time varying natural mortality rate with splines              //
// TODO:    add cubic spline interpolation for time varying M               //
// CHANGED  Fix the type 6 selectivity implementation. not working.          //
// TODO:  fix cubic spline selectivity for only years when data avail        //
// CHANGED: fixed a bug in the simulation model log_ft_pars goes out         //
//        of bounds.                                                         //
// TODO: write a projection routine and verify equilibrium calcs             //
// TODO: add DIC calculation for MCMC routines (in -mcveal phase)            //
// TODO: add SOK fishery a) egg fishing mort 2) bycatch for closed ponds     //
//                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
// ------------------------------------------------------------------------- //
//-- CHANGE LOG:                                                           --//
//--  Nov 30, 2010 -modified survey biomass by the fraction of total       --//
//--                mortality that occurred during the time of the         --//
//--                survey. User specifies this fraction (0-1) in the      --//
//--                data file as the last column of the relative           --//
//--                abundance index.                                       --//
//--                                                                       --//
//--  Dec 6, 2010 -modified the code to allow for empiracle weight-        --//
//--               at-age data to be used.                                 --//
//--              -rescaled catch and relative abundance /1000, this       --//
//--               should be done in the data file and not here.           --//
//--                                                                       --//
//--  Dec 20, 2010-added prior to survey q's in control file               --//
//--                                                                       --//
//--  Dec 24, 2010-added random walk for natural mortality.                --//
//--                                                                       --//
//--  Jan 23, 2011-in Penticton Hospital with my mom in ICU, adopting      --//
//--               the naming conventions in The Elements of C++           --//
//--               style to keep my mind busy.                             --//
//--                                                                       --//
//-- May 5, 2011- added logistic selectcitivty as a fucntion of            --//
//--              mean body weight.  3 parameter logistic.                 --//
//--              NOT WORKING YET                                          --//
//--                                                                       --//
//-- May 6, 2011- added pre-processor commands to determin PLATFORM        --//
//--              either "Windows" or "Linux"                              --//
//--                                                                       --//
//-- use -mcmult 1.5 for MCMC with log_m_nodes with SOG herrning           --//
//--                                                                       --//
//--                                                                       --//
//-- Dec 11, 2011- added halibut branch to local git repository aim is to  --//
//--               add gender dimension and stock dimension.               --//
//--               This was created on the "twosex" branch in git merged   --//
//--                                                                       --//
//-- Dec 30, 2011- working on length-based selectivity for halibut.        --//
//--                                                                       --//
//-- Feb 5, 2012- working on two sex model for halibut.                    --//
//--                                                                       --//
//--                                                                       --//
//--                                                                       --//
//--                                                                       --//
//-- Jan 5, 2012 - adding spawn on kelp fishery as catch_type ivector      --//
//--             - modified the following routines:                        --//
//--             - calcFisheryObservations                                 --//
//--             - calcTotalMortality                                      --//
//-- TODO: add catch_type to equilibrium calculations for reference points --//
//--                                                                       --//
//-- Feb 14, 2012.  Approval for work on the halibut simulation model      --//
//--                                                                       --//
//-- TODO: add sex-based selectivity paraemters. Currently assuming same   --//
//--       selectivity curves for both sexes.                              --//
//--                                                                       --//
//--                                                                       --//
//--                                                                       --//
//--                                                                       --//
//--                                                                       --//
//--                                                                       --//
//--                                                                       --//
//--                                                                       --//
//--                                                                       --//
// ------------------------------------------------------------------------- //


DATA_SECTION
	!! cout<<"iSCAM has detected that you are on a "<<PLATFORM<<" box"<<endl;
	

	init_adstring DataFile;
	init_adstring ControlFile;
	
	
	!! BaseFileName=stripExtension(ControlFile);
	!! cout<<BaseFileName<<endl;
	!! ReportFileName = BaseFileName + adstring(".rep");
	
	!! ad_comm::change_datafile_name(DataFile);
	
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
			cout<<"______________________________________________________\n"<<endl;
			cout<<"    **Implementing Retrospective analysis** "<<endl;
			cout<<"    **Number of retrospective years = "<<retro_yrs<<endl;
			cout<<"______________________________________________________"<<endl;
		}
	END_CALCS
	
	//Read in objects from data file using the init_ prefix
	init_int syr;
	init_int nyr;
	!! cout<<"syr\t"<<syr<<endl;
	!! cout<<"nyr\t"<<nyr<<endl;
	
	init_int sage;
	init_int nage;
	!! cout<<"sage\t"<<sage<<endl;
	!! cout<<"nage\t"<<nage<<endl;
	vector age(sage,nage);
	!! age.fill_seqadd(sage,1);
	
	init_int ngear;		//number of gear types with unique selectivities
	init_int nsex;		//number of sexes
	!! cout<<"ngear\t"<<ngear<<endl;
	!! cout<<"nsex\t"<<nsex<<endl;
	init_vector allocation(1,ngear);
	init_ivector catch_type(1,ngear);
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
	
	init_vector fixed_m(1,nsex);		//FIXME: depricate this from data files
	init_vector linf(1,nsex);
	init_vector vonbk(1,nsex);
	init_vector to(1,nsex);
	init_vector a(1,nsex);
	init_vector b(1,nsex);
	init_vector ah(1,nsex);
	init_vector gh(1,nsex);
	
	matrix la(1,nsex,sage,nage);		//length-at-age
	matrix wa(1,nsex,sage,nage);		//weight-at-age
	LOC_CALCS
	  cout<<"linf\t"<<linf<<endl;
	  int h;
	  for(h=1;h<=nsex;h++)
	  {
	  	la(h) = linf(h)*( 1.-exp(-vonbk(h)*(age-to(h))) );
	  	wa(h) = a(h)*pow(la(h),b(h));
	  }
	  //la=linf*(1.-exp(-vonbk*(age-to)));
	  //wa=a*pow(la,b);
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
	END_CALCS
	
	
	//Age comps for all gears.
	init_int na_gears;	//total number of aging observations
	init_ivector na_nobs(1,na_gears);
	init_ivector a_sage(1,na_gears);	//youngest age in the ageing matrix
	init_ivector a_nage(1,na_gears);	//oldest age in the ageing matrix

	init_3darray A(1,na_gears,1,na_nobs,a_sage-2,a_nage);
	
	//Mean weight-at-age data (units are kg) (if exists)(by sex)
	init_int n_wt_nobs;
	init_3darray tmp_wt_obs(1,nsex,1,n_wt_nobs,sage-1,nage);
	
	3darray wt_obs(1,nsex,syr,nyr+1,sage,nage);		//weight-at-age by sex
	3darray wt_dev(1,nsex,syr,nyr+1,sage,nage);		//standardized deviations in weight-at-age
	3darray fec(1,nsex,syr,nyr+1,sage,nage);		//fecundity-at-age
	matrix avg_fec(1,nsex,sage,nage);				//average fecundity-at-age
	LOC_CALCS
		int j,jyr;
		avg_fec.initialize();
		for(h=1;h<=nsex;h++)
		{
			for(j=syr;j<=nyr+1;j++)
			{
				wt_obs(h)(j)=wa(h);			
				fec(h)(j)=elem_prod(plogis(age,ah(h),gh(h)),wt_obs(h)(j));
			}
		}
		//if empiracle weight-at-age data exist, the overwrite wt_obs & fec.
		for(h=1;h<=nsex;h++)
		{
			for(j=1;j<=n_wt_nobs;j++)
			{
				jyr=tmp_wt_obs(h)(j,sage-1);  //index for year
				wt_obs(h)(jyr)=tmp_wt_obs(h)(j)(sage,nage);
				fec(h)(jyr)=elem_prod(plogis(age,ah(h),gh(h)),wt_obs(h)(jyr));
			}
		}
		//Average fecundity
		for(h=1;h<=nsex;h++)
		{
			int nfec = fec(h).rowmax()-fec(h).rowmin()+1;
			avg_fec(h)=colsum(fec(h))/nfec;
			//from Jake Schweigert: use mean-weight-at-age data
			//from the last 5 years for the projected mean wt.
			dvector tmp=colsum(wt_obs(h).sub(nyr-5,nyr))/6.;
			wt_obs(h)(nyr+1) = tmp;
			//July 14, 2011 Error handler to ensure non-zero wt_obs
			//in the data. Otherwise causes and error with weight-based
			//selectivities.
			if(min(wt_obs(h))==0)
			{
				cout<<"Cannont have a observed 0 mean weight at age\n";
				cout<<"in the data file.  Please fix.\n Aborting program!"<<endl;
				exit(2);
			}
		}
		cout<<"n_wt_nobs\t"<<n_wt_nobs<<endl;
		cout<<"Ok after empiracle weight-at-age data"<<endl;
		
		/*	
			FEB 14. Changed documentation regarding selectivity based
			on deviations in mean weight-at-age.
			
			Trying a compromize where estimating an additional parameter
			that attemps to explain residual variation in age-comps
			via changes in standardized mean weights at age. This is 
			implemented as:
			
			log_sel = log(plogis(age,s1,s2)) + s3*wa_dev
			where wa_dev is a matrix of standardized deviations
			(mu=0, std=1) of weights at age.  delta is calculated
			based on the wt_dev matrix above.
		*/
		wt_dev.initialize();
		for(h=1;h<=nsex;h++)
		{
			dmatrix mtmp = trans(wt_obs(h));
			for(i=sage;i<=nage;i++)
			{
				dvector wa_dev = (mtmp(i)-mean(mtmp(i)))/sqrt(var(mtmp(i)));
				mtmp(i) = wa_dev;
			}
			wt_dev(h) = trans(mtmp);	//each column has mean=0 sd=1
		}
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

	
	
	
	number fmsy;					//Fishing mortality rate at Fmsy
	number msy;						//Maximum sustainable yield
	number bmsy;					//Spawning biomass at MSY
	vector age_tau2(1,na_gears);	//MLE estimate of the variance for the age comps
	//catch-age for simulation model (could be declared locally 3d_array)
	3darray d3C(1,ngear,syr,nyr,sage,nage);		
	
	
	
	
	
	
	
	
	
	
	// ***************************************************
	// ** Read parameter controls from ControlFile
	// ***************************************************
	!! ad_comm::change_datafile_name(ControlFile);
	
	init_int npar;
	init_matrix theta_control(1,npar,1,7);
	!! cout<<theta_control<<endl;
	
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
	// ** Feb 14, 2012 added index for sex.
	// ** With nsex>1 separate table for female & male
	// ***************************************************
	// type 1 = logistic (2pars)
	// type 2 = selcoffs (A-1 pars)
	// type 3 = cubic spline (age_nodes)
	// type 4 = time varying cubic spline (age_nodes for each year)
	// type 5 = bicubic spline with age_nodes adn yr_nodes
	// type 6 = fixed logistic by turning sel_phz to (-ve)
	// type 7 = logistic (3pars) as a function of body weight.
	init_ivector isel_type(1,ngear);  	//Switch for selectivity
	ivector isel_npar(1,ngear);			//ivector for # of parameters for each gear.
	ivector jsel_npar(1,ngear);			//ivector for the number of rows for time-varying selectivity.
	init_vector ahat(1,ngear);			//age-at-50% vulnerbality for logistic function
	init_vector ghat(1,ngear);			//std at 50% age of vulnerability for logistic funciton
	init_vector age_nodes(1,ngear);		//No. of age-nodes for bicubic spline.
	init_vector yr_nodes(1,ngear);		//No. of year-nodes for bicubic spline.
	init_ivector sel_phz(1,ngear);		//Phase for estimating selectivity parameters.
	init_vector sel_2nd_diff_wt(1,ngear);	//Penalty weight for 2nd difference in selectivity.
	init_vector sel_dome_wt(1,ngear);		//Penalty weight for dome-shaped selectivity.
	!! cout<<isel_type<<endl;
	
	LOC_CALCS
		//cout up the number of selectivity parameters
		//depending on the value of isel_type
		isel_npar.initialize();
		for(i=1;i<=ngear;i++)
		{	
			jsel_npar(i)=1;
			switch(isel_type(i))
			{
				case 1:
					// logistic selectivity
					isel_npar(i) = 2; 
					break;
					
				case 2:
					// age-specific coefficients
					isel_npar(i) = (nage-sage); 
					break;
					
				case 3:
				 	// cubic spline 
					isel_npar(i) = age_nodes(i);
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
					break;
					
				case 8:
					// Alternative logistic selectivity with wt_dev coefficients.
					isel_npar(i) = 3;
					break;
					
				case 11:
					// Logistic length-based selectivity.
					isel_npar(i) = 2;
					break;
					
				case 12:
					// Length-based selectivity coeffs with cubic spline interpolation
					isel_npar(i) = (nage-sage);
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
	
	
	
	//Miscellaneous controls
	// 1 -> verbose
	// 2 -> recruitment model (1=beverton-holt, 2=rickers)
	// 3 -> std in catch first phase
	// 4 -> std in catch in last phase
	// 5 -> assumed unfished in first year (0=FALSE, 1=TRUE)
	// 6 -> minimum proportion at age to consider in the dmvlogistic likelihood
	// 7 -> mean fishing mortality rate to regularize the solution
	// 8 -> standard deviation of mean F penalty in first phases
	// 9 -> standard deviation of mean F penalty in last phase.
	// 10-> phase for estimating deviations in natural mortality.
	// 11-> std in natural mortality deviations.
	// 12-> number of estimated nodes for deviations in natural mortality
	// 13-> fraction of total mortality that takes place prior to spawning
	// 14-> switch for age-composition likelihood (1=dmvlogistic,2=dmultinom)
	// FIXME: document cntrl(14).
	init_vector cntrl(1,14);
	int verbose;
	
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
	
	ivector ilvec(1,6);
	!!ilvec=ngear;			//number of fisheries
	!!ilvec(2)=nit;			//number of surveys
	!!ilvec(3)=na_gears;	//number of age-comps
	!!ilvec(4)=1;
	
	// SM Oct 31, 2010.  Implementing retrospective analysis.
	//Dont read in any more data below the retrospective reset of nyr
	!! nyr = nyr - retro_yrs;
	
	
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
	!! for(int i=1;i<=npar;i++) theta(i)=theta_ival(i);
	
	//Selectivity parameters (A very complicated ragged array)
	//Not sure how to handle this for sex dimension.
	//Could mulitiply jsel_npar by nsex or isel_npar?.
	init_bounded_matrix_vector sel_par(1,ngear,1,jsel_npar,1,isel_npar,-15.,15.,sel_phz);
	LOC_CALCS
		//initial values for logistic selectivity parameters
		//set phase to -1 for fixed selectivity.
		for(int k=1;k<=ngear;k++)
		{
			if( isel_type(k)==1 || isel_type(k)==6 || isel_type(k)>=7 )
			{
				sel_par(k,1,1) = log(ahat(k));
				sel_par(k,1,2) = log(ghat(k));
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
	number bo;					//unfished spawning stock biomass
	number kappa;				//Goodyear compensation ratio
	number log_avgrec;			//log of average recruitment.
	number log_recinit;			//log of initial recruitment in syr.
	//number log_avg_f;			//log of average fishing mortality DEPRICATED
	number rho;					//proportion of the observation error
	number varphi				//total precision in the CPUE & Rec anomalies.
	number so;
	number beta;
	
	vector m(1,nsex);			//initial natural mortality rate
	vector m_bar(1,nsex);		//average natural mortality rate
	vector log_rt(syr-nage+sage,nyr);
	vector vax(sage,nage);		//survey selectivity coefficients
	vector q(1,nit);			//survey catchability coefficients
	
	vector sbt(syr,nyr+1);		//spawning stock biomass
	vector rt(syr+sage,nyr); 	//predicted sage recruits from S-R curve
	
	vector delta(syr+sage,nyr);	//residuals for stock recruitment
	vector avg_log_sel(1,ngear);//conditional penalty for objective function
	vector log_m_devs(syr+1,nyr);// log deviations in natural mortality
	
	matrix nlvec(1,6,1,ilvec);	//matrix for negative loglikelihoods
	
	//matrix jlog_sel(1,ngear,sage,nage);		//selectivity coefficients for each gear type.
	//matrix log_sur_sel(syr,nyr,sage,nage);	//selectivity coefficients for survey.
	matrix ft(1,ngear,syr,nyr);					//Gear specific fishing mortality rates
	
	3darray N(1,nsex,syr,nyr+1,sage,nage);		//Numbers at age
	3darray M_tot(1,nsex,syr,nyr,sage,nage);	//Age-specific natural mortality
	3darray log_ft(1,nsex,1,ngear,syr,nyr);		//Gear specific log fishing mortlity rates
	3darray F(1,nsex,syr,nyr,sage,nage);		//Age-specific fishing mortality
	3darray Z(1,nsex,syr,nyr,sage,nage);
	3darray S(1,nsex,syr,nyr,sage,nage);
	matrix ct(1,ngear,syr,nyr);				//predicted catch biomass
	matrix epsilon(1,nit,1,nit_nobs);		//residuals for survey abundance index
	matrix pit(1,nit,1,nit_nobs);			//predicted relative abundance index
	matrix qt(1,nit,1,nit_nobs);			//catchability coefficients (time-varying)
	

	3darray Ahat(1,na_gears,1,na_nobs,a_sage-2,a_nage);		//predicted age proportions by gear & year
	3darray A_nu(1,na_gears,1,na_nobs,a_sage-2,a_nage);		//residuals for age proportions by gear & year
	
	4darray log_sel(1,nsex,1,ngear,syr,nyr,sage,nage);		//selectivity coefficients for each gear type.
	4darray Chat(1,nsex,1,ngear,syr,nyr,sage,nage);			//predicted catch-at-age
	
	sdreport_number sd_depletion;
	
	
PRELIMINARY_CALCS_SECTION
  //Run the model with input parameters to simulate real data.
  nf=0;
  if(SimFlag) 
  {
    initParameters();
    calcSelectivities();
    calcTotalMortality();
    //simulation_model(rseed);
  }

RUNTIME_SECTION
    maximum_function_evaluations 100,200,500,25000,25000
    convergence_criteria 0.01,0.01,1.e-5,1.e-5


PROCEDURE_SECTION
	initParameters();

	calcSelectivities();
	
	calcTotalMortality();
	
	calcNumbersAtAge();
	
	calcFisheryObservations();
	
	calcAgeProportions();
	
	calcSurveyObservations();
	
	calcStockRecruitment();
	
	calcObjectiveFunction();

	sd_depletion=sbt(nyr)/bo;
	
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
	
	Need to have a system for leading parameters to have sex based m.
	*/
	
	ro          = mfexp(theta(1));
	dvariable h = theta(2);
	m           = mfexp(theta(3));  //FIXME: for sex-based m
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


FUNCTION calcSelectivities
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
		
		CHANGED: add an option for length-based selectivity.  Use inverse of
		allometric relationship w_a = a*l_a^b; to get mean length-at-age from
		empirical weight-at-age data, then calculate selectivity based on 
		mean length. 
		
		TODO: Add sex based selectivity for 1,nsex dimension.
		
		FEB 15, 2012.  Adding sex based selectivity option.
		log_sel is now a 4darray(1,nsex,1,ngear,syr,nyr,sage,nage)
		
	
	*/
	int h,i,j,k;
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
	
	for(h=1;j<=nsex;h++)
	{
		for(j=1;j<=ngear;j++)
		{
			tmp.initialize(); tmp2.initialize();
			dvector iy(1,yr_nodes(j));
			dvector ia(1,age_nodes(j));
		
			switch(isel_type(j))
			{
				case 1:
					// logistic selectivity for case 1 or 6
					p1 = mfexp(sel_par(j,1,1));
					p2 = mfexp(sel_par(j,1,2));
					for(i=syr; i<=nyr; i++)
					{
						log_sel(h)(j)(i) = log( plogis(age,p1,p2)+tiny );
					}
					break;
			
				case 6:
					// logistic selectivity for case 1 or 6
					p1 = mfexp(sel_par(j,1,1));
					p2 = mfexp(sel_par(j,1,2));
					for(i=syr; i<=nyr; i++)
					{
						log_sel(h)(j)(i) = log( plogis(age,p1,p2) );
					}
					break;
				
				case 2:		
					// age-specific selectivity coefficients
					for(i=syr; i<=nyr; i++)
					{
						for(k=sage;k<=nage-1;k++)
						log_sel(h)(j)(i)(k) = sel_par(j)(1)(k-sage+1);
						log_sel(h)(j)(i,nage) = log_sel(h)(j)(i,nage-1);
					}
					break;
				
				case 3:		
					// cubic spline
					log_sel(h)(j)(syr)=cubic_spline( sel_par(j)(1) );
					for(i=syr; i<nyr; i++)
					{
						log_sel(h)(j)(i+1) = log_sel(h)(j)(i);
					}
					break;
				
				case 4:		
					// time-varying cubic spline every year
					for(i=syr; i<=nyr; i++)
					{
						log_sel(h)(j)(i) = cubic_spline(sel_par(j)(i-syr+1));
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
					log_sel(h)(j) = tmp2; 
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
						dvar_vector tmpwt=log(wt_obs(h)(i)*1000)/mean(log(wt_obs(h)*1000.));
						log_sel(h)(j)(i) = log( plogis(tmpwt,p1,p2)+tiny );
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
						tmp2(i) = p3*wt_dev(h)(i);
						log_sel(h)(j)(i) = log( plogis(age,p1,p2)+tiny ) + tmp2(i);
					}
					break;
				
				case 11:
					//logistic selectivity based on mean length-at-age
					p1 = mfexp(sel_par(j,1,1));
					p2 = mfexp(sel_par(j,1,2));
				
					for(i=syr; i<=nyr; i++)
					{
						dvector tmp = wt_obs(h)(i) / a(h);
						dvector len = pow( tmp,1./b(h) );
						log_sel(h)(j)(i) = log( plogis(len,p1,p2) );
					}
					break;
				
				case 12:
					//length-specific selectivity coefficients
					//based on cubic spline interpolation
					for(i=syr; i<=nyr; i++)
					{
						dvector tmp = wt_obs(h)(i) / a(h);
						dvector len = pow( tmp,1./b(h) );
						log_sel(h)(j)(i)=cubic_spline( sel_par(j)(1), len );
					}
					break;
				
				
				default:
					log_sel(h)(j)=0;
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
				log_sel(h)(j)(i) -= log( mean(mfexp(log_sel(h)(j)(i)))+tiny );
				//log_sel(j)(i) -= log(max(mfexp(log_sel(j)(i))));
			
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
	occur on the adult component.  Added if(catch_type(k)!=3) //exclude roe fisheries.
	
	Feb 15, 2012 Added sex based fishing mortality.
	The F, Z M_tot S arrays are now 3darray(1,nsex,syr,nyr,sage,nage)
	
	For now (Feb 15) assuming Fishign mort rate is the same for both sexs
	and differences in Fa is due to difference in selectivity.  THis must be
	fixed when there are sex-based selectivity parameters.  
	Should we also be estimating fishing mortality by sex?
	
	FIXME: watch for obs_ct by sex in this routine. May need to be modified
	if there is sex based fishing mortality rates.
	*/
	int h,j,k,ki;
	dvariable ftmp;
	F.initialize();
	ft.initialize();
	log_ft.initialize();
	
	//Fishing mortality
	for(h=1;h<=nsex;h++)
	{
		ki=1;
		for(k=1;k<=ngear;k++)
		{
		
			for(j=syr;j<=nyr;j++)
			{	
				ftmp=0;
				if( obs_ct(k,j)>0  )
					ftmp = mfexp(log_ft_pars(ki++));
			
				ft(k,j)=ftmp;
			
				if(catch_type(k)!=3){	//exclude roe fisheries
					F(h)(j)+=ftmp*mfexp(log_sel(h)(k)(j));
					//cout<<obs_ct(k,i)<<endl;
				}
			}
		}
		//Natural mortality (year and age specific)
		//M_tot(syr,nyr,sage,nage);
		M_tot(h) = m(h);
	}


	// Cubic spline to interpolate log_m_devs (log_m_nodes)
	/* Assume deviations in natural mortality rates are the 
	   same for both sexes at this time.
	*/
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
	if(active(log_m_nodes)&&i>syr)
	{
		for(h=1;h<=nsex;h++)
		{
			for(j=syr+1;j<=nyr;j++)
			{
				M_tot(h)(j)=M_tot(h)(j-1)*exp(log_m_devs(j));
			}
			m_bar(h) = mean(M_tot(h));
	
			Z(h)=M_tot(h)+F(h);
			S(h)=mfexp(-Z(h));
		}
	}
	if(verbose) cout<<"**** OK after calcTotalMortality ****"<<endl;
	
  }
	
	
FUNCTION calcNumbersAtAge
  {
	/*
		TODO Need to check the difference between the initialization 
		of the numbers at age here at the margins in comparison to the
		simulation model.
		
		Feb 15, 2012:  Added sex calculations N is now a 3darray
		for 3darry N(1,nsex,syr,nyr+1,sage,nage).
		
		ro is now the total number of sage-recruits (female + male)
	*/
	
	int h,i,j;
	N.initialize();
	
	
	if(cntrl(5)){	//If initializing in at unfished conditions
		log_rt(syr) = log(ro);
		for(h=1;h<=nsex;h++)
		{
			for(j=sage;j<=nage;j++)
			{
				N(h)(syr,j)=ro/nsex * exp(-m_bar(h)*(j-1.));
			}
			N(h)(syr,nage)/=(1.-exp(-m_bar(h)));
		}
	}
	else{			//If starting at fished conditions
		log_rt(syr) = log_avgrec+log_rec_devs(syr);
		for(h=1;h<=nsex;h++)
		{
			N(h)(syr,sage)=mfexp(log_rt(syr))/nsex;
			for(j=sage+1;j<=nage;j++)
			{
				dvariable tmp_rt = mfexp(log_recinit+init_log_rec_devs(j));
				N(h)(syr,j)      = tmp_rt/nsex * exp(-m_bar(h)*(j-sage));
			}
			N(h)(syr,nage)/=(1.-exp(-m_bar(h)));
		}
	}
	
	
	
	//initial number of sage recruits from year syr+1, nyr;
	for(h=1;h<=nsex;h++)
	{
		for(i=syr+1;i<=nyr;i++)
		{
			log_rt(i)=log_avgrec+log_rec_devs(i);
			N(h)(i,sage)=mfexp(log_rt(i))/nsex;
		}
		N(h)(nyr+1,sage)=mfexp(log_avgrec)/nsex;
	
		/* Dynamic state variables */
		for(i=syr;i<=nyr;i++)
		{
			N(h)(i+1)(sage+1,nage)=++elem_prod(N(h)(i)(sage,nage-1),S(h)(i)(sage,nage-1));
			N(h)(i+1,nage)+=N(h)(i,nage)*S(h)(i,nage);
		}
	}
	if(verbose)cout<<"**** Ok after calcNumbersAtAge ****"<<endl;
	
  }

FUNCTION calcAgeProportions
  {
	/*
	This function loops over each gear and year
	and calculates the predicted proportions at age
	sampled based on the selectivity of that gear and
	the numbers-at-age in the population.
	
	Feb 16, 2012 added nsex calculations 
	
	**
	The index for sex (ih) need to be added to the age-comps
	in the data file.  The following should be used:
	 - 0 for both sexes combined
	 - 1 for female age comps
	 - 2 for male age comps
	**
	*/
	
	int h,i,k,iyr,ig,ih;
	for(k=1;k<=na_gears;k++)
	{
		for(i=1;i<=na_nobs(k);i++)
		{
			iyr = A(k,i,a_sage(k)-2);	//index for year
			ig  = A(k,i,a_sage(k)-1);	//index for gear
			ih  = 0; //A(k,i,a_sage(k)-3); //index for sex once in data file.
			if(iyr>nyr)break;		//trap for retrospective analysis
			
			A_nu(k,i,a_sage(k)-2) = iyr;
			A_nu(k,i,a_sage(k)-1) = ig;
			Ahat(k,i,a_sage(k)-2) = iyr;
			Ahat(k,i,a_sage(k)-1) = ig;
			dvar_vector ctmp(a_sage(k),a_nage(k));
			ctmp.initialize();
			switch(ih)
			{
				case 0:	//add both sexes
					for(h=1;h<=nsex;h++)
					{
						ctmp += Chat(h)(k)(iyr)(a_sage(k),a_nage(k));
					}
				break;
				default:
					ctmp = Chat(ih)(k)(iyr)(a_sage(k),a_nage(k));
				break;
			}
			Ahat(k)(i)(a_sage(k),a_nage(k)) = ctmp / sum(ctmp);
			//Ahat(k)(i)(a_sage(k),a_nage(k))=Chat(k)(iyr)(a_sage(k),a_nage(k))
			//							/sum(Chat(k)(iyr)(a_sage(k),a_nage(k)));
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
	
	Feb 12, 2012	Added nsex calculations 
	*/
	
	/*
		FIXED Reconcile the difference between the predicted catch 
		here and in the simulation model.
	*/
	int h,i,k;
	ct.initialize();
	for(h=1;h<=nsex;h++)
	{
		for(i=syr;i<=nyr;i++)
		{
			for(k=1;k<=ngear;k++)
			{
				dvar_vector log_va=log_sel(h)(k)(i);
				
				//SJDM Jan 16, 2011 Modification as noted above.
				//SJDM Jan 06, 2012 Modification as noted above.
				if(obs_ct(k,i)>0)
				{	/*
					If there is a commercial fishery, then calculate the
					catch-at-age (in numbers) and total catch (in weight).
				
					Note that fa is nsex due to log_va, bt ft(k,i) is the same
					for both sexes
					*/
					dvar_vector fa = ft(k,i)*mfexp(log_va);
					dvar_vector d1 = elem_div(fa,Z(h)(i));
					Chat(h)(k,i)   = elem_prod(elem_prod(d1,1.-S(h)(i)),N(h)(i));
					switch(catch_type(k))
					{
						case 1:	//catch in weight
							ct(k,i) += Chat(h)(k,i)*wt_obs(h)(i);
						break;
						case 2:	//catch in numbers
							ct(k,i) += sum(Chat(h)(k,i));
						break;
						case 3:	//catch in roe that does not contribute to SSB
							dvariable ssb = elem_prod(N(h)(i),exp(-Z(h)(i)*cntrl(13)))*fec(h)(i);
							ct(k,i) += ( 1.-mfexp(-ft(k,i)) )*ssb;
						break;
					}
				}
				else
				{
					/*
					If there is no commercial fishery the set Chat equal to 
					the expected proportions at age.
					*/
					dvar_vector fa = mfexp(log_va);
					dvar_vector d1 = elem_div(fa,Z(h)(i));
					Chat(h)(k,i)=elem_prod(elem_prod(d1,1.-S(h)(i)),N(h)(i));
				}
			}
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
	
	Feb 12, 2012.  Added nsex calculations to this routine.
	FIXME: check roe fishery calculations when nsex>1
	*/
	/*
		CHANGED add capability to accomodate priors for survey q's.
		DONE
	*/
	
	int h,i,j,ii,k,kk;
	
	
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
			// Trap for retrospective analysis.
			if(ii>nyr) break;
			
			
			
			// Get fishing mortality rate on spawn.
			dvariable ftmp = 0;
			for(kk=1;kk<=ngear;kk++)
				if(catch_type(kk)==3)
					ftmp += ft(kk,ii);
			//
			
			for(h=1;h<=nsex;h++)
			{
				
				// Adjust survey biomass by the fraction of the mortality that 
				// occurred during the time of the survey.
				dvar_vector Np = elem_prod(N(h)(ii),exp( -Z(h)(ii)*it_timing(i,j) ));

				dvar_vector log_va=log_sel(h)(k)(ii);
				switch(survey_type(i))
				{
					case 1:
						V(j) += elem_prod(Np,mfexp(log_va));
					break;
					case 2:
						V(j) += elem_prod(elem_prod(Np,mfexp(log_va)),wt_obs(h)(ii));
					break;
					case 3:
						//SJDM Jan 6, 2012 Modified for roe fishery
						// If nsex > 1 the fec for males should equal 0.
						V(j) += elem_prod(Np,fec(h)(ii))*exp(-ftmp);
					break;
				}
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
	
	Jan 6, 2012.  Need to adjust stock-recruitment curvey for reductions 
	in fecundity associated with removal of roe from a spawn on kelp fishery.
	
	**
	Feb 16,2012.  Added nsex calculations to this routine.
	
	*/ 
	int h,i,j,k;
	
	// Process error variance
	dvariable tau = sqrt(1.-rho)/varphi;
	
	/*
	Need to calculate the parameters of the stock recruitment 
	function based on unfished conditions.  In the case where
	nsex > 1, base the stock recruitment relationship on mature
	female biomass only.
	
	In the case of time-varying natural mortality rate, the 
	average natural mortality rate is used to assume unfished
	conditions.
	*/
	dvar_vector tmp_rt(syr+sage,nyr);
	dvar_vector lx(sage,nage); 
	lx = 1;
	h  = 1;		// females only if nsex > 1
	for(j=sage+1;j<=nage;j++)
	{
		lx(j)=lx(j-1)*exp(-m_bar(h));	
	} 
	lx(nage)/=(1.-exp(-m_bar(h)));
	
	// Average fecundity of females (h=1) only.
	dvariable beta;
	
	// Female spawning biomass per recruit	
	dvariable phib = (lx*exp(-m_bar(h)*cntrl(13))) * avg_fec(h);
	dvariable so   = kappa/phib;		//max recruits per spawner
	bo             = ro*phib;  			//unfished female spawning biomass
	
	
	//SJDM Jan 6, 2012 Need to adjust sbt to reflect roe fishery 
	//in the sbt calculation below.
	for(i=syr;i<=nyr;i++)
	{
		sbt(i) = elem_prod(N(h)(i),exp(-Z(h)(i)*cntrl(13)))*fec(h)(i);
		
		//Adjustment to female spawning biomass for roe fisheries
		for(k=1;k<=ngear;k++)
			if(catch_type(k)==3)
				sbt(i) *= mfexp(-ft(k,i));
	}
	sbt(nyr+1) = N(h)(nyr+1)*fec(h)(nyr+1);
	//cout<<"sbt\n"<<sbt<<endl;
	//exit(1);
	
	dvar_vector tmp_st=sbt(syr,nyr-sage).shift(syr+sage);
	
	switch(int(cntrl(2)))
	{
		case 1:
			//Beverton-Holt model
			beta = (kappa-1.)/bo;
			tmp_rt=elem_div(so*tmp_st,1.+beta*tmp_st);
			break;
		case 2:
			//Ricker model
			//dvariable so=kappa/phib;
			beta=log(kappa)/bo;
			tmp_rt=elem_prod(so*tmp_st,exp(-beta*tmp_st));
		break;
	}
	
	//residuals in stock-recruitment curve
	rt=exp(log_rt(syr+sage,nyr));//trans(N)(1)(syr+1,nyr);
	delta = log(rt)-log(tmp_rt)+0.5*tau*tau;
	
	if(verbose)cout<<"**** Ok after calcStockRecruitment ****"<<endl;
	
  }
	
FUNCTION calcObjectiveFunction
  {
	/*
	Dec 20, 2010.  SJDM added prior to survey qs.
	
	q_prior is an ivector with current options of 0 & 1.
	0 is a uniform density (ignored) and 1 is a normal
	prior density applied to log(q).
	
	Feb 16, 2012.  Added nsex calculations for the penalties on
	log_sel.
	
	*/
	/*
	There are several components to the objective function
	Likelihoods:
		-1) likelihood of the catch data
		-2) likelihood of the survey abundance index
		-3) likelihood of the survey age comps
		-4) likelihood of the fishery age comps
		-5) likelihood for stock-recruitment relationship
		-6) likelihood for fishery selectivities
	*/
	int h,i,j,k;
	double o=1.e-10;

	dvar_vector lvec(1,6); lvec.initialize();
	nlvec.initialize();
	
	//1) likelihood of the catch data (retro)
	double sig_c =cntrl(3);
	if(last_phase())sig_c=cntrl(4);
	if(active(log_ft_pars))
	for(k=1;k<=ngear;k++){
		for(i=syr;i<=nyr;i++)
		{
			if(obs_ct(k,i)!=0)
				nlvec(1,k)+=dnorm(log(ct(k,i)),log(obs_ct(k,i)),sig_c);
		}
		//if(active(log_ft_pars))
		//	nlvec(1,k)=dnorm(log(obs_ct(k).sub(syr,nyr)+o)-log(ct(k).sub(syr,nyr)+o),sig_c);
	}
	
	
	//2) likelihood of the survey abundance index (retro)
	for(k=1;k<=nit;k++)
	{
		dvar_vector sig = (sqrt(rho)/varphi)/it_wt(k);
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
	dvariable tau = sqrt(1.-rho)/varphi;
	if(active(theta(1)))
		nlvec(4,1)=dnorm(delta,tau);
	
	
	
	//5-6) likelihood for selectivity paramters
	for(k=1;k<=ngear;k++)
	{
		if(active(sel_par(k))){
			//if not using logistic selectivity then
			//CHANGED from || to &&  May 18, 2011 Vivian
			if( isel_type(k)!=1 && 
				isel_type(k)!=7 && 
				isel_type(k)!=8 &&
				isel_type(k)!=11 )  
			{
				for(h=1;h<=nsex;h++)
				{
					for(i=syr;i<=nyr;i++)
					{
						//curvature in selectivity parameters
						dvar_vector df2=first_difference(first_difference(log_sel(h)(k)(i)));
						nlvec(5,k)+=sel_2nd_diff_wt(k)/(nage-sage+1)*norm2(df2);
				
						//penalty for dome-shapeness
						for(j=sage;j<=nage-1;j++)
							if(log_sel(h,k,i,j)>log_sel(h,k,i,j+1))
								nlvec(6,k)+=sel_dome_wt(k)
											*square(log_sel(h,k,i,j)-log_sel(h,k,i,j+1));
					}
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
					lvec(1)+=1000.*s*s;
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
					lvec(1)+=1000.*s*s;
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
		cout<<"nlvec\t"<<nlvec<<endl;
		cout<<"lvec\t"<<lvec<<endl;
		cout<<"priors\t"<<priors<<endl;
		cout<<"penalties\t"<<pvec<<endl;
	}
	f=sum(nlvec)+sum(lvec)+sum(priors)+sum(pvec)+sum(qvec);
	//cout<<f<<endl;
	nf++;
	if(verbose)cout<<"**** Ok after calcObjectiveFunction ****"<<endl;
	
  }

//FUNCTION void equilibrium(const double& fe, const double& ro, const double& kap, const dvector& m, const dvector& age, const dvector& wa, const dvector& fa, d3array& va, const dvector& allocation, double& re, double& ye, double& be, double& phiq, double& dphiq_df, double& dre_df)
//  {
//	/*
//	This is the equilibrium age-structured model that is 
//	conditioned on fe (the steady state fishing mortality rate).
//	
//	In the case of multiple fisheries, fe is to be considered as the
//	total fishing mortality rate and each fleet is given a specified
//	allocation based on its selectivity curve.  The allocation to 
//	each fleet must be specified a priori.
//	
//	args:
//	fe	-steady state fishing mortality
//	ro	-unfished sage recruits
//	kap	-recruitment compensation ration
//	m	-instantaneous natural mortality rate
//	age	-vector of ages
//	wa	-mean weight at age
//	fa	-mean fecundity at age
//	va	-mean vulnerablity at age for fe gear.
//	
//	Modified args:
//	re	-steady state recruitment
//	ye	-steady state yield
//	be	-steady state spawning biomass
//	phiq		-per recruit yield
//	dre_df		-partial of recruitment wrt fe
//	dphiq_df	-partial of per recruit yield wrt fe
//	
//	FIXME add Ricker model to reference points calculations.
//	FIXME partial derivatives for dphif_df need to be fixed when cntrl(13)>0.
//	
//	FEB 20, 2012.  Adding nsex calculations to this routine.
//	
//	- With more than 1 gear and an allocation to each of the gears, for a
//	  given fishing rate, need to figure out how to adjust fe_k to meet the
//	  allocation.  For now, just set fe_k = fe*allocation(k);
//	
//	*/
//	int h,i,j,k;
//	
//	int nage = max(age);
//	int sage = min(age);
//	int nk   = allocation.indexmax();//number of gears involved (or from allocation)
//	dvector fe_k(1,nk);
//	fe_k = fe*allocation;
//	dmatrix lx(1,nsex,sage,nage);
//	dmatrix lz(1,nsex,sage,nage);
//	dmatrix za(1,nsex,sage,nage);
//	dmatrix sa(1,nsex,sage,nage);
//	3darray qa(1,nsex,1,nk,sage,nage);
//	
//	for(h=1;h<=nsex;h++)
//	{
//		lx(h)      = pow(exp(-m(h)),age-double(sage));
//		lx(h,nage) = lx(h,nage)/(1.-exp(-m(h)));
//		
//		lz(h)      = lx(h);
//		za(h)      = m(h);
//		for(k=1;k<=nk;k++)
//		{
//			za(h) += fe_k(k)*va(h)(k);
//		}
//		sa(h)      = (1.-exp(-zh(h)));
//		for(k=1;k<=ng;k++)
//		{
//			qa(h)(k) = elem_prod(elem_div(va(h)(k),za(h)),sa(h));
//		}
//		
//	}
//	//dvector lx=pow(exp(-m),age-double(sage));
//	//lx(nage)/=(1.-exp(-m));
//	//dmatrix lz=lx;
//	//dvector za=m+fe*va;
//	//dvector sa=1.-exp(-za);
//	//dvector qa=elem_prod(elem_div(va,za),sa);
//	
//	double phie = lx*fa;		//eggs per recruit
//	double so = kap/phie;
//	double beta = (kap-1.)/(ro*phie);
//	
//	
//	double dlz_df = 0, dphif_df = 0;
//	dphiq_df=0; dre_df=0;
//	for(i=sage; i<=nage; i++)
//	{
//		if(i>sage) lz[i]=lz[i-1]*exp(-za[i-1]);
//		if(i>sage) dlz_df=dlz_df*exp(-za[i-1]) - lz[i-1]*va[i-1]*exp(-za[i-1]);
//		if(i==nage){ //6/11/2007 added plus group.
//					lz[i]/=(1.-mfexp(-za[i]));
//					
//					dlz_df=dlz_df/(1.-mfexp(-za[i]))
//							-lz[i-1]*mfexp(-za[i-1])*va[i]*mfexp(-za[i])
//					/((1.-mfexp(-za[i]))*(1.-mfexp(-za[i])));
//				}
//		dphif_df=dphif_df+fa(i)*dlz_df;
//		dphiq_df=dphiq_df+wa(i)*qa(i)*dlz_df+(lz(i)*wa(i)*va(i)*va(i))/za(i)*(exp(-za[i])-sa(i)/za(i));
//	}
//	//CHANGED need to account for fraction of mortality that occurs
//	//before the spawning season in the recruitment calculation.
//	//cout<<"lz\t"<<elem_prod(lz,exp(-za*cntrl(13)))<<endl;
//	//exit(1);
//	//double phif = lz*fa;
//	double phif = elem_prod(lz,exp(-za*cntrl(13)))*fa;
//	phiq=sum(elem_prod(elem_prod(lz,wa),qa));
//	re=ro*(kap-phie/phif)/(kap-1.);
//	//cout<<fe<<" spr ="<<phif/phie<<endl;
//	if(re<=0) re=0;
//	dre_df=(ro/(kap-1.))*phie/square(phif)*dphif_df;
//	ye=fe*re*phiq;
//	be=re*phif;	//spawning biomass
//	
//	//cout<<"Equilibrium\n"<<ro<<"\n"<<re<<"\n"<<ye<<endl;
//	
//  }
//	
//FUNCTION void calcReferencePoints()
//  {
//	/**
//	\file iscam.tpl
//	\author Steven Martell
//	Uses Newton_Raphson method to determine Fmsy and MSY 
//	based reference points.  
//	
//	Code check: appears to find the correct value of MSY
//	in terms of maximizing ye.  Check to ensure rec-devs
//	need a bias correction term to get this right.
//	
//	Modification for multiple fleets:
//	-	Need to pass a weighted average vector of selectivities
//		to the equilibrium routine, where the weights for each
//		selectivity is based on the allocation to each fleet.
//		
//	-	Perhaps as a default, assign an equal allocation to each
//		fleet.  Eventually,user must specify allocation in 
//		control file.
//		
//	-	Use selectivity in the terminal year to calculate reference
//		points.
//	
//	Feb 16, 2012.  Starting to add nsex calculations to this routine 
//	& its dependencies.  If nsex > 1, now have to calculate catch of
//	each sex.  This routine needs to be over-hauled to properly account
//	for allocations among gear types.
//	
//	PSEUDO CODE
//	1) set initial guess for fe = 1.5*mean(m)
//	2) call equilibrium routine passing (theta,fe,allocation,avg_log_sel(1,ngear))
//	   where theta is a vector of biological parameters. The equilibrium routine returns
//	   the total yield and derivatives of the catch equation.
//	
//	*/
//	int h,i,j;
//	double re,ye,be,phiq,dphiq_df,dre_df,fe;
//	double dye_df,ddye_df,spr;
//	fe = 1.5*value(mean(m_bar));
//	
//	
//	
//	
//	/* DEPRECATE THE CODE BELOW */
//	
//	/*Calculate average vulnerability*/
//	dvector va_bar(sage,nage);
//	va_bar.initialize();
//	/*CHANGED user now specifies allocation for MSY based reference points 
//	in the data file.  Used to be fsh_flag, but now is an allocation for gear k*/
//	//dvector allocation(1,ngear);
//	//allocation = dvector(fsh_flag/sum(fsh_flag));
//	
//	
//	/*CHANGED Allow for user to specify allocation among gear types.*/
//	/*FIXME:  this allocation should be on the catch on the vulnerabilities*/
//	for(j=1;j<=ngear;j++)
//	{
//		va_bar+=allocation(j)*value(exp(log_sel(j)(nyr)));
//		/*cout<<exp(log_sel(j)(nyr))<<endl;*/
//	}
//	
//	/*CHANGED Changed equilibrium calculations based on average m */
//	//FIXME: change Bmsy calculations to be based on average fecundity & weight at age. 
//	for(i=1;i<=20;i++)
//	{
//		//equilibrium(fe,value(ro),value(kappa),value(m),age,wa,fa,value(exp(log_sel(1)(nyr))),re,ye,be,phiq,dphiq_df,dre_df);
//		equilibrium(fe,value(ro),value(kappa),value(m_bar),age,wt_obs(nyr),
//					fec(nyr),va_bar,re,ye,be,phiq,dphiq_df,dre_df);
//		dye_df = re*phiq+fe*phiq*dre_df+fe*re*dphiq_df;
//		ddye_df = phiq*dre_df + re*dphiq_df;
//		fe = fe - dye_df/ddye_df;
//		if(verbose) cout<<"fe\t"<<fe<<"\t"<<dye_df<<"\t"<<ye<<endl;
//		if(sfabs(dye_df)<1.e-5)break;
//	}
//	fmsy=fe;
//	equilibrium(fmsy,value(ro),value(kappa),value(m_bar),age,wt_obs(nyr),
//				fec(nyr),va_bar,re,ye,be,phiq,dphiq_df,dre_df);
//	msy=ye;
//	bmsy=be;
//	
//	/* DEPRECATE THE CODE ABOVE */
//	
//	/*TODO print this to the REPORT file for plotting.*/
//	/*SM Loop over discrete value of fe and ensure above code is 
//	finding the correct value of msy.*/
//	
//	
//	ofstream report_file("iscam.eql");
//	
//	if(report_file.is_open())
//	{
//		report_file<<"index\t fe \t ye \t be \t re \t spr\n";
//		
//		fe = 0; i=0;
//		while(i < 1500)
//		{
//			equilibrium(fe,value(ro),value(kappa),value(m_bar),age,wt_obs(nyr),
//						fec(nyr),va_bar,re,ye,be,phiq,dphiq_df,dre_df);
//			if(re<=0)break;
//			
//			double spr = value(-ro/((kappa-1)*re-ro*kappa));
//			report_file<<i++<<"\t"<<fe<<"\t"<<ye<<"\t"<<be<<"\t";
//			report_file<<re<<"\t"<<spr<<endl;
//			
//			fe += 0.01;
//		}
//	}
//	//exit(1);
//	
//	if(verbose)cout<<"**** Ok after calcReferencePoints ****"<<endl;
//  }
//	
//FUNCTION void simulation_model(const long& seed)
//  {
//	/*
//	Call this routine to simulate data for simulation testing.
//	The random number seed can be used to repeat the same 
//	sequence of random number for simulation testing.
//	
//	Implemented using the "-SimFlag 99" command line option where
//	99 is the random number seed.
//	
//	-SimFlag 99 is a special case used for the manuscript for case 1.
//	-SimFlag 000 is a special case with 0 error (exact data)
//	
//	-This routine will over-write the observations in memory
//	with simulated data, where the true parameter values are
//	the initial values.  Change the standard deviations of the 
//	random number vectors epsilon (observation error) or 
//	recruitment devs wt (process error).
//	*/
//	
//	
//	cout<<"___________________________________________________\n"<<endl;
//	cout<<"  **Implementing Simulation--Estimation trial**    "<<endl;
//	cout<<"___________________________________________________"<<endl;
//	cout<<"\tRandom Seed No.:\t"<< rseed<<endl;
//	cout<<"___________________________________________________\n"<<endl;
//	
//	
//	//Indexes:
//	int i,j,k,ii,ki;
//
//	//3darray Chat(1,ngear,syr,nyr,sage,nage);
//	//C.initialize();
//	
//	
//	
//	/*----------------------------------*/
//	/*	-- Generate random numbers --	*/
//	/*----------------------------------*/
//	random_number_generator rng(seed);
//	dvector wt(syr-nage-1,nyr);			//recruitment anomalies
//	dmatrix epsilon(1,nit,1,nit_nobs);  //observation errors in survey
//	
//	double sig = value(sqrt(rho)/varphi);
//	double tau = value(sqrt(1.-rho)/varphi);
//	
//	if(seed==000)
//	{
//		cout<<"No Error\n";
//		sig=0;
//		tau=0;
//	}
//	wt.fill_randn(rng); wt *= tau;
//	epsilon.fill_randn(rng); 
//	
//	//now loop over surveys and scale the observation errors
//	for(k=1;k<=nit;k++)
//	{
//		for(j=1;j<=nit_nobs(k);j++)
//			epsilon(k,j) *= sig/it_wt(k,j);
//	}
//	
//	cout<<"	OK after random numbers\n";
//	/*----------------------------------*/
//	
//	
//	
//	
//	
//	/*----------------------------------*/
//    /*		--Initialize model--		*/
//	/*CHANGED now calculating phie based on m_bar and avg_fec*/
//	/*----------------------------------*/
//	dvector lx=pow(exp(-value(m_bar)),age-min(age));
//	lx(nage)/=(1.-exp(-value(m_bar)));
//	double phie=(lx*exp(-value(m_bar)*cntrl(13)))*avg_fec;//fec(syr);
//	so=kappa/phie;
//	
//	
//	if(cntrl(2)==1) beta=(kappa-1.)/(ro*phie);
//	if(cntrl(2)==2) beta=log(kappa)/(ro*phie);
//	
//	//Initial numbers-at-age with recruitment devs
//	/*for(i=syr;i < syr+sage;i++)
//			N(i,sage)=exp(log_avgrec+wt(i));
//			
//		for(j=sage+1;j<=nage;j++)
//			N(syr,j)=exp(log_avgrec+wt(syr-j))*lx(j);
//		*/
//	
//	N.initialize();
//	if(cntrl(5)){	//If initializing in at unfished conditions
//		log_rt(syr) = log(ro);
//		for(h=1;h<=nsex;h++)
//		{
//			for(j=sage;j<=nage;j++)
//			{
//				N(h)(syr,j)=ro/nsex*exp(-m_bar*(j-1.));
//			}
//		}
//	}
//	else{			//If starting at unfished conditions
//		log_rt(syr) = log_avgrec;
//		N(syr,sage)=mfexp(log_rt(syr));
//		for(j=sage+1;j<=nage;j++)
//		{
//			N(syr,j)=mfexp(log_recinit+init_log_rec_devs(j))*exp(-m_bar*(j-sage));
//		}
//	}
//	N(syr,nage)/=(1.-exp(-m_bar));
//	
//	//log_rt=log_avgrec+log_rec_devs;
//	//log_rt(syr) = log(ro);
//
//	for(i=syr+1;i<=nyr;i++){
//		log_rt(i)=log_avgrec+log_rec_devs(i);
//		N(i,sage)=mfexp(log_rt(i));
//	}
//	N(nyr+1,sage)=mfexp(log_avgrec);
//
//	/*
//	for(j=sage;j<=nage;j++) 
//	{
//		if(cntrl(5))  //if starting at unfished state
//		{
//			N(syr,j)=ro*exp(-m_bar*(j-1));
//		}
//		else{
//			log_rt(syr-j+sage)=log_avgrec+log_rec_devs(syr-j+sage);
//			N(syr,j)=mfexp(log_rt(syr-j+sage))*exp(-m_bar*(j-sage));
//		}
//	}
//
//	N(syr,nage)/=(1.-exp(-m_bar));
//	*/
//	cout<<"	Ok after initialize model\n";
//	/*----------------------------------*/
//	
//	
//	
//	/*----------------------------------*/
//    /*		--    Selectivity   --		*/
//	/*----------------------------------*/
//	/*
//		-Based on values in the control file.
//		-Add autocorrelated random numbers
//		for time varying or something to that
//		effect.
//		
//		-If seed==99 then set up a special case
//		for the cubic spline manunscript using
//		the eplogistic function where g goes from
//		strongly domed to asymptotic, e.g.,
//		g = 0.2 * (nyr-i)/(nyr-syr);
//		
//	*/
//
//	/*CHANGED May 15, 2011 calcSelectivities gets called from PRELIMINARY_CALCS*/
//	dmatrix va(1,ngear,sage,nage);			//fishery selectivity
//	d3_array dlog_sel(1,ngear,syr,nyr,sage,nage);
//	dlog_sel=value(log_sel);
//	/*
//	for(k=1;k<=ngear;k++)
//		for(i=syr;i<=nyr;i++)
//		{
//			//sel(k)(i)=plogis(age,ahat(k),ghat(k));
//			log_sel(k)(i)=log(plogis(age,ahat(k),ghat(k)));
//			log_sel(k)(i) -= log(mean(exp(log_sel(k)(i))));
//		}
//		//log_sel(j)(i) -= log(mean(mfexp(log_sel(j)(i))));
//	*/
//	cout<<"	Ok after selectivity\n";
//
//	/*----------------------------------*/
//	
//	
//	/*----------------------------------*/
//    /*	--  Population dynamics  --		*/
//	/*----------------------------------*/
//	
//	dmatrix zt(syr,nyr,sage,nage);			//total mortality
//	zt.initialize();
//	dmatrix ft(syr,nyr,1,ngear);
//	ft.initialize();
//	dvector sbt(syr,nyr+1);
//	sbt.initialize();
//	
//	
//	for(i=syr;i<=nyr;i++)
//	{   
//		
//		//total biomass at age
//		//dvector bt = elem_prod(value(N(i)),wa);
//		dvector bt = elem_prod(value(N(i)),wt_obs(i));
//
//		/*calculate instantaneous fishing mortalities
//		based on Baranov's catch equation and the 
//		observed catch from each fleet.*/
//		dvector oct = trans(obs_ct)(i);
//		
//		
//		for(k=1;k<=ngear;k++)
//			va(k)=exp(dlog_sel(k)(i));
//		
//		//get_ft is defined in the Baranov.cxx file
//		//CHANGED these ft are based on biomass at age, should be numbers at age
//		//ft(i) = get_ft(oct,value(m),va,bt);
//		ft(i) = get_ft(oct,value(m),va,value(N(i)),wt_obs(i));
//		//cout<<trans(obs_ct)(i)<<"\t"<<oct<<endl;
//		
//		// overwrite observed catch incase it was modified by get_ft
//		for(k=1;k<=ngear;k++)
//			obs_ct(k,i)=oct(k);
//		
//		//total age-specific mortality
//		//dvector zt(sage,nage);
//		zt(i)=value(m);
//		for(k=1;k<=ngear;k++){
//			zt(i)+= ft(i,k)*exp(dlog_sel(k)(i));
//		}
//		
//		
//		//CHANGED definition of spawning biomass based on ctrl(12)
//		sbt(i) = value(elem_prod(N(i),exp(-zt(i)*cntrl(13)))*fec(i));
//		
//		//Update numbers at age
//		if(i>=syr+sage-1)
//		{
//			double rt;
//			//double et=value(N(i-sage+1))*fec(i-sage+1);
//			double et=sbt(i-sage+1);
//			if(cntrl(2)==1)rt=value(so*et/(1.+beta*et));
//			if(cntrl(2)==2)rt=value(so*et*exp(-beta*et));
//			N(i+1,sage)=rt*exp(wt(i)-0.5*tau*tau);
//			
//			/*CHANGED The recruitment calculation above is incosistent
//			  with the assessment model.  Below recruitment is based on
//			  rt=exp(log_avgrec + wt + rt_dev), where the rt_dev calculation
//			is based on the BH or Ricker model.*/
//			//double rt_dev = log(rt)-value(log_avgrec);
//			//N(i+1,sage)=exp(log_avgrec+wt(i));
//			
//		}
//
//
//		N(i+1)(sage+1,nage)=++elem_prod(N(i)(sage,nage-1),exp(-zt(i)(sage,nage-1)));
//		N(i+1,nage)+=N(i,nage)*exp(-zt(i,nage));
//		
//		
//		//Catch & Catch-at-age
//		for(k=1;k<=ngear;k++)
//		{
//			if(ft(i,k)>0)
//			{
//				dvector sel = exp(dlog_sel(k)(i));
//				d3C(k)(i)=elem_prod(elem_div(ft(i,k)*sel,zt(i)),elem_prod(1.-exp(-zt(i)),value(N(i))));
//				obs_ct(k,i)=d3C(k)(i)*wt_obs(i);
//			}
//			else	//if this is a survey
//			{
//				dvector sel = exp(dlog_sel(k)(i));
//				d3C(k)(i)=elem_prod(elem_div(sel,zt(i)),elem_prod(1.-exp(-zt(i)),value(N(i))));
//			}
//		}
//		
//	}
//	
//	
//	//initial values of log_ft_pars set to true values
//	ki=1;
//	for(k=1;k<=ngear;k++)
//		for(i=syr;i<=nyr;i++)
//			if(obs_ct(k,i)>0){
//				log_ft_pars(ki++)=log(ft(i,k));
//			}
//	
//	// Error handler to inform user population went extinct.
//	if(min(sbt(syr,nyr))<=1.e-5)
//	{
//		cout<<"---------------------------------------\n";
//		cout<<"Simulated population went extinct, try\n";
//		cout<<"increasing steepness, Ro and Rbar\n";
//		cout<<sbt<<endl;
//		cout<<"Minimum spawning biomass="<<min(sbt(syr,nyr))<<endl;
//		cout<<"---------------------------------------\n";
//		exit(1);
//	}
//	
//	//Average recruitment calculation
//	
//	
//	cout<<"	log(mean(column(N,sage))) = "<<mean(log(column(N,sage)))<<endl;
//	cout<<"	log_avgrec = "<<log_avgrec<<endl;
//	cout<<"	Ok after population dynamics\n";
//	/*----------------------------------*/
//	
//	/*----------------------------------*/
//    /*	--  Observation models  --		*/
//	/*----------------------------------*/
//	
//	//Simulated Age-compositions
//	int ig;
//	for(k=1;k<=na_gears;k++)
//	{
//		for(i=1;i<=na_nobs(k);i++)
//		{
//			ii=A(k,i,a_sage(k)-2);	//index for year
//			ig=A(k,i,a_sage(k)-1);	//index for gear
//			dvector pa = d3C(ig)(ii);	//
//			pa/=sum(pa);
//			
//			
//			dvector t1=pa(a_sage(k),a_nage(k));
//			t1/=sum(t1);
//			A(k)(i)(a_sage(k),a_nage(k))=rmvlogistic(t1,0.3,i+seed);
//			if(seed==000)
//			{
//				A(k)(i)(a_sage(k),a_nage(k))=t1;
//			}
//			//cout<<iyr<<"\t"<<k<<endl;
//		}
//	}
//
//	//cout<<Ahat<<endl;
//	
//	//Relative abundance indices
//	//CHANGED fixed this to reflect survey timing etc & survey_type
//	for(k=1;k<=nit;k++)
//	{   
//		for(i=1;i<=nit_nobs(k);i++)
//		{
//			ii=iyr(k,i);
//			ig=igr(k,i);
//			dvector sel = exp(dlog_sel(ig)(ii));
//			dvector Np = value(elem_prod(N(ii),exp(-zt(ii)*it_timing(k,i))));
//			switch(survey_type(k))
//			{
//				case 1: //survey based on numbers
//					Np = elem_prod(Np,sel);
//				break;
//				case 2: //survey based on biomass
//					Np = elem_prod(elem_prod(Np,sel),wt_obs(ii));
//				break;
//				case 3: //survey based on spawning biomass
//					Np = elem_prod(Np,fec(ii));
//				break;
//			}
//			it(k,i) = sum(Np) * exp(epsilon(k,i));
//		}
//	}
//	
//	
//
//	cout<<"	OK after observation models\n";
//	/*----------------------------------*/
//	
//	//CHANGED Fixed bug in reference points calc call from simulation model,
//	//had to calculate m_bar before running this routine.
//	
//	calcReferencePoints();
//	//cout<<"	OK after reference points\n"<<fmsy<<endl;
//	//exit(1);
//	//	REPORT(fmsy);
//	//	REPORT(msy);
//	//	REPORT(bmsy);
//	
//	
//	cout<<"___________________________________________________"<<endl;
//	ofstream ofs("iscam.sim");
//	ofs<<"fmsy\n"<<fmsy<<endl;
//	ofs<<"msy\n"<<msy<<endl;
//	ofs<<"bmsy\n"<<bmsy<<endl;
//	ofs<<"va\n"<<va<<endl;
//	ofs<<"sbt\n"<<sbt<<endl;//<<rowsum(elem_prod(N,fec))<<endl;
//	ofs<<"rt\n"<<rt<<endl;
//	ofs<<"ct\n"<<obs_ct<<endl;
//	ofs<<"ft\n"<<trans(ft)<<endl;
//	ofs<<"ut\n"<<elem_div(colsum(obs_ct),N.sub(syr,nyr)*wa)<<endl;
//	ofs<<"iyr\n"<<iyr<<endl;
//	ofs<<"it\n"<<it<<endl;
//	ofs<<"N\n"<<N<<endl;
//	ofs<<"A\n"<<A<<endl;
//	ofs<<"dlog_sel\n"<<dlog_sel<<endl;
//	cout<<"  -- Simuation results written to iscam.sim --\n";
//	cout<<"___________________________________________________"<<endl;
//	
//	//cout<<N<<endl;
//	//exit(1);
//  }
//	
//FUNCTION dvector cis(const dvector& na)
//  {
//	//Cohort Influenced Selectivity
//	//This function returns a vector of residuals from a
//	//linear regression of log(pa)= a+b*age+res that can be
//	//used to modify age-based selectivity according to relative
//	//cohort strengths.
//	
//	//SM  Currently not used at all in iscam and should be deprecated.
//	dvector y = log(na);
//	dvector x = age;
//	double b = sum(elem_prod(x-mean(x),y-mean(y)))/sum(square(x-mean(x)));
//	double a = mean(y)-b*mean(x);
//	dvector res = y - (a+b*x);
//	return(res);
//  }

REPORT_SECTION
  {
	if(verbose)cout<<"Start of Report Section..."<<endl;
	int i,j,k;
	REPORT(ControlFile);
	REPORT(f);
	REPORT(nlvec);
	REPORT(ro);
	double rbar=value(exp(log_avgrec));
	REPORT(rbar);
	double rinit=value(exp(log_recinit));
	REPORT(rinit);
	REPORT(bo);
	REPORT(kappa);
	double steepness=value(theta(2));
	REPORT(steepness);
	REPORT(m);
	double tau = value(sqrt(1.-rho)/varphi);
	double sig = value(sqrt(rho)/varphi);
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
	REPORT(ft);
	/*FIXED small problem here with array bounds if using -retro option*/
	//report<<"ut\n"<<elem_div(colsum(obs_ct)(syr,nyr),N.sub(syr,nyr)*wa)<<endl;
	//report<<"bt\n"<<rowsum(elem_prod(N,wt_obs))<<endl;
	//report<<"sbt\n"<<sbt<<endl;

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

//	if(last_phase())
//	{	calcReferencePoints();
//		REPORT(fmsy);
//		REPORT(msy);
//		REPORT(bmsy);
//	}
		
	//Parameter controls
	dmatrix ctrl=theta_control;
	REPORT(ctrl);
	
	
	//if(last_phase()) projection_model(0);
	
	dvector rt3(1,3);
//	if(last_phase())
//	{
//		dvector rt3 = age3_recruitment(value(column(N,3)),wt_obs(nyr+1,3),value(M_tot(nyr,3)));
//		REPORT(rt3);
//	}
	
//	dvector future_bt = value(elem_prod(elem_prod(N(nyr+1),exp(-M_tot(nyr))),wt_obs(nyr+1)));
//	REPORT(future_bt);
//	double future_bt4 = sum(future_bt(4,nage));
//	REPORT(future_bt4);
//	

	
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
	
	
  }
	
FUNCTION mcmc_output
  {
	if(nf==1){
		adstring str_q;
		str_q="lnq";
		ofstream ofs("iscam.mcmc");
		ofs<<"log.ro\t h\t log.m\t log.rbar\t log.rinit\t rho\t vartheta\t";
		ofs<<"bo\t bmsy\t msy\t fmsy\t";
		ofs<<"SSB\t Age-4\t Poor\t Average\t Good\t";
		for(int i=1;i<=nit;i++)ofs<<str_q<<i<<"\t";
		ofs<<"f\t"<<endl;
		
		ofstream of1("sbt.mcmc");
		ofstream of2("rt.mcmc");
		
	}
	
	// leading parameters & reference points
//	calcReferencePoints();
	// decision table output
//	dvector future_bt = value(elem_prod(elem_prod(N(nyr+1),exp(-M_tot(nyr))),wt_obs(nyr+1)));
//	double future_bt4 = sum(future_bt(4,nage));
//	dvector rt3 = age3_recruitment(value(column(N,3)),wt_obs(nyr+1,3),value(M_tot(nyr,3)));	
//	ofstream ofs("iscam.mcmc",ios::app);
//	ofs<<theta;
//	ofs<<" "<<bo<<" "<<bmsy<<" "<<msy<<" "<<fmsy<<"\t\t";
//	ofs<<sbt(nyr)<<" "<<future_bt4<<" "<<future_bt4+rt3<<"\t\t";
//	ofs<<log(q)<<" "<<f<<endl;
	
	// output spawning stock biomass
	ofstream of1("sbt.mcmc",ios::app);
	of1<<sbt(syr,nyr)<<endl;
	
	// output age-1 recruits
	ofstream of2("rt.mcmc",ios::app);
	of2<<rt<<endl;
	
	// Projection model.
//	for(int i=0;i<=10;i++)
//	{
//		double tac = double(i)/10. * 1.5*msy;
//		projection_model(tac);
//	}
	
	// Deviance Information Criterion
	/*
		DIC = pd + Dbar, where:
		pd is the effective number of parameters,
		Dbar is the expected (average deviance)
		
		Dbar = mean(D(theta))
		pd = Dbar - D(mean(theta))
		Agorithm:
			-for each sample compute Dtotal += 2*f;
			-compute a running Dbar = Dtotal/nf;
			
			
	
	//cout<<initial_params::nvarcalc()<<endl;
	int nvar = initial_params::nvarcalc();
	dvector y(1,nvar);
	static dvector sum_y(1,nvar);
	initial_params::xinit(y);
	sum_y = sum_y + y;
	cout<<y(1,3)<<endl;
	cout<<sum_y(1,3)<<endl;
	get_monte_carlo_value(nvar,y);
	if(nf==2)exit(1);
	*/
	
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

//FUNCTION void projection_model(const double& tac);
//  {
//	/*
//	This routine conducts population projections based on 
//	the estimated values of theta.  Note that all variables
//	in this routine are data type variables.
//	
//	Arguments:
//	tac is the total allowable catch that must be allocated 
//	to each gear type based on allocation(k)
//	
//	theta(1) = log_ro
//	theta(2) = h
//	theta(3) = log_m
//	theta(4) = log_avgrec
//	theta(5) = log_recinit
//	theta(6) = rho
//	theta(7) = vartheta
//	
//	*/
//	
//	int i,j,k;
//	int pyr = nyr+2;	//projection year.
//	
//	// --derive stock recruitment parameters
//	dvector lx(sage,nage); lx=1;
//	for(i=sage+1;i<=nage;i++) lx(i)=lx(i-1)*exp(-value(m_bar));
//	lx(nage)/=(1.-exp(-value(m_bar)));
//	
//	double phib = lx*fec(nyr);//(lx*exp(-value(m_bar))) * fec(nyr);//avg_fec;  
//	double so = value(kappa)/phib;		//max recruits per spawner
//	double beta;
//	double bo = value(ro)*phib;  				//unfished spawning biomass	
//	switch(int(cntrl(2)))
//	{
//		case 1:
//			beta = (value(kappa)-1.)/bo;
//		break;
//		case 2:
//			beta = log(value(kappa))/bo;
//		break;
//	}
//	
//	
//	dvector p_sbt(syr,pyr);
//	dvector p_ct(1,ngear);
//	dmatrix p_ft(nyr+1,pyr);
//	dmatrix p_N(syr,pyr+1,sage,nage);
//	dmatrix p_Z(syr,pyr,sage,nage);
//	p_N.initialize();
//	p_N.sub(syr,nyr+1) = value(N.sub(syr,nyr+1));
//	p_sbt(syr,nyr)=value(sbt(syr,nyr));
//	p_Z.sub(syr,nyr) = value(Z.sub(syr,nyr));
//	
//	//selecticity
//	/*CHANGED User to specifies allocation among gear types in data file.*/
//	dmatrix va_bar(1,ngear,sage,nage);
//	for(k=1;k<=ngear;k++)
//	{
//		//va_bar+=allocation(k)*value(exp(log_sel(k)(nyr)));
//		/*cout<<exp(log_sel(j)(nyr))<<endl;*/
//		p_ct(k)   = allocation(k)*tac;
//		va_bar(k) = exp(value(log_sel(k)(nyr)));
//	}
//		
//	
//	for(i = nyr+1; i<=pyr; i++)
//	{
//		
//		//get_ft is defined in the Baranov.cxx file
//		//(wt_obs(nyr+1) is the average wt at age in the last 5 years)
//		p_ft(i) = get_ft(p_ct,value(m_bar),va_bar,p_N(i),wt_obs(nyr+1));
//		
//		//Calculate mortality
//		p_Z(i) = value(m_bar);
//		for(k=1;k<=ngear;k++)
//		{
//			p_Z(i)+=p_ft(i,k)*va_bar(k);
//		}
//		
//		
//		//Spawning biomass
//		p_sbt(i) = elem_prod(p_N(i),exp(-p_Z(i)*cntrl(13)))*fec(nyr);//avg_fec;
//		
//		//Age-sage recruits
//		double tau = value(sqrt(1.-rho)/varphi); 
//		double xx = randn(int(tac)+i)*tau;
//		
//		if(i>=syr+sage-1)
//		{
//			double rt;
//			double et=p_sbt(i-sage+1);
//			if(cntrl(2)==1)rt=(so*et/(1.+beta*et));
//			if(cntrl(2)==2)rt=(so*et*exp(-beta*et));
//			p_N(i+1,sage)=rt*exp(xx-0.5*tau*tau); 
//		}
//		
//		//Update numbers at age
//		p_N(i+1)(sage+1,nage)=++elem_prod(p_N(i)(sage,nage-1),exp(-p_Z(i)(sage,nage-1)));
//		p_N(i+1,nage)+=p_N(i,nage)*exp(-p_Z(i,nage));
//		
//	}
//	if(nf==1)
//	{
//		ofstream ofs(BaseFileName + ".proj");
//		ofs<<"TAC\t PSC\t PSS\t Ut"<<endl;
//	}
//	ofstream ofs(BaseFileName + ".proj",ios::app);
//	ofs<<tac<<"\t"
//	<<0.25*bo/p_sbt(pyr)<<"\t"
//	<<p_sbt(pyr-1)/p_sbt(pyr)<<"\t"
//	<<(tac/(p_N(pyr-1)(3,nage)*wt_obs(nyr+1)(3,nage)))<<endl;
//	
//  }

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

	#if defined(WIN32) && !defined(__linux__)
		const char* PLATFORM = "Windows";
	#else
		const char* PLATFORM = "Linux";
	#endif

	#include <admodel.h>
	#include <time.h>
	#include <string.h>
	#include "stats.cxx"
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
		
			cout<<"Copied MCMC Files"<<endl;
		}
	}

	if( last_phase() && PLATFORM =="Linux" && retro_yrs )
	{
		//copy report file with .ret# extension for retrospective analysis
		adstring bscmd = "cp iscam.rep " + BaseFileName + ".ret" + str(retro_yrs);
		system(bscmd);
	}


