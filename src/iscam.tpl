// ------------------------------------------------------------------- //
//         integrated Statistical Catch Age Model (iSCAM)              //
//                                                                     //
//                           VERSION 1.0                               //
//                            2010-11-27                               //
//                                                                     //
//                                                                     //
//           Created by Steven Martell on 2010-04-09                   //
//           Copyright (c) 2010. All rights reserved.                  //
//                                                                     //
// AUTHORS: SJDM Steven Martell                                        //
//                                                                     //
// CONVENTIONS: Formatting conventions are based on the The            //
//               Elements of C++ Style (Misfeldt et al. 2004)          //
//                                                                     //
// NAMING CONVENTIONS:                                                 //
//             Macros       -> UPPERCASE                               //
//             Constants    -> UpperCamelCase                          //
//             Functions    -> lowerCamelCase                          //
//             Variables    -> lowercase                               //
//                                                                     //
// CHANGED add option for using empirical weight-at-age data           //
// TODO:    add gtg options for length based fisheries                //
// CHANGED add time varying natural mortality rate with splines        //
// TODO:    add cubic spline interpolation for time varying M         //
// CHANGED  Fix the type 6 selectivity implementation. not working.    //
// TODO:  fix cubic spline selectivity for only years when data avail  //
// CHANGED: fixed a bug in the simulation model log_ft_pars goes out   //
//        of bounds.                                                   //
//                                                                     //
//                                                                     //
//                                                                     //
// ------------------------------------------------------------------- //
//-- CHANGE LOG:                                                     --//
//--  Nov 30, 2010 -modified survey biomass by the fraction of total --//
//--                mortality that occurred during the time of the   --//
//--                survey. User specifies this fraction (0-1) in the--//
//--                data file as the last column of the relative     --//
//--                abundance index.                                 --//
//--                                                                 --//
//--  Dec 6, 2010 -modified the code to allow for empiracle weight-  --//
//--               at-age data to be used.                           --//
//--              -rescaled catch and relative abundance /1000, this --//
//--               should be done in the data file and not here.     --//
//--                                                                 --//
//--  Dec 20, 2010-added prior to survey q's in control file         --//
//--                                                                 --//
//--  Dec 24, 2010-added random walk for natural mortality.          --//
//--                                                                 --//
//--  Jan 23, 2011-in Penticton Hospital with my mom in ICU, adopting--//
//--               the naming conventions in The Elements of C++     --//
//--               style to keep my mind busy.                       --//
//--                                                                 --//
//-- May 5, 2011- added logistic selectcitivty as a fucntion of      --//
//--              mean body weight.  3 parameter logistic.           --//
//--              NOT WORKING YET                                    --//
//--                                                                 --//
//-- May 6, 2011- added pre-processor commands to determin PLATFORM  --//
//--              either "Windows" or "Linux"                        --//
//--                                                                 --//
//--                                                                 --//
//--                                                                 --//
// ------------------------------------------------------------------- //


DATA_SECTION
	!! cout<<"iSCAM has detected that you are on a "<<PLATFORM<<" box"<<endl;
	

	init_adstring DataFile;
	init_adstring ControlFile;
	
	
	!! BaseFileName=stripExtension(DataFile);
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
	vector age(sage,nage);
	!! age.fill_seqadd(sage,1);
	
	init_int ngear;				//number of gear types with unique selectivities
	!! cout<<"ngear\t"<<ngear<<endl;
	init_ivector fsh_flag(1,ngear);
	ivector ft_phz(1,ngear);
	LOC_CALCS
		int k;
		ft_phz=1;
		for(k=1;k<=ngear;k++)
		{
			if(!fsh_flag(k))
				ft_phz(k)=-1;
		}
	END_CALCS
	
	init_number fixed_m;
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
	  cout<<"linf\t"<<linf<<endl;
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
	LOC_CALCS
		ft_count=0;
		for(k=1;k<=ngear;k++)
			obs_ct(k)=column(catch_data,k+1);
		int i;
		for(k=1;k<=ngear;k++)
		{	
			for(i=syr;i<=nyr;i++)
				if(obs_ct(k,i)>0) ft_count++;
		}
		cout<<"ft_count\n"<<ft_count<<endl;
		cout<<"last row of catch \n"<<catch_data(nyr)<<endl;
		cout<<"Ok after catch extraction"<<endl;
	END_CALCS
	
	init_int nit;
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
		cout<<"OK after relative abundance index"<<endl;
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
	vector avg_fec(sage,nage);				//average fecundity-at-age
	LOC_CALCS
		int iyr;
		avg_fec.initialize();
		for(i=syr;i<=nyr+1;i++)
		{
			wt_obs(i)=wa;			
			fec(i)=elem_prod(plogis(age,ah,gh),wt_obs(i));
		}
		//if empiracle weight-at-age data exist, the overwrite wt_obs & fec.
		for(i=1;i<=n_wt_nobs;i++)
		{
			iyr=tmp_wt_obs(i,sage-1);  //index for year
			wt_obs(iyr)=tmp_wt_obs(i)(sage,nage);
			fec(iyr)=elem_prod(plogis(age,ah,gh),wt_obs(iyr));
		}
		//CHANGED average fecundity
		int nfec = fec.rowmax()-fec.rowmin()+1;
		avg_fec=colsum(fec)/nfec;
		
		
		//from Jake Schweigert: use mean-weight-at-age data
		//from the last 5 years for the projected mean wt.
		dvector tmp=colsum(wt_obs.sub(nyr-5,nyr))/6.;
		wt_obs(nyr+1) = tmp;
		cout<<"n_wt_nobs\t"<<n_wt_nobs<<endl;
		cout<<"Ok after empiracle weight-at-age data"<<endl;
		
		//May 5, 2011 SJDM: Calculating standardized deviates
		//in mean weights-at-age from empiracle data to use
		//with selectivity as a function of mean weight at age.
		//Idea borrowed from Vivian Haist in HCAM model.
		/*
			FIXME: selectivity function based on average weight.
			Based on the logistic function;
			1) calculate a matrix of standardized deviates
			wt_dev = (wt_obs-mean(wt_obs))/sd(wt_obs);
		*/
		wt_dev.initialize();
		dmatrix mtmp = trans(wt_obs);
		for(i=sage;i<=nage;i++)
		{
			mtmp(i) = (mtmp(i)-mean(mtmp(i)))/sqrt(var(mtmp(i)));
		}
		wt_dev = trans(mtmp);
		//cout<<colsum(wt_dev)<<endl;
		//exit(1);
		
		
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
	3darray d3C(1,ngear,syr,nyr,sage,nage);		//catch-age for simulation model (could be declared locally 3d_array)
	
	
	
	
	
	
	
	
	
	
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
		for(int i=1;i<=ngear;i++)
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
					// FIXME: This is not working. Need to figure this out.
					// logistic (3 parameters) with mean body 
					// weight deviations. 
					isel_npar(i) = 2;
					break;
					
				default: break;
			}
		}
		//cout<<"Number of estimated selectivity parameters\n"<<isel_npar<<endl;
	END_CALCS
	
	//Controls for prior on survey q.
	init_int nits;
	init_ivector q_prior(1,nits);
	init_vector q_mu(1,nits);
	init_vector q_sd(1,nits);
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
	// 12-> fraction of total mortality that takes place prior to spawning
	
	init_vector cntrl(1,12);
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
	!! nyr = nyr - retro_yrs;
	
	
PARAMETER_SECTION
	//Leading parameters
	//theta[1]		log_ro, or log_msy
	//theta[2]		steepness(h), or log_fmsy
	//theta[3]		log_m
	//theta[4]		log_avgrec
	//theta[5]		rho
	//theta[6]		kappa
	
	
	init_bounded_number_vector theta(1,npar,theta_lb,theta_ub,theta_phz);
	!! for(int i=1;i<=npar;i++) theta(i)=theta_ival(i);
	
	//Selectivity parameters (A very complicated ragged array)
	init_bounded_matrix_vector sel_par(1,ngear,1,jsel_npar,1,isel_npar,-15.,15.,sel_phz);
	LOC_CALCS
		//initial values for logistic selectivity parameters
		//set phase to -1 for fixed selectivity.
		for(int k=1;k<=ngear;k++)
		{
			if( isel_type(k)==1 || isel_type(k)==6 || isel_type(k)==7 )
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
	
	
	//FIXME Vivian Haist: possible bias in Bo due to dev_vector for log_rec_devs
	//		Try and compare estimates using init_bounded_vector version.
	
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
	init_bounded_vector log_m_devs(syr+1,nyr,-5.0,5.0,m_dev_phz);
	
	objective_function_value f;
    
	number ro;					//unfished age-1 recruits
	number bo;					//unfished spawning stock biomass
	number kappa;
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
	
	matrix nlvec(1,6,1,ilvec);	//matrix for negative loglikelihoods
	
	matrix jlog_sel(1,ngear,sage,nage);		//selectivity coefficients for each gear type.
	//matrix log_sur_sel(syr,nyr,sage,nage);	//selectivity coefficients for survey.
	
	matrix N(syr,nyr+1,sage,nage);			//Numbers at age
	matrix F(syr,nyr,sage,nage);			//Age-specific fishing mortality
	matrix M_tot(syr,nyr,sage,nage);		//Age-specific natural mortality
	matrix ft(1,ngear,syr,nyr);				//Gear specific fishing mortality rates
	matrix log_ft(1,ngear,syr,nyr);			//Gear specific log fishing mortlity rates
	matrix Z(syr,nyr,sage,nage);
	matrix S(syr,nyr,sage,nage);
	matrix ct(1,ngear,syr,nyr);				//predicted catch biomass
	matrix epsilon(1,nit,1,nit_nobs);		//residuals for survey abundance index
	matrix pit(1,nit,1,nit_nobs);			//predicted relative abundance index

	3darray Ahat(1,na_gears,1,na_nobs,a_sage-2,a_nage);		//predicted age proportions by gear & year
	3darray A_nu(1,na_gears,1,na_nobs,a_sage-2,a_nage);		//residuals for age proportions by gear & year
	
	3darray log_sel(1,ngear,syr,nyr,sage,nage);		//selectivity coefficients for each gear type.
	3darray Chat(1,ngear,syr,nyr,sage,nage);		//predicted catch-at-age
	
	sdreport_number sd_depletion;
	
PRELIMINARY_CALCS_SECTION
  //Run the model with input parameters to simulate real data.
  nf=0;
  if(SimFlag) 
  {
    initParameters();
    calcSelectivities();
    calcTotalMortality();
    simulation_model(rseed);
  }

RUNTIME_SECTION
    maximum_function_evaluations 100,100,500,25000,25000
    convergence_criteria 0.01,0.01,1.e-4,1.e-4


PROCEDURE_SECTION
	initParameters();

	calcSelectivities();
	
	calcTotalMortality();
	
	calcNumbersAtAge();
	
	calcFisheryObservations();
	
	calcAgeProportions();
	
	calc_survey_observations();
	
	calc_stock_recruitment();
	
	calc_objective_function();

	sd_depletion=sbt(nyr)/bo;
	
	if(mceval_phase()) mcmc_output();
	
	//The following causes a linker error
	//duplicate symbol in libdf1b2o.a
	//dvariable a=3.0;
	//cout<<"testing gammln(dvariable)"<<gammln(a)<<endl;
	
FUNCTION initParameters
	/*
	This function is used to extract the specific parameter values
	from the init_bounded_number_vector to the specific variables
	used in the code.
	
	Note that you must call this routine before runnning the 
	simulation model to generate fake data.
	*/
	
	ro = mfexp(theta(1));
	dvariable h = theta(2);
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
	
	
	m = mfexp(theta(3));
	log_avgrec = theta(4);
	log_recinit = theta(5);
	
	rho=theta(6);
	varphi=theta(7);
	
	if(verbose)cout<<"**** Ok after initParameters ****"<<endl;
	
FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs)
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

FUNCTION dvar_matrix cubic_spline_matrix(const dvar_matrix& spline_coffs)
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
	


FUNCTION calcSelectivities
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
		7) logistic selectivity with 3 parameters. 3rd parameter incorporates
		the age-specific deviations in mean weight at age; if 0 then no effect.
		
		Following the initializatoin of the selectivity curves, time-varying 
		considerations are implemented.
		
		TODO: Add penality (10.*square(avg_log_sel)) to objective function 
		in cases where estimating sel_coffs to mimic an init_bounded_dev vector.
		
	
	*/
	int i,j;
	dvariable p1,p2;//,p3;
	dvar_vector age_dev=age;
	dvar_matrix t1;
	dvar_matrix tmp(syr,nyr-1,sage,nage);
	dvar_matrix tmp2(syr,nyr,sage,nage);
	dvar_matrix ttmp2(sage,nage,syr,nyr);
	jlog_sel.initialize();
	log_sel.initialize();
	avg_log_sel.initialize();
	
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
				jlog_sel(j) = log( plogis(age,p1,p2) );
				break;
			
			case 6:
				// logistic selectivity for case 1 or 6
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				jlog_sel(j) = log( plogis(age,p1,p2) );
				break;
				
			case 2:		
				// age-specific selectivity coefficients
				jlog_sel(j)(sage,nage-1) = sel_par(j)(sage);
				jlog_sel(j,nage) = jlog_sel(j,nage-1);
				break;
				
			case 3:		
				// cubic spline
				jlog_sel(j)=cubic_spline( sel_par(j)(1) );
				break;
				
			case 4:		
				// time-varying cubic spline every year
				jlog_sel(j) = cubic_spline(sel_par(j)(1));
				t1 = cubic_spline_matrix(sel_par(j).sub(2,jsel_npar(j)));
				for(i = t1.indexmin(); i <= t1.indexmax(); i++)
				{
					tmp( syr+(i-t1.indexmin()) ) = t1(i);
				}
				break;
				
			case 5:		
				// time-varying bicubic spline
				ia.fill_seqadd( 0,1./(age_nodes(j)-1) );
				iy.fill_seqadd( 0,1./(yr_nodes(j)-1) );
				// bicubic_spline function is located in stats.cxx library
				bicubic_spline( iy,ia,sel_par(j),tmp2 );  
				break;
				
			case 7:
				// time-varying selectivity based on deviations in weight-at-age
				// FIXME This is not working and should not be used. (May 5, 2011)
				// SJDM:  I was not able to get this to run very well.
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				
				jlog_sel(j) = 0.;
				for(i = syr; i<=nyr; i++)
				{
					dvar_vector tmpwt=log(wt_obs(i)*1000)/mean(log(wt_obs*1000.));
					tmp2(i) = log( plogis(tmpwt,p1,p2) );
					
					//tmp2(i) = log(1./(1.+exp(p1-p2*wt_obs(i))));
				}	 
				break;
				
			default:
				jlog_sel(j)=0;
				break;
				
		}  // switch
		

		log_sel(j)(syr) = jlog_sel(j)+tmp2(syr);
		for(int i=syr;i<nyr;i++)
		{
			log_sel(j)(i+1)(sage,nage) = log_sel(j)(i)(sage,nage)+tmp(i)+tmp2(i+1);			
		}
		
		//subtract mean to ensure sum(log_sel)==0
		for(int i=syr;i<=nyr;i++)
			log_sel(j)(i) -= log(mean(mfexp(log_sel(j)(i))));
			
		//TODO 
			
		
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
	
	
	

FUNCTION calcTotalMortality
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
		//for(i=f_syr(k);i<=f_nyr(k);i++)
		{
			/*if(i==0) break;
						if(active(log_avg_f(k)))
							log_ft(k,i)=log_avg_f(k)+log_ft_devs(k,i);
						else log_ft(k,i)=-70.;*/	
			
			ftmp=0;
			if(obs_ct(k,i)>0)
				ftmp = mfexp(log_ft_pars(ki++));
			
			ft(k,i)=ftmp;
			
			//F(i)+=mfexp(log_ft(k,i)+log_sel(k)(i));
			F(i)+=ftmp*mfexp(log_sel(k)(i));			
		}
	}
	
	//Natural mortality (year and age specific)
	//M_tot(syr,nyr,sage,nage);
	M_tot = m;
	
	//Random walk in natural mortality.
	for(i=syr;i<=nyr;i++)
	{
		if(active(log_m_devs)&&i>syr)
		{
			M_tot(i)=M_tot(i-1)*exp(log_m_devs(i));
		}
		
	}
	m_bar = mean(M_tot);
	
	
	Z=M_tot+F;
	S=mfexp(-Z);
	if(verbose) cout<<"**** OK after calcTotalMortality ****"<<endl;
	
FUNCTION calcNumbersAtAge
	/*
		TODO Need to check the difference between the initialization 
		of the numbers at age here at the margins in comparison to the
		simulation model.
	*/
	
	int i,j;
	N.initialize();
	
	
	if(cntrl(5)){	//If initializing in at unfished conditions
		log_rt(syr) = log(ro);
		for(j=sage;j<=nage;j++)
		{
			N(syr,j)=ro*exp(-m_bar*(j-1.));
		}
	}
	else{			//If starting at unfished conditions
		log_rt(syr) = log_avgrec+log_rec_devs(syr);
		N(syr,sage)=mfexp(log_rt(syr));
		for(j=sage+1;j<=nage;j++)
		{
			N(syr,j)=mfexp(log_recinit+init_log_rec_devs(j))*exp(-m_bar*(j-sage));
		}
	}
	N(syr,nage)/=(1.-exp(-m_bar));
	
	
	//initial number of sage recruits from year syr+1, nyr;
	for(i=syr+1;i<=nyr;i++){
		log_rt(i)=log_avgrec+log_rec_devs(i);
		N(i,sage)=mfexp(log_rt(i));
	}
	N(nyr+1,sage)=mfexp(log_avgrec);
	
	/*
	for(j=sage;j<=nage;j++) 
	{
		if(cntrl(5))  //if starting at unfished state
		{
			N(syr,j)=ro*exp(-m_bar*(j-1));
		}
		else{
			log_rt(syr-j+sage)=log_avgrec+log_rec_devs(syr-j+sage);
			N(syr,j)=mfexp(log_rt(syr-j+sage))*exp(-m_bar*(j-sage));
		}
	}*/
	
	for(i=syr;i<=nyr;i++)
	{
		N(i+1)(sage+1,nage)=++elem_prod(N(i)(sage,nage-1),S(i)(sage,nage-1));
		N(i+1,nage)+=N(i,nage)*S(i,nage);
	}
	if(verbose)cout<<"**** Ok after calcNumbersAtAge ****"<<endl;

FUNCTION calcAgeProportions
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
	
FUNCTION calcFisheryObservations
	/*
	Dec 6, 2010.  Modified ct calculations to include
				  empirical weight at age data (wt_obs);
				
	Jan 16, 2011. 	modified this code to get age-comps from surveys, rather than 
					computing the age-comps in calc_fisheries_observations
	*/
	
	/*
		FIXED Reconcile the difference between the predicted catch here and in the simulation model.
	*/
	int i,k;
	ct.initialize();
	for(i=syr;i<=nyr;i++)
	{
		for(k=1;k<=ngear;k++)
		{
			
			dvar_vector log_va=log_sel(k)(i);// - log(mean(mfexp(log_sel(k)(i))));
			//dvar_vector fa=mfexp(log_ft(k,i)+log_va);
			
			
			//SJDM Jan 16, 2011 Modification as noted above.
			dvar_vector fa=ft(k,i)*mfexp(log_va);
			if(obs_ct(k,i)>0)
			{/*If there is a commercial fishery, then calculate the
			   catch-at-age (in numbers) and total catch (in weight)*/
				Chat(k,i)=elem_prod(elem_prod(elem_div(fa,Z(i)),1.-S(i)),N(i));
				ct(k,i) = Chat(k,i)*wt_obs(i);
			}
			else
			{/*If there is no commercial fishery the set Chat equal to 
			   the expected proportions at age.*/
				fa = mfexp(log_va);
				Chat(k,i)=elem_prod(elem_prod(elem_div(fa,Z(i)),1.-S(i)),N(i));
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
	
	
FUNCTION calc_survey_observations
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
	
	*/
	/*
		CHANGED add capability to accomodate priors for survey q's.
		DONE
	*/
	
	int i,j,ii,k;
	
	
	//survey abudance index residuals
	epsilon.initialize();
	pit.initialize();
	
	for(i=1;i<=nit;i++)
	{	
		int nx=0;		//counter for retrospective analysis
		dvar_matrix V(1,nit_nobs(i),sage,nage);  //vulnerable units for survey comparison
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
			//V(j)=elem_prod(N(ii),mfexp(log_va));
			switch(survey_type(i))
			{
				case 1:
					V(j)=elem_prod(Np,mfexp(log_va));
				break;
				case 2:
					V(j)=elem_prod(elem_prod(Np,mfexp(log_va)),wt_obs(ii));
				break;
				case 3:
					V(j)=elem_prod(Np,fec(ii));
				break;
			}
			
			//VB(j)=elem_prod(V(j),wt_obs(ii));		//SM Dec 6, 2010
			
			//If the survey is a spawn index, then need to account
			//for changes in fecundity.
			
			if(iyr(i,j)<=nyr) nx++;
		}
		dvar_vector t1 = rowsum(V);//V*wa;
		//cout<<"V\n"<<t1<<endl;
		
		dvar_vector zt=log(it(i).sub(1,nx))-log(t1(1,nx));
		epsilon(i).sub(1,nx) = zt-mean(zt);
		q(i) = exp(mean(zt));
		//pit(i).sub(1,nx)=(V*wa)*q(i);	//predicted index
		pit(i).sub(1,nx)=t1(1,nx)*q(i);	//predicted index
		
	}
	
	
	
	
	if(verbose)cout<<"**** Ok after calc_survey_observations ****"<<endl;
	
FUNCTION calc_stock_recruitment
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
	*/ 
	int i;
	dvariable tau = (1.-rho)/varphi;
	dvar_vector tmp_rt(syr+sage,nyr);
	dvar_vector lx(sage,nage); lx=1;
	for(i=sage+1;i<=nage;i++) lx(i)=lx(i-1)*exp(-m_bar);
	lx(nage)/=(1.-exp(-m_bar));
	
	dvariable phib = (lx*exp(-m_bar*cntrl(12))) * avg_fec;  //fec(syr);		//SM Dec 6, 2010
	dvariable so = kappa/phib;		//max recruits per spawner
	dvariable beta;
	bo = ro*phib;  					//unfished spawning biomass
	
	//sbt=rowsum(elem_prod(N,fec));			//SM Dec 6, 2010
	//CHANGED adjusted spawning biomass downward by ctrl(12)
	for(i=syr;i<=nyr;i++)
	{
		sbt(i) = elem_prod(N(i),exp(-Z(i)*cntrl(12)))*fec(i);
	}
	sbt(nyr+1) = N(nyr+1)*fec(nyr+1);
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
	
	if(verbose)cout<<"**** Ok after calc_stock_recruitment ****"<<endl;
	
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
		-3) likelihood of the survey age comps
		-4) likelihood of the fishery age comps
		-5) likelihood for stock-recruitment relationship
		-6) likelihood for fishery selectivities
	*/
	int i,j,k;
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
		
			nlvec(3,k)=2*dmvlogistic(O,P,nu,age_tau2(k),cntrl(6));
		
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
		if(active(sel_par(k))){
			//if not using logistic selectivity then
			if( isel_type(k)!=1 || isel_type(k)!=7 )  
			{
				for(i=syr;i<=nyr;i++)
				{
					//curvature in selectivity parameters
					dvar_vector df2=first_difference(first_difference(log_sel(k)(i)));
					nlvec(5,k)+=sel_2nd_diff_wt(k)/(nage-sage+1)*norm2(df2);
				
					//penalty for dome-shapeness
					for(j=sage;j<=nage-1;j++)
						if(log_sel(k,i,j)>log_sel(k,i,j+1))
							nlvec(6,k)+=sel_dome_wt(k)*square(log_sel(k,i,j)-log_sel(k,i,j+1));
				}
			}
		}
	}
	
	
	// CONSTRAINT FOR SELECTIVITY DEV VECTORS
	// Ensure vector of sel_par sums to 0. (i.e., a dev_vector)
	// TODO for isel_type==2 ensure mean 0 as well (ie. a dev_vector)
	for(k=1;k<=ngear;k++)
	{
		if( active(sel_par(k)) && isel_type(k)!=1 && isel_type(k)!=7 )
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
			if(isel_type(k)==4 || isel_type(k)==3)
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
	
	/*if(active(spline_coffs)||active(spline2_coffs))
		{   
			//lvec(6)=1.*norm2(spline2_coffs);
			for(i=syr;i<=nyr;i++)
			{
				//curvature in selectivity parameters
				dvar_vector df2=first_difference(first_difference(log_sel(i)));
				lvec(6)+=50.0/(nyr-syr+1)*norm2(df2);
				
				//penalty for dome-shapeness
				for(j=sage;j<=nage-2;j++)
					if(log_sel(i,j)>log_sel(i,j+1))
						lvec(6)+=3.125*square(log_sel(i,j)-log_sel(i,j+1));
			}
			//penalty for random walk in age-changes
			for(j=sage;j<=nage;j++){
				lvec(6)+=0.*norm2(first_difference(trans(log_sel)(j)));
			}
		}
		else
		{
			dvar_vector df2=first_difference(first_difference(log_sel(syr)));
			lvec(6)=50./(nyr-syr+1)*norm2(df2);
		}*/
	 
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
				
			case 2:		//lognormal
				ptmp=dlnorm(mfexp(theta(i)),theta_control(i,6),theta_control(i,7));
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
		if(q_prior(i)==1) qvec(i)=dnorm(log(q(i)),q_mu(i),q_sd(i));
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
	if(active(log_m_devs))
	{
		double std_mdev = cntrl(11);
		dvar_vector fd_mdevs=first_difference(log_m_devs);
		//pvec(2) = dnorm(log_m_devs,std_mdev);
		pvec(2) = dnorm(fd_mdevs,std_mdev);
	}
	
	
	if(verbose)
	{
		cout<<"nlvec\t"<<nlvec<<endl;
		cout<<"lvec\t"<<lvec<<endl;
		cout<<"priors\t"<<priors<<endl;
		cout<<"penalties\t"<<pvec<<endl;
	}
	f=sum(nlvec)+sum(lvec)+sum(priors)+sum(pvec)+sum(qvec);
	nf++;
	if(verbose)cout<<"**** Ok after initParameters ****"<<endl;
	}


FUNCTION void equilibrium(const double& fe,const double& ro, const double& kap, const double& m, const dvector& age, const dvector& wa, const dvector& fa, const dvector& va,double& re,double& ye,double& be,double& phiq,double& dphiq_df, double& dre_df)
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
	
	Modified args:
	re	-steady state recruitment
	ye	-steady state yield
	be	-steady state spawning biomass
	phiq		-per recruit yield
	dre_df		-partial of recruitment wrt fe
	dphiq_df	-partial of per recruit yield wrt fe
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
					//dlz_df=dlz_df*mfexp(-za[i-1]) - lz[i-1]*va[i-1]*mfexp(-za[i-1])/(1.-mfexp(-za[i]))
					dlz_df=dlz_df/(1.-mfexp(-za[i]))
							-lz[i-1]*mfexp(-za[i-1])*va[i]*mfexp(-za[i])
					/((1.-mfexp(-za[i]))*(1.-mfexp(-za[i])));
				}
		dphif_df=dphif_df+fa(i)*dlz_df;
		dphiq_df=dphiq_df+wa(i)*qa(i)*dlz_df+(lz(i)*wa(i)*va(i)*va(i))/za(i)*(exp(-za[i])-sa(i)/za(i));
	}
	double phif = lz*fa;
	phiq=sum(elem_prod(elem_prod(lz,wa),qa));
	re=ro*(kap-phie/phif)/(kap-1.);
	//re<0?re=0:NULL;
	if(re<=0) re=0;
	dre_df=(ro/(kap-1.))*phie/square(phif)*dphif_df;
	ye=fe*re*phiq;
	be=re*phif;	//spawning biomass
	
	//cout<<"Equilibrium\n"<<ro<<"\n"<<re<<"\n"<<ye<<endl;
	

FUNCTION void calc_reference_points()
	{
	/*
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
		control file.
		
		Use selectivity in the terminal year to calculate reference
		points.
	
	*/
	int i,j;
	double re,ye,be,phiq,dphiq_df,dre_df,fe;
	double dye_df,ddye_df;
	fe = 1.5*value(m_bar);
	
	/*Calculate average vulnerability*/
	dvector va_bar(sage,nage);
	va_bar.initialize();
	/*TODO user should specify allocation for MSY based reference points*/
	dvector allocation(1,ngear);
	allocation = dvector(fsh_flag/sum(fsh_flag));
	
	
	/*FIXME Allow for user to specify allocation among gear types.*/
	for(j=1;j<=ngear;j++)
	{
		va_bar+=allocation(j)*value(exp(log_sel(j)(nyr)));
		/*cout<<exp(log_sel(j)(nyr))<<endl;*/
	}
	
	/*CHANGED Changed equilibrium calculations based on average m */
	for(i=1;i<=20;i++)
	{
		//equilibrium(fe,value(ro),value(kappa),value(m),age,wa,fa,value(exp(log_sel(1)(nyr))),re,ye,be,phiq,dphiq_df,dre_df);
		equilibrium(fe,value(ro),value(kappa),value(m_bar),age,wt_obs(nyr),fec(nyr),va_bar,re,ye,be,phiq,dphiq_df,dre_df);
		dye_df = re*phiq+fe*phiq*dre_df+fe*re*dphiq_df;
		ddye_df = phiq*dre_df + re*dphiq_df;
		fe = fe - dye_df/ddye_df;
		//cout<<"fe\t"<<fe<<"\t"<<dye_df<<"\t"<<ye<<endl;
		if(sfabs(dye_df)<1.e-5)break;
	}
	fmsy=fe;
	equilibrium(fmsy,value(ro),value(kappa),value(m_bar),age,wt_obs(nyr),fec(nyr),va_bar,re,ye,be,phiq,dphiq_df,dre_df);
	msy=ye;
	bmsy=be;
	if(verbose)cout<<"**** Ok after calc_reference_points ****"<<endl;
	}

FUNCTION void simulation_model(const long& seed)
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
	*/
	
	
	cout<<"___________________________________________________\n"<<endl;
	cout<<"  **Implementing Simulation--Estimation trial**    "<<endl;
	cout<<"___________________________________________________"<<endl;
	cout<<"\tRandom Seed No.:\t"<< rseed<<endl;
	cout<<"___________________________________________________\n"<<endl;
	
	
	//Indexes:
	int i,j,k,ii,ki;

	//3darray Chat(1,ngear,syr,nyr,sage,nage);
	//C.initialize();
	
	
	
	/*----------------------------------*/
	/*	-- Generate random numbers --	*/
	/*----------------------------------*/
	random_number_generator rng(seed);
	dvector wt(syr-nage-1,nyr);			//recruitment anomalies
	dmatrix epsilon(1,nit,1,nit_nobs);  //observation errors in survey
	
	double sig = value(rho/varphi);
	double tau = value((1.-rho)/varphi);
	
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
			epsilon(k,j) *= sig/it_wt(k,j);
	}
	
	cout<<"	OK after random numbers\n";
	/*----------------------------------*/
	
	
	
	
	
	/*----------------------------------*/
    /*		--Initialize model--		*/
	/*CHANGED now calculating phie based on m_bar and avg_fec*/
	/*----------------------------------*/
	dvector lx=pow(exp(-value(m_bar)),age-min(age));
	lx(nage)/=(1.-exp(-value(m_bar)));
	double phie=(lx*exp(-value(m_bar)*cntrl(12)))*avg_fec;//fec(syr);
	so=kappa/phie;
	
	
	if(cntrl(2)==1) beta=(kappa-1.)/(ro*phie);
	if(cntrl(2)==2) beta=log(kappa)/(ro*phie);
	
	//Initial numbers-at-age with recruitment devs
	/*for(i=syr;i < syr+sage;i++)
			N(i,sage)=exp(log_avgrec+wt(i));
			
		for(j=sage+1;j<=nage;j++)
			N(syr,j)=exp(log_avgrec+wt(syr-j))*lx(j);
		*/
	
	N.initialize();
	if(cntrl(5)){	//If initializing in at unfished conditions
		log_rt(syr) = log(ro);
		for(j=sage;j<=nage;j++)
		{
			N(syr,j)=ro*exp(-m_bar*(j-1.));
		}
	}
	else{			//If starting at unfished conditions
		log_rt(syr) = log_avgrec;
		N(syr,sage)=mfexp(log_rt(syr));
		for(j=sage+1;j<=nage;j++)
		{
			N(syr,j)=mfexp(log_recinit+init_log_rec_devs(j))*exp(-m_bar*(j-sage));
		}
	}
	N(syr,nage)/=(1.-exp(-m_bar));
	
	//log_rt=log_avgrec+log_rec_devs;
	//log_rt(syr) = log(ro);

	for(i=syr+1;i<=nyr;i++){
		log_rt(i)=log_avgrec+log_rec_devs(i);
		N(i,sage)=mfexp(log_rt(i));
	}
	N(nyr+1,sage)=mfexp(log_avgrec);

	/*
	for(j=sage;j<=nage;j++) 
	{
		if(cntrl(5))  //if starting at unfished state
		{
			N(syr,j)=ro*exp(-m_bar*(j-1));
		}
		else{
			log_rt(syr-j+sage)=log_avgrec+log_rec_devs(syr-j+sage);
			N(syr,j)=mfexp(log_rt(syr-j+sage))*exp(-m_bar*(j-sage));
		}
	}

	N(syr,nage)/=(1.-exp(-m_bar));
	*/
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
		
	*/

	/*CHANGED May 15, 2011 calcSelectivities gets called from PRELIMINARY_CALCS*/
	dmatrix va(1,ngear,sage,nage);			//fishery selectivity
	d3_array dlog_sel(1,ngear,syr,nyr,sage,nage);
	dlog_sel=value(log_sel);
	/*
	for(k=1;k<=ngear;k++)
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
		
		
		for(k=1;k<=ngear;k++)
			va(k)=exp(dlog_sel(k)(i));
		
		//get_ft is defined in the Baranov.cxx file
		//CHANGED these ft are based on biomass at age, should be numbers at age
		//ft(i) = get_ft(oct,value(m),va,bt);
		ft(i) = get_ft(oct,value(m),va,value(N(i)),wt_obs(i));
		//cout<<trans(obs_ct)(i)<<"\t"<<oct<<endl;
		
		// overwrite observed catch incase it was modified by get_ft
		for(k=1;k<=ngear;k++)
			obs_ct(k,i)=oct(k);
		
		//total age-specific mortality
		//dvector zt(sage,nage);
		zt(i)=value(m);
		for(k=1;k<=ngear;k++){
			zt(i)+= ft(i,k)*exp(dlog_sel(k)(i));
		}
		
		
		//CHANGED definition of spawning biomass based on ctrl(12)
		sbt(i) = value(elem_prod(N(i),exp(-zt(i)*cntrl(12)))*fec(i));
		
		//Update numbers at age
		if(i>=syr+sage-1)
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
			A(k)(i)(a_sage(k),a_nage(k))=rmvlogistic(t1,0.3,i+seed);
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
	
FUNCTION dvector cis(const dvector& na)
	//Cohort Influenced Selectivity
	//This function returns a vector of residuals from a
	//linear regression of log(pa)= a+b*age+res that can be
	//used to modify age-based selectivity according to relative
	//cohort strengths.
	dvector y = log(na);
	dvector x = age;
	double b = sum(elem_prod(x-mean(x),y-mean(y)))/sum(square(x-mean(x)));
	double a = mean(y)-b*mean(x);
	dvector res = y - (a+b*x);
	return(res);

REPORT_SECTION
	if(verbose)cout<<"Start of Report Section..."<<endl;
	int i,j,k;
	REPORT(ControlFile);
	REPORT(f);
	REPORT(nlvec);
	REPORT(ro);
	double rbar=value(exp(log_avgrec));
	REPORT(rbar);
	REPORT(bo);
	REPORT(kappa);
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
	REPORT(ft);
	report<<"ut\n"<<elem_div(colsum(obs_ct),N.sub(syr,nyr)*wa)<<endl;
	report<<"bt\n"<<rowsum(elem_prod(N,wt_obs))<<endl;
	report<<"sbt\n"<<sbt<<endl;
	int rectype=int(cntrl(2));
	REPORT(rectype);
	REPORT(rt);
	dvector ln_rt=value(log_rt(syr,nyr));
	REPORT(ln_rt);
	REPORT(delta);
	REPORT(q);
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
	
	if(last_phase())
	{	calc_reference_points();
		REPORT(fmsy);
		REPORT(msy);
		REPORT(bmsy);
	}
		
	//Parameter controls
	dmatrix ctrl=theta_control;
	REPORT(ctrl);
	
	if(verbose)cout<<"END of Report Section..."<<endl;
	
	//Make copies of the report file using the ReportFileName
	//to ensure the results are saved to the same directory 
	//that the data file is in.
	if(last_phase() && PLATFORM =="Linux")
	{
		adstring copyrep = "cp iscam.rep " +ReportFileName;
		system(copyrep);
	}
	
	/*IN the following, I'm renaming the report file
	in the case where retrospective analysis is occurring*/
	if(retro_yrs && last_phase() && PLATFORM =="Linux")
	{
		//adstring rep="iscam.ret"+str(retro_yrs);
		//rename("iscam.rep",rep);
		adstring copyrep = "cp iscam.rep iscam.ret"+str(retro_yrs);
		system(copyrep);
	}
	
	

FUNCTION mcmc_output
	if(nf==1){
		ofstream ofs("iscam.mcmc");
		ofs<<"log.ro\t h\t log.m\t log.rbar\t rho\t kappa\t";
		ofs<<"bo\t bmsy\t msy\t fmsy\t"<<endl;
		
		ofstream of1("sbt.mcmc");
		ofstream of2("rt.mcmc");
	}
	
	// leading parameters & reference points
	calc_reference_points();
	ofstream ofs("iscam.mcmc",ios::app);
	ofs<<theta;
	ofs<<" "<<bo<<" "<<bmsy<<" "<<msy<<" "<<fmsy<<endl;
	
	// output spawning stock biomass
	ofstream of1("sbt.mcmc",ios::app);
	of1<<sbt<<endl;
	
	// output age-1 recruits
	ofstream of2("rt.mcmc",ios::app);
	of2<<rt<<endl;

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

