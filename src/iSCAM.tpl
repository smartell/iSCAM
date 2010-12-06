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
// TO DO: -add option for using empiracle weight-at-age data           //
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
//--                                                                 --//
//--                                                                 --//
//--                                                                 --//
//--                                                                 --//
//--                                                                 --//
//--                                                                 --//
//--                                                                 --//
// ------------------------------------------------------------------- //


DATA_SECTION
	init_adstring data_file;
	init_adstring control_file;


	!! ad_comm::change_datafile_name(data_file);
	
	int sim;
	int rseed;
	int retro_yrs;
	LOC_CALCS
		sim=0;
		rseed=999;
		int on,opt;
		//the following line checks for the "-sim" command line option
		//if it exists the if statement retreives the random number seed
		//that is required for the simulation model
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			sim=1;
			rseed=atoi(ad_comm::argv[on+1]);
			cout<<"______________________________________________________\n"<<endl;
			cout<<"    **Implementing Simulation--Estimation trial** "<<endl;
			cout<<"______________________________________________________"<<endl;
			cout<<"\tRandom Seed No.:\t"<< rseed<<endl;
			cout<<"______________________________________________________\n"<<endl;
			//if(sim)exit(1);
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
	
	init_int sage;
	init_int nage;
	vector age(sage,nage);
	!! age.fill_seqadd(sage,1);
	
	init_int ngear;				//number of gear types with unique selectivities
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
	vector fa(sage,nage);		//fecundity-at-age
	vector la(sage,nage);		//length-at-age
	vector wa(sage,nage);		//weight-at-age
	LOC_CALCS
	  la=linf*(1.-exp(-vonbk*(age-to)));
	  wa=a*pow(la,b);
	  fa=elem_prod(plogis(age,ah,gh),wa);
	  cout<<setprecision(2);		//2 decimal places for output
	  cout<<"la\n"<<la<<endl;
	  cout<<"wa\n"<<wa<<endl;
	  cout<<"fa\n"<<fa<<endl;
	  cout<<setprecision(5);
	END_CALCS
	
	//Time series data
	init_matrix catch_data(syr,nyr,1,ngear+1);
	matrix obs_ct(1,ngear,syr,nyr);
	int ft_count;
	LOC_CALCS
		for(k=1;k<=ngear;k++)
			obs_ct(k)=column(catch_data,k+1);
		int i;
		for(k=1;k<=ngear;k++)
		{	
			for(i=syr;i<=nyr;i++)
				obs_ct(k,i)>0 ? ft_count++ : NULL;
		}
		cout<<"ft_count\n"<<ft_count<<endl;
		cout<<"Ok after catch extraction"<<endl;
	END_CALCS
	
	init_int nit;
	init_ivector nit_nobs(1,nit);
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
		
	
	//End of data file
	init_int eof;	
	LOC_CALCS
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
	// ** Read parameter controls from control_file
	// ***************************************************
	!! ad_comm::change_datafile_name(control_file);
	
	init_int npar;
	init_matrix theta_control(1,npar,1,7);
	
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
	
	
	LOC_CALCS
		//cout up the number of selectivity parameters
		//depending on the value of isel_type
		isel_npar.initialize();
		for(int i=1;i<=ngear;i++)
		{	
			jsel_npar(i)=1;
			switch(isel_type(i))
			{
				case 1: isel_npar(i)=2; break;
				case 2: isel_npar(i)=(nage-sage); break;
				case 3: 
					isel_npar(i)=age_nodes(i);
					break;
				case 4:	 //annual cubic splines
					isel_npar(i)=age_nodes(i);
					jsel_npar(i)=(nyr-syr-retro_yrs)+1;
					break;
				case 5:  //bicubic spline
					jsel_npar(i)=age_nodes(i);
					isel_npar(i)=yr_nodes(i);
					break;
				default: break;
			}
		}
		//cout<<"Number of estimated selectivity parameters\n"<<isel_npar<<endl;
	END_CALCS
	
	
	//Miscellaneous controls
	// 1 -> verbose
	// 2 -> recruitment model (1=beverton-holt, 2=rickers)
	// 3 -> std in catch first phase
	// 4 -> std in catch in last phase
	// 5 -> assumed unfished in first year (0=FALSE, 1=TRUE)
	
	init_vector cntrl(1,5);
	int verbose;
	
	init_int eofc;
	LOC_CALCS
		verbose = cntrl(1);
		if(eofc==999){
			cout<<"\n -- END OF CONTROL FILE -- \n"<<endl;
		}else{
			cout<<"\n *** ERROR CONTROL FILE *** \n"<<endl; exit(1);
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
	//theta[5]		log_avg_f
	
	init_bounded_number_vector theta(1,npar,theta_lb,theta_ub,theta_phz);
	!! for(int i=1;i<=npar;i++) theta(i)=theta_ival(i);
	
	//Selectivity parameters (A very complicated ragged array)
	init_bounded_matrix_vector sel_par(1,ngear,1,jsel_npar,1,isel_npar,-5.,5.,sel_phz);
	LOC_CALCS
		//initial values for logistic selectivity parameters
		//set phase to -1 for fixed selectivity.
		for(int k=1;k<=ngear;k++)
		{
			if(isel_type(k)==1)
			{
				sel_par(k,1,1) = log(ahat(k));
				sel_par(k,1,2) = log(ghat(k));
			}
		}
	END_CALCS
	
	//Fishing mortality rate parameters
	//init_bounded_vector log_avg_f(1,ngear,-15.0,15.0,1);
	//init_bounded_matrix log_ft_devs(1,ngear,syr,nyr,-15.,15.,3);
	//init_bounded_number_vector log_avg_f(1,ngear,-15.,15.0,-1/*f_phz*/);
	//init_bounded_vector_vector log_ft_devs(1,ngear,f_syr,f_nyr,-5,5,-1/*f_phz*/);
	init_bounded_vector log_ft_pars(1,ft_count,-30.,3.0,1);
	
	LOC_CALCS
		log_ft_pars = log(0.1);
	END_CALCS
	
	
	//Annual recruitment deviations
	!! int ii;
	!! ii=syr-nage+sage;
	!! if(cntrl(5)) ii = syr;  //if initializing with ro
	init_bounded_dev_vector log_rec_devs(ii,nyr,-15.,15.,2);
	
	
	objective_function_value f;
    
	number ro;					//unfished age-1 recruits
	number bo;					//unfished spawning stock biomass
	number kappa;
	number m;
	number log_avgrec;			//log of average recruitment.
	//number log_avg_f;			//log of average fishing mortality
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
	
	
	matrix nlvec(1,6,1,ilvec);	//matrix for negative loglikelihoods
	
	matrix jlog_sel(1,ngear,sage,nage);		//selectivity coefficients for each gear type.
	//matrix log_sur_sel(syr,nyr,sage,nage);	//selectivity coefficients for survey.
	
	matrix N(syr,nyr+1,sage,nage);			//Numbers at age
	matrix F(syr,nyr,sage,nage);			//Age-specific fishing mortality
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
  if(sim) {
	initialize_parameters();
	simulation_model(rseed);
	}

RUNTIME_SECTION
    maximum_function_evaluations 100,100,500,5000,5000
    convergence_criteria 0.01,0.01,1.e-4,1.e-4


PROCEDURE_SECTION
	initialize_parameters();

	calc_selectivities();
	
	calc_mortality();
	
	calc_numbers_at_age();
	
	calc_fishery_observations();
	
	calc_age_proportions();
	
	calc_survey_observations();
	
	calc_stock_recruitment();
	
	calc_objective_function();

	sd_depletion=sbt(nyr)/bo;
	
	if(mceval_phase()) mcmc_output();
	
FUNCTION initialize_parameters
	/*
	This function is used to extract the specific parameter values
	from the init_bounded_number_vector to the specific variables
	used in the code.
	
	Note that you must call this routine before runnning the 
	simulation model to generate fake data.
	*/
	
	ro = mfexp(theta(1));
	dvariable h = theta(2);
	kappa = (4.*h/(1.-h));
	
	//Alternative parameterization using MSY and FMSY as leading parameters
	
	
	m = mfexp(theta(3));
	log_avgrec = theta(4);
	
	rho=theta(5);
	varphi=theta(6);
	
	//I'm not sure about the following.
	//dvariable tau = (1.-rho)/varphi;
	//ro *= exp(-0.5*tau*tau);
	if(verbose)cout<<"**** Ok after initialize_parameters ****"<<endl;
	
FUNCTION dvar_vector cubic_spline(const dvar_vector& spline_coffs)
	RETURN_ARRAYS_INCREMENT();
	int nodes=size_count(spline_coffs);
	dvector ia(1,nodes);
	dvector fa(sage,nage);
	ia.fill_seqadd(0,1./(nodes-1));
	//fa.fill_seqadd(sage,1);			//Thanks Dr. Ahrens for helping me with this dumb bug.
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
	int nodes=spline_coffs.colmax()-spline_coffs.colmin()+1;
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
	


FUNCTION calc_selectivities
	/*
		This function loops over each ngear and calculates the corresponding
		selectivity for that gear type. It first uses a switch statement 
		to calculate selectivity curves based on isel_type where:
		1) logistic selectivity with 2 parameters
		2) age-specific selectivity coefficients with (nage-sage) parameters
		   and the last two age-classes are assumed to have the same selectivity.
		3) a reduced age-specific parameter set based on a bicubic spline.
		
		Following the initializatoin of the selectivity curves, time-varying 
		considerations are implemented.
		
		Having trouble getting a positive definate hessian with the cubic
		spline and bicubic spline interpolations.  May have to revisit these options.
	
	*/
	int i,j;
	dvariable p1,p2;
	dvar_matrix t1;
	dvar_matrix tmp(syr,nyr-1,sage,nage);
	dvar_matrix tmp2(syr,nyr,sage,nage);
	dvar_matrix ttmp2(sage,nage,syr,nyr);
	jlog_sel.initialize();
	log_sel.initialize();
	
	for(j=1;j<=ngear;j++)
	{
		tmp.initialize(); tmp2.initialize();
		dvector iy(1,yr_nodes(j));
		dvector ia(1,age_nodes(j));
		
		switch(isel_type(j))
		{
			case 1:		//logistic selectivity
				p1 = mfexp(sel_par(j,1,1));
				p2 = mfexp(sel_par(j,1,2));
				jlog_sel(j)=log(plogis(age,p1,p2));
				break;
			case 2:		//age-specific selectivity coefficients
				jlog_sel(j)(sage,nage-1)=sel_par(j)(sage);
				jlog_sel(j,nage)=jlog_sel(j,nage-1);
				break;
			case 3:		//cubic spline
				jlog_sel(j)=cubic_spline(sel_par(j)(1));
				break;
			case 4:		//time-varying cubic spline every year
				jlog_sel(j) = cubic_spline(sel_par(j)(1));
				t1 = cubic_spline_matrix(sel_par(j).sub(2,jsel_npar(j)));
				for(i=t1.indexmin();i<=t1.indexmax();i++) tmp(syr+(i-t1.indexmin()))=t1(i);
				break;
			case 5:		//time-varying bicubic spline
				ia.fill_seqadd(0,1./(age_nodes(j)-1));
				iy.fill_seqadd(0,1./(yr_nodes(j)-1));
				bicubic_spline(iy,ia,sel_par(j),tmp2);  //function located in stats.cxx
				break;
			default:
				jlog_sel(j)=0;
		}
		

		log_sel(j)(syr) = jlog_sel(j)+tmp2(syr);
		for(int i=syr;i<nyr;i++)
		{
			log_sel(j)(i+1)(sage,nage) = log_sel(j)(i)(sage,nage)+tmp(i)+tmp2(i+1);			
		}
		
		//subtract mean to ensure sum(log_sel)==0
		for(int i=syr;i<=nyr;i++)
			log_sel(j)(i) -= log(mean(mfexp(log_sel(j)(i))));
			
			
		
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
	
	if(verbose)cout<<"**** Ok after calc_selectivities ****"<<endl;
	
	
	

FUNCTION calc_mortality
	/*
	This routing calculates fishing mortality, total mortality
	and annaul survival rates (exp(-Z)) for each age in each
	year.
	
	There is a complication in that if there is a survey gear
	then the fishing mortality rate is set to an extremely small
	value: exp(-70.)~3.975e-31, which is effectively 0.
	
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
			
			//if(obs_ct(k,i)>0)
			obs_ct(k,i)>0? ftmp = mfexp(log_ft_pars(ki++)): ftmp=0;
			ft(k,i)=ftmp;
			
			//F(i)+=mfexp(log_ft(k,i)+log_sel(k)(i));
			F(i)+=ftmp*mfexp(log_sel(k)(i));			
		}
	}
	
	Z=m+F;
	S=mfexp(-Z);
	if(verbose) cout<<"**** OK after calc_mortality ****"<<endl;
	
FUNCTION calc_numbers_at_age
	int i,j;
	N.initialize();
	
	//log_rt=log_avgrec+log_rec_devs;
	
	for(i=syr+1;i<=nyr;i++){
		log_rt(i)=log_avgrec+log_rec_devs(i);
		N(i,sage)=mfexp(log_rt(i));
	}
	
	for(j=sage;j<=nage;j++) 
	{
		if(cntrl(5))  //if starting at unfished state
		{
			N(syr,j)=ro*exp(-m*(j-1));
		}
		else{
			log_rt(syr-j+sage)=log_avgrec+log_rec_devs(syr-j+sage);
			N(syr,j)=mfexp(log_rt(syr-j+sage))*exp(-m*(j-1));
		}
	}
	N(syr,nage)/=(1.-exp(-m));
	
	for(i=syr;i<=nyr;i++)
	{
		N(i+1)(sage+1,nage)=++elem_prod(N(i)(sage,nage-1),S(i)(sage,nage-1));
		N(i+1,nage)+=N(i,nage)*S(i,nage);
	}
	if(verbose)cout<<"**** Ok after calc_numbers_at_age ****"<<endl;

FUNCTION calc_age_proportions
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
	if(verbose)cout<<"**** Ok after calc_age_proportions ****"<<endl;
	
FUNCTION calc_fishery_observations
	int i,k;
	for(i=syr;i<=nyr;i++)
	{
		for(k=1;k<=ngear;k++)
		{
			
			dvar_vector log_va=log_sel(k)(i);// - log(mean(mfexp(log_sel(k)(i))));
			//dvar_vector fa=mfexp(log_ft(k,i)+log_va);
			
			
			/*need to add a tiny constant to deal with catch-age data 
			for non-extractive survey age comps*/
			dvar_vector fa=ft(k,i)*mfexp(log_va)+1.e-30;
			
			//Catch-at-age by gear
			Chat(k,i)=elem_prod(elem_prod(elem_div(fa,Z(i)),1.-S(i)),N(i));
			
			//Catch weight by gear
			ct(k,i)=Chat(k,i)*wa;  
		}
	}
	if(verbose)cout<<"**** Ok after calc_fishery_observations ****"<<endl;
	
	
FUNCTION calc_survey_observations
	/*This code needs to be modified to accomodate
	multiple surveys or block changes in survey q.
	
	Oct 31, 2010, added retrospective counter.
	
	Nov 22, 2010, adding multiple surveys. Still need to check with retrospective option
	
	Nov 30, 2010, adjust the suvery biomass by the fraction of Z that has occurred 
	when the survey was conducted. For herring spawning biomass this would be after the 
	fishery has taken place.
	*/
	
	int i,j,ii,k;
	
	
	//survey abudance index residuals
	epsilon.initialize();
	
	for(i=1;i<=nit;i++)
	{	
		int nx=0;		//counter for retrospective analysis
		dvar_matrix V(1,nit_nobs(i),sage,nage);  //vulnerable numbers 
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
			V(j)=elem_prod(Np,mfexp(log_va));
		
			if(iyr(i,j)<=nyr) nx++;
		}
		dvar_vector t1 = V*wa;
		dvar_vector zt=log(it(i).sub(1,nx))-log(t1(1,nx));
		epsilon(i).sub(1,nx) = zt-mean(zt);
		q(i) = exp(mean(zt));
		pit(i).sub(1,nx)=(V*wa)*q(i);	//predicted index
		
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
	*/ 
	int i;
	dvariable tau = (1.-rho)/varphi;
	dvar_vector tmp_rt(syr+sage,nyr);
	dvar_vector lx(sage,nage); lx=1;
	for(i=sage+1;i<=nage;i++) lx(i)=lx(i-1)*exp(-m);
	lx(nage)/=(1.-exp(-m));
	dvariable phib = lx * fa;
	dvariable so = kappa/phib;		//max recruits per spawner
	dvariable beta;
	bo = ro*phib;  
	sbt=(N*fa);
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
	for(k=1;k<=ngear;k++){
		if(active(log_ft_pars))
			nlvec(1,k)=dnorm(log(obs_ct(k).sub(syr,nyr)+o)-log(ct(k).sub(syr,nyr)+o),sig_c);
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
		
			nlvec(3,k)=dmvlogistic(O,P,nu,age_tau2(k),0.02);
		
			for(i=1;i<=na_nobs(k);i++)
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
			if(isel_type(k)!=1)  //if not using logistic selectivity then
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
	for(k=1;k<=ngear;k++)
	{
		if(active(sel_par(k))&&isel_type(k)!=1)
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
		default:	//uniform density
			ptmp=1./(theta_control(i,3)-theta_control(i,2));
			break;
		}
		priors(i)=ptmp;
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
		-2) normal prior on fishing mortality deviates with
			a large standard deviation of 50.
		-3) normal prior for log rec devs with std=50.
		-4) penalty for the spline coefficients, which
			should keep the random walk in check.  Might
			want to try first differences.
	*/
	
	dvar_vector pvec(1,7);
	pvec.initialize();
	
	//Penalties to regularize the solution for fishing mortality rates
	dvariable log_fbar = mean(log_ft_pars);
	if(last_phase())
	{
		/*for(k=1;k<=ngear;k++){
					if(active(log_avg_f(k)))
					{
						pvec(1)+= 0.01*square(log_avg_f(k)-log(0.2));
						dvariable s=mean(log_ft_devs(k));
						pvec(2)+= 100000.0*s*s;
					}
				}*/
		pvec(1) = 0.01*square(log_fbar-log(0.2));
		//Penalty for log_rec_devs (large variance here)
		pvec(4) = dnorm(log_rec_devs,5.);
	}
	else
	{
		/*for(k=1;k<=ngear;k++){
					if(active(log_avg_f(k)))
					{
						pvec(1)+=500000.*square(log_avg_f(k)-log(0.2));
						dvariable s=mean(log_ft_devs(k));
						pvec(2)+= 100000.0*s*s;
					}
				}*/
		pvec(1) = 50000.*square(log_fbar-log(0.2));
		//Penalty for log_rec_devs (CV ~ 0.0707) in early phases
		pvec(4)=100.*norm2(log_rec_devs);
	}
	
	
	
	if(verbose)
	{
		cout<<"nlvec\t"<<nlvec<<endl;
		cout<<"lvec\t"<<lvec<<endl;
		cout<<"priors\t"<<priors<<endl;
		cout<<"penalties\t"<<pvec<<endl;
	}
	f=sum(nlvec)+sum(lvec)+sum(priors)+sum(pvec);
	nf++;
	if(verbose)cout<<"**** Ok after initialize_parameters ****"<<endl;



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
	re<0?re=0:NULL;
	dre_df=(ro/(kap-1.))*phie/square(phif)*dphif_df;
	ye=fe*re*phiq;
	be=re*phif;	//spawning biomass
	
	//cout<<"Equilibrium\n"<<ro<<"\n"<<re<<"\n"<<ye<<endl;
	

FUNCTION void calc_reference_points()
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
	fe = 1.5*value(m);
	
	/*Calculate average vulnerability*/
	dvector va_bar(sage,nage);
	va_bar.initialize();
	dvector allocation(1,ngear);
	allocation = dvector(fsh_flag/sum(fsh_flag));
	for(j=1;j<=ngear;j++)
	{
		va_bar+=allocation(j)*value(exp(log_sel(j)(nyr)));
		/*cout<<exp(log_sel(j)(nyr))<<endl;*/
	}


	for(i=1;i<=20;i++)
	{
		//equilibrium(fe,value(ro),value(kappa),value(m),age,wa,fa,value(exp(log_sel(1)(nyr))),re,ye,be,phiq,dphiq_df,dre_df);
		equilibrium(fe,value(ro),value(kappa),value(m),age,wa,fa,va_bar,re,ye,be,phiq,dphiq_df,dre_df);
		dye_df = re*phiq+fe*phiq*dre_df+fe*re*dphiq_df;
		ddye_df = phiq*dre_df + re*dphiq_df;
		fe = fe - dye_df/ddye_df;
		//cout<<"fe\t"<<fe<<"\t"<<dye_df<<"\t"<<ye<<endl;
		if(sfabs(dye_df)<1.e-5)break;
	}
	fmsy=fe;
	equilibrium(fmsy,value(ro),value(kappa),value(m),age,wa,fa,va_bar,re,ye,be,phiq,dphiq_df,dre_df);
	msy=ye;
	bmsy=be;
	if(verbose)cout<<"**** Ok after calc_reference_points ****"<<endl;
	

FUNCTION void simulation_model(const long& seed)
	/*
	Call this routine to simulate data for simulation testing.
	The random number seed can be used to repeat the same 
	sequence of random number for simulation testing.
	
	Implemented using the "-sim 99" command line option where
	99 is the random number seed.
	
	-sim 99 is a special case used for the manuscript for case 1.
	-sim 000 is a special case with 0 error (exact data)
	
	-This routine will over-write the observations in memory
	with simulated data, where the true parameter values are
	the initial values.  Change the standard deviations of the 
	random number vectors epsilon (observation error) or 
	recruitment devs wt (process error).
	*/
	int i,j,k,ii;

	//3darray Chat(1,ngear,syr,nyr,sage,nage);
	//C.initialize();
	
	//Random number generator
	random_number_generator rng(seed);
	dvector wt(syr-nage-1,nyr);		//recruitment anomalies
	dvector epsilon(1,nit);			//observation errors in survey
	double sig = value(rho/varphi);
	double tau = value((1.-rho)/varphi);
	if(seed==000)
	{
		sig=0;
		tau=0;
	}
	wt.fill_randn(rng); wt *= tau;
	epsilon.fill_randn(rng); epsilon *= sig;
	

    //Initialize model
	dvector lx=pow(exp(-value(m)),age-1.);
	lx(nage)/=(1.-exp(-value(m)));
	double phie=lx*fa;		//this is the sum of products for two vectors
	so=kappa/phie;
	if(cntrl(2)==1) beta=(kappa-1.)/(ro*phie);
	if(cntrl(2)==2) beta=log(kappa)/(ro*phie);
	
	//Initial numbers-at-age with recruitment devs
	N(syr,1)=exp(log_avgrec+wt(syr));
	for(j=2;j<=nage;j++)N(syr,j)=exp(log_avgrec+wt(syr-j))*lx(j);
	
	
	dmatrix va(1,ngear,sage,nage);				//fishery selectivity
	d3_array vai(1,ngear,syr,nyr,sage,nage);	//time-varying selectivity
	for(k=1;k<=ngear;k++)
	{
		va(k)=plogis(age,0.33*(nage-sage),1.5);
		for(i=syr;i<=nyr;i++)
			vai(k)(i)=va(k);
	}

	//Special case for manunscript
	if(seed==99)
	{
		double g=0.3;
		for(k=1;k<=ngear;k++)
		{
			va(k)=plogis(age,5-k,0.5);
			for(i=syr;i<=nyr;i++)
			{
				g = 0.2 * (nyr-i)/(nyr-syr);  //goes from strongly domed to asymptotic.
				if(fsh_flag(k))
					vai(k)(i)=eplogis(age,1.5,5-k,g);
				else
					vai(k)(i)=va(k);
			}
				
		}
	}
	//dvector va=plogis(age,4,1.5);				//fishery selectivity
	//dmatrix vat(syr,nyr,sage,nage);				//matrix for time-varying selectivity 
	//dvector va=eplogis(age,1/1.5,4,0.1);		//fishery selectivity domed
	//dvector vax=plogis(age,4.5,1.8);			//survey selectivity
	
	
	for(i=syr;i<=nyr;i++)
	{   
		//Approximate fishing mortality with Popes approximation
		//double btmp = value(N(i))*exp(-m/2.)*elem_prod(va,wa);
		//double ftmp = obs_ct(i)/btmp;
		//dvector va_dev = cis(value(N(i)));
		//va = elem_prod(plogis(age,4,1.5),exp(0.0*va_dev));
		//vat(i) = va;
		//va = va/max(va);
		//cout<<setprecision(3)<<"va\t"<<va<<endl;
		
		
		dvector bt = elem_prod(value(N(i)),wa);
		//double ft = get_ft(obs_ct(1,i),value(m),va,bt);
		dvector oct = trans(obs_ct)(i);
		for(k=1;k<=ngear;k++)
			va(k)=vai(k)(i);
		dvector ft = get_ft(oct,value(m),va,bt);
		// overwrite observed catch incase it was modified by get_ft
		for(k=1;k<=ngear;k++)
			obs_ct(k,i)=oct(k);
		
		dvector zt(sage,nage);
		zt=value(m);
		for(k=1;k<=ngear;k++){
			zt+= ft(k)*va(k);
		}
		 
		
		//Update numbers at age
		double et=value(N(i))*fa;
		double rt;
		if(cntrl(2)==1)rt=value(so*et/(1.+beta*et));
		if(cntrl(2)==2)rt=value(so*et*exp(-beta*et));
		N(i+1,1)=rt*exp(wt(i)-0.5*tau*tau);
		N(i+1)(2,nage)=++elem_prod(N(i)(1,nage-1),exp(-zt(1,nage-1)));
		N(i+1,nage)+=N(i,nage)*exp(-zt(nage));
		
		//Catch & Catch-at-age
		for(k=1;k<=ngear;k++)
		{
			if(ft(k)>0)
			{
				d3C(k)(i)=elem_prod(elem_div(ft(k)*va(k),zt),elem_prod(1.-exp(-zt),value(N(i))));
				obs_ct(k,i)=d3C(k)(i)*wa;
			}
			else	//if this is a survey
			{
				d3C(k)(i)=elem_prod(elem_div(va(k),zt),elem_prod(1.-exp(-zt),value(N(i))));
			}
		}
		
	}
	
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
			
			//Ahat(k,i,a_sage(k)-2)=ii;
			//Ahat(k,i,a_sage(k)-1)=ig;
			dvector t1=pa(a_sage(k),a_nage(k));
			t1/=sum(t1);
			A(k)(i)(a_sage(k),a_nage(k))=rmvlogistic(t1,0.3,i+seed);
			if(seed==000)
				A(k)(i)(a_sage(k),a_nage(k))=t1;
			//cout<<iyr<<"\t"<<k<<endl;
		}
	}
	//cout<<Ahat<<endl;
	
	/*
		Calculations for survey time series
	*/
	for(i=1;i<=nit;i++)
	{   
		//ii=iyr(i);
		//it(i) = (value(N(ii))*elem_prod(wa,va(igr(i))))*exp(epsilon(i));
		
	}
	//cout<<"Simulated exploitation rate\n"<<elem_div(ct,(N.sub(syr,nyr)*elem_prod(wa,va)))<<endl;
	
	cout<<"_________________________________"<<endl;
	ofstream ofs("iscam.sim");
	ofs<<"va\n"<<va<<endl;
	ofs<<"sbt\n"<<N*fa<<endl;
	ofs<<"ct\n"<<obs_ct<<endl;
	ofs<<"ut\n"<<elem_div(colsum(obs_ct),N.sub(syr,nyr)*wa)<<endl;
	ofs<<"iyr\n"<<iyr<<endl;
	ofs<<"it\n"<<it<<endl;
	ofs<<"N\n"<<N<<endl;
	ofs<<"A\n"<<A<<endl;
	
	
	
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
	REPORT(control_file);
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
	REPORT(fa);
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
	report<<"bt\n"<<N*wa<<endl;
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

	REPORT(a_sage);
	REPORT(a_nage);
	REPORT(A); 
	REPORT(Ahat);
	REPORT(A_nu);
	REPORT(N);
	
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
	
	
	/*IN the following, I'm renaming the report file
	in the case where retrospective analysis is occurring*/
	if(retro_yrs && last_phase())
	{
		adstring rep="iscam.ret"+str(retro_yrs);
		rename("iscam.rep",rep);
	}

FUNCTION mcmc_output
	if(nf==1){
		ofstream ofs("iscam.mcmc");
		ofs<<"log.ro\t h\t log.m\t log.rbar\t rho\t kappa\t";
		ofs<<"bo\t bmsy\t msy\t fmsy\t"<<endl;
		
		ofstream of1("sbt.mcmc");
	}
	
	//leading parameters & reference points
	calc_reference_points();
	ofstream ofs("iscam.mcmc",ios::app);
	ofs<<theta;
	ofs<<" "<<bo<<" "<<bmsy<<" "<<msy<<" "<<fmsy<<endl;
	
	//output spawning stock biomass
	ofstream of1("sbt.mcmc",ios::app);
	of1<<sbt<<endl;

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

	#include <admodel.h>
	#include <time.h>
	#include <stats.cxx>
	#include <baranov.cxx>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
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
	cout<<"*******************************************"<<endl;

