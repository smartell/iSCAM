//  -----------------------------------------------------------------------  //
//  test: A simple program for testing the MSY algorithm.
//  
//  
//  Created by Martell on 2012-07-29.
//  Copyright (c) 2012. All rights reserved.
//
//  Comments:
//  This program calculates MSY-based reference points. 
//  
//  INSTRUCTIONS:
//    -1) Define a Msy class object using: Msy cMSY(ro,h,m,wa,fa,V);
//        where: 
//              ro is unfished recruitment
//              h  is steepness
//              m  is natural mortality rate
//              wa is a vector of mean weights-at-age
//              fa is a vector of mean mature biomass at age (Spawning bio)
//              V  is a matrix of selctivities (row = gear, col = age)
//    -2) Use the class methods to obtain information you desire.
//        Methods:
//  -----------------------------------------------------------------------  //
DATA_SECTION
	int A;
	!! A = 20;
	int ngear;  //number of fishing fleets
	!! ngear = 3;
	vector age(1,A);
	vector ah(1,ngear);
	vector gh(1,ngear);
	vector wa(1,A);
	vector fa(1,A);
	matrix V(1,ngear,1,A);
	LOC_CALCS
		age.fill_seqadd(1,1);
		ah = 6.;
		gh = 0.5 * ah;
		wa = pow((1.-exp(-0.2*age)),3);
		fa = wa;
		fa(1,4) = 0;
		int k;
		// Selectivities for each gear.
		// Chang the ah gh values here to not how they impact MSY calcs.
		for(k=1;k<=ngear;k++)
		{
			V(k) = log(plogis(age,ah(k)+double(k),gh(k)));  // progressive increase in age-50% harvest
		}
		for(k=1;k<=ngear;k++)
		{
			V(k) = exp(V(k) - log(mean(exp(V(k)))) );
		}
		cout<<"fa\n"<<fa<<endl;
	END_CALCS
	

PARAMETER_SECTION
	init_number dum;
	objective_function_value f;

PROCEDURE_SECTION
	
	
FINAL_SECTION
	double ro    = 500.0;
	double h     = 0.75;
	double m     = 0.3;
	double rho   = 0.00;
	int i;
	dvector fe(1,ngear);
	
	Msy cMSY(ro,h,m,rho,wa,fa,V);
	
	/* Lets have a look at the unfished spawning biomass per recruit */
	double phie = cMSY.getPhie();
	cout<<"\tUnfished spawning biomass per recruit = "<<phie<<endl;
	
	/* Now specify an initial guess for Fmsy and get Fmsy from Msy class */
	dvector fmsy(1,ngear);
	fmsy = m*h/0.8;  // initial guess for fmsy
	cMSY.get_fmsy(fmsy);
	//cMSY.get_fmsy_safe(fmsy,1e-3,5.0);
	
	
	/* Now print out other associated reference points */
	cout<<"\t---------------------------------------"<<endl;
	cout<<setprecision(5)<<endl;
	cout<<"\t M = 0.3 and h = 0.75"<<endl;
	cout<<"\tUnfished spawning biomass             = "<<cMSY.getBo()<<endl;
	cout<<"\tF_MSY values for each fishing fleet   ="<<fmsy<<endl;
	cout<<"\tMSY values for each fishing fleet     ="<<cMSY.getMsy()<<endl;
	cout<<"\tSpawning biomass at MSY (Bmsy)        = "<<cMSY.getBmsy()<<endl;
	cout<<"\tRecruitment at MSY (Rmsy)             = "<<cMSY.getRmsy()<<endl;
	cout<<"\tSpawning potential ratio at MSY       = "<<cMSY.getSprMsy()<<endl;
	exit(1);
	
	/* Now change M and steepness and recalculate reference poinots */
	h = 0.80;
	m = 0.20;
	cMSY.set_h(h);
	cMSY.set_m(m);
	cMSY.calc_bo(m,fa);  //necessary to update phie due to change in m.
	fmsy = 0.6*m/ngear;  // initial guess for fmsy
	cMSY.get_fmsy(fmsy);
	/* Now print out other associated reference points */
	cout<<"\t---------------------------------------"<<endl;
	cout<<setprecision(4)<<endl;
	cout<<"\t M = 0.2 and h = 0.80"<<endl;
	cout<<"\tUnfished spawning biomass             = "<<cMSY.getBo()<<endl;
	cout<<"\tF_MSY values for each fishing fleet   ="<<fmsy<<endl;
	cout<<"\tMSY values for each fishing fleet     ="<<cMSY.getMsy()<<endl;
	cout<<"\tSpawning biomass at MSY (Bmsy)        = "<<cMSY.getBmsy()<<endl;
	cout<<"\tRecruitment at MSY (Rmsy)             = "<<cMSY.getRmsy()<<endl;
	cout<<"\tSpawning potential ratio at MSY       = "<<cMSY.getSprMsy()<<endl;
	
	/* Now run the model with a fixed allocation among gears */
	dvector allocation(1,ngear);
	allocation.fill_seqadd(1,1);
	allocation/=sum(allocation);
	cout<<"\t---------------------------------------"<<endl;
	cMSY.get_fmsy(fmsy,allocation);
	cout<<setprecision(4)<<endl;
	cout<<"\t M = 0.2 and h = 0.80, & allocations  "<<endl;
	cout<<"\tAllocaton for each gear               ="<<allocation<<endl;
	cout<<"\tUnfished spawning biomass             = "<<cMSY.getBo()<<endl;
	cout<<"\tFishing rate for each fishing fleet   ="<<fmsy<<endl;
	cout<<"\tYield values for each fishing fleet   ="<<cMSY.getYe()<<endl;
	cout<<"\tSpawning biomass at given allocation  = "<<cMSY.getBe()<<endl;
	cout<<"\tRecruitment at given allocaiton       = "<<cMSY.getRe()<<endl;
	cout<<"\tSpawning potential ratio              = "<<cMSY.getSpr()<<endl;
	cout<<"\t---------------------------------------"<<endl;
    
	
	cout<<"fe\t ye \t dye \t re"<<endl;
	for(i=1;i<=100;i++)
	{
		fe = double(i-1.)/100;
		cMSY.calc_equilibrium(fe);
		cout<<fe<<"\t"<<cMSY.getYe()<<cMSY.getdYe()<<"\t"<<cMSY.getRe()<<endl;
	}
	
	
	cout<<"\t---------------------------------------"<<endl;
	ro     =  3598.1;
	h      =  0.62303;
	m      =  0.3;
	
	/* The data below are from the Flak Lake Model (catage demo in ADMB)*/
	dvector waa  =   ("{1.0121, 1.3858, 1.725, 2.0174, 2.2609, 2.4589, 2.6171}");
	dvector faa  =   ("{0.072924, 0.66935, 1.5841, 2.0026, 2.2595, 2.4588, 2.6171}");
	dvector tV   =   ("{0.0041512, 0.084884, 0.60249, 1.8459, 1.7443, 1.3592, 1.3592}");
	dmatrix VV(1,ngear,1,7);
	VV(1) = tV;
	cout<<"Flak Lake Trout"<<endl;
	Msy cFLK(ro,h,m,rho,waa,faa,VV);
	rho = 0.5;
	Msy cFLK5(ro,h,m,rho,waa,faa,VV);
	fe = 0.1;
	
	cout<<"fe\t ye1 \t ye2\t dye1 \t dye2\t re"<<endl;
	for(i=1;i<=100;i++)
	{
		fe = double(i-1.)/100;
		cFLK.calc_equilibrium(fe);;
		cFLK5.calc_equilibrium(fe);
		cout<<fe<<"\t"<<cFLK.getYe()<<cFLK5.getYe()<<cFLK.getdYe()<<cFLK5.getdYe()<<"\t"<<cFLK.getRe()<<endl;
	}
	
	
	
REPORT_SECTION

GLOBALS_SECTION
	#include <contrib.h>
	#include <msy.cpp>
