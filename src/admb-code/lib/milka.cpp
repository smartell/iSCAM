// milka.cpp
/**
 * @Milka Source code for Operating Model
 * @author Steven Martell & Catarina Wor
 * @details The default constuctort uses model_data as the base
 * class for the OperatingModel class.  
 * 
 * The OperatingModel class has a major function called 
 * runScenario:
 * 
 * STATUS LEGEND
 *  - : partially implemented
 *  + : implemented & testing
 *   : Good to go! 
 * 
 * runScenario:                                STATUS
 * 		|- readMSEcontrols                     [-]			
 * 		|- initParameters                      [-]			
 * 			|- surveyQ			               [ ]
 * 			|- stock-recruitment parameters	   [ ]	
 * 		|- initMemberVariables			       [-]
 * 		|- conditionReferenceModel			   [-]
 * 		|- setRandomVariables			       [-]
 * 		|- | getReferencePointsAndStockStatus  [-]		
 * 		   | calculateTAC                      [-]
 * 		   | allocateTAC                       [-]			
 * 		   | implementFisheries				   [-]
 * 		   		|- calcSelectivity			   [ ]
 * 		   		|- calcRetentionDiscards	   [ ]		
 * 		   		|- calcTotalMortality		   [ ]	
 * 		   | updateReferenceModel			   [ ]
 * 		   | writeDataFile					   [-]
 * 		   | runStockAssessment				   [ ]
 * 		|- |			
 * 		|- writeSimulationVariables			
 * 		|- calculatePerformanceMetrics			
 * 					
 */

#include <admodel.h>
#include "milka.h"
#include "baranov.h"
#include <contrib.h>

// Destructor
OperatingModel::~OperatingModel(){}

// Constructor
OperatingModel::OperatingModel(ModelVariables _mv,int argc,char * argv[])
:model_data(argc,argv), mv(_mv)
{
	cout<<"Inheritance version using model_data as base class"<<endl;
	cout<<"Ngroup "<<ngroup<<endl;
	cout<<"Catch Data\n"<<dCatchData<<endl;
	cout<<"d3 Survey Data\n"<<d3_survey_data<<endl;
	cout<<"eof "<<eof<<endl;

}


void OperatingModel::runScenario(const int &seed)
{
	readMSEcontrols();

	initParameters();

	initMemberVariables();

	conditionReferenceModel();

	setRandomVariables(seed);

	for(int i = nyr+1; i <= m_nPyr; i++ )
	{
		 
		getReferencePointsAndStockStatus();

		calculateTAC();
		
		allocateTAC(i);

		implementFisheries(i);

		updateReferenceModel();

		writeDataFile();

		runStockAssessment();
	}

	cout<<m_dCatchData<<endl;
}

/**
 * @brief Read control file for Management Strategy Evaluation.
 * @details Use cifstream to read in controls for MSE related options.
 * 
 */
void OperatingModel::readMSEcontrols()
{
	if(verbose) cout<<"MSE Control file\n"<<ProjControlFile<<endl;

	cifstream ifs(ProjControlFile);
	ifs>>m_nPyr;
	ifs>>m_nHCR;

	m_nGearIndex.allocate(1,ngear);
	m_nCSex.allocate(1,ngear);
	m_nASex.allocate(1,ngear);
	m_nAGopen.allocate(1,ngear,1,narea);
	
	// Controls for sexing catch and comps and fishing in given areas.
	imatrix tmp(1,ngear,-2,narea);
	ifs >> tmp;
	m_nGearIndex = column(tmp,-2);
	m_nCSex = column(tmp,-1);
	m_nASex = column(tmp,0);
	m_nAGopen = trans(trans(tmp).sub(1,narea));
	

}

/**
 * @brief Initialize model parameters based on model variable struct.
 * @details [long description]
 */
void OperatingModel::initParameters()
{
	
	// Initializing data members
	m_nNyr = nyr; // needs to be updated for each year inside the mse loop do we need this here??
	m_irow = nCtNobs; // counter for current number of rows in the catch table.

	// needs to be updated for each year in the mse loop
	int nn = 0;
	for( k = 1; k <= ngear; k++ )
	{
		nn += sum(m_nAGopen(k));
		nn += m_nCSex(k)*nn;
	}
	cout<<n_ags<< " \t"<<nn<<endl;
	m_nCtNobs = nCtNobs + (m_nPyr - nyr)*nn;
	
	m_dCatchData.allocate(1,m_nCtNobs,1,7);
	m_dCatchData.initialize();
	m_dCatchData.sub(1,nCtNobs) = dCatchData;

	//m_nItNobs = nItNobs;
	m_n_it_nobs.allocate(1,nItNobs);
	m_n_it_nobs.initialize();
	
	m_d3SurveyData.allocate(1,nItNobs,1,m_n_it_nobs,1,8);
	m_d3SurveyData.initialize();

	for(int k=1;k<=nItNobs;k++)
	{
		m_n_it_nobs(k) = n_it_nobs(k) + (m_nPyr-nyr);
		m_d3SurveyData(k).sub(1,n_it_nobs(k)) = d3_survey_data(k);	
	}
		
	m_n_A_nobs.allocate(1,nAgears);
	m_n_A_nobs.initialize();
	
	m_d3_A.allocate(1,nAgears,1,m_n_A_nobs,n_A_sage-5,n_A_nage);
	m_d3_A.initialize();
	
	for(int k=1;k<=nAgears;k++)
	{
		m_n_A_nobs(k) = n_A_nobs(k) + (m_nPyr-nyr);
		m_d3_A(k).sub(1,n_A_nobs(k)) = d3_A(k);	
	}
	 
	m_nWtNobs.allocate(1,nWtTab);
	m_nWtNobs.initialize();

	m_d3_inp_wt_avg.allocate(1,nWtTab,1,m_nWtNobs,sage-5,nage);
	m_d3_inp_wt_avg.initialize();

	for(int k=1;k<=nWtTab;k++)
	{
		m_nWtNobs(k)= nWtNobs(k)+(m_nPyr-nyr);
		m_d3_inp_wt_avg(k).sub(1,nWtNobs(k)) = d3_inp_wt_avg(k);	
	}
	

	// initializing population parameters
	m_dRo        = exp(mv.log_ro);
	m_dSteepness = mv.steepness;
	m_dM         = exp(mv.m);
	m_dRho       = mv.rho;
	m_dVarphi    = sqrt(1.0/mv.varphi);
	m_dSigma     = sqrt(m_dRho) * m_dVarphi;
	m_dTau       = sqrt(1.0-m_dRho)*m_dVarphi;

	for(int ih = 1; ih <= n_ag; ih++ )
	{
		m_dRbar  = exp(mv.log_rbar(ih));
		m_dRinit = exp(mv.log_rinit(ih));
	}

	switch(int(d_iscamCntrl(2)))
	{
		case 1:
			//Beverton-Holt model
			m_dKappa = elem_div(4.*m_dSteepness,(1.-m_dSteepness));
			break;
		case 2:
			//Ricker model
			m_dKappa = pow((5.*m_dSteepness),1.25);
		break;
	}
}


void OperatingModel::initMemberVariables()
{
	m_N.allocate(1,n_ags,syr,m_nPyr,sage,nage); m_N.initialize();
	m_M.allocate(1,n_ags,syr,m_nPyr,sage,nage); m_M.initialize();
	m_F.allocate(1,n_ags,syr,m_nPyr,sage,nage); m_F.initialize();
	m_Z.allocate(1,n_ags,syr,m_nPyr,sage,nage); m_Z.initialize();
	m_S.allocate(1,n_ags,syr,m_nPyr,sage,nage); m_S.initialize();
	m_ft.allocate(1,n_ags,1,ngear,syr,m_nPyr);  m_ft.initialize();
	m_d3_wt_avg.allocate(1,n_ags,syr,m_nPyr+1,sage,nage); m_d3_wt_avg.initialize();

	m_log_rt.allocate(1,n_ag,syr-nage+sage,nyr); m_log_rt.initialize();
	
	m_est_bo.allocate(1,ngroup);
	m_est_bmsy.allocate(1,ngroup);
	m_est_sbtt.allocate(1,ngroup);
	m_est_btt.allocate(1,ngroup);
	m_est_fmsy.allocate(1,ngroup,1,nfleet);
	m_est_msy.allocate(1,ngroup,1,nfleet);

	m_dTAC.allocate(1,ngroup,1,nfleet);

	// Initialize Mortality arrays from ModelVariables (mv)
	for(int ig = 1; ig <= n_ags; ig++ )
	{
		m_M(ig).sub(syr,nyr) = (*mv.d3_M)(ig);
		m_F(ig).sub(syr,nyr) = (*mv.d3_F)(ig);
		m_Z(ig).sub(syr,nyr) = m_M(ig).sub(syr,nyr) + m_F(ig).sub(syr,nyr);
		m_S(ig).sub(syr,nyr) = exp(-m_Z(ig).sub(syr,nyr));
		m_d3_wt_avg(ig).sub(syr,nyr+1) = d3_wt_avg(ig).sub(syr,nyr+1);

		// Temporary extend natural mortality out to m_nPyr
		for( i = nyr+1; i <= m_nPyr; i++ )
		{
			m_M(ig)(i) = m_M(ig)(nyr);
			m_d3_wt_avg(ig)(i+1) = d3_wt_avg(ig)(nyr+1);
		}
	}

	// Selectivity
	d4_logSel.allocate(1,ngear,1,n_ags,syr,m_nPyr,sage,nage);
	d4_logSel.initialize();
	for( k = 1; k <= ngear; k++ )
	{
		for(int ig = 1; ig <= n_ags; ig++ )
		{
			d4_logSel(k)(ig).sub(syr,nyr) = (*mv.d4_logSel)(k)(ig);

			// Temporarily extend selectivity out to m_nPyr
			for( i = nyr+1; i <= m_nPyr; i++ )
			{
				d4_logSel(k)(ig)(i) = d4_logSel(k)(ig)(nyr);
			}
		}
	}
}

void OperatingModel::conditionReferenceModel()
{
	int ig,ih;

	for( ig = 1; ig <= n_ags; ig++ )
	{
		f  = n_area(ig);
		g  = n_group(ig);
		ih = pntr_ag(f,g);

		dvector lx(sage,nage);
		dvector tr(sage,nage);
		lx(sage) = 1.0;
		for(j=sage;j< nage;j++)
		{
			lx(j+1) = lx(j) * exp( -m_M(ig)(syr)(j) );
		}
		lx(nage) /= (1.-exp(-m_M(ig)(syr,nage)));
		
		if( d_iscamCntrl(5) ) // initialize at unfished conditions.
		{
			tr =  log( m_dRo(g) ) + log(lx);
		}
		else if ( !d_iscamCntrl(5) )
		{
			tr(sage)        = ( mv.log_rbar(ih)+mv.log_rec_devs(ih)(syr));
			tr(sage+1,nage) = (mv.log_rinit(ih)+mv.init_log_rec_devs(ih));
			tr(sage+1,nage) = tr(sage+1,nage)+log(lx(sage+1,nage));
		}
		m_N(ig)(syr)(sage,nage) = 1./nsex * mfexp(tr);
		m_log_rt(ih)(syr-nage+sage,syr) = tr.shift(syr-nage+sage);

		for(i=syr;i<=nyr;i++)
		{
			if( i>syr )
			{
				m_log_rt(ih)(i) = (mv.log_rbar(ih)+mv.log_rec_devs(ih)(i));
				m_N(ig)(i,sage) = 1./nsex * mfexp( m_log_rt(ih)(i) );				
			}

			m_N(ig)(i+1)(sage+1,nage) =++elem_prod(m_N(ig)(i)(sage,nage-1)
			                                     ,m_S(ig)(i)(sage,nage-1));
			m_N(ig)(i+1,nage)        +=  m_N(ig)(i,nage)*m_S(ig)(i,nage);

			// average biomass for group in year i
			//bt(g)(i) += N(ig)(i) * d3_wt_avg(ig)(i);
		}
		m_N(ig)(nyr+1,sage) = 1./nsex * mfexp( mv.log_rbar(ih));
	}

}

void OperatingModel::setRandomVariables(const int& seed)
{
	m_nSeed = seed;
	random_number_generator rng(m_nSeed);

}


void OperatingModel::getReferencePointsAndStockStatus()
{
	// read iscam.res file to get this information.
	cifstream ifs("iSCAM.res");
	ifs >> m_est_bo;
	ifs >> m_est_fmsy;
	ifs >> m_est_msy;
	ifs >> m_est_bmsy;
	ifs >> m_est_sbtt;
	ifs >> m_est_btt;

	//cout<<m_est_btt<<endl;
}

/**
 * @brief Calculate the Total Allowable Catch
 * @details Total Allowable Catch is based on the estimates of current
 * stock statuts, reference points and the harvest control rule.
 * 
 * Uses a switch statement for HCR.  The HCR is set in the pcf file.
 */
void OperatingModel::calculateTAC()
{
	for( g = 1; g <= ngroup; g++ )
	{
		switch( int(m_nHCR) )
		{
			case 1: // Constant harvest rate
				 m_dTAC(g)  = (1.0-exp(-m_est_fmsy(g))) * m_est_btt(g);
				//m_dTAC(g) = 1.0;
			break; 
		}
	}
}


void OperatingModel::allocateTAC(const int& iyr)
{
	static int irow = nCtNobs;
	//m_dCatchdata(year,gear,area,group,sex,type,value)
	int h;
	for( k = 1; k <= nfleet; k++ )
	{
		h = m_nCSex(k);
		for( f = 1; f <= narea; f++ )
		{
			if(m_nAGopen(k,f))
			{ 
			for( g = 1; g <= ngroup; g++ )
			{
				if(!h)
				{
					irow ++;
					m_dCatchData(irow,1) = iyr;
					m_dCatchData(irow,2) = nFleetIndex(k);
					m_dCatchData(irow,3) = f;
					m_dCatchData(irow,4) = g;
					m_dCatchData(irow,5) = h;
					m_dCatchData(irow,6) = 1;  //TODO: Fix this catch type
					m_dCatchData(irow,7) = m_dTAC(g)(k);  // TODO: call a manager!
				}
				if(h)
				{	
					for( h = 1; h <= nsex; h++ )
					{
						irow ++;
						m_dCatchData(irow,1) = iyr;
						m_dCatchData(irow,2) = nFleetIndex(k);
						m_dCatchData(irow,3) = f;
						m_dCatchData(irow,4) = g;
						m_dCatchData(irow,5) = h;
						m_dCatchData(irow,6) = 1;  //TODO: Fix this
						m_dCatchData(irow,7) = m_dTAC(g)(k);  // TODO: call a manager!
					}
				}
			}
			}
		}
	}

	

	
}

/**
 * @brief Implement spatially explicity fishery.
 * @details Implement the spatially epxlicity fishery using the Baranov catch equation
 * to determine the instantaneous fishing mortality rate in each area by each gear. This
 * routine uses the BaranovCatchEquation class object to do this.
 * 
 * Notes:
 * 	m_dTAC is a vector of allocated catches assiged to each fleet.
 * 	
 * 	Algorithm:
 * 	|- Apportion m_dTAC by area (f) for each stock (g)
 * 	|- Loop over each area and allocate catch in area (f) to gear (k),
 * 	|- Add implementation error to each gear and catch.
 * 	|- Assemble arguments for BaranovCatchEquation Class. 
 * 	     -> .getFishingMortality(ct,ma,&Va,na,wa,_hCt)
 * 	|- Calculate Fishing mortality rates on reference population.
 * 	|- Calculate Total landed catch.
 * 	|- Calculate total discards based on size-limits.
 * 	|- Calculate total discards from non-retention fisheries.
 * 	
 */
void OperatingModel::implementFisheries(const int &iyr)
{
	dvector tac(1,narea);
	dvector  ct(1,nfleet);
	dmatrix  ma(1,nsex,sage,nage);
	dmatrix  na(1,nsex,sage,nage);
	dmatrix  wa(1,nsex,sage,nage);
	dmatrix  d_allocation(1,narea,1,nfleet);
	dmatrix  _hCt(1,nsex,1,nfleet);
	d3_array d3_Va(1,nsex,1,nfleet,sage,nage);
	tac.initialize();
	na.initialize();

	BaranovCatchEquation cBCE;

	for(int f = 1; f <= narea; f++ )
	{
		for(int g = 1; g <= ngroup; g++ )
		{
			ct = m_dTAC(g);  // Catch for each fleet.
			for(int h = 1; h <= nsex; h++ )
			{
				int ig = pntr_ags(f,g,h);
				ma(h) = m_M(ig)(iyr);			// natural mortality
				na(h) = m_N(ig)(iyr);			// numbers-at-age
				wa(h) = m_d3_wt_avg(ig)(iyr);	// weight-at-age
				for(int k = 1; k <= nfleet; k++ )
				{
					int kk = nFleetIndex(k);
					d3_Va(h)(k) = exp(d4_logSel(kk)(ig)(iyr));
				}
			}  // nsex
			cout<<"Start"<<endl;
			cout<<d3_Va(1)<<endl;
			// Calculate instantaneous fishing mortality rates.
			dvector ft = cBCE.getFishingMortality(ct,ma,&d3_Va,na,wa,_hCt);
			cout<<"Ft =\t"<<ft<<endl;

			// Fill m_dCatchData array with actual catches taken by each fleet.
			for(int k = 1; k <= nfleet; k++ )
			{
				if( ft(k) > 0 )
				{
					int kk = nFleetIndex(k);
					int hh = m_nCSex(k);   // flag for sex
					for( h = 1; h <= hh+1; h++ )
					{
						m_irow ++;
						m_dCatchData(m_irow,1) = iyr;
						m_dCatchData(m_irow,2) = kk;
						m_dCatchData(m_irow,3) = f;
						m_dCatchData(m_irow,4) = g;
						m_dCatchData(m_irow,5) = hh>0?h:0;
						m_dCatchData(m_irow,6) = 1;  //TODO: Fix this
						m_dCatchData(m_irow,7) = hh>0?_hCt(h,k):colsum(_hCt)(k);
					}
				}
			}

		}  // ngroup g
	} // narea f
	cout<<m_dCatchData<<endl;
	cout<<"END"<<endl;
	ad_exit(1);

}

void OperatingModel::updateReferenceModel()
{

 
}

void OperatingModel::writeDataFile()
{

		adstring sim_datafile_name = "Simulated_Data_"+str(rseed)+".dat";
	  	ofstream dfs(sim_datafile_name);
	  	dfs<<"#Model dimensions"<<endl;
	  	dfs<< narea 		<<endl;
	  	dfs<< ngroup		<<endl;
	  	dfs<< nsex			<<endl;
	  	dfs<< syr   		<<endl;
	  	dfs<< i   			<<endl;
	  	dfs<< sage  		<<endl;
	  	dfs<< nage  		<<endl;
	  	dfs<< ngear 		<<endl;
	 
	  	dfs<<"#Allocation"	<<endl;
	  	dfs<< dAllocation 	<<endl;
	  	
	  	dfs<<"#Age-schedule and population parameters"<<endl;
	  	dfs<< d_linf  			<<endl;
	  	dfs<< d_vonbk  			<<endl;
	  	dfs<< d_to  			<<endl;
	  	dfs<< d_a  				<<endl;
	  	dfs<< d_b  				<<endl;
	  	dfs<< d_ah  			<<endl;
	  	dfs<< d_gh  			<<endl;
	  	dfs<< n_MAT				<<endl;
		dfs<< d_maturityVector  <<endl;
	
	  	dfs<<"#Observed catch data"<<endl;
	  	
	  		int nn = 0; // calculating nn again make initParameters return it???? 
			for( k = 1; k <= ngear; k++ )
			{
				nn += sum(m_nAGopen(k));
				nn += m_nCSex(k)*nn;
			}

	  	dfs<< m_nCtNobs + (i-nyr)*nn  		<<endl; 
	  	dfs<< m_dCatchData.sub(1,nCtNobs+(i-nyr)*nn)    <<endl;
	
	  	dfs<<"#Abundance indices"	<<endl;
		
	  	ivector tmp_n_it_nobs(1,nItNobs);
	  	d3_array tmp_d3SurveyData(1,nItNobs,1,tmp_n_it_nobs,1,8);

	  		for(int k=1;k<=nItNobs;k++)
			{
				tmp_n_it_nobs(k) = n_it_nobs(k) + (i-nyr);
				tmp_d3SurveyData(k) = m_d3SurveyData(k).sub(1,tmp_n_it_nobs(k));
			}
	  	
	  	dfs<< nItNobs 					<<endl;
	  	dfs<< tmp_n_it_nobs 				<<endl;
	  	dfs<< n_survey_type 			<<endl;
	  	dfs<< tmp_d3SurveyData		<<endl;
	
	  	dfs<<"#Age composition"		<<endl;

	  		ivector tmp_n_A_nobs(1,nAgears);
	  		d3_array tmp_d3_A(1,nAgears,1,tmp_n_A_nobs,n_A_sage-5,n_A_nage);
	  			
	  		for(int k=1;k<=nAgears;k++)
			{
				tmp_n_A_nobs(k) = n_A_nobs(k) + (i-nyr);
				tmp_d3_A(k) = m_d3_A(k).sub(1,tmp_n_A_nobs(k));	
			}

	  	dfs<< nAgears				<<endl;
	  	dfs<< tmp_n_A_nobs			<<endl;
	  	dfs<< n_A_sage				<<endl;
	  	dfs<< n_A_nage				<<endl;
	  	dfs<< inp_nscaler 			<<endl;
	  	dfs<< tmp_d3_A				<<endl;
	
	  	dfs<<"#Empirical weight-at-age data"	<<endl;

	  	ivector tmp_nWtNobs(1,nAgears);
	  	d3_array tmp_d3_inp_wt_avg(1,nWtTab,1,tmp_nWtNobs,sage-5,nage);

	  		for(int k=1;k<=nWtTab;k++)
			{
				tmp_nWtNobs(k)= nWtNobs(k)+(i-nyr);
				tmp_d3_inp_wt_avg(k)= m_d3_inp_wt_avg(k).sub(1,tmp_nWtNobs(k)) ;	
			}

	  	dfs<< nWtTab 				<<endl;
	  	dfs<< tmp_nWtNobs				<<endl;
		dfs<< tmp_d3_inp_wt_avg			<<endl; 
	
		dfs<<"#EOF"	<<endl;
		dfs<< 999	<<endl;

}

void OperatingModel::runStockAssessment()
{

}
