// |---------------------------------------------------------------------------------|
// | GENERAL OPERATING MODEL FOR USE WITH iSCAM.                                     |
// |---------------------------------------------------------------------------------|



#include <admodel.h>
#include "OpMod.h"
#include "Selex.h"
#include "msy.h"
#include "baranov.h"

double getSurveyQ(const dvector& bt, const dvector& it);
/** \brief Destructor
	\author Steven Martell 
**/
OperatingModel::~OperatingModel()
{

};

/** \brief Constructor for operating model

	The Operating Model class is derived from model_data,
	this is done intentionally such that if new data are added to 
	iscam, then these data are readily available for the operating
	model.

	\author Steven Martell
**/
OperatingModel::OperatingModel(	const s_iSCAMdata&  mse_data, 
                               	const s_iSCAMvariables& mse_vars,
								int argc,
								char * argv[] )
								: model_data(argc,argv)
{
	initializeConstants(mse_data);  // deprecate, 
	initializeVariables(mse_vars);

	cout<<" testing the inheritance of model_data"<<endl;
	cout<<"narea\t"<<narea<<endl;
	cout<<"pntr_ags\n"<<pntr_ags<<endl;
	// cout<<"catch_array\n"<<catch_array<<endl;
	// exit(1);
}

/** \brief Initialize constants in Operating model

	Initialize private member variables based on scenario class 
	
	\author Steven Martell
**/
void OperatingModel::initializeConstants(const s_iSCAMdata& cS)
{
	// |------------------|
	// | Model dimensions |
	// |------------------|
	nStock = ngroup;       //cS.nStock;  
	nArea  = narea;        //cS.nArea;
	nSex   = nsex;         //cS.nSex;
	nSyr   = syr;          //cS.nSyr;
	nNyr   = nyr;          //cS.nNyr;
	nPyr   = mse_cntrl(1); //cS.nPyr;
	nSage  = sage;         //cS.nSage;
	nNage  = nage;         //cS.nNage;
	nGear  = ngear;        //cS.nGear;
	nFleet = nfleet;       //cS.nFleet;
	dAge   = age;

	// |----------------------------------------------------------------------|
	// | Allocate arrays for operating model output and assessment data files |
	// |----------------------------------------------------------------------|
	// | Catch data
	int nCount = (nPyr-nNyr+1);
	m_nCtNobs_counter = nCtNobs;

	m_nCtNobs  = nCtNobs + nCount*nSex*nGear;  /**< Length of catch array */

	m_dCatchData.allocate(1,m_nCtNobs,1,7);
	m_dCatchData.initialize();
	m_dCatchData.sub(1,nCtNobs) = dCatchData;
	m_catch_sex_composition = catch_sex_composition;
	m_catch_type = catch_type;

	// | Survey data
	m_n_it_nobs.allocate(1,nItNobs);
	m_n_it_nobs = n_it_nobs + nCount*nItNobs;
	m_n_it_counter = n_it_nobs;

	m_d3_survey_data.allocate(1,nItNobs,1,m_n_it_nobs,1,8);
	m_d3_survey_data.initialize();
	for(int k = 1; k <= nItNobs; k++ )
	{
		m_d3_survey_data(k).sub(1,n_it_nobs(k)) = d3_survey_data(k);
	}

	// | Age composition
	m_n_A_nobs.allocate(1,nAgears);
	m_n_A_nobs = n_A_nobs + nCount*nAgears;

	m_d3_A.allocate(1,nAgears,1,m_n_A_nobs,n_A_sage-5,n_A_nage);
	m_d3_A.initialize();
	for(int k = 1; k <= nAgears; k++ )
	{
		m_d3_A(k).sub(1,n_A_nobs(k)) = d3_A(k);
	}

	// |Empirical weight-at-age data (-99 = NA)
	m_nWtNobs = nWtNobs + nCount*nSex;

	m_imp_wt_avg.allocate(1,m_nWtNobs,sage-5,nage);
	m_imp_wt_avg.sub(1,nWtNobs) = inp_wt_avg;

	m_d3_wt_avg.allocate(1,n_ags,syr,nPyr+1,sage,nage);
	m_d3_wt_mat.allocate(1,n_ags,syr,nPyr+1,sage,nage);
	m_d3_wt_avg.initialize();
	m_d3_wt_mat.initialize();
	for(int ig = 1; ig <= n_ags; ig++ )
	{
		m_d3_wt_avg(ig).sub(syr,nyr+1) = d3_wt_avg(ig);
		m_d3_wt_mat(ig).sub(syr,nyr+1) = d3_wt_mat(ig);
	}

	
	cout<<"initializeConstants"<<endl;
	// cout<<"pntr_ags\n"<<pntr_ags<<endl;

}


void OperatingModel::initializeVariables(const s_iSCAMvariables& cS)
{
	// Leading parameters
	cout<<"initializeVariables"<<endl;
	dRo                = exp(cS.d_log_ro);
	dSteepness         = cS.d_steepness;
	dM                 = exp(cS.d_log_m);
	dAvgRec            = exp(cS.d_log_rbar);
	dInitRec           = exp(cS.d_log_rinit);
	dRho               = cS.d_rho;
	dVartheta          = sqrt(1.0/cS.d_varphi);
	dLog_m_devs        = cS.dLog_M_devs;
	dLog_rbar_devs     = cS.dLog_Rbar_devs;
	dLog_init_rec_devs = cS.dLog_Rinit_devs;
	
	nSel_type  = cS.nSel_type;
	nSel_block = cS.nSel_block;
	d3_selPars.allocate(*cS.dSelPars);
	d3_selPars = *cS.dSelPars;
	
	d4_log_sel.allocate(*cS.d4_log_sel);
	d4_log_sel = *cS.d4_log_sel;

	
	d3_Mt.allocate(*cS.d3_Mt);
	d3_Mt = *cS.d3_Mt;

	d3_St.allocate(*cS.d3_St);
	d3_St = *cS.d3_St;

	d3_Nt.allocate(1,n_ags,nSyr,nPyr+1,nSage,nNage);
	d3_Nt.initialize();

	d3_Zt.allocate(1,n_ags,nSyr,nPyr,nSage,nNage);
	d3_Zt.initialize();

	// Fishing mortality rates
	d3_Ft.allocate(*cS.d3_Ft);
	d3_Ft = *cS.d3_Ft;
	

	// cout<<dFt<<endl;
	// recruitment compensation
	m_kappa.allocate(1,nStock);
	m_so.allocate(1,nStock);
	m_beta.allocate(1,nStock);
	m_dSbo.allocate(1,nStock);
	m_kappa = elem_div(4.0 * dSteepness, 1. - dSteepness);

	// survey catchability
	m_survey_q.allocate(1,nItNobs);
	m_survey_q.initialize();
}

void OperatingModel::runScenario(const int &seed)
{
	cout<<"Top of runScenario"<<endl;
	/*
	PSUEDOCODE:
		- declare local variables --> done in the constructor.
		- initialize stock-recruitment parameters
		- condition operating model on historical assessment
		|-- calculate reference points
		| - calculate harvest control rule
		| - implement harvest on reference population
		| - update reference population
		| - generate data for annual assessment model
		| - run stock-assessment procedures
		|-- repeat:
		- write simulation variables to output file.
		- write performance statistics to output file.
	*/

	/*
	- Local variables
	*/
	m_nSeed = seed;

	/*
	- Initialize stock-recruitment parameters for each stock.
	- Survivorship and fecundity is based on average mortality & weight-at-age
	- (so, beta, sbo, ro)
	*/
	calcStockRecruitment();
	cout<<"Ok apres calcStockRecruitment.   \t tout bon"<<endl;

	/*
	- Condition operating model on historical assessment
	- Run model from syr-nyr based on input parameters.
	*/
	conditionReferenceModel();
	cout<<"Ok apres conditionReferenceModel. \t pas fini"<<endl;

	/*
	- Calculate survey catchability coefficients based on data & operating model
	*/
	calcSurveyCatchability();
	cout<<"Ok apres calcSurveyCatchability.   \t pas fini"<<endl;

	cout<<nNyr<<"\t"<<nPyr<<endl;
	for(int i = nNyr + 1; i <= nPyr; i++ )
	{
		cout<<"year = "<<i<<endl;
		cout<<d3_Nt(1)(i)<<endl;
		/*
		- Calculate reference points that are required for the harvest control rule.
		*/
		calcReferencePoints();
		cout<<"Ok apres clacReferencePoints.      \t pas fini"<<endl;

		/*
		- Use Harvest Control Rule to calculate TAC
		*/
		calcTAC();
		cout<<"Tac = "<<m_dTac<<endl;
		cout<<"OK apres calcTAC.              	\t pas fini"<<endl;

		/*
		- Implement harvest on reference population.
		- Ensure that actual catch does not exceed available biomass
		- Implementation errors occur here.
		*/
		implementFisheries(i);
		cout<<"Ok apres implementFisheries       \t presque fini"<<endl;

		/*
		- Run Fisheries independent survey.
		*/
		calcRelativeAbundance(i);
		cout<<"Ok apres calcRelativeAbundance     \t pas fini"<<endl;

		/*
		- Update reference population
		*/	
		updateReferencePopulation(i);
		cout<<"Ok apres updateReferencePopulation \t pas fini"<<endl;


		/*
		- Growth
		*/
		calcGrowth(i);
		cout<<"Ok apres calcGrowth                \t pas fini"<<endl;

		/*
		- Movement
		*/

		/*
		- Generate data for stock assessment.
		*/
		generateStockAssessmentData(i);
		cout<<"Ok apres generateStockAssessmentData \t pas fini"<<endl;


		/*
		- Run stock assessment procedures.
		*/
		runStockAssessment();
		cout<<"Ok apres runStockAssessment          \t pas fini"<<endl;

	}
}

void OperatingModel::runStockAssessment()
{

}

/** \brief generateStockAssessmentData
	
		This routine writes the data file for iSCAM.dat based on the input variables
		upto iyr.
	
	\author  Steve Martell
	\date Sept 4, 2013
	\param  iyr terminal year of assessment.
	
	\return null
	\sa
**/
void OperatingModel::generateStockAssessmentData(const int& iyr)
{
	adstring sim_datafile_name = "Simulated_Data_"+str(m_nSeed)+".dat";
  	ofstream dfs(sim_datafile_name);
  	dfs<<"#Model dimensions"<<endl;
  	dfs<< nArea 		    <<endl;
  	dfs<< nStock		    <<endl;
  	dfs<< nSex			    <<endl;
  	dfs<< nSyr   		    <<endl;
  	dfs<< iyr   		    <<endl;
  	dfs<< nSage  		    <<endl;
  	dfs<< nNage  		    <<endl;
  	dfs<< nGear 		    <<endl;
 
  	dfs<<"#Allocation"	<<endl;
  	dfs<< dAllocation 	<<endl;
  	

  	dfs<<"#Age-schedule and population parameters"<<endl;
  	dfs<< d_linf  		                          <<endl;
	dfs<< d_vonbk  		                          <<endl;
	dfs<< d_to  		                          <<endl;
	dfs<< d_a  			                          <<endl;
	dfs<< d_b  			                          <<endl;
	dfs<< d_ah  		                          <<endl;
	dfs<< d_gh  		                          <<endl;

  	dfs<<"#Observed catch data"                 <<endl;
  	dfs<< m_nCtNobs_counter                     <<endl;
  	dfs<< m_dCatchData.sub(1,m_nCtNobs_counter) <<endl;  

  	dfs<<"#Abundance indices"	<<endl;
  	dfs<< nItNobs 				<<endl;
  	dfs<< n_it_nobs 			<<endl;
  	dfs<< n_survey_type 		<<endl;
  	dfs<< d3_survey_data 			<<endl;

 //  	dfs<<"#Age composition"		<<endl;
 //  	dfs<< na_gears				<<endl;
 //  	dfs<< na_nobs				<<endl;
 //  	dfs<< a_sage				<<endl;
 //  	dfs<< a_nage				<<endl;
 //  	dfs<< A						<<endl;

 //  	dfs<<"#Empirical weight-at-age data"	<<endl;
 //  	dfs<< n_wt_nobs				<<endl;
	// dfs<< inp_wt_avg			<<endl;

	// dfs<<"#EOF"	<<endl;
	// dfs<< 999	<<endl;
	
	// | END OF WRITING SIMULATED DATAFILE.

}

/** \brief Calculate survey catchabiltiy for each relative abundance index.
	
		Use the observed abundance indicies and the conditioned operating model
		to determine the appropriate scaler.
	
	\author  Steve Martell
	\date `date +%Y-%m-%d`
	\param  description of parameter
	\param  description of parameter
	\return description of return value
	\sa calcRelativeAbundance
**/
void OperatingModel::calcSurveyCatchability()
{
	int i,k,f,g,h,ig;
	double di;

	for(int kk = 1; kk <= nItNobs; kk++ )
	{
		// | Vulnerable numbers-at-age to survey
		dmatrix V(1,n_it_nobs(kk),sage,nage);
		V.initialize();

		for(int ii = 1; ii <= n_it_nobs(kk); ii++ )
		{
			i    = d3_survey_data(kk)(ii)(1);  // year
			k    = d3_survey_data(kk)(ii)(3);  // gear
			f    = d3_survey_data(kk)(ii)(4);  // area
			g    = d3_survey_data(kk)(ii)(5);  // stock
			h    = d3_survey_data(kk)(ii)(6);  // sex
			di   = d3_survey_data(kk)(ii)(8);  // timing of survey

			dvector na(sage,nage);
			dvector sa(sage,nage);
			dvector va(sage,nage);
			dvector wa(sage,nage);
			dvector ma(sage,nage);
			na.initialize();
			for( h = 1; h <= nsex; h++ )
			{
				ig = pntr_ags(f,g,h);
				va = exp( d4_log_sel(k)(ig)(i) );
				sa = mfexp( -d3_Zt(ig)(i)*di );
				na = elem_prod(d3_Nt(ig)(i),sa);
				wa = d3_wt_avg(ig)(i);
				ma = d3_wt_mat(ig)(i);
				switch(n_survey_type(kk))
				{
					case 1:
						V(ii) += elem_prod(na,va);
					break;

					case 2:
						V(ii) += elem_prod(elem_prod(na,va),wa);
					break;

					case 3:
						V(ii) += elem_prod(elem_prod(na,va),ma);
					break;
				}
			}
		}
		dvector it  = trans(d3_survey_data(kk))(2); //relative abundance index
		dvector wt  = trans(d3_survey_data(kk))(7); //relative weight (multiplier)
		        wt  = wt / sum(wt);
		dvector t1  = rowsum(V);
		dvector zt  = log(it) - log(t1);
		double zbar = zt * wt;
		m_survey_q(kk) = exp(zbar);
		
	}  // end of nItNobs
	
}


/** \brief Calculate and append relative abundance index.
	
		This routine requires calcSurveyCatchabililty to determine the 
		appropriate scaler for each of the relativen abundance indices.
	
	\author  Steve Martell
	\date `date +%Y-%m-%d`
	\param  description of parameter
	\param  description of parameter
	\return description of return value
	\sa calcSurveyCatchability
**/
void OperatingModel::calcRelativeAbundance(const int& iyr)
{
	// Update m_d3_survey_data
	// m_d3_survey data is a ragged object based on m_n_it_nobs

	
	m_n_it_counter = m_n_it_counter + 1;
	for(int kk = 1; kk <= nItNobs; kk++ )
	{
		// m_d3_survey_data
	}
	// cout<<"q = "<<q<<endl;
	
	

}


void OperatingModel::calcGrowth(const int& iyr)
{
	// Implement alternative growth models.
	for(int  ig = 1; ig <= n_ags; ig++ )
	{
		m_d3_wt_avg(ig)(iyr+1) = d3_wt_avg(ig)(nNyr+1);
		m_d3_wt_mat(ig)(iyr+1) = d3_wt_mat(ig)(nNyr+1);
	}
}

/**
 * @brief Update the reference population
 *
 	\todo movement transition matrix needs to be implemented here.
 * 
 * @param iyr index for the current assesment/projection year
 *
 * @return null
 */
void OperatingModel::updateReferencePopulation(const int& iyr)
{
	int f,g;//,h,i,j,k;
	int ig,ih;

	
	for( ig = 1; ig <= n_ags; ig++ )
	{	
		f  = n_area(ig);
		g  = n_group(ig);
		ih = pntr_ag(f,g);

		// cout<<dAvgRec(ih)<<endl;
		// Current recruitment is based on historical average.
		// Need to add recruitment variation here, environmental effects, 
		// and stock-recruitment relationship.

		d3_Nt(ig)(iyr+1,nSage) = dAvgRec(ih);
		

		// Update numbers-at-age.
		dvector st = exp(-d3_Zt(ig)(iyr)(nSage,nNage-1));
		d3_Nt(ig)(iyr+1)(nSage+1,nNage) =++ elem_prod(d3_Nt(ig)(iyr)(nSage,nNage-1),st);
		d3_Nt(ig)(iyr+1,nNage)     += d3_Nt(ig)(iyr,nNage) * exp(-d3_Zt(ig)(iyr,nNage));
	}
	cout<<"Got Here"<<endl;
}

/**
 @brief Calculate the approprate area based TAC based on 
        apportioned biomass, and application of harvest rates.

 @author Steve Martell
 @date   August 13, 2013
 @location  St Paul Island 

Procedure currently used by the IPHC
	- 1. Estimate coast-wide biomass.
	- 2. Apportion biomass in each regulatory area based on survey.
	- 3. Apply area-specific harvest rates to apportioned biomass.
		- a. If stock status > 30% Bo use full Frate
		- b. Frate ramps down to 0 at limit reference point (20% Bo).
		- c. Frate=0 at or below limit reference point.
	- 4. Subtract off estimated wastage & bycatch.
	- 5. For area 2B use catch allocation

*/
void OperatingModel::calcTAC()
{
	cout<<"Fmsy\t"<<m_dFmsy<<endl;	

	/* Get stock assessment results -> Ian Stewart */
	m_est_bo.allocate(1,nStock);
	m_est_fmsy.allocate(1,nFleet);
	m_est_msy.allocate(1,nFleet);
	m_est_bmsy.allocate(1,nStock);
	m_est_sbt.allocate(1,nStock);
	m_est_bt.allocate(1,nStock);

	ifstream ifs("iSCAM.res");
	ifs >> m_est_bo;
	ifs >> m_est_fmsy;
	ifs >> m_est_msy ;
	ifs >> m_est_bmsy;
	ifs >> m_est_sbt ;
	ifs >> m_est_bt;

	cout<<nStock<<endl;
	cout<<nFleet<<endl;
	cout<<m_est_bo<<endl;
	cout<<m_est_bt<<endl;
	/* HARVEST CONTROL RULE */
	/* Compute apportionment schedule -> Ray Webster */

	// SJDM Simple 20% hr for now.
	m_dTac.allocate(1,nStock);
	m_dTac.initialize();
	int g;
	for( g = 1; g <= nStock; g++ )
	{
		m_dTac(g) = 0.2 * m_est_bt(g);
	}


}



void OperatingModel::implementFisheries(const int& iyr)
{
	/*
		Author Steve Martell
		Notes:
			- This is going to be a bit tricky due to the number of factors at play.
			- Essentially this routine must calculate the F-at-age for each gear that
			  is given a TAC > 0.  Also concerned with bycatch fisheries (the effort
			  in bycatch fisheries is likely to be independent of halibut stock size).

			- For the directed fishery, the TAC value includes wastage, so must also 
			  keep track of landed fish, U32 discards, and lost gear.  Use the same model
			  that is described by Lew Coggins, and the work I did with Bill Pine where
			  there is a joint probability of capture * (probability of retention + 
			  (1-probability of retention)*discard mortality) that is a function of size-
			  at-age and the legal size limit.

			- Aug 27, for now keep it simple  

			- Sep 4, update the catch data array for writing to simulated data file.

			- Sep 10, 2013. Changed the logic in this routine.  Implemented as follows:
				- 1) Apportion m_dTac by area (f) from each stock (g)
				- 2) Loop over each area and allocate catch in area (f) to gear (k),
				- 3) Assemble arguments for BaranovCatchEquation class
					-> catch by gear,
					-> natural mortality rate by sex(row) and age(col)
					-> Selectivity array sex,gear,age
					-> Numbers by sex(row) and age (col)
					-> Weight by sex(row) and age (col)
				- 4) Calculate Fishing Mortality rates using class BaranovCatchEquation
				- 5) Calculate Total mortality rate (d3_Zt) 
	*/
	 /* 
	 Given a tac in year iyr, figure out what the F is.
	 getFishingMortality in Baranov.h arguments are:
	 	- ct for each fleet
	 	- Natural mortality by age
	 	- selectivity for each fleet
	 	- numbers -at-age 
	 	- [weight -at-age if catch is in biomass units]
	 */

	int f,g,h,k;
	int ig,kk,hh;

	
	dvector tac(1,nArea);
	dvector  ct(1,nFleet);
	dmatrix  ma(1,nSex,nSage,nNage);
	dmatrix  na(1,nSex,nSage,nNage);
	dmatrix  wa(1,nSex,nSage,nNage);
	d3_array d_Va(1,nSex,1,nFleet,nSage,nNage);
	dmatrix d_alloc(1,nArea,1,nFleet);  // Allocation for each fleet in each area.
	tac.initialize();
	na.initialize();
	for( f = 1; f <= nArea; f++ )
	{
		// -1) Apportion catch to each area (f) from each stock (g)
		for( g = 1; g <= nStock; g++ )
		{
		 	tac(f) += m_dTac(g);
		}

		// -2) Allocate catch in each area (f) to gear (k).
		//     -[ ] TODO: Will need to add area specific rules here.
		//     -[ ] TODO: add implementation error here.
		for( k = 1; k <= nFleet; k++ )
		{
			 d_alloc(f,k) = dAllocation(k);
			 ct(k)        = d_alloc(f,k)*tac(f);
		}

		// -3) Assemble arguments for BarnovCatchEquation class.
		// [ ] TODO: allow for time-varying M in future
		// [ ] TODO: allow for Selectivity to change in future.
		// [ ]
		m_dFt.allocate(1,nStock,1,nFleet);
		dmatrix _hCt(1,nSex,1,nFleet);
		BaranovCatchEquation cBCE;
		for( g = 1; g <= nStock; g++ )
		{
			for( h = 1; h <= nSex; h++ )
			{
				ig    = pntr_ags(f,g,h);
				ma(h) = d3_Mt(ig)(nNyr);
				for( k = 1; k <= nFleet; k++ )
				{
					kk         = nFleetIndex(k);
					d_Va(h)(k) = exp(d4_log_sel(kk)(ig)(nNyr));
				}
				na(h)+= d3_Nt(ig)(iyr);
				wa(h) = m_d3_wt_avg(ig)(iyr);
			}

			// 4) - calculate fishing mortality rate based on baranov catch eqn.
			// Potential issue here if nStock > 1, what wa ma vector should be used?
			dvector ft = cBCE.getFishingMortality(ct,ma,&d_Va,na,wa,_hCt);
			cout<<"ft = "<<ft<<endl;
			cout<<"hCt = "<<_hCt<<endl;
			m_dFt(g) = ft;

			// -4) Fill catch data array
			for( k = 1; k <= nFleet; k++ )
			{
				if( ft(k)>0 )
				{
					kk = nFleetIndex(k);
					// determine if the fleet has sex-specific or aggregated catches
			 		hh = m_catch_sex_composition(k);
			 		int nn = hh>0?hh:1;
			 		for(h=1; h<=nn; h++)
			 		{
						m_nCtNobs_counter ++;
						m_dCatchData(m_nCtNobs_counter)(1) = iyr;
				 		m_dCatchData(m_nCtNobs_counter)(2) = kk;
				 		m_dCatchData(m_nCtNobs_counter)(3) = f;
				 		m_dCatchData(m_nCtNobs_counter)(4) = g;
				 		m_dCatchData(m_nCtNobs_counter)(5) = hh>0?h:0;	
				 		m_dCatchData(m_nCtNobs_counter)(6) = m_catch_type(k);
				 		m_dCatchData(m_nCtNobs_counter)(7) = hh>0?_hCt(h,k):colsum(_hCt)(k);
			 		}
				}
			}
		}  // end of stock loop
	}



	/* DEPRECATE THE CODE BELOW, See Sept 10 Note.*/

	/* Get selectivities and allocations for each fleet. */
	// dvector d_ak(1,nFleet);
	// d3_array d_V(1,n_ags,1,nFleet,nSage,nNage);
	// for( k = 1; k <= nFleet; k++ )
	// {
	// 	kk = nFleetIndex(k);
	// 	d_ak(k) = dAllocation(k);
	// 	for( ig = 1; ig <= n_ags; ig++ )
	// 	{
	// 		d_V(ig)(k) = exp(d4_log_sel(kk)(ig)(nNyr));
	// 	}

	// }

	/* Calculate fleet specific fishing mortality rates by group */
	// TO DO
	// [ ] - add joint capture probability for size-based selectivity and 32" size limit.
	// [ ] - add time-varying natural mortality rates.
	// [ ] - add implementation error.
	// BaranovCatchEquation cBaranov;
	// for( ig = 1; ig <= n_ags; ig++ )
	// {

	// 	// dvector na = d3_Nt(ig)(iyr);
	// 	// dvector wa = d3_wt_avg(ig)(iyr);
	// 	// dvector ma = d3_Mt(ig)(nNyr);
	// 	// dmatrix va = d_V(ig);

	// 	// // Add implementation error here.
	// 	// // m_dTac comes from the harvest control rule.
	// 	// // Also need to record the actual catch in m_dCatchData for data file.
	// 	// dvector ct = m_dTac;

	// 	// // SM Working here end of Sept 4. 2013


	// 	// // cout<<"Na\t"<<na<<endl;
	// 	// // cout<<"wa\t"<<wa<<endl;
	// 	// // cout<<"ma\t"<<ma<<endl;
	// 	// // cout<<"va\t"<<va<<endl;
	// 	// // cout<<"ct\t"<<ct<<endl;

	// 	// dvector ft = cBaranov.getFishingMortality(ct,ma,va,na,wa);
	// 	// m_dFt = ft;

	// 	// Total mortality
	// 	d3_Zt(ig)(iyr) = ma;
	// 	for( k = 1; k <= nFleet; k++ )
	// 	{
	// 		d3_Zt(ig)(iyr) += m_dFt(k) * va(k);
	// 	}
	// 	// cout<<"Zt\t"<<d3_Zt(1)(iyr)<<endl;
	// }

}

/**
* @brief calculate reference points.
* @author Steve Martell
*/
void OperatingModel::calcReferencePoints()
{
	/*
	TO DO:
		1) add an const int iyr argument to determine average weights
		at age and maturity at age upto the current assessment year.
	*/

	int g,k,kk;
	int ig;
	double d_rho = d_iscamCntrl(13);
	dvector d_ak(1,nFleet);
	d3_array d_V(1,n_ags,1,nFleet,nSage,nNage);
	for( k = 1; k <= nFleet; k++ )
	{
		kk = nFleetIndex(k);
		d_ak(k)  = dAllocation(k);
		for( ig = 1; ig <= n_ags; ig++ )
		{
			d_V(ig)(k) = (exp(d4_log_sel(kk)(ig)(nNyr)));
		}
	}

	// average weight-at-age in referecne points.
	dmatrix wt_bar(1,n_ags,nSage,nNage);
	dmatrix fa_bar(1,n_ags,nSage,nNage);
	dmatrix M_bar(1,n_ags,nSage,nNage);
	// cout<<"Average weight"<<endl;
	for( ig = 1; ig <= n_ags; ig++ )
	{
		wt_bar(ig) = colsum(d3_wt_avg(ig))/(nNyr-nSyr+1);
		fa_bar(ig) = colsum(d3_wt_mat(ig))/(nNyr-nSyr+1);
		M_bar(ig)  = colsum(d3_Mt(ig))/(nNyr-nSyr+1);
	}
	
	
	m_dFmsy.allocate(1,nFleet);
	m_dFmsy=0.1 / nFleet;

	/* Instantiate the MSY class for each stock */
	for( g = 1; g <= nStock; g++ )
	{
		// Working here.
		Msy cMSY(dRo(g),dSteepness(g),M_bar,d_rho,wt_bar,fa_bar,&d_V);
		cMSY.get_fmsy(m_dFmsy);

	}
}


/**
* @brief Condition the reference model for numbers-at-age between syr and nyr
* @author Steve Martell
* 
*/
void OperatingModel::conditionReferenceModel()
{

	int f,g,i,j;
	int ig,ih;
	d3_Nt.initialize();
	for( ig = 1; ig <= n_ags; ig++ )
	{
		f  = n_area(ig);
		g  = n_group(ig);
		ih = pntr_ag(f,g);

		dvector lx(nSage,nNage);
		dvector tr(nSage,nNage);
		lx(nSage) = 1.0;
		for( j = nSage; j <nNage; j++ )
		{
			lx(j+1) = lx(j) * exp(-d3_Mt(ig)(nSyr)(j));
		}
		lx(nNage) /= 1.0-exp(-d3_Mt(ig)(nSyr,nNage));
		
		// initialize unfished conditions.
		if( d_iscamCntrl(5) )
		{
			tr = log(dRo(g)) + log(lx);
		}
		else if (! d_iscamCntrl(5) )
		{
			tr(nSage)         = log(dAvgRec(ih)) + dLog_rbar_devs(ih)(nSyr);
			tr(nSage+1,nNage) = log(dInitRec(ih)) + dLog_init_rec_devs(ih);
			tr                += log(lx); 
		}

		// fill numbers-at-age arrays
		d3_Nt(ig)(nSyr) = 1./nSex * mfexp(tr);

		for( i = nSyr; i <= nNyr; i++ )
		{
			d3_Zt(ig)(i) = -1.0*log(d3_St(ig)(i));
			if( i > nSyr )
			{
				d3_Nt(ig)(i,nSage) = 1./nSex*dAvgRec(ih)*mfexp( dLog_rbar_devs(ih)(i) );
			}
			d3_Nt(ig)(i+1)(nSage+1,nNage) =++ elem_prod(d3_Nt(ig)(i)(nSage,nNage-1)
			                                            ,d3_St(ig)(i)(nSage,nNage-1));
			d3_Nt(ig)(i+1,nNage) += d3_Nt(ig)(i,nNage)*d3_St(ig)(i,nNage);

		}
		d3_Nt(ig)(nNyr+1,nSage) = 1./nSex * dAvgRec(ih);
		
		// cout<<d3_Nt(ig)<<endl; // ce code cest bon.
	}
}

/** \brief Determine stock-recruitment parameters
	
	Determine Beverton-Holt SR parameters given Ro and steepness.

		\f$ R = \frac{s_o S}{(1+ b S)} \f$

	\author  Steve Martell
	
	\sa
**/
void OperatingModel::calcStockRecruitment()
{
	int f,g,h,j;
	int ig;
	double  dt = d_iscamCntrl(13);      // get this from cntrl(13) in the control file.
	double  phib;
	dvector ma(nSage,nNage);
	dvector fa(nSage,nNage);
	dvector lx(nSage,nNage);
	dvector lw(nSage,nNage);

	m_so.allocate(1,nStock);
	m_beta.allocate(1,nStock);
	m_dSbo.allocate(1,nStock);


	for( g = 1; g <= nStock; g++ )
	{
		// Survivorship
		lx.initialize();
		lw.initialize();
		lx(nSage) = 1.0;
		lw(nSage) = 1.0;
		phib = 0;
		for( f = 1; f <= nArea; f++ )
		{
			for( h = 1; h <= nSex; h++ )
			{
				ig = pntr_ags(f,g,h);	
				for(j=nSage; j<=nNage; j++)
				{
					ma(j) = mean(trans(d3_Mt(ig))(j));
					fa(j) = mean(trans(d3_wt_mat(ig))(j));
			
					if(j>nSage)
					{
						lx(j) = lx(j-1) * mfexp(-ma(j-1));
					}
					lw(j) = lx(j) * mfexp(-ma(j)*dt);
				}
				lx(nNage) /= 1.0 - mfexp(-ma(nNage));
				lw(nNage) /= 1.0 - mfexp(-ma(nNage));

				// Average spawning biomass per recruit.
				phib += 1./(nArea * nSex) * lw*fa;
			}
		}
		// Stock-recruitment parameters
		m_so(g) = m_kappa(g)/phib;
		m_dSbo(g) = dRo(g) * phib;
		m_beta(g) = (m_kappa(g)-1.0) / m_dSbo(g);
	}
 cout<<"OpMod\n"<<m_beta<<endl; 
}

/* calculate Selectivity coefficients */
// Aug 12. Possibly Deprecate this function.
void OperatingModel::calcSelectivities()
{
	/*
		Calculate selectivity arrays for each gear,stock,sex,year,age
		based on selectivty type (logistic, coefficients etc.) and 
		selectivity parameters defined in the Scenario class.
	*/

	// int f;
	int g;
	int h;
	// int i;
	// int j;
	int k;
	d5_logSel.allocate(1,nGear,1,nStock,1,nSex,nSyr,nPyr,nSage,nNage);

	cout<<"In calcSelectivities"<<endl;
	cout<<nSel_type<<endl;
	Selex cAgeSelex(dAge);
	Selex cLenSelex;

	// logistic_selectivity clogisticSelex(dAge);
	for( k = 1; k <= nGear; k++ )
	{
		for( g = 1; g <= nStock; g++ )
		{
			for( h = 1; h <= nSex; h++ )
			{
				switch(nSel_type(k))
				{
					default:
						d5_logSel(k,g,h) = 0.0;
					break;

					// case 1: logistic curve
					case 1:
						cAgeSelex.logistic(d3_selPars(k),nSel_block(k),d5_logSel(k,g,h));
					break;

					// case 2: selectivity coefficients
					case 2:
						cAgeSelex.selcoeff(d3_selPars(k),nSel_block(k),d5_logSel(k,g,h));
					break;

					// case 11: length-based logistic curve
					case 11:
					// Bug here in that there is no length-info in the projection years.
						cout<<"sex "<<h<<endl;
						cout<<d3_wt_avg(h)<<endl;
						cout<<"SuCK EGGS"<<endl;
						dmatrix len = pow( d3_wt_avg(h)/d_a(h), 1./d_b(h) );
						cLenSelex.logistic(d3_selPars(k),nSel_block(k),len,d5_logSel(k,g,h));
					break;
				}
				cout<<d5_logSel(k,g,h)<<endl;
			}
		}
	}
	
}


