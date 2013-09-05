// |---------------------------------------------------------------------------------|
// | GENERAL OPERATING MODEL FOR USE WITH iSCAM.                                     |
// |---------------------------------------------------------------------------------|



#include <admodel.h>
#include "OpMod.h"
#include "Selex.h"
#include "msy.h"
#include "baranov.h"



/** \brief Destructor
	\author Steven Martell 
**/
OperatingModel::~OperatingModel()
{

};

/** \brief Constructor 
	\author Steven Martell
**/
OperatingModel::OperatingModel(const s_iSCAMdata&  mse_data, const s_iSCAMvariables& mse_vars)
{
	initializeConstants(mse_data);
	initializeVariables(mse_vars);
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
	nStock = cS.nStock;  
	nArea  = cS.nArea;
	nSex   = cS.nSex;
	nSyr   = cS.nSyr;
	nNyr   = cS.nNyr;
	nPyr   = cS.nPyr;
	nSage  = cS.nSage;
	nNage  = cS.nNage;
	nGear  = cS.nGear;
	nFleet = cS.nFleet;


	// | links for array indexing
	n_ags = nArea * nStock * nSex;
	n_ag  = nArea * nStock;
	n_gs  = nStock * nSex;
	n_area.allocate(1,n_ags);
	n_group.allocate(1,n_ags);
	n_sex.allocate(1,n_ags);
	pntr_ag.allocate(1,nArea,1,nStock);
	pntr_gs.allocate(1,nStock,1,nSex);
	pntr_ags.allocate(1,nArea,1,nStock,1,nSex);

	int f,g,h,i,j,k;
	int ig = 0;
	int ih = 0;
	int is = 0;
	for( f = 1; f <= nArea; f++ )
	{
		for( g = 1; g <= nStock; g++ )
		{
			ih ++;
			pntr_ag(f,g) = ih;
			for( h = 1; h <= nSex; h++ )
			{
				ig ++;
				n_area(ig)  = f;
				n_group(ig) = g;
				n_sex(ig)   = h;
				pntr_ags(f,g,h) = ig;
				if ( f==1 )
				{
					is ++;
					pntr_gs(g,h) = is;
				}
			}
		}
	}
	// cout<<"Houston, we have a problem"<<endl;
	nFleetIndex.allocate(1,nFleet);
	nFleetIndex = cS.nFleetIndex;

	dAllocation  = cS.dAllocation;
	dAge.allocate(nSage,nNage);
	dAge.fill_seqadd(nSage,1);

	// Growth & Maturity
	d_linf  = cS.d_linf;
	d_vonbk = cS.d_vonbk;
	d_to   	= cS.d_to;
	d_a     = cS.d_a;
	d_b 	= cS.d_b;
	d_ah 	= cS.d_ah;
	d_gh 	= cS.d_gh;

	// Catch Data
	nCtNobs 	= cS.nCtNobs;
	dCatchData 	= cS.dCatchData;
	d3_Ct.allocate(1,n_ags,nSyr,nPyr,1,nGear);
	d3_Ct.initialize();
	for(int ii=1;ii<=nCtNobs;ii++)
	{
		i = dCatchData(ii)(1);
		k = dCatchData(ii)(2);
		f = dCatchData(ii)(3);
		g = dCatchData(ii)(4);
		h = dCatchData(ii)(5);
		int ncut;
		h == 0? ncut=nSex: ncut=1;
	
		for(int hh=1;hh<=nSex;hh++)
		{
			ig = pntr_ags(f,g,hh);
			d3_Ct(ig)(i)(k) = 1./ncut*dCatchData(ii)(7);
		}
			
	}
	// Protected member catch array required for writing data file
	int nCount = nCtNobs + (nPyr-nNyr+1)*nSex*nGear;
	m_nCtNobs  = nCtNobs;  /**< Initialize counter */
	m_dCatchData.allocate(1,nCount,1,7);
	m_dCatchData.initialize();
	m_dCatchData.sub(1,nCtNobs) = dCatchData;
	



	// Survey data
	nIt         = cS.nIt;
	nItNobs     = cS.nItNobs;
	nSurveyType = cS.nSurveyType;
	dSurveyData.allocate(*cS.dSurveyData);
	dSurveyData = *cS.dSurveyData;

	// Composition data
	nAgears = cS.nAgears;
	nAnobs  = cS.nAnobs;
	nAsage  = cS.nAsage;
	nAnage  = cS.nAnage;
	dA.allocate(*cS.dA);
	dA      = *cS.dA;

	// Empirical weight-at-age
	nWtNobs = cS.nWtNobs;
	dWt_avg.allocate(*cS.dWt_avg);
	dWt_mat.allocate(*cS.dWt_mat);
	dWt_avg = *cS.dWt_avg;
	dWt_mat = *cS.dWt_mat;
	dWt_bar.allocate(1,n_ags,nSage,nNage);
	dEt_bar.allocate(1,n_ags,nSage,nNage);
	d3_wt_avg.allocate(1,n_ags,nSyr,nPyr,nSage,nNage);
	d3_wt_mat.allocate(1,n_ags,nSyr,nPyr,nSage,nNage);
	dWt_bar.initialize();
	dEt_bar.initialize();
	d3_wt_avg.initialize();
	d3_wt_mat.initialize();
	for( ig = 1; ig <= n_ags; ig++ )
	{
		d3_wt_avg(ig).sub(nSyr,nNyr+1) = dWt_avg(ig);
		d3_wt_mat(ig).sub(nSyr,nNyr+1) = dWt_mat(ig);
		
		dWt_bar(ig) = d3_wt_avg(ig)(nNyr+1);
		dEt_bar(ig) = d3_wt_mat(ig)(nNyr+1);

		//  For now assume average weight in future 
		//  Could assume density-dependent growth etc. & overwrite these variables
		for( i = nNyr+2; i <= nPyr; i++ )
		{
			d3_wt_avg(ig)(i) = dWt_bar(ig);
			d3_wt_mat(ig)(i) = dEt_bar(ig);			
		}
		
		// cout<<"d3_wt_avg\n"<<d3_wt_avg(ig)<<endl;
		// cout<<"dWt_bar\n"<<dWt_bar(ig)<<endl;
		// cout<<"dEt_bar\n"<<dEt_bar(ig)<<endl;
	}
	// cntrl vector
	nCntrl = cS.cntrl;
	cout<<"initializeConstants"<<endl;
	cout<<"pntr_ags\n"<<pntr_ags<<endl;

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

	d3_Nt.allocate(1,n_ags,nSyr,nPyr,nSage,nNage);
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
}

void OperatingModel::runScenario(const int &seed)
{
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
		cout<<"Ok apres implementFisheries       \t pas fini"<<endl;
		
		/*
		- Update reference population
		*/	
		updateReferencePopulation(i);
		cout<<"Ok apres updateReferencePopulation \t pas fini"<<endl;

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
  	dfs<< nArea 		<<endl;
  	dfs<< nStock		<<endl;
  	dfs<< nSex			<<endl;
  	dfs<< nSyr   		<<endl;
  	dfs<< iyr   		<<endl;
  	dfs<< nSage  		<<endl;
  	dfs<< nNage  		<<endl;
  	dfs<< nGear 		<<endl;
 
  	dfs<<"#Allocation"	<<endl;
  	dfs<< dAllocation 	<<endl;
  	

  	dfs<<"#Age-schedule and population parameters"<<endl;
  	dfs<< d_linf  		<<endl;
	dfs<< d_vonbk  		<<endl;
	dfs<< d_to  		<<endl;
	dfs<< d_a  			<<endl;
	dfs<< d_b  			<<endl;
	dfs<< d_ah  		<<endl;
	dfs<< d_gh  		<<endl;

  	dfs<<"#Observed catch data"<<endl;
  	dfs<< nCtNobs 		<<endl;
  	dfs<< dCatchData    <<endl;  

 //  	dfs<<"#Abundance indices"	<<endl;
 //  	dfs<< nit 					<<endl;
 //  	dfs<< nit_nobs 				<<endl;
 //  	dfs<< survey_type 			<<endl;
 //  	dfs<< survey_data 			<<endl;

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
	int f,g,h,i,j,k;
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
	m_dTac.allocate(1,nFleet);
	m_dTac.initialize();
	int g;
	for( g = 1; g <= nStock; g++ )
	{
		m_dTac += 0.2/nFleet * m_est_bt(g);
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

	int f,g,h,i,j,k;
	int ig,kk;

	/* Get selectivities and allocations for each fleet. */
	dvector d_ak(1,nFleet);  // Allocation for each fleet.
	d3_array d_V(1,n_ags,1,nFleet,nSage,nNage);
	for( k = 1; k <= nFleet; k++ )
	{
		kk = nFleetIndex(k);
		d_ak(k) = dAllocation(k);
		for( ig = 1; ig <= n_ags; ig++ )
		{
			d_V(ig)(k) = exp(d4_log_sel(kk)(ig)(nNyr));
		}

	}

	/* Calculate fleet specific fishing mortality rates by group */
	// TO DO
	// [ ] - add joint capture probability for size-based selectivity and 32" size limit.
	// [ ] - add time-varying natural mortality rates.
	// [ ] - add implementation error.
	BaranovCatchEquation cBaranov;
	for( ig = 1; ig <= n_ags; ig++ )
	{

		dvector na = d3_Nt(ig)(iyr);
		dvector wa = d3_wt_avg(ig)(iyr);
		dvector ma = d3_Mt(ig)(nNyr);
		dmatrix va = d_V(ig);

		// Add implementation error here.
		// m_dTac comes from the harvest control rule.
		// Also need to record the actual catch in m_dCatchData for data file.
		dvector ct = m_dTac;

		// SM Working here end of Sept 4. 2013


		// cout<<"Na\t"<<na<<endl;
		// cout<<"wa\t"<<wa<<endl;
		// cout<<"ma\t"<<ma<<endl;
		// cout<<"va\t"<<va<<endl;
		// cout<<"ct\t"<<ct<<endl;

		dvector ft = cBaranov.getFishingMortality(ct,ma,va,na,wa);
		m_dFt = ft;

		// Total mortality
		d3_Zt(ig)(iyr) = ma;
		for( k = 1; k <= nFleet; k++ )
		{
			d3_Zt(ig)(iyr) += m_dFt(k) * va(k);
		}
		// cout<<"Zt\t"<<d3_Zt(1)(iyr)<<endl;
	}

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

	int f,g,h,i,j,k,kk;
	int ig;
	double d_rho = nCntrl(13);
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
		wt_bar(ig) = colsum(dWt_avg(ig))/(nNyr-nSyr+1);
		fa_bar(ig) = colsum(dWt_mat(ig))/(nNyr-nSyr+1);
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

	int f,g,h,i,j,k;
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
		if( nCntrl(5) )
		{
			tr = log(dRo(g)) + log(lx);
		}
		else if (! nCntrl(5) )
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
	int f,g,h,i,j,k;
	int ig;
	double  dt = nCntrl(13);      // get this from cntrl(13) in the control file.
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
					fa(j) = mean(trans(dWt_mat(ig))(j));
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

	int f;
	int g;
	int h;
	int i;
	int j;
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
					// Bug here in that there is now length-info in the projection years.
						cout<<"sex "<<h<<endl;
						cout<<dWt_avg(h)<<endl;
						cout<<"SuCK EGGS"<<endl;
						dmatrix len = pow( dWt_avg(h)/d_a(h), 1./d_b(h) );
						cLenSelex.logistic(d3_selPars(k),nSel_block(k),len,d5_logSel(k,g,h));
					break;
				}
				cout<<d5_logSel(k,g,h)<<endl;
			}
		}
	}
	
}


