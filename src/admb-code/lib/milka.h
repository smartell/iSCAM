// milka.h

#ifndef MILKA_H
#define MILKA_H

#include <admodel.h>
#include "../iscam.htp"

/**
 * @defgroup Milka Operating model for iSCAM
 * @Milka Classes and functions for the operating model
 * 
 * @author  Steven Martell
 * @deprecated  Feb 23, 2014
 * 
 * @details The namespace is mse, short for management strategy evaluation.
 * 
 * Steve  April 1, 2014.
 * There should be two structs for passing datastructures into this class:
 * ModelData:       -> struct that reflects the DATA_SECTION
 * ModelVariables:  -> struct that reflects the PARAMETER_SECTION
 * 
 * CLASS Objects:
 * Operating model:
 * 
 */

// namespace mse {

	struct ModelVariables
	{
		dvector log_ro;
		dvector steepness;
		dvector m;
		dvector log_rbar;
		dvector log_rinit;
		dvector rho;
		dvector varphi;


		// Selectivity parameters
		d3_array *d3_log_sel_par;
	};

	class OperatingModel: public model_data
	{
	private:

		int m_nNyr;   
		
		int m_nCtNobs;
		dmatrix m_dCatchData;
	
		ivector m_n_it_nobs;
		d3_array m_d3SurveyData;

		ivector m_n_A_nobs;
		d3_array m_d3_A;

		ivector m_nWtNobs;
		d3_array m_d3_inp_wt_avg;

		int m_nPyr;				/// Terminal year for Operating Model.

		dvector m_dRo;
		dvector m_dSteepness;
		dvector m_dM;
		dvector m_dRbar;
		dvector m_dRinit;
		dvector m_dRho;
		dvector m_dVarphi;
		dvector m_dSigma;
		dvector m_dTau;
		dvector m_dKappa;

		d4_array d4_logSel;
		d3_array m_log_sel_par;

		ModelVariables mv;		// Structure for model variables.

	public:
		OperatingModel(ModelVariables _mv,int argc,char * argv[]);		
		~OperatingModel();
	
		void runScenario();

	protected:
		void readMSEcontrols();
		void initParameters();
		void conditionReferenceModel();
		void setRandomVariables();
		void getReferencePointsAndStockStatus();
		void calculateTAC();
		void allocateTAC();
		void implementFisheries();

		void updateReferenceModel();
		void writeDataFile();
		void runStockAssessment();
		

		
	};


dvector cubic_spline(const dvector& spline_coffs, const dvector& la);

// } // mse namespace


#endif