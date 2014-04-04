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

		dmatrix log_rec_devs;
		dmatrix init_log_rec_devs;

		// Selectivity parameters
		d3_array *d3_log_sel_par;
		d4_array *d4_logSel;

		// Mortality
		d3_array *d3_M;
		d3_array *d3_F;

	};

	class OperatingModel: public model_data
	{
	private:

		int m_nNyr; 
		ivector m_nGearIndex;
		ivector m_nCSex;
		ivector m_nASex;
		imatrix m_nAGopen;
		
		int m_nCtNobs;
		dmatrix m_dCatchData;
	
		ivector m_n_it_nobs;
		d3_array m_d3SurveyData;

		ivector m_n_A_nobs;
		d3_array m_d3_A;

		ivector m_nWtNobs;
		d3_array m_d3_inp_wt_avg;

		int m_nPyr;				/// Terminal year for Operating Model.
		int m_nSeed;			/// random number seed

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

		// Assessment model results
		dvector m_est_bo;
		dmatrix m_est_fmsy;
		dmatrix m_est_msy;
		dvector m_est_bmsy;
		dvector m_est_sbtt;
		dvector m_est_btt;

		dmatrix  m_dTAC;
		int     m_nHCR;

		dmatrix  m_log_rt;

		d3_array m_N;
		d3_array m_M;
		d3_array m_F;
		d3_array m_Z;
		d3_array m_S;
		d3_array m_d3_wt_avg;
		d3_array m_ft;
		d3_array m_log_sel_par;

		d4_array d4_logSel;

		ModelVariables mv;		// Structure for model variables.

	public:
		OperatingModel(ModelVariables _mv,int argc,char * argv[]);		
		~OperatingModel();
	
		void runScenario(const int &seed);

	protected:
		void readMSEcontrols();
		void initParameters();
		void initMemberVariables();
		void conditionReferenceModel();
		void setRandomVariables(const int &seed);
		void getReferencePointsAndStockStatus();
		void calculateTAC();
		void allocateTAC(const int& iyr);
		void implementFisheries(const int& iyr);

		void updateReferenceModel();
		void writeDataFile();
		void runStockAssessment();
		

		
	};


dvector cubic_spline(const dvector& spline_coffs, const dvector& la);

// } // mse namespace


#endif