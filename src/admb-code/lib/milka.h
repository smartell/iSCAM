// milka.h

#ifndef MILKA_H
#define MILKA_H

#undef REPORT
#define REPORT(object) report << #object "\n" << object << endl;

#include <admodel.h>
#include <contrib.h>
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
		dvector q;
		dmatrix sbt;

		dmatrix log_rec_devs;
		dmatrix init_log_rec_devs;

		// Selectivity parameters
		d3_array *d3_log_sel_par;
		d4_array *d4_logSel;

		// Mortality
		d3_array *d3_M;
		d3_array *d3_F;
		d3_array *d3_ft;

		//recruitment
		dvector sbo;
		dvector so;
	};

	class OperatingModel: public model_data
	{
	private:

		int m_nNyr; 
		int m_irow;
		int m_nyrs;  // number of simulation years

		ivector m_nGearIndex;
		ivector m_nCSex;
		ivector m_nASex;
		dvector m_nATau;
		ivector m_nWSex;
		dvector m_dLslim;
		dvector m_dUslim;
		dvector m_dDiscMortRate;
		imatrix m_nAGopen;

		//random variables
		dmatrix m_epsilon;  /// observation errors
		dmatrix m_delta;    /// recruitment deviations (process errors)
		double  m_gamma_r;  /// recruitment autocorrelation.
		
		// catch arrays
		int      m_nCtNobs;
		dmatrix  m_dCatchData;
		dmatrix  m_dSubLegalData;
		d3_array m_d3_Ct;
		
		// survey arrays
		ivector  m_n_it_nobs;
		d3_array m_d3SurveyData;

		// composition arrays
		ivector m_A_irow;
		ivector  m_n_A_nobs;
		d3_array m_d3_A;

		//weight at age arrays
		ivector m_W_irow;
		ivector m_nWtNobs;
		d3_array m_d3_inp_wt_avg;

		// MSE controls
		int m_nPyr;				/// Terminal year for Operating Model.
		int m_nSeed;			/// random number seed
		int m_nRecType;
		double m_dBthreshold;
		double m_dBlimit;
		double m_maxf;

		dmatrix m_dispersal; 

		adstring MseCtlFile;
		adstring MsePfcFile;
		
		int m_nn;
		
		dvector m_dRo;
		dvector m_dBo;
		dvector m_dSteepness;
		dvector m_dM;
		dvector m_dRbar;
		dvector m_dRinit;
		dvector m_dRho;
		dvector m_dVarphi;
		dvector m_dSigma;
		dvector m_dTau;
		dvector m_dKappa;
		dvector m_dbeta;
		dvector m_q;



		// Assessment model results
		dvector m_est_bo;
		dmatrix m_est_fmsy;
		dmatrix m_est_msy;
		dvector m_est_bmsy;
		dvector m_est_sbtt;
		dvector m_est_btt;
		dmatrix m_est_N;
		dmatrix m_est_wa;
		dmatrix m_est_M;
		dmatrix m_est_log_sel;

		dmatrix m_sbt;

		dmatrix  m_dTAC;
		int     m_nHCR;

		dmatrix  m_log_rt;

		d3_array m_N;
		d3_array m_M;
		d3_array m_F;
		d3_array m_Z;
		d3_array m_S;
		d3_array m_d3_wt_avg;
		d3_array m_d3_wt_mat;
		d3_array m_ft;
		//d3_array m_log_sel_par;

		d4_array d4_logSel;

		//PerformanceVariables
		dmatrix m_AAV;

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
		void calcTotalMortality(const int& iyr);
		void calcRelativeAbundance(const int& iyr);
		void calcCompositionData(const int& iyr);
		void calcEmpiricalWeightAtAge(const int& iyr);
		void updateReferenceModel(const int& iyr);
		void writeDataFile(const int& iyr);
		void runStockAssessment();
		void writeSimulationVariables();

		
	};


dvector cubic_spline(const dvector& spline_coffs, const dvector& la);

template<typename T, typename T1>
T1 retention_probability(const T &lsl, const T&usl, const T1 &mu, const T1 &sd)
{
	// Use the cumd_norm function to determine proportion of mu > lsl
	int x1 = mu.indexmin();
	int x2 = mu.indexmax();
	T1 p(x1,x2);
	for(int i = x1; i <= x2; i++ )
	{
		const T zl = (mu(i) - lsl)/sd(i);
		const T zu = (mu(i) - usl)/sd(i);
		p(i)=cumd_norm(zl) - cumd_norm(zu);
	}
	return (p);
}

// } // mse namespace


#endif