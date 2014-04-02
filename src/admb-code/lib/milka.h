// milka.h

#ifndef MILKA_H
#define MILKA_H

#include <admodel.h>

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

namespace mse {

	/**
	 * @brief Data structure for storing model data
	 * @details Data variables treated as Constants in the operating model.
	 * 
	 */
	struct ModelData
	{
		// Model dimensions
		int nStock;
		int nArea;
		int nSex;
		int nSyr;
		int nNyr;
		int nSage;
		int nNage;
		int nGear;
		int nFleet;
		dvector dAllocation;

		// Growth And Maturity Parameters
		dvector d_linf;
		dvector d_vonbk;
		dvector d_to;
		dvector d_a;
		dvector d_b;
		dvector d_ah;
		dvector d_gh;
		dvector d_maturityVector;
		int n_MAT;

		// Catch array
		int nCtNobs;
		d3_array *d3_Ct;

		// Abundance Indices
		int nItNobs;
		ivector n_it_nobs;
		ivector n_survey_type;
		d3_array *d3_survey_data;

		// Composition Data
		int nAgears;
		ivector  n_A_nobs;
		ivector  n_A_sage;
		ivector  n_A_nage;
		d3_array *d3_A;

		// Weight-At-Age Data
		const d3_array *d3_wt_avg;
		const d3_array *d3_wt_mat;
		const d3_array *d3_len_age;
	};

	struct ModelVariables
	{
		dvector log_ro;
		dvector steepness;
		dvector m;
		dvector log_rbar;
		dvector log_rinit;
		dvector rho;
		dvector varphi;
	};

	class OperatingModel
	{
	private:
		ModelData      md;		// Structure for model data.
		ModelVariables mv;		// Structure for model variables.

	public:
		OperatingModel(const ModelData &_md, const ModelVariables &_mv);
		~OperatingModel();
		
	};


} // mse namespace


#endif