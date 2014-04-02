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


	/**
	 * @brief Class for storing data for operating model
	 * @details [long description]
	 * @return NULL
	 */
	class omData{
	private:
		// |------------------|
		// | Model dimensions |
		// |------------------|
		int m_nStock;
		int m_nArea;
		int m_nSex;
		int m_nSyr;
		int m_nNyr;
		int m_nSage;
		int m_nNage;
		int m_nGear;
		int m_nFleet;
		dvector m_dAllocation;

		ModelData s_md;
		
		d3_array m_d3_ct;

	    //friend class OperatingModel;
	public:
		omData();
		~omData();
		omData(const d3_array &ct)
		:m_d3_ct(ct)
		{
			cout<<"Hi Catarina"<<endl;
			cout<<m_d3_ct<<endl;
		}


		omData(const ModelData &_md)
		:s_md(_md)
		{
			cout<<"THis is fucking cool"<<endl;
			cout<<s_md.nNage<<endl;
			cout<<s_md.d_linf<<endl;
			cout<<s_md.nCtNobs<<endl;
			cout<<s_md.nItNobs<<endl;
			cout<<"Average weight\n"<<endl;
			cout<<*s_md.d3_wt_avg<<endl;			
		}

		// |---------|
		// | Setters |
		// |---------|
		void set_nStock(const int &n) { m_nStock = n; cout<<n<<endl;}
		void set_nArea (const int &n) { m_nArea  = n; }
		void set_nSex  (const int &n) { m_nSex   = n; }
		void set_nSyr  (const int &n) { m_nSyr   = n; }
		void set_nNyr  (const int &n) { m_nNyr   = n; }
		void set_nSage (const int &n) { m_nSage  = n; }
		void set_nNage (const int &n) { m_nNage  = n; }
		void set_nGear (const int &n) { m_nGear  = n; }
		void set_nFleet(const int &n) { m_nFleet = n; }

		void set_dAllocation (const dvector &d) {m_dAllocation = d;}
		

	};

	/**
	 * @brief Class for storing variables for operating model
	 * @details [long description]
	 * 
	 */
	class omVariables{
	private:
		dvector m_log_Ro;
		friend class OperatingModel;
	public:
		omVariables();
		~omVariables();

		void set_log_Ro(const dvector &n) { m_log_Ro = n; }
	};

	class OperatingModel
	{
	private:
		int          m_nSeed;
		omData       m_data;
		omVariables  m_vars;
		
	protected:
		
	public:
		~OperatingModel();

		OperatingModel(const omData& md, const omVariables& mv, const int &seed);

		void runScenario(const int &seed);

	};
} // mse namespace


#endif