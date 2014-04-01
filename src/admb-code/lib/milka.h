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
 * <br> The following class objects are: <br><br>
 * <br> Class
 */

namespace mse {

	struct ModelDimensions
	{
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
	};


	struct GrowthAndMaturityParameters
	{
		dvector d_linf;
		dvector d_vonbk;
		dvector d_to;
		dvector d_a;
		dvector d_b;
		dvector d_ah;
		dvector d_gh;
		dvector d_maturityVector;
		int n_MAT;
	};

	struct TimeSeriesData
	{
		int nCtNobs;
		d3_array *d3_Ct;
	};

	struct AbundanceIndices
	{

		int nItNobs;
		ivector n_it_nobs;
		ivector n_survey_type;
		d3_array *d3_survey_data;
	};

	struct CompositionData
	{
		int nAgears;
		ivector  n_A_nobs;
		ivector  n_A_sage;
		ivector  n_A_nage;
		d3_array *d3_A;

	};
	struct WeightAtAgeData
	{
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

		ModelDimensions s_md;
		GrowthAndMaturityParameters s_gamp;
		TimeSeriesData s_tsd;
		AbundanceIndices s_ai;
		WeightAtAgeData s_waad;
		CompositionData s_cd;


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


		omData(const ModelDimensions &_md, 
			   const GrowthAndMaturityParameters &_gamp,
			   const TimeSeriesData &_tsd,
			   const AbundanceIndices &_ai,
		       const CompositionData &_cd,
		       const WeightAtAgeData &_waad)

		:s_md(_md), s_gamp(_gamp), s_tsd(_tsd), s_ai(_ai),s_waad(_waad),s_cd(_cd)
		{
			cout<<"THis is fucking cool"<<endl;
			cout<<s_md.nNage<<endl;
			cout<<s_gamp.d_linf<<endl;
			cout<<s_tsd.nCtNobs<<endl;
			cout<<s_ai.nItNobs<<endl;
			cout<<"Average weight\n"<<endl;
			cout<<*s_waad.d3_wt_avg<<endl;			
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