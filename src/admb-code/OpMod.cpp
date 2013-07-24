#include <admodel.h>
#include "OpMod.h"
#include "Selex.h"
// |---------------------------------------------------------------------------------|
// | ScenarioData
// |---------------------------------------------------------------------------------|
// |
ScenarioData::~ScenarioData(){};
ScenarioData::ScenarioData(
			const int&      _nStock,
			const int&      _nArea,
			const int&      _nSex,
			const int&      _nSyr,
			const int&      _nNyr,
			const int&      _nPyr,
			const int&      _nSage,
			const int&      _nNage,
			const int&      _nGear,
			const dvector&  _dAllocation,
			const dvector&  _d_linf,
			const dvector&  _d_vonbk,
			const dvector&  _d_to,
			const dvector&  _d_a,
			const dvector&  _d_b,
			const dvector&  _d_ah,
			const dvector&  _d_gh,
			const int&      _nCtNobs,
			const dmatrix&  _dCatchData,
			const int&      _nIt,
			const ivector&  _nItNobs,
			const ivector&  _nSurveyType,
			const d3_array& _dSurveyData,
			const int&      _nWtNobs,
			const d3_array& _dWt_avg,
			const d3_array& _dWt_mat
                           )
:
m_nStock(_nStock),
m_nArea(_nArea),
m_nSex(_nSex),
m_nSyr(_nSyr),
m_nNyr(_nNyr),
m_nPyr(_nPyr),
m_nSage(_nSage),
m_nNage(_nNage),
m_nGear(_nGear),
m_dAllocation(_dAllocation),
m_d_linf(_d_linf),
m_d_vonbk(_d_vonbk),
m_d_to(_d_to),
m_d_a(_d_a),
m_d_b(_d_b),
m_d_ah(_d_ah),
m_d_gh(_d_gh),
m_nCtNobs(_nCtNobs),
m_dCatchData(_dCatchData),
m_nIt(_nIt),
m_nItNobs(_nItNobs),
m_nSurveyType(_nSurveyType),
m_dSurveyData(_dSurveyData),
m_nWtNobs(_nWtNobs),
m_dWt_avg(_dWt_avg),
m_dWt_mat(_dWt_mat)
{

};


// |---------------------------------------------------------------------------------|
// | SCENARIO MEMBER FUNCTIONS, CONSTRUCTOR & DESTRUCTOR
// |---------------------------------------------------------------------------------|
// |
Scenario::~Scenario(){};

Scenario::Scenario(
            const int nStock, 
         	const int nArea, 
         	const int nSex, 
         	const int nSyr, 
           	const int nNyr, 
           	const int nPyr, 
           	const int nSage, 
           	const int nNage,
           	const int nGear,
           	const dvector log_ro, 
           	const dvector steepness,
           	const dvector log_m,
           	const dvector log_avgrec,
           	const dvector log_initrec,
           	const dvector rho,
           	const dvector vartheta,
           	const d3_array& selpar,
           	const ivector & sel_type)
  			: 
  			ModelData(nStock,nArea,nSex,nSyr,nNyr,nPyr,nSage,nNage,nGear),
   			ModelParams(log_ro,steepness,log_m, log_avgrec, log_initrec, rho, vartheta,
   			            selpar), m_sel_type(sel_type)

  {
    cout<<"In Constructor"<<endl;
    cout<<m_nStock<<endl;
    cout<<m_selpars<<endl;
    cout<<"Leaving constructor"<<endl;
  }

ModelData::~ModelData(){};

ModelParams::~ModelParams(){};

// |---------------------------------------------------------------------------------|
// | OPERTATING MODEL MEMBER FUNCTIONS, CONSTRUCTORS, AND DESTRUCTOR
// |---------------------------------------------------------------------------------|
// |

/* Destructor */
OperatingModel::~OperatingModel(){};

/* Constructor */
OperatingModel::OperatingModel(Scenario &cScenario)
: m_cScenario(cScenario)
{
	initializeVariables(cScenario);
	calcSelectivities();
}

/* Initialize private member variables based on scenario class */
void OperatingModel::initializeVariables(Scenario& cS)
{
	// Model dimensions
	nStock = cS.m_nStock;
	nArea  = cS.m_nArea;
	nSex   = cS.m_nSex;
	nSyr   = cS.m_nSyr;
	nNyr   = cS.m_nNyr;
	nPyr   = cS.m_nPyr;
	nSage  = cS.m_nSage;
	nNage  = cS.m_nNage;
	nGear  = cS.m_nGear;

	dAge.allocate(nSage,nNage);
	// dAge.fill_seqadd(nSage,1);


	// Leading parameters
	dRo        = mfexp(cS.m_log_ro);
	dSteepness = cS.m_steepness;
	dM         = mfexp(cS.m_log_m);
	dAvgRec    = mfexp(cS.m_log_avgrec);
	dInitRec   = mfexp(cS.m_log_initrec);
	dRho       = mfexp(cS.m_rho);
	dVartheta  = mfexp(cS.m_vartheta);

	dKappa     = elem_div(4. * dSteepness, 1.-dSteepness);

	// Selectivity parameters
	nSel_type  = cS.m_sel_type;
	d3_selPars.allocate(cS.m_selpars);
	d3_selPars = cS.m_selpars;

}

/* calculate Selectivity coefficients */
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
	Selex cSelex();
	logistic_selectivity clogisticSelex(dAge);
	for( k = 1; k <= nGear; k++ )
	{
		for( g = 1; g <= nStock; g++ )
		{
			for( h = 1; h <= nSex; h++ )
			{
				for( i = nSyr; i <= nPyr; i++ )
				{
					d5_logSel(k,g,h,i) = clogisticSelex(d3_selPars(k,g));
					
				}
			}
		}
	}
	
}


