#include <admodel.h>
#include "OpMod.h"
#include "Selex.h"

test::test(const s_iSCAMdata& data)
{
	cout<<"Im here in the test constructor"<<endl;
	cout<<data.nStock<<endl;
	d3_array junk = *data.dSurveyData;
	cout<<junk<<endl;
	cout<<"Leaving test constructor"<<endl;
};

// |---------------------------------------------------------------------------------|
// | ScenarioParameters
// |---------------------------------------------------------------------------------|
// |
ScenarioParameters::~ScenarioParameters(){};
// ScenarioParameters::ScenarioParameters(){};
ScenarioParameters::ScenarioParameters(
        	const dvector&  _dBo,
			const dvector&  _dSteepness,
			const dvector&  _dM,
			const dvector&  _dRbar,
			const dvector&  _dRinit,
			const dvector&  _dRho,
			const dvector&  _dVarphi,
			const dvector&  _dLog_M_devs,
			const dmatrix&  _dLog_Rbar_devs,
			const dmatrix&  _dLog_Rinit_devs,
			const d3_array& _dSelPars,
			const dmatrix&  _dFt
                                       )
:
m_dBo(_dBo),
m_dSteepness(_dSteepness),
m_dM(_dM),
m_dRbar(_dRbar),
m_dRinit(_dRinit),
m_dRho(_dRho),
m_dVarphi(_dVarphi),
m_dLog_M_devs(_dLog_M_devs),
m_dLog_Rbar_devs(_dLog_Rbar_devs),
m_dLog_Rinit_devs(_dLog_Rinit_devs),
m_dSelPars(_dSelPars),
m_dFt(_dFt)
{
	cout<<"O=k to here dude"<<endl;
};

// |---------------------------------------------------------------------------------|
// | SCENARIO MEMBER FUNCTIONS, CONSTRUCTOR & DESTRUCTOR
// |---------------------------------------------------------------------------------|
// |
Scenario::~Scenario(){};

// Scenario::Scenario(const ScenarioData& data, const ScenarioParameters& params)
// {
// 	cout<<"Scenario is all set"<<endl;
// 	cout<<m_nPyr<<endl;
// }

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


