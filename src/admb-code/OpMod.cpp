#include <admodel.h>
#include "OpMod.h"

// |---------------------------------------------------------------------------------|
// | SCENARIO MEMBER FUNCTIONS, CONSTRUCTOR & DESTRUCTOR
// |---------------------------------------------------------------------------------|
// |
Scenario::~Scenario(){};

Scenario::Scenario(const int nStock, 
         	const int nArea, 
         	const int nSex, 
         	const int nSyr, 
           	const int nNyr, 
           	const int nPyr, 
           	const int nSage, 
           	const int nNage,
           	const dvector log_ro, 
           	const dvector steepness,
           	const dvector log_m,
           	const dvector log_avgrec,
           	const dvector log_initrec,
           	const dvector rho,
           	const dvector vartheta,
           	d3_array selpar)
  			: 
  			ModelData(nStock,nArea,nSex,nSyr,nNyr,nPyr,nSage,nNage),
   			ModelParams(log_ro,steepness,log_m, log_avgrec, log_initrec, rho, vartheta,
   			            selpar)
  {
    cout<<"In Constructor"<<endl;
    cout<<m_nStock<<endl;
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

	// Leading parameters
	dRo        = mfexp(cS.m_log_ro);
	dSteepness = cS.m_steepness;
	dM         = mfexp(cS.m_log_m);
	dAvgRec    = mfexp(cS.m_log_avgrec);
	dInitRec   = mfexp(cS.m_log_initrec);
	dRho       = mfexp(cS.m_rho);
	dVartheta  = mfexp(cS.m_vartheta);

	dKappa     = elem_div(4. * dSteepness, 1.-dSteepness);
}