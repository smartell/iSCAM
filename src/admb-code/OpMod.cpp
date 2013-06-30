#include <admodel.h>
#include "OpMod.h"
#include "Selex.h"

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
	int f;
	int g;
	int h;
	int i;
	int j;
	int k;
	d5_logSel.allocate(1,nGear,1,nStock,1,nSex,nSyr,nPyr,nSage,nNage);

	cout<<"In calcSelectivities"<<endl;
	
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


