// |---------------------------------------------------------------------------------|
// | OPERATING MODEL FOR USE WITH ISCAM
// |---------------------------------------------------------------------------------|
// |

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
OperatingModel::~OperatingModel()
{

};

/* Constructor */
OperatingModel::OperatingModel(const s_iSCAMdata&  mse_data, const s_iSCAMvariables& mse_vars)
{
	initializeConstants(mse_data);
	initializeVariables(mse_vars);
}

/* Initialize private member variables based on scenario class */
void OperatingModel::initializeConstants(const s_iSCAMdata& cS)
{
	// Model dimensions
	nStock = cS.nStock;
	nArea  = cS.nArea;
	nSex   = cS.nSex;
	nSyr   = cS.nSyr;
	nNyr   = cS.nNyr;
	nPyr   = cS.nPyr;
	nSage  = cS.nSage;
	nNage  = cS.nNage;
	nGear  = cS.nGear;

	// links for array indexing
	n_ags = nArea * nStock * nSex;
	n_ag  = nArea * nStock;
	n_gs  = nStock * nSex;
	n_area.allocate(1,n_ags);
	n_group.allocate(1,n_ags);
	n_sex.allocate(1,n_ags);
	pntr_ag.allocate(1,nArea,1,nStock);
	pntr_gs.allocate(1,nStock,1,nSex);
	pntr_ags.allocate(1,nArea,1,nStock,1,nSex);
	int f,g,h;
	int ig = 0;
	int ih = 0;
	int is = 0;
	for( f = 1; f <= nArea; f++ )
	{
		for( g = 1; g <= nStock; g++ )
		{
			ih ++;
			pntr_ag(f,g) = ih;
			for( h = 1; h <= nSex; h++ )
			{
				ig ++;
				n_area(ig)  = f;
				n_group(ig) = g;
				n_sex(ig)   = h;
				pntr_ags(f,g,h) = ig;
				if ( f==1 )
				{
					is ++;
					pntr_gs(g,h) = is;
				}
			}
		}
	}


	dAge.allocate(nSage,nNage);
	dAge.fill_seqadd(nSage,1);

	// Growth & Maturity
	d_linf  = cS.d_linf;
	d_vonbk = cS.d_vonbk;
	d_to   	= cS.d_to;
	d_a     = cS.d_a;
	d_b 	= cS.d_b;
	d_ah 	= cS.d_ah;
	d_gh 	= cS.d_gh;

	// Catch Data
	nCtNobs 	= cS.nCtNobs;
	dCatchData 	= cS.dCatchData;

	// Survey data
	nIt         = cS.nIt;
	nItNobs     = cS.nItNobs;
	nSurveyType = cS.nSurveyType;
	dSurveyData.allocate(*cS.dSurveyData);
	dSurveyData = *cS.dSurveyData;

	// Composition data
	nAgears = cS.nAgears;
	nAnobs  = cS.nAnobs;
	nAsage  = cS.nAsage;
	nAnage  = cS.nAnage;
	dA.allocate(*cS.dA);
	dA      = *cS.dA;

	// Empirical weight-at-age
	nWtNobs = cS.nWtNobs;
	dWt_avg.allocate(*cS.dWt_avg);
	dWt_mat.allocate(*cS.dWt_mat);
	dWt_avg = *cS.dWt_avg;
	dWt_mat = *cS.dWt_mat;

	cout<<"initializeConstants"<<endl;
	cout<<pntr_ags<<endl;

}


void OperatingModel::initializeVariables(const s_iSCAMvariables& cS)
{
	// Leading parameters
	cout<<"initializeVariables"<<endl;
	dRo                = exp(cS.d_log_ro);
	dSteepness         = cS.d_steepness;
	dM                 = exp(cS.d_log_m);
	dAvgRec            = exp(cS.d_log_rbar);
	dInitRec           = exp(cS.d_log_rinit);
	dRho               = cS.d_rho;
	dVartheta          = sqrt(1.0/cS.d_varphi);
	dLog_m_devs        = cS.dLog_M_devs;
	dLog_rbar_devs     = cS.dLog_Rbar_devs;
	dLog_init_rec_devs = cS.dLog_Rinit_devs;
	
	nSel_type  = cS.nSel_type;
	nSel_block = cS.nSel_block;
	d3_selPars.allocate(*cS.dSelPars);
	d3_selPars = *cS.dSelPars;
	
	dFt                = cS.dFt;
	
	d3_Mt.allocate(*cS.d3_Mt);
	d3_Mt = *cS.d3_Mt;

	// cout<<dFt<<endl;
}

void OperatingModel::runScenario()
{
	/*
	PSUEDOCODE:
		- declare local variables.
		- initialize stock-recruitment parameters
		- condition operating model on historical assessment
		|-- calculate reference points
		| - calculate harvest control rule
		| - implement harvest on reference population
		| - update reference population
		| - generate data for annual assessment model
		| - run stock-assessment procedures
		|-- repeat:
		- write simulation variables to output file.
		- write performance statistics to output file.
	*/

	/*
	- Local variables
	*/
	int f,g,h,i,j,k;

	/*
	- Initialize stock-recruitment parameters for each stock.
	- Survivorship and fecundity is based on average mortality & weight-at-age
	- (so, beta, sbo, ro)
	*/

	// calcStockRecruitment
	{
		double  dt = 0.5;      // get this from cntrl(13) in the control file.
		double  phib;
		dvector ma(nSage,nNage);
		dvector fa(nSage,nNage);
		dvector lx(nSage,nNage);
		dvector lw(nSage,nNage);

		m_so.allocate(1,nStock);
		m_beta.allocate(1,nStock);
		m_dSbo.allocate(1,nStock);
		for( g = 1; g <= nStock; g++ )
		{
			for( h = 1; h <= nSex; h++ )
			{
				
				// Survivorship
				lx.initialize();
				lw.initialize();
				lx(nSage) = 1.0;
				lw(nSage) = 1.0;
				for(j=nSage; j<=nNage; j++)
				{
					ma(j) = mean(trans(d3_Mt(g))(j));
					fa(j) = mean(trans(dWt_mat(g))(j));
					if(j>nSage)
					{
						lx(j) = lx(j-1) * mfexp(-ma(j-1));
					}
					lw(j) = lx(j) * mfexp(-ma(j)*dt);
				}
				lx(nNage) /= 1.0 - mfexp(-ma(nNage));
				lw(nNage) = lx(nNage);

				// Average spawning biomass per recruit.
				phib += lw*fa;
			}
		}
	cout<<ma<<endl;
	}
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
	cout<<nSel_type<<endl;
	Selex cAgeSelex(dAge);
	Selex cLenSelex;

	// logistic_selectivity clogisticSelex(dAge);
	for( k = 1; k <= nGear; k++ )
	{
		for( g = 1; g <= nStock; g++ )
		{
			for( h = 1; h <= nSex; h++ )
			{
				switch(nSel_type(k))
				{
					default:
						d5_logSel(k,g,h) = 0.0;
					break;

					// case 1: logistic curve
					case 1:
						cAgeSelex.logistic(d3_selPars(k),nSel_block(k),d5_logSel(k,g,h));
					break;

					// case 2: selectivity coefficients
					case 2:
						cAgeSelex.selcoeff(d3_selPars(k),nSel_block(k),d5_logSel(k,g,h));
					break;

					// case 11: length-based logistic curve
					case 11:
					// Bug here in that there is now length-info in the projection years.
						cout<<"sex "<<h<<endl;
						cout<<dWt_avg(h)<<endl;
						cout<<"SuCK EGGS"<<endl;
						dmatrix len = pow( dWt_avg(h)/d_a(h), 1./d_b(h) );
						cLenSelex.logistic(d3_selPars(k),nSel_block(k),len,d5_logSel(k,g,h));
					break;
				}
				cout<<d5_logSel(k,g,h)<<endl;
			}
		}
	}
	
}


