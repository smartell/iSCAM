#include <admodel.h>
#include "OpMod.h"

// |---------------------------------------------------------------------------------|
// | ModelDimensions: CONSTRUCTORS & GETTERS & SETTERS 
// |---------------------------------------------------------------------------------|
// |
ModelDimension::ModelDimension()
{}

ModelDimension::ModelDimension(
                   const int nArea,
	               const int nGroup,
	               const int nSex,
	               const int nSyr,
	               const int nNyr,
	               const int nSage,
	               const int nNage,
	               const int nGear)
:	m_nArea(nArea),
	m_nGroup(nGroup),
	m_nSex(nSex),
	m_nSyr(nSyr),
	m_nNyr(nNyr),
	m_nSage(nSage),
	m_nNage(nNage),
	m_nGear(nGear)
{
	m_nAGS = nArea  * nGroup * nSex;
	m_nAG  = nArea  * nGroup;
	m_nGS  = nGroup * nSex;

	m_pAGS.allocate(1,m_nArea,1,m_nGroup,1,m_nSex);
	int f,g,h,ig;
	ig = 0;
	for( f = 1; f <= m_nArea; f++ )
	{
		for( g = 1; g <= m_nGroup; g++ )
		{
			for( h = 1; h <= m_nSex; h++ )
			{
				ig ++;
				m_pAGS(f,g,h) = ig;				
			}
		}
	}
	cout<<m_pAGS(1)<<endl;
}


ModelDimension::~ModelDimension()
{

}

// |---------------------------------------------------------------------------------|
// | StockParameters: 
// |---------------------------------------------------------------------------------|
// |
StockParameters::StockParameters(const int f,
				                 const int g,
				                 const int h,
				                 const int syr,
				                 const int nyr,
				                 const int sage,
				                 const int nage,
				                 const int ngear,
				                 const dvector ro,
				                 const dvector steepness,
				                 const dvector m,
				                 const dvector rbar,
				                 const dvector rinit,
				                 const double rho,
				                 const double vartheta)
: 	ModelDimension(f,g,h,syr,nyr,sage,nage,ngear),
	m_dRo(ro),
	m_dSteepness(steepness),
	m_dNaturalMortality(m),
	m_dRbar(rbar),
	m_dRinit(rinit),
	m_dRho(rho),
	m_dVartheta(vartheta)
{

};

StockParameters::~StockParameters()
{};

