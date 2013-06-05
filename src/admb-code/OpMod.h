/**
 * \file OpMod.h
 * \author Steve Martell
 * \date Jun 4, 2013
**/

#include <admodel.h>

#ifndef _MODEL_DIMENSION_H
#define _MODEL_DIMENSION_H
/**
 * \brief A class that defines the dimensions of the Operating Model.
 * \author Steve Martell
**/
class ModelDimension
{
private:

public:
	int m_nArea;
	int m_nGroup;
	int m_nSex;
	int m_nSyr;
	int m_nNyr;
	int m_nSage;
	int m_nNage;
	int m_nGear;

	int m_nAGS;
	int m_nAG;
	int m_nGS;

	d3_array m_pAGS;

	~ModelDimension();
	ModelDimension();
	ModelDimension(const int nArea,
	               const int nGroup,
	               const int nSex,
	               const int nSyr,
	               const int nNyr,
	               const int nSage,
	               const int nNage,
	               const int nGear);
};

#endif

#ifndef _STOCK_PARAMETERS_H
#define _STOCK_PARAMETERS_H
/**
 * \brief A class object for the stock parameters
 * \author Steven Martell
 * \
**/
class StockParameters: public ModelDimension
{
private:

public:
	dvector m_dRo;
	dvector m_dSteepness;
	dvector m_dNaturalMortality;
	dvector m_dRbar;
	dvector m_dRinit;
	double m_dRho;
	double m_dVartheta;

	~StockParameters();
	StockParameters(const int f,
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
	                const double vartheta);

	/* data */
};
#endif
