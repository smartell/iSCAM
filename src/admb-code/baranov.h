#include <admodel.h>

#ifndef BARANOV_CATCH_EQUATION_H
#define BARANOV_CATCH_EQUATION_H

#define MAXITS 50
#define MAXF   5.0
#ifndef TOL 
#define TOL    1.e-9
#endif


/** \brief  Baranov Catch Equation
	
	A class for iteratively solving the Baranov Catch Equation
	when conditioned on catch.
	
© Copyright `2013`  - . All Rights Reserved.

	\author  Steve Martell
	\version 1.0
	\sa
**/
class BaranovCatchEquation
{
	dmatrix m_hCt;
public:
	//! Constructor
	BaranovCatchEquation();

	//! Destructor
	/**
	 * No pointer objects are used in this class.
	**/
	~BaranovCatchEquation();

	double get_ft(const double& ct, const double& m, const dvector& va, const dvector& ba);
	dvector getFishingMortality(const dvector &ct, const double &m, const dmatrix &V, const dvector &na);
	dvector getFishingMortality(const dvector &ct, const double &m, const dmatrix &V, const dvector &na, const dvector &wa);

	dvector getFishingMortality(const dvector &ct, const dvector &ma, const dmatrix &V, const dvector &na);
	dvector getFishingMortality(const dvector &ct, const dvector &ma, const dmatrix &V, const dvector &na, const dvector &wa);
	
	// SM Added 2sex version Sept. 10, 2013.  See Baranov2Sex.R for example.
	dvector getFishingMortality(const dvector &ct, const dmatrix &ma, const d3_array *_V, const dmatrix &na);
	dvector getFishingMortality(const dvector &ct, const dmatrix &ma, const d3_array *_V, const dmatrix &na, const dmatrix &wa);

	// SM 2sex version which modifes a catch by sex matrix
	dvector getFishingMortality(const dvector &ct, const dmatrix &ma, const d3_array *p_V, const dmatrix &na, dmatrix &_hCt);
	dvector getFishingMortality(const dvector &ct, const dmatrix &ma, const d3_array *p_V, const dmatrix &na, const dmatrix &wa, dmatrix &_hCt);
	
};


// dvector get_ft(dvector& ct,const double& m, const dmatrix& V,const dvector& na, const dvector& wa);
// dvector get_ft(dvector& ct,const double& m, const dmatrix& V,const dvector& ba);




#endif