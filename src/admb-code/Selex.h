#include <admodel.h>
#include "iscam.htp"


/**
 * @file Selex.h
 * @author Steve Martell
 * @Description
 * Selex is the main aggregation class.  User instantiates Selex, this will then initialize
 * class members for various selectivity functions.
 * 
 * This class calculated selectivity coefficients based on parametric and
 * nonparametric functions.  Idea here is to create a Selex class object then
 * call the appropriate function based on the selectivity-type.
 * 
 * Returns a log_sel matrix for a given fishery that may be time-invariate, blocked, 
 * or continuous changes over year.
**/



/** \brief  Class object for selectivity functions
	
	A series of alternative selectivity functions based on:
		- Logistic
		- coefficients
		- linear approximation
	
© Copyright `date +%Y`  - . All Rights Reserved.

	\author  
	\author $LastChangedBy$
	\date `date +%Y-%m-%d`
	\date $LastChangedDate$
	\version $Rev$	\sa
**/
#ifndef SELEX_H
#define SELEX_H

class Selex
{
private:
	dvector     m_x;
	

	
public:
	Selex();
	Selex(const dvector& _x);
	~Selex();

	void set_x(const dvector& _x){m_x = _x;} 	//!< set independent variables
	dvector get_x(){return(m_x);}				//!< return independent variables

	/* logistic prototypes */
	dvector     logistic( const dvector& x,const double& mu, const double& sd );
	dvar_vector logistic( const dvector& x,const dvariable& mu, const dvariable& sd );
	dvar_vector logistic( const dvar_vector& log_selpar );
	void        logistic( const dmatrix& theta, const ivector& blk, dmatrix& log_sel);
	void        logistic( const dmatrix& theta, const ivector& blk, const dmatrix& len, dmatrix& log_sel);

	/* selectivity coefficients */
	void selcoeff( const dmatrix& theta, const ivector& blk, dmatrix& log_sel);

	/* eplogistic */
	// dvector     eplogis( const dvector& x, const double& x1, 
                        // const double& x2, const double& gamma );
	// dvar_vector eplogis( const dvector& x, const dvariable& x1, 
                            // const dvariable& x2, const dvariable& gamma );


	/* Linear approximation */
	dvector     linapprox(const dvector& x, const dvector& y, const dvector& xout);
	dvar_vector linapprox(const dvector& x, const dvar_vector& y, const dvector& xout);
};

#endif

