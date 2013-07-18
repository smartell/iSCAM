

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

#ifndef LOGISTIC_SELECTIVITY_H
#define LOGISTIC_SELECTIVITY_H
class logistic_selectivity
{
private: 
	dvector x;
	dvector dy;
	dvar_vector y;

public:
	logistic_selectivity(){};
	logistic_selectivity(const dvector & _x);
	~logistic_selectivity();

	dvector     operator () (const dvector & log_selpar);
	dvar_vector operator () (const dvar_vector & log_selpar);


};
#endif

#ifndef SELEX_H
#define SELEX_H

class Selex
{
private:
	int         m_selType;
	dvector     m_x;
	

	dmatrix     m_dlogsel;
	dvar_matrix m_dvar_logsel;
	dvar_matrix m_dvar_selPar;

	// Selectivity classes
	logistic_selectivity m_plogis;

public:
	Selex();
	Selex(const int& _selType, const dvector& _x, const dvar_matrix& _dvar_selPar);
	~Selex();

	void fill_selex_array(dvar_matrix& log_sel);

	/* logistic prototypes */
	dvector     logistic( const dvector& x,const double& mu, const double& sd );
	dvar_vector logistic( const dvector& x,const dvariable& mu, const dvariable& sd );
	dvar_vector logistic( const dvar_vector& log_selpar );

	/* eplogistic */
	dvector     eplogis( const dvector& x, const double& x1, 
                        const double& x2, const double& gamma );
	dvar_vector eplogis( const dvector& x, const dvariable& x1, 
                            const dvariable& x2, const dvariable& gamma );


	/* Linear approximation */
	dvector     linapprox(const dvector& x, const dvector& y, const dvector& xout);
	dvar_vector linapprox(const dvector& x, const dvar_vector& y, const dvector& xout);
};

#endif

