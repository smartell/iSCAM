

/**
 * @file Selex.h
 * @author Steve Martell
 * @Description
 * This class calculated selectivity coefficients based on parametric and
 * nonparametric functions.  Idea here is to create a Selex class object then
 * call the appropriate function based on the selectivity-type.
 * 
 * Returns a log_sel matrix for a given fishery that may be time-invariate, blocked, 
 * or continuous changes over year.
**/
#ifndef SELEX_H
#define SELEX_H
enum eSelType
{
	LOGISTIC = 1,
	EPLOGISTIC,
	LINEARAPPROX
};

class Selex
{
private:
	int     m_selType;
	dvector m_x;
	dmatrix m_dlogsel;
	dvar_matrix m_dvar_logsel;

public:
	Selex();
	Selex(const dvector& x);
	~Selex();

	/* logistic prototypes */
	dvector     logistic( const dvector& x,const double& mu, const double& sd );
	dvar_vector logistic( const dvector& x,const dvariable& mu, const dvariable& sd );
	dvar_vector logistic( const dvar_vector& selpar );

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