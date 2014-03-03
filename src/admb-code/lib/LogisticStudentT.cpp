#include <admodel.h>
#include "LogisticNormal.h"

// Constructor
logistic_student_t::logistic_student_t(const dmatrix& _O,const dvar_matrix& _E,
	                				   const double _minp,const double _eps)
:logistic_normal(_O,_E,_minp,_eps)
{

}

dvariable logistic_student_t::operator () ()
{
	
	m_nll = 0;

	// Get correlation vector rho
	get_rho();

	// Construct convariance (m_V)
	compute_correlation_array();

	// Compute weighted sumofsquares
	compute_weighted_sumofsquares();

	// Compute negative loglikelihood
	m_nll = negative_log_likelihood();

	return m_nll;
}


dvariable logistic_student_t::negative_log_likelihood()
{
	// 7) Compute nll using student t-distrbution.
	// v is the degrees of freedom.
	
	cout<<"Student T"<<endl;
	RETURN_ARRAYS_INCREMENT();
	double v = 5.0e6;
	double p;
  	dvariable nll = 0.0;
  	for(int i = m_y1; i <= m_y2; i++ )
  	{
  		p    = size_count(m_Op(i));
  		nll += -1.0 * gammln(0.5*(v + p));
  		nll += gammln(0.5*v) + 0.5*p*log(v) + 0.5*p*log(PI);
  		nll += 0.5*log(det(m_V(i)));
  	}
	nll += 0.5*(p+v) * log(1.0 + m_wss/v);

	RETURN_ARRAYS_DECREMENT();
	return nll;
}
