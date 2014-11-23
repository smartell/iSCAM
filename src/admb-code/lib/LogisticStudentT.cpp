#include <admodel.h>
#include "../include/LogisticNormal.h"

// Constructor
logistic_student_t::logistic_student_t(const dmatrix& _O,const dvar_matrix& _E,
	                				   const double _minp,const double _eps)
:logistic_normal(_O,_E,_minp,_eps)
{
	m_v = 300000.0;
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

	// mle of the student T variance
	// v/(v-2)∑, where ∑ is the covariance matrix
	// mle of the variance
	double n = size_count(m_Op) * (m_y2-m_y1+1);
	m_sigma2 = m_wss / n; 
	m_sigma  = sqrt(m_sigma2);

	// Compute negative loglikelihood
	m_nll = negative_log_likelihood();

	return m_nll;
}

dvariable logistic_student_t::operator () (const dvariable& _df)
{
	
	m_nll = 0;

	// Get correlation vector rho
	get_rho();

	// Construct convariance (m_V)
	compute_correlation_array();

	// Compute weighted sumofsquares
	compute_weighted_sumofsquares();

	// mle of the student T variance
	// v/(v-2)∑, where ∑ is the covariance matrix
	m_v = _df;
	//m_sigma2 = m_v/(m_v-2.0);

	// mle of the variance
	double n = size_count(m_Op) * (m_y2-m_y1+1);
	m_sigma2 = m_wss / n; 
	m_sigma  = sqrt(m_sigma2);

	// Compute negative loglikelihood
	m_nll = negative_log_likelihood();

	return m_nll;
}

dvariable logistic_student_t::negative_log_likelihood()
{
	// 7) Compute nll using student t-distrbution.
	// v is the degrees of freedom.
	
	cout<<"Student T, df = "<<m_v<<endl;
	RETURN_ARRAYS_INCREMENT();
	dvariable v = m_v;
	int p;
  	dvariable nll = 0.0;
  	for(int i = m_y1; i <= m_y2; i++ )
  	{
  		p    = size_count(m_Op(i)) - 1;
  		const double lppi2 = 0.5*p*log(PI);
  		nll += -1.0 * gammln(0.5*(v + p));
  		nll += gammln(0.5*v) + 0.5*p*log(v) + lppi2;
  		nll += 0.5*log(det(m_V(i)));
  	}
	nll += 0.5*(p+v) * log(1.0 + m_wss/v);

	RETURN_ARRAYS_DECREMENT();
	return nll;
}
