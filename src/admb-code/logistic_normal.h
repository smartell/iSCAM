/**
	This is a class for implementing the logistic normal negative
	loglikelihoood for composition data.
*/

#include <admodel.h>

#ifndef LOGISTIC_NORMAL_H
#define LOGISTIC_NORMAL_H

class logistic_normal
{
private:
	int         m_b1;  				///> index for first bin (col)
	int         m_b2;				///> inded for last bin
	
	int         m_y1; 				///> index for first year (row)
	int         m_y2;				///> index for last year 

	double m_eps;                   ///> small constant to add to proportions.
	double m_dMinimumProportion;	///> smallest proportion for aggregation.


	dvariable   m_nll;				///> negative loglikelihood
	dvariable   m_sig2;				///> variance weighting parameter.
	dvariable   m_sig;              ///> stdev weighting parameter.

	ivector     m_nB2;				///> integer vector for # of rows in ragged O & E.
	dvector     m_dWy;				///> relative weight for each year.

	imatrix     m_nAgeIndex;		///> Index for aggregated residuals
	dmatrix     m_O;				///> Raw Data
	dmatrix     m_Op;				///> Observed proportions

	dvar_matrix m_E;				///> Raw expected data
	dvar_matrix m_Ep;				///> Expected proportions

	dmatrix m_residual;			    ///> residuals.
	dvar_matrix m_w;        		///> Residual relative to O_{By}
	dvar3_array m_V;				///> Covariance matrix for likelihood
	
	adtimer m_gprof;				///> Timer class for profiling.

	logistic_normal();

	dvariable calc_sigma_mle();

	void add_constant(const double& eps);
	void aggregate_and_compress_arrays();
	
	void compute_relative_weights();
	void compute_residual_arrays();
	void compute_covariance_matrix();
	void compute_covariance_matrix(const dvariable& phi);
	void compute_covariance_matrix(const dvariable& phi,const dvariable& psi);
	
	void compute_mle_sigma(const dvar3_array& dCor);
	void compute_negative_loglikelihood();
	void compute_standardized_residuals();
	

public:
	~logistic_normal();

	// default constructor
	logistic_normal(const dmatrix *_O,const dvar_matrix *_E,
                    const double _minProportion,const double eps=0);

	/* getters */
	dvariable& operator ()();
	dvariable& operator ()(const dvariable& phi);
	dvariable& operator ()(const dvariable& phi, const dvariable& psi);

	dvariable get_nll()       { return(m_nll);         }
	dmatrix   get_residuals() { compute_standardized_residuals(); return(m_residual); }
	double get_sig2()         { return(value(m_sig2)); }
	double get_sig()          { return(value(m_sig));  }
};

/**
 * Template fucntion to return the product of an array.
**/
template <typename T1, typename T2>
T2 prod(T1 x, T2 count)
{	
	RETURN_ARRAYS_INCREMENT();
	int lb = x.indexmin();
	T2 p = x(lb);
	for(int i = lb+1; i <= count; i++ )
	{
		p *= x(i);
	}
	RETURN_ARRAYS_DECREMENT();
	return(p);
}

#endif