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
	int         m_B;    			///> total number of bins (cols)

	int         m_y1; 				///> index for first year (row)
	int         m_y2;				///> index for last year 
	int         m_Y;    			///> total number of years (rows)

	double m_eps;                   ///> small constant to add to proportions.
	double m_dMinimumProportion;	///> smallest proportion for aggregation.


	dvariable   m_nll;				///> negative loglikelihood
	dvariable   m_sig2;				///> variance weighting parameter.
	dvariable   m_sig;              ///> stdev weighting parameter.
	dvariable   m_rho;				///> correlation coefficient.

	ivector     m_nB2;				///> integer vector for # of rows in ragged O & E.
	ivector     m_nNminp;			///> index for aggregated proportions
	dvector     m_dWy;				///> relative weight for each year.

	imatrix     m_nAgeIndex;		///> Index for aggregated residuals
	dmatrix     m_O;				///> Raw Data
	dmatrix     m_Op;				///> Observed proportions
	dmatrix     m_Ox;				///> Logistic transform using Schnute's approach.
	dmatrix     m_Oy;				///> Logistic transform using Aitichson's approach.
	dmatrix     m_Oz;         		///> Lostitic transform for residual calculations.
	dmatrix     m_Oa;				///> Aggregated matrix for tail compression and zeros

	dvar_matrix m_E;				///> Raw expected data
	dvar_matrix m_Ep;				///> Expected proportions
	dvar_matrix m_Ex;				///> Logistic transform using Schnute's approach.
	dvar_matrix m_Ey;				///> Logistic transform using Aitichson's approach.
	dvar_matrix m_Ez;         		///> Lostitic transform for residual calculations.
	dvar_matrix m_Ea;				///> Aggregated matrix for tail compression and zeros

	dmatrix m_residual;			    ///> residuals.
	dvar_matrix m_w;        		///> Residual relative to O_{By}
	dvar_matrix m_s;				///> Residual based on multivariate normal.
	dvar_matrix m_covar;
	dvar_matrix m_Vmat;
	dvar_matrix m_Vinv;
	adtimer m_gprof;
	dvar3_array m_V;				///> Covariance matrix for likelihood
	dvar3_array m_C;				///> Correlation matrix
	

	logistic_normal();

	dvariable calc_sigma_mle();
	void add_constant(const double& eps);
	void aggregate_and_compress_arrays();
	
	void compute_relative_weights();
	void compute_residual_arrays();
	void compute_correlation_matrix();
	
	void compute_mle_sigma();
	void compute_negative_loglikelihood();
	void compute_standardized_residuals();
	// void aggregate_arrays();
	// void correlation_matrix();
	// dvar_matrix calc_covmat(const dvariable& sig2,const int& n);
	// dvar_matrix calc_KCK(const dvar_matrix& V);

public:
	~logistic_normal();
	// logistic_normal(const dmatrix& _O, const dvar_matrix& _E);
	logistic_normal(const dmatrix _O,const dvar_matrix _E,
                    const double _minProportion,const double eps=0);

	/* data */
	// dvariable nll(const dvariable& sig2);
	// dvariable negative_loglikelihood(const dvariable& tau2);
	// dvariable negative_loglikelihood();
	// dvar_matrix standardized_residuals();


	/* setters */
	// void set_MinimumProportion(double &p) {m_dMinimumProportion = p;}

	/* getters */
	dvariable get_nll()       { compute_negative_loglikelihood(); return(m_nll);      }
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