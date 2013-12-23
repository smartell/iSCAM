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

	double m_dMinimumProportion;	///> smallest proportion for aggregation.


	dvariable   m_nll;				///> negative loglikelihood
	dvariable   m_sig2;				///> variance weighting parameter
	dvariable   m_rho;				///> correlation coefficient.

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

	dvar_matrix m_residual;			///> residuals.
	dvar_matrix m_w;        		///> Residual relative to O_{By}
	dvar_matrix m_covar;
	dvar_matrix m_Vmat;
	dvar_matrix m_Vinv;

	dvar3_array m_V;				///> Covariance matrix
	dvar3_array m_C;				///> Correlation matrix

	dvariable calc_sigma_mle();
	void aggregate_arrays();
	void correlation_matrix();
	dvar_matrix calc_covmat(const dvariable& sig2,const int& n);
	dvar_matrix calc_KCK(const dvar_matrix& V);

public:
	~logistic_normal();
	logistic_normal();
	logistic_normal(const dmatrix& _O, const dvar_matrix& _E);

	/* data */
	dvariable nll(const dvariable& sig2);
	dvariable negative_loglikelihood(const dvariable& tau2);
	dvariable negative_loglikelihood();
	dvar_matrix standardized_residuals();


	/* setters */
	void set_MinimumProportion(double &p) {m_dMinimumProportion = p;}
};


#endif