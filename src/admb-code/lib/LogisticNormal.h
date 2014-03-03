/**
 * \brief Class for implementing the logistic normal negative log-likelihood.
 * \author Steven Martell
 * \addtogroup Likelihoods [group-title]
 */


/// Jimmy can ride a bike

/**
Logistic-normal likelihood for size-age composition data.
Much of this code is based on a paper written by Chris Francis

Author: Steve Martell
Date  : Nov 22, 2013.

Descrition:
	Let  Y = number of years of composition information (rows)
	let  B = number of age/size bins in the composition matrix (cols)

	---- CONSTRUCTOR ----
	The constructor requires and observed (O) and expected (E) matrix
	of composition data.  These data are immediately transformed into
	proportions that sum to 1 in each year.

	There are 3 different arrays that are constructed from (O) and (E).
	O_x is the logistic transform O_x = exp(X_b)/sum_b exp(X_b)
	O_y is based on the Aitchison (2003) approach where vector of 
	observations was transformed as:
		O_y = exp(Y_b)/(1+sum_b^{B-1} exp(Y_b)) , for b = 1, ..., B-1
		O_y = 1 - sum_b^{B-1} O_b               , for b = B
	O_z is a multivariate normal vector constructed from (O):
		O_z = exp(Z_b)/sum_b exp(Z_b) - mean[exp(Z_b)/sum_b exp(Z_b)]



	Then I use the aggregate arrays function to almalgamate 0s and perform
	tail compression.  An member array of residuals m_w is created in this 
	function where: 

		m_w = ln(O_b/O_B) - ln(E_b/E_B)

	Note that each m_w vector (year) is of length (B-1).

	---- LIKELIHOODS ----
	There are three separate likelihoods where each differes by the number of
	weighting parameters used to define the covariance matrix.
		1) LN1: based on sig2
		2) LN2: based on sig2, phi_1
		3) LN3: based on sig2, phi_1, and phi_2.

	In call three cases differences in among year sample sizes are weighted by:
		W_y = [mean(N)/N_y]^0.5
	where N_y is sample size each year.
	
	-- LN1 --
	For the simplist likelihood with no correlation and an estimated variance 
	tau2, the likelihood for each year is given by:
	
	nll = 0.5(B-1)*ln(2*\pi) + sum(ln(O_b)) + 0.5*ln(|V|) + 
	      (B-1)*ln(W_y) + 0.5* m_w^2/(V*W_y^2).

	where V is the covariance matrix which is just a simple scaler sig2.
		V = sig2*diag(B)

	---- METHODS FOR SUPRESSING ZEROS ----

	---- CONDITIONAL MLE FOR THE VARIANCE SCALER ----
	The maximum likelihood estimate for sigma can be calculated 
	conditional on the observed and expected values, as well as the
	autocorrelation parameters (correlation matrix).

	The variance is given by the correlated weighted sum of squares.

	sig2 = {sum_y (w_y V_y^{-1} w_y)/W_y^2} / sum_y(B_y-1)

	where V_y^{-1} = K C K', where C is the correlation matrix.
	

	—————————————————————————————————————————————————————————————————————————————————————
	PSUEDOCODE FOR IMPLEMENTATION OF THE NEGATIVE LOG-LIKELIHOOD
	—————————————————————————————————————————————————————————————————————————————————————
	1). In the constructor, suppress 0s using:
		a) add a constant and renormalize the rows of _O and _E.
		b) aggregate & compress tails based on minimum proportions.
	
	3). Compute annual relative weights for variance scaling.
	4). Compute residual arrays w_y & s_y
		w_y = log(O_yb/O_yB) - log(E_yb/E_yB);             //for sum of squares
		s_y = log(O_yb-mean(O_yb)) - log(E-yb-mean(E_yb)); //for residuals
	5). Compute a vector of covariance arrays for each year:
		V_y = K C K';             // for negloglikelihood
		G_y = F'H^{-1} VH^{-1} F; // for std residuals
	6). Compute conditional mle estimate of sigma
		sigma = [{sum_y (w_y' V^{-1} w_y)/ W_y^2} / sum_y (B_y-1) ]^0.5
	7). Compute negative loglikelihood:
		nll  = t1 + t2 + t3 + t4 + t5 + t6; where:
		bym1 = sum_y(B_y-1);
		t1   = 0.5*log(2pi) * bym1;
		t2   = sum_{by} log(O_{by});
		t3   = log(sigma) * bym1;
		t4   = 0.5*sum_y log(det(V_y))
		t5   = sum_y (By-1)*log(W_y)
		t6   = 0.5*sigma^{-2}* sum_y (w'_y V_y^{-1} w_y)/W_y^2
	8). Compute standardized residuals.
		res{_by} = [log(O_{by}/\tilde(O)) - log(E_{by}/\tilde(E))] / (W_y G_y^{0.5})
	—————————————————————————————————————————————————————————————————————————————————————


**/

#include <admodel.h>
#ifndef __LOGISTIC_NORMAL_H
#define __LOGISTIC_NORMAL_H

class logistic_normal
{
private:
	int m_y1;
	int m_y2;
	int m_b1;
	int m_b2;
	double minp;
	double eps;
	double m_bm1;
	ivector m_nb1;
	ivector m_nb2;
	dvector m_Wy;

	dvariable m_wss;	///> weighted sum of squares  
	dvariable m_nll;
	dvariable m_sigma;
	dvariable m_sigma2;

	imatrix m_nAgeIndex;

	dmatrix m_O;
	dmatrix m_Op;
	dmatrix m_Oz;
	dmatrix m_nAidx;
	dmatrix m_std_residual;

	dvar_vector m_rho;

	dvar_matrix m_E;
	dvar_matrix m_Ep;
	dvar_matrix m_Ez;
	dvar_matrix m_ww;
	dvar3_array m_V;

	logistic_normal();
	dvariable negative_log_likelihood();
	dvector compute_relative_weights(const dmatrix &O);
	
	void get_rho();
	void get_rho(const dvariable &phi);
	void get_rho(const dvariable &phi, const dvariable &psi);

	void std_residuals();

	void compute_correlation_array();
	void compute_likelihood_residuals();
	void compute_weighted_sumofsquares();
	void aggregate_and_compress_arrays();

	friend class logistic_student_t;

public:
	// Constructor
	logistic_normal(const dmatrix& _O,const dvar_matrix& _E,
	                const double _minp,const double _eps=0);

	// Four alternative methods for calculating the nll.
	dvariable operator() ();
	dvariable operator() (const dvariable &sigma2);
	dvariable operator() (const dvariable &sigma2 ,const dvariable &phi);
	dvariable operator() (const dvariable &sigma2 ,const dvariable &phi,
	                      const dvariable &psi);


	// Return the estimated (or mle) of the variance
	double      get_sigma () { return value(m_sigma ); }
	double      get_sigma2() { return value(m_sigma2); }
	dmatrix     get_standardized_residuals() { std_residuals(); return m_std_residual; }
	

};


class logistic_student_t: public logistic_normal
{
private:
	dvariable m_v;  // degrees of freedom.
	dvariable negative_log_likelihood();

public:
	logistic_student_t(const dmatrix& _O,const dvar_matrix& _E,
	                   const double _minp,const double _eps=0);

	dvariable operator() ();

};

// void aggregate(dmatrix Op, dvar_matrix Ep, const double &minp);
// dmatrix get_tail_compressed_index(const dmatrix &O, const double &minp);
// dvector compute_relative_weights(const dmatrix &O);
// dvar3_array compute_correlation_matrix(const dmatrix &n_Age);
// dvar3_array compute_correlation_matrix(const dmatrix &n_Age,const dvariable &phi);
// dvar_matrix compute_residual_difference(const dmatrix &O, const dvar_matrix &E);
// dvariable compute_weighted_sumofsquares(const dvector &Wy, const dvar_matrix &wwy,
//                                         const dvar3_array &V);
// dvariable nll_logistic_normal(const dmatrix &O, const dvar_matrix &E, 
//                               const double &minp, const double &eps,
//                               double &age_tau2);
// dvariable nll_logistic_normal(const dmatrix &O, const dvar_matrix &E, 
//                               const double &minp, const double &eps,
//                               double &age_tau2, const dvariable &phi);
// dvariable nll_logistic_normal(const dmatrix &O, const dvar_matrix &E, 
//                               const double &minp, const double &eps,
//                               const dvariable &theta, const dvariable &phi);

template <typename T>
void add_constant_normalize(T M, const double &eps)
{
	int i,y1,y2;
	y1 = M.rowmin();
	y2 = M.rowmax();
	for( i = y1; i <= y2; i++ )
	{
		M(i) = M(i) / sum(M(i)) + eps;
		M(i) = M(i) / sum(M(i));
	}
}

/**
 * Template fucntion to return the product of an array.
**/
template <typename T1,typename T>
T1 geomean(T x)
{	
	int lb = x.indexmin();
	int ub = x.indexmax();
	int n  = ub - lb + 1;
	T1 p = x(lb);
	for(int i = lb+1; i <= ub; i++ )
	{
		p *= x(i);
	}
	return(pow(p,1./n));
}

// template <typename T>
// T tail_compress(const T &O,const dmatrix &n_Age)
// {
// 	int i,y1,y2;
// 	int b1,b2;

// 	y1 = O.rowmin();
// 	y2 = O.rowmax();
// 	b1 = O.colmin();
// 	b2 = O.colmax();

// 	T P;
// 	P.allocate(n_Age);
// 	P.initialize();
// 	T px;
// 	px.allocate(O);
// 	px.initialize();
	
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		ivector ix = n_Age(i);
// 		px(i) = O(i)/sum(O(i));
// 		P(i)(min(ix),max(ix)) = px(i)(ix);
// 		if( min(ix) > b1 )
// 		{
// 			P(i)(min(ix)) = sum(px(i)(b1,min(ix)));
// 		}
// 		if( max(ix) < b2 )
// 		{
// 			P(i)(max(ix)) = sum(px(i)(max(ix),b2));
// 		}
// 	}
	
// 	return(P);
// }

#endif