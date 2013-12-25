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
	PSUEDOCODE FOR IMPLEMENTATION
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
#include "logistic_normal.h"

logistic_normal::~logistic_normal()
{

}

logistic_normal::logistic_normal()
{}

logistic_normal::logistic_normal(const dmatrix& _O,const dvar_matrix _E,
                                 const double _minProportion,const double& eps)
: m_O(_O), m_E(_E), m_eps(eps), m_dMinimumProportion(_minProportion)
{
	m_y1 = m_O.rowmin();
	m_y2 = m_O.rowmax();
	m_Y  = m_y2 - m_y1 + 1;

	m_b1 = m_O.colmin();
	m_b2 = m_O.colmax();


	// 1). Suppress zeros, if eps==0 use aggregate_and_compress_arrays
	if( m_eps > 0.0 )
	{
		// add constant
		add_constant(eps);
	}
	else if( eps == 0.0 )
	{
		// use aggregate_and_compress_arrays
		aggregate_and_compress_arrays();
	}

	

	// 3). Compute relative weights
	compute_relative_weights();

	// 4). Compute residual arrays
	compute_residual_arrays();

	// 5). Compute vector of covariance arrays for each year.
	compute_covariance_arrays();

	// 6). Compute the conditional mle of sigma.
	compute_mle_sigma();

	// 7). Compute negative loglikelihood.
	compute_negative_loglikelihood();

	// 8). Compute standardized residuals.
	compute_standardized_residuals();
	// cout<<"Here I am "<<endl;
}


template <typename T1, typename T2>
T2 prod(T1 x, T2 count)
{	
	T2 p = x(1);
	for(int i = 2; i <= count; i++ )
	{
		p *= x(i);
	}
	return(p);
}

/**
 * Compute standardized residuals.
 * res{_by} = [log(O_{by}/\tilde(O)) - log(E_{by}/\tilde(E))] / (W_y G_y^{0.5})
**/
void logistic_normal::compute_standardized_residuals()
{
	m_residual.allocate(m_O);
	m_residual.initialize();

	int i,j,k;
	double n;
	
	for( i = m_y1; i <= m_y2; i++ )
	{
		n             = m_nB2(i);
		double gOmean = pow(prod(m_Op(i),n),1./n);
		dvector   t1  = log(m_Op(i)/gOmean);

		double gEmean = pow(prod(value(m_Ep(i)),n),1./n);
		dvector   t2  = log(value(m_Ep(i))/gEmean);

		dvector sd    = sqrt(diagonal(m_S(i)));
		for( j = m_b1; j <= m_nB2(i); j++ )
		{
			k = m_nAgeIndex(i)(j);
			m_residual(i)(k) = (t1(j)-t2(j)) / (sd(j)*m_dWy(i));
		}
	}
}

/**
 * Compute negative loglikelihood:
 *      nll  = t1 + t2 + t3 + t4 + t5 + t6; where:
 *      bym1 = sum_y(B_y-1);
 *      t1   = 0.5*log(2pi) * bym1;
 *      t2   = sum_{by} log(O_{by});
 *      t3   = log(sigma) * bym1;
 *      t4   = 0.5*sum_y log(det(V_y))
 *      t5   = sum_y (By-1)*log(W_y)
 *      t6   = 0.5*sigma^{-2}* sum_y (w'_y V_y^{-1} w_y)/W_y^2
**/
void logistic_normal::compute_negative_loglikelihood()
{
	m_nll.initialize();

	double bym1 = sum(dvector(m_nB2) - 1.);
	double t1   = 0.5 * log(2. * PI) * bym1;
	double t2   = sum( log(m_Op) );
	dvariable t3   = log(m_sig) * bym1;
	dvariable t4   = 0;
	dvariable t5   = 0;
	dvariable t6   = 0;
	for(int i = m_y1; i <= m_y2; i++ )
	{
		t4 += 0.5*log( det(m_V(i)) );
		t5 += (m_nB2(i)-1) * log(m_dWy(i));
		t6 += 0.5/m_sig2 * (m_w(i) * inv(m_V(i)) * m_w(i)) / (m_dWy(i)*m_dWy(i));
	}

	m_nll = t1 + t2 + t3 + t4 + t5 + t6;
}



/**
 * Conditional maximum likelihood estimate of the variance scaler:
 *  sigma = [{sum_y (w_y' V^{-1} w_y)/ W_y^2} / sum_y (B_y-1) ]^0.5
 *
 * Conditional on the covariance matrix V and the weights W_y
**/
void logistic_normal::compute_mle_sigma()
{
	int i;
	//cout<<m_dWy<<endl;
	dvariable SS = 0;
	for( i = m_y1; i <= m_y2; i++ )
	{
		double wt = m_dWy(i) * m_dWy(i);
		SS       += ( m_w(i) * inv(m_V(i)) * m_w(i) ) / wt;
	}
	m_sig2 = SS/sum(m_nB2-1);
	m_sig  = pow(m_sig2,0.5);
}


/**
 * Compute the covariance arrays.
 * There are two arrays of interests here, one for the likelihoods
 * and the second for calculating the standardized residuals.

 This function should have a correlation matrix as an arg (C)

 * Computes the V_y and S_y arrays for likleihood and residuals.
**/
void logistic_normal::compute_covariance_arrays()
{
	// V_y = K C K'   (Eq. A3 in Francis paper)
	m_V.allocate(m_y1,m_y2,m_b1,m_nB2-1,m_b1,m_nB2-1);
	m_V.initialize();

	int i,j,k,nb;
	for( i = m_y1; i <= m_y2; i++ )
	{
		nb = m_nB2(i);
		dmatrix tK(1,nb,1,nb-1);
		dmatrix I      = identity_matrix(1,nb-1);
		tK.sub(1,nb-1) = I;
		tK(nb)         = -1;
		dmatrix K      = trans(tK);
		dmatrix C      = identity_matrix(1,nb);
		
		m_V(i) = K * C * tK;         	// Eq. A3 		
	}

	// G_y = F'H^{-1} VH^{-1} F; (Eq. A7 in Francis paper)
	m_S.allocate(m_y1,m_y2,m_b1,m_nB2,m_b1,m_nB2);
	m_S.initialize();

	for( i = m_y1; i <= m_y2; i++ )
	{
		nb            = m_nB2(i);
		dmatrix tF(1,nb,1,nb-1);
		dmatrix J(1,nb,1,nb);
		dmatrix I     = identity_matrix(1,nb-1);
		J             = 1;
		tF.sub(1,nb-1)= I;
		tF(nb)        = 1;
		
		dmatrix Hinv  = inv(I + 1);
		dmatrix FHinv = tF * Hinv;
		m_S(i)        = FHinv * value(m_V(i)) * trans(FHinv); 
	}

}


/**
 * Compute relative weights based on sample size in m_O;
**/
void logistic_normal::compute_relative_weights()
{
	m_dWy.allocate(m_y1,m_y2);
	dvector dNy = rowsum(m_O);	
	double dN   = mean(dNy);
	m_dWy       = sqrt( dN / dNy );
}



/**
 * Compute residuals between observed and expected proportions for likleihoods & residuals
**/
void logistic_normal::compute_residual_arrays()
{
	m_w.allocate(m_y1,m_y2,1,m_nB2-1);
	m_s.allocate(m_y1,m_y2,1,m_nB2);

	m_w.initialize();
	m_s.initialize();
	int i;
	for( i = m_y1; i <= m_y2; i++ )
	{
		// matrix of log_residuals for likelihood.
		int nB = m_nB2(i);
		m_w(i) = log(m_Op(i)(m_b1,nB-1)/m_Op(i,nB))
		        -log(m_Ep(i)(m_b1,nB-1)/m_Ep(i,nB));

		// matrix of log_residuals for standardized residuals.
		// based on multivariate normal
		m_s(i) = log(m_Op(i)) - log(m_Ep(i));
		m_s(i) = m_s(i) - mean(m_s(i));
	}
}



/**
 * Add small constant (eps) and renormalize arrays
**/
void logistic_normal::add_constant(const double& eps)
{
	int i;
	m_Op.allocate(m_O);
	m_Ep.allocate(m_E);
	m_nB2.allocate(m_y1,m_y2);
	m_nB2 = m_Op.colmax();

	for( i = m_y1; i <= m_y2; i++ )
	{
		m_O(i)  = m_O(i) + eps;
		m_Op(i) = m_O(i) / sum(m_O(i));

		m_E(i)  = m_E(i) + eps;
		m_Ep(i) = m_E(i) / sum(m_E(i));
	}
}


/**
 * Aggregate 0 observations into adjacent smaller bin and compress +group tail.
**/
void logistic_normal::aggregate_and_compress_arrays()
{
	
	m_nB2.allocate(m_y1,m_y2);
	int i;
	for( i = m_y1; i <= m_y2; i++ )
	{
		int n = 0;  // number of non-zero observations in each year
		for(int j = m_b1; j <= m_b2; j++ )
		{
			double op = m_O(i,j)/sum(m_O(i));
			if( op > m_dMinimumProportion )
			{
				n ++;
			}
		}
		m_nB2(i) = n;
	}
	m_nAgeIndex.allocate(m_y1,m_y2,m_b1,m_nB2);

	m_Op.allocate(m_y1,m_y2,m_b1,m_nB2);
	m_Ep.allocate(m_y1,m_y2,m_b1,m_nB2);
	m_Op.initialize();
	m_Ep.initialize();
	
	for( i = m_y1; i <= m_y2; i++ )
	{
		dvector     oo = m_O(i) / sum(m_O(i));
		dvar_vector ee = m_E(i) / sum(m_E(i));
		int k=m_b1;
		for(int j = m_b1; j <= m_b2; j++ )
		{
			if( oo(j) <= m_dMinimumProportion )  // Check this < gives -Inf need to compress left tails as well.
			{
				m_Op(i)(k) += oo(j);
				m_Ep(i)(k) += ee(j);
			}
			else
			{
				m_Op(i)(k) += oo(j);
				m_Ep(i)(k) += ee(j);
				if( k <=m_nB2(i) ) m_nAgeIndex(i,k) = k;
				if( k < m_nB2(i) ) k++;
			}
		}
	}
}





















logistic_normal::logistic_normal(const dmatrix& _O, const dvar_matrix& _E)
: m_O(_O), m_E(_E)
{
	/*
	O - observed numbers-at-age or proportions-at-age.
	E - expected numbers-at-age or porportions-at-age.

	Matrix rows correspond to years, cols age or length bins.
	*/
	// cout<<"Entered the constructor"<<endl;
	m_b1 = m_O.colmin();
	m_b2 = m_O.colmax();
	m_B  = m_b2 - m_b1 + 1;

	m_y1 = m_O.rowmin();
	m_y2 = m_O.rowmax();
	m_Y  = m_y2 - m_y1 + 1;

	m_dWy.allocate(m_y1,m_y2);
	m_dWy.initialize();

	// Observed matrixes
	m_Op.allocate(m_y1,m_y2,m_b1,m_b2);
	m_Ox.allocate(m_y1,m_y2,m_b1,m_b2);
	m_Oy.allocate(m_y1,m_y2,m_b1,m_b2);
	m_Oz.allocate(m_y1,m_y2,m_b1,m_b2);

	// Expected matrixes
	m_Ep.allocate(m_y1,m_y2,m_b1,m_b2);
	m_Ex.allocate(m_y1,m_y2,m_b1,m_b2);
	m_Ey.allocate(m_y1,m_y2,m_b1,m_b2);
	m_Ez.allocate(m_y1,m_y2,m_b1,m_b2);

	int i;
	int l = m_b1;
	int bm1 = m_b2-1;
	for( i = m_y1; i <= m_y2; i++ )
	{
		// Ensure proportions-at-age/size sum to 1.
		m_Op(i)        = m_O(i) / sum( m_O(i) );
		m_Ep(i)        = m_E(i) / sum( m_E(i) );

		
		// Logistic transformation based on Schnute's approach (Ox)
		m_Ox(i)        = exp(m_Op(i)) / sum( m_Op(i) );
		m_Ex(i)        = exp(m_Ep(i)) / sum( m_Ep(i) );

		
		// Logistic transformation using Aitichson's approach (Oy)
		m_Oy(i)(l,bm1) = exp(m_Op(i)(l,bm1)) / (1. + sum( exp(m_Op(i)(l,bm1)) ));
		m_Oy(i)(m_b2)  = 1. - sum( m_Oy(i)(l,bm1) );
		m_Ey(i)(l,bm1) = exp(m_Ep(i)(l,bm1)) / (1. + sum( exp(m_Ep(i)(l,bm1)) ));
		m_Ey(i)(m_b2)  = 1. - sum( m_Ey(i)(l,bm1) );

		
		// Logistic transformation using multivariate normal vector (Oz)
		double mu_O    = mean(m_Op(i));
		dvector   Z    = m_Op(i) - mu_O;
		m_Oz(i)        = exp(Z) / sum( exp(Z) );

		dvariable mu_E = mean(m_Ep(i));
		dvar_vector Zv = m_Ep(i) - mu_E;
		m_Ez(i)        = exp(Zv) / sum( exp(Zv) );
	}
	
	// Calculate mean weighting parameters for each year.
	dvector dNy = rowsum(m_O);	
	double dN   = mean(dNy);
	m_dWy       = sqrt( dN / dNy );

	
	// Minimum proportion to pool into adjacent cohort.
	m_dMinimumProportion = 0;
	m_rho = 0.1;
	
	//cout<<"Ok at the end of the constructor"<<endl;

}

dvariable logistic_normal::nll(const dvariable& sig2)
{
	/*
		Returns the negative loglikelihood based on the logistic normal
		based on the covariance matrix.

		-PSEUDOCODE:
		 1) suppress 0's or aggregate adjacent bins.
		 2) compute the covariance matrix
		 3) compute the negative log-likelihood
	*/

	// 1) suppress 0's or aggregate adjacent bins.
	aggregate_arrays();

	int i;
	m_nll = 0;
	for( i = m_y1; i <= m_y2; i++ )
	{
		// 2) compute covariance matrix
		int nB            = m_nNminp(i);
		dvar_matrix covar = calc_covmat(sig2,nB);
		dvar_matrix V     = calc_KCK(covar);
		dvar_matrix inv_V = inv(V);


		// 3) compute negative log-likelihood
		m_nll += 0.5 * (nB-1.) * log(2.*PI);
		m_nll += sum(log(m_Oa(i)));
		m_nll += 0.5 * log( det(V) );
		m_nll += (nB-1.) * log(m_dWy(i));
		m_nll += (0.5 / (m_dWy(i)*m_dWy(i))) * m_w(i) * inv_V * m_w(i);
		
	}
	return(m_nll);
}

/**
	Calculate conditional maximum likelihood estimate of the 
	variance parameter (sig2) based on the correlation matrix (m_C)
*/
dvariable logistic_normal::calc_sigma_mle()
{
	int i;
	for( i = m_y1; i <= m_y2; i++ )
	{
		
	}
	return(m_sig2);
}

dvar_matrix logistic_normal::calc_covmat(const dvariable& sig2,const int& n)
{
	/*
		Calculate the covariance matrix give sig2
	
	*/
	int i;
	dvar_matrix Vmat(1,n,1,n);
	Vmat.initialize();
	for( i = 1; i <= n; i++ )
	{
		 Vmat(i,i) = sig2;
	}
	return(Vmat);
}

dvar_matrix logistic_normal::calc_KCK(const dvar_matrix& V)
{
	int j;
	int n = V.colmax() - V.colmin() + 1;
	dmatrix I=identity_matrix(1,n-1);
	dmatrix K(1,n-1,1,n);
	for( j = 1; j <= n-1; j++ )
	{
		 K.colfill(j, extract_column(I,j) );
	}
	dvector tmp(1,n-1);
	tmp = -1;
	K.colfill(n, tmp);

	dvar_matrix KCK = K * V * trans(K);
	
	return(KCK); 
}


// dvariable logistic_normal::negative_loglikelihood(const dvariable& tau2)
// {

// 	aggregate_arrays();

// 	m_nll = 0;

// 	int i;
// 	m_sig2=tau2;

// 	for( i = m_y1; i <= m_y2; i++ )
// 	{	
// 		int    nB     = m_nNminp(i);
// 		double t1     = 0.5 * (nB - 1.0) * log(2.0*PI);
// 		double t3     = sum( log(m_Oa(i)) );
// 		dvariable det = pow(nB * m_sig2, 2.*(nB-1.));
// 		dvariable t5  = 0.5 * log(det);
// 		double t7     = (nB - 1.0) * log(m_dWy(i));
// 		dvariable t9  = 0.5 * sum( square(m_w(i)) / (m_sig2 * square(m_dWy(i))) );
		
// 		m_nll        += t1 + t3 + t5 + t7 + t9;
// 	}
		
// 	return(m_nll);
// }

// dvariable logistic_normal::negative_loglikelihood()
// {
// 	// negative loglikelihood evaluated at the MLE of the variance tau2.
// 	m_nll = 0;
// 	cout<<"Integrate over the variance"<<endl;
// 	int i;
// 	correlation_matrix();
// 	m_sig2 = 0;
// 	dvariable sws = 0;
// 	for( i = m_y1; i <= m_y2; i++ )
// 	{
// 		cout<<m_w(i)<<endl<<endl;
// 		cout<<m_V(i)<<endl<<endl;
// 		cout<<m_w(i)*m_V(i)<<endl<<endl;
// 		cout<<(m_w(i)*m_V(i))*m_w(i)<<endl;
// 		sws += (m_w(i) * m_V(i)) * m_w(i);		///> equation A11
// 		exit(1);
// 	}

// 	double t1 = 0.5 * log(2.0*PI) * sum( m_nNminp-1 );
// 	double t2 = sum( log(m_Oa) );

// 	return(m_nll);
// }

dvar_matrix logistic_normal::standardized_residuals()
{
	/*
		Standardized residuals
		nu = [log(O/~O) - log(E/~E)]/(Wy m_covar2^0.5)
		where O & E are the observed and expected proportion vectors
		~O and ~E is the geometric means of each proportion vector
		Wy is the relative weight each year.

		m_covar is
		
	*/
	int i,j;
	double    tilde_O;
	dvariable tilde_E;
	m_residual.allocate(m_Op);
	m_residual.initialize();
	for( i = m_y1; i <= m_y2; i++ )
	{
		dvector     oo = m_Oa(i);
		dvar_vector pp = m_Ea(i);
		tilde_O = exp( mean(log(oo)) );
		tilde_E = exp( mean(log(pp)) );

		dvar_vector	w = log(oo/tilde_O) - log(pp/tilde_E);
		w            /= m_sig2/m_dWy(i);
		
		for( j = 1; j <= m_nNminp(i); j++ )
		{
			int idx = m_nAgeIndex(i,j);
			m_residual(i)(idx) = value(w(j));
		}
		
	}
	
	return(m_residual);
}


void logistic_normal::aggregate_arrays()
{
	/*
		- Aggregate minimum proportions in to adjacent cohorts (tail compression)
		- This routine populates m_Oa and m_Ea as ragged objects to be used for
		  likelihood calculations.
	*/
	m_nNminp.allocate(m_y1,m_y2);
	

	int k;
	int nB;
	for(int i = m_y1; i <= m_y2; i++ )
	{
		int n = 0;  // number of non-zero observations in each year
		for(int j = m_b1; j <= m_b2; j++ )
		{
			if( m_Op(i,j) > m_dMinimumProportion )
			{
				n ++;
			}
		}
		m_nNminp(i) = n;
	}

	// Ragged arrays with tail compression and zeros omitted.
	m_nAgeIndex.allocate(m_y1,m_y2,1,m_nNminp);
	m_Oa.allocate(m_y1,m_y2,1,m_nNminp);
	m_Ea.allocate(m_y1,m_y2,1,m_nNminp);
	m_w.allocate(m_y1,m_y2,1,m_nNminp-1);
	
	m_nAgeIndex.initialize();
	m_Oa.initialize();
	m_Ea.initialize();
	m_w.initialize();

	// cout<<"Minimum proportion = "<<m_dMinimumProportion<<endl;
	for(int i = m_y1; i <= m_y2; i++ )
	{
		dvector     oo = m_Oy(i);
		dvar_vector pp = m_Ey(i);
		k = 1;
		for(int j = m_b1; j <= m_b2; j++ )
		{
			if( oo(j) <= m_dMinimumProportion )
			{
				m_Oa(i,k) += oo(j);
				m_Ea(i,k) += pp(j);
			}
			else
			{
				m_Oa(i,k) += oo(j);
				m_Ea(i,k) += pp(j);

				if(k <=m_nNminp(i)) m_nAgeIndex(i,k) = k;
				if(k < m_nNminp(i)) k++;
			}
		}

		// compressed matrix of log-residuals
		nB = m_nNminp(i);
		m_w(i) = log(m_Oa(i)(1,nB-1)/m_Oa(i,nB))
		        -log(m_Ea(i)(1,nB-1)/m_Ea(i,nB));
		
	}
}


/**
	Calculate the correlation matrix for each year based on the tail compressed 
	arrays.
*/
void logistic_normal::correlation_matrix()
{
	// Calculate the covariance matrix for each year, based on ragged arrays m_w;

	int i,j,k,n;
	m_V.allocate(m_y1,m_y2,1,m_nNminp,1,m_nNminp);
	m_C.allocate(m_y1,m_y2,1,m_nNminp,1,m_nNminp);
	m_V.initialize();
	m_C.initialize();


	for( i = m_y1; i <= m_y2; i++ )
	{
		n = m_nNminp(i);
		dmatrix I=identity_matrix(1,n-1);
		dmatrix K(1,n-1,1,n);
		for( j = 1; j <= n-1; j++ )
		{
			 K.colfill(j, extract_column(I,j) );
		}
		dvector tmp(1,n-1);
		tmp = -1;
		K.colfill(n, tmp);
		 //dmatrix K = identity_matrix(1,m_nNminp(i));
		 //dvar_matrix C(1,m_nNminp(i),1,m_nNminp(i));
		 for( j = 1; j <= m_nNminp(i); j++ )
		 {
		 	 for( k = 1; k <= m_nNminp(i); k++ )
		 	 {
		 	 	m_C(i)(j,k) = pow(m_rho,abs(double(j)-double(k)));
		 	 }
		 }
		 m_V(i) = K * m_C(i) * trans(K);
	}	

}













