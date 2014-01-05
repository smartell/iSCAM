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
#include "logistic_normal.h"

logistic_normal::~logistic_normal()
{

}

logistic_normal::logistic_normal()
{}

logistic_normal::logistic_normal(const dmatrix *_O,const dvar_matrix *_E,
                                 const double _minProportion,const double _eps)
: m_O(*_O), m_E(*_E), m_eps(_eps), m_dMinimumProportion(_minProportion)
{
	m_y1 = m_O.rowmin();
	m_y2 = m_O.rowmax();
	m_Y  = m_y2 - m_y1 + 1;

	m_b1 = m_O.colmin();
	m_b2 = m_O.colmax();

	if(m_E.colmin() != m_b1 || m_E.colmax() != m_b2 )
	{
		cerr<<"Observed and Predicted composition columns do not match"<<endl;
		ad_exit(1);
	}
	if(m_E.rowmin() != m_y1 || m_E.rowmax() != m_y2 )
	{
		cerr<<"Observed and Predicted composition rows do not match"<<endl;
		ad_exit(1);
	}

	// 1). Suppress zeros, if eps==0 use aggregate_and_compress_arrays
	if( m_eps > 0.0 )
	{
		// add constant
		add_constant(m_eps);

	}
	else if( m_eps == 0.0 && m_dMinimumProportion >= 0 )
	{
		// use aggregate_and_compress_arrays
		aggregate_and_compress_arrays();
	}

	// 3). Compute relative weights (m_dWy)
	compute_relative_weights();
	
}

/**
 * No autocorrelation LN1
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
	// 4). Compute residual arrays (specific only to the likelihood)
	compute_residual_arrays();

	// 5). Compute vector of covariance arrays for each year. (V_y)
	compute_correlation_matrix();

	// 6). Compute the conditional mle of sigma.
	compute_mle_sigma();
	// 7). Compute negative loglikelihood.
	m_nll.initialize();
	
	double bym1 = sum( dvector(m_nB2-m_b1) - 1.);
	double t1   = 0.5 * log(2. * PI) * bym1;
	double t2   = sum( log(m_Op) );
	dvariable t3   = log(m_sig) * bym1;
	dvariable t4   = 0;
	dvariable t5   = 0;
	dvariable t6   = 0;
	for(int i = m_y1; i <= m_y2; i++ )
	{
		// scale covariance matrix
		for(int j = m_b1; j < m_nB2(i); j++ )
		{
			m_V(i)(j,j) *= m_sig2;
		}
		dvar_matrix Vinv = inv(m_V(i));

		t4 += 0.5*log( det(m_V(i)) );
		t5 += (m_nB2(i)-m_b1-1) * log(m_dWy(i));
		t6 += (0.5/(m_dWy(i)*m_dWy(i))) * m_w(i) * Vinv * m_w(i);
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

	dvariable SS = 0;
	double wt2;
	for( i = m_y1; i <= m_y2; i++ )
	{
		wt2 = m_dWy(i) * m_dWy(i);
		SS += ( m_w(i) * inv(m_V(i)) * m_w(i) ) / wt2;
	}
	m_sig2 = SS/sum(m_nB2-1);
	m_sig  = pow(m_sig2,0.5);
}


/**
 * Compute covariance array for no autocorrelation case
 * This just sets m_V to an identity matrix + 1.
**/
void logistic_normal::compute_correlation_matrix()
{
	m_V.allocate(m_y1,m_y2,m_b1,m_nB2-1,m_b1,m_nB2-1);
	m_V.initialize();

	int i,nb;
	for( i = m_y1; i <= m_y2; i++ )
	{
		nb = m_nB2(i);
		dmatrix I = identity_matrix(m_b1,nb-1);
		m_V(i) = 1 + I;
	}
}


/**
 * Compute standardized residuals.
 * res{_by} = [log(O_{by}/\tilde(O)) - log(E_{by}/\tilde(E))] / (W_y G_y^{0.5})
**/
void logistic_normal::compute_standardized_residuals()
{

	// Assumed the covariance matrix (m_V) has already been calculated.
	// when the likelihood was called.
	

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

		dmatrix  I = identity_matrix(m_b1,n-1);

		dmatrix tF(m_b1,n,m_b1,n-1);
		tF.sub(m_b1,n-1) = I;
		tF(n)         = 1;

		dmatrix  Hinv    = inv(I + 1);
		dmatrix FHinv    = tF * Hinv;
		dmatrix     G    = FHinv * value(m_V(i)) * trans(FHinv);
		
		dvector sd    = sqrt(diagonal(G));
		for( j = m_b1; j <= m_nB2(i); j++ )
		{
			k = m_nAgeIndex(i)(j);
			m_residual(i)(k) = (t1(j)-t2(j)) / (sd(j)*m_dWy(i));
		}
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
	m_w.allocate(m_y1,m_y2,m_b1,m_nB2-1);
	// m_s.allocate(m_y1,m_y2,1,m_nB2);

	m_w.initialize();
	// m_s.initialize();
	int i;
	for( i = m_y1; i <= m_y2; i++ )
	{
		// matrix of log_residuals for likelihood.
		int nB = m_nB2(i);
		
		m_w(i) = log(m_Op(i)(m_b1,nB-1)/m_Op(i,nB))
		        -log(m_Ep(i)(m_b1,nB-1)/m_Ep(i,nB));
		
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
	m_nAgeIndex.allocate(m_y1,m_y2,m_b1,m_nB2);

	for( i = m_y1; i <= m_y2; i++ )
	{
		m_O(i)  = m_O(i) + eps;
		m_Op(i) = m_O(i) / sum(m_O(i));

		m_E(i)  = m_E(i) + eps;
		m_Ep(i) = m_E(i) / sum(m_E(i));

		m_nAgeIndex(i).fill_seqadd(m_b1,1);
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
		int n = m_b1-1;  // index for maximum column each year for ragged objects
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








