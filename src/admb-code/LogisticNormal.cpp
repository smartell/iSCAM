#include <admodel.h>
#include "LogisticNormal.h"

/*
	Implementation of the logistic normal negative loglikelihood.

	The overloaded operator () can have the following args:
	()             -no args, and compute the MLE of the variance
	(sig)          -estimate the variance
	(sig,phi1)     -AR1 process
	(sig,phi1,phi2)-AR2 process
*/


/**
 * Psuedocode:
 * 1) (aggregate || add constant) && compress tails
 * 2) Compute relative weights for each year W_y
 * 3) Compute covariance matrix V_y
 * 4) Compute residual differences (w_y)
 * 5) Compute weighted sum of squares (wSS)
 * 6) Compute MLE of variance
 * 7) Compute nll_logistic_normal
**/

// Constructor
logistic_normal::logistic_normal(const dmatrix& _O,const dvar_matrix& _E,
	                			const double _minp,const double _eps)
: minp(_minp),eps(_eps),m_O(_O),m_E(_E)
{
	m_y1 = m_O.rowmin();
	m_y2 = m_O.rowmax();
	m_b1 = m_O.colmin();
	m_b2 = m_O.colmax();

	m_nAidx = get_tail_compressed_index(m_O,minp);
	m_Op    = tail_compress(m_O,m_nAidx);
	m_Ep    = tail_compress(m_E,m_nAidx);
	if( eps )
	{
		add_constant_normalize(m_Op,eps);
		add_constant_normalize(m_Ep,eps);
	}
	else if ( eps == 0 )
	{
		aggregate(m_Op,m_Ep,minp);
	}

	m_nb1.allocate(m_y1,m_y2);
	m_nb2.allocate(m_y1,m_y2);
	for(int i = m_y1; i <= m_y2; i++ )
	{
		m_nb1(i) = min(m_nAidx(i));
		m_nb2(i) = max(m_nAidx(i));
	}
	m_V.allocate(m_y1,m_y2,m_nb1,m_nb2-1,m_nb1,m_nb2-1);
	m_V.initialize();

	m_bm1 = size_count(m_Op) - (m_y2-m_y1+1.0);

	m_Wy = compute_relative_weights(m_O);
	
	compute_likelihood_residuals();
	

}

void logistic_normal::compute_likelihood_residuals()
{
	int i;
	m_ww.allocate(m_y1,m_y2,m_nb1,m_nb2-1);
	m_ww.initialize();
	for( i = m_y1; i <= m_y2; i++ )
	{
		int l = m_nb1(i);
		int u = m_nb2(i);
		m_ww(i) = (log(m_Op(i)(l,u-1)) - log(m_Op(i,u)))
				- (log(m_Ep(i)(l,u-1)) - log(m_Ep(i,u)));
	}
}


/**
 * Returns the negative loglikelihood for the logistic normal distribution
 * assuming no autocorrelation in the composition data.
**/
dvariable logistic_normal::operator() ()
{
	m_nll = 0;

	// Construct covariance (m_V)
	compute_correlation_array();

	// Compute weighted sum of squares
	compute_weighted_sumofsquares();

	// mle of the variance
	m_sigma2 = m_wss / m_bm1; 
	m_sigma  = sqrt(m_sigma2);

	// compute negative loglikelihood
	m_nll = negative_log_likelihood();
	return m_nll;
}

/**
 * Returns the negative loglikelihood for the logistic normal distribution
 * assuming AR1 autocorrelation in the composition data.
**/
dvariable logistic_normal::operator() (const prevariable phi)
{
	m_nll = 0;

	// Construct covariance (m_V)
	compute_correlation_array(phi);

	// Compute weighted sum of squares
	compute_weighted_sumofsquares();

	// mle of the variance
	m_sigma2 = m_wss / m_bm1; 
	m_sigma  = sqrt(m_sigma2);

	// compute negative loglikelihood
	m_nll = negative_log_likelihood();
	return m_nll;
}





dvariable logistic_normal::negative_log_likelihood()
{
	// 7) Compute nll_logistic_normal
	RETURN_ARRAYS_INCREMENT();
	dvariable nll;
	nll  = 0.5 * log(2.0 * PI) * m_bm1;
	nll += sum( log(m_Op) );
	nll += log(m_sigma) * m_bm1;

	for(int i = m_y1; i <= m_y2; i++ )
	{
		nll += 0.5 * log(det(m_V(i)));
		nll += (size_count(m_Op(i))-1) * log(m_Wy(i));
	}
	nll += 0.5 / m_sigma2 * m_wss;
	RETURN_ARRAYS_DECREMENT();
	return nll;
}

void logistic_normal::compute_weighted_sumofsquares()
{
	int i;
	m_wss=0;
	for( i = m_y1; i <= m_y2; i++ )
	{
		dvar_matrix Vinv = inv(m_V(i));
		m_wss += (m_ww(i) * Vinv * m_ww(i)) / (m_Wy(i) * m_Wy(i));
	}
}

void logistic_normal::compute_correlation_array()
{
	int i;
	for( i = m_y1; i <= m_y2; i++ )
	{
		m_V(i) = 1 + identity_matrix(m_nb1(i),m_nb2(i)-1);
	}
}

void logistic_normal::compute_correlation_array(const prevariable phi)
{
	int i,j,k;
	RETURN_ARRAYS_INCREMENT();
	for( i = m_y1; i <= m_y2; i++ )
	{
		int l = m_nb1(i);
		int u = m_nb2(i);
		
		dvar_vector rho(l,u);
		for( j = l, k = 1; j <= u; j++, k++ )
		{
			rho(j) = pow(phi,k);
			// rho(j) = mfexp( k * log(*phi) );
		}

		dvar_matrix C = identity_matrix(l,u);
		for( j = l; j <= u; j++ )
		{
			for( k = l; k <= u; k++ )
			{
				if(j != k) C(j,k) = rho(l-1+abs(j-k));
			}
		}

		dmatrix I = identity_matrix(l,u-1);
		dmatrix tK(l,u,l,u-1);
		tK.sub(l,u-1) = I;
		tK(u)         = -1;
		dmatrix K = trans(tK);
		m_V(i) = K * C * tK;
	}
	RETURN_ARRAYS_DECREMENT();	
}











dvariable nll_logistic_normal(const dmatrix &O, const dvar_matrix &E, 
                              const double &minp, const double &eps,
                              double &age_tau2)
{
	int i,y1,y2;
	RETURN_ARRAYS_INCREMENT();
	y1 = O.rowmin();
	y2 = O.rowmax();

	// 1) (aggregate || add constant) && compress tails
	dmatrix n_Age  = get_tail_compressed_index(O,minp);
	dmatrix Op     = tail_compress(O,n_Age);
	dvar_matrix Ep = tail_compress(E,n_Age);
	
	
	if( eps )
	{
		// cout<<"adding constant"<<endl;
		add_constant_normalize(Op,eps);
		add_constant_normalize(Ep,eps);
	}
	else
	{
		// cout<<"aggregating cohorts"<<endl;
		aggregate(Op,Ep,minp);
	}

	// 2) Compute relative weights for each year W_y
	dvector Wy = compute_relative_weights(O);
	
	// 3) Compute covariance matrix V_y = K C K'
	dvar3_array Vy = compute_correlation_matrix(n_Age);
	
	// 4) Compute residual differences (w_y)
	dvar_matrix wwy = compute_residual_difference(Op,Ep);

	// 5) Compute weighted sum of squares (wSS)
	dvariable ssw = compute_weighted_sumofsquares(Wy,wwy,Vy);

	// 6) Compute MLE of variance
	double bm1 = size_count(Op) - (y2-y1+1.0);
	dvariable sigma2 = ssw / bm1;
	dvariable sigma  = sqrt(sigma2);
	age_tau2         = value(sigma2);
	
	// 7) Compute nll_logistic_normal
	dvariable nll;
	nll  = 0.5 * log(2.0 * PI) * bm1;
	nll += sum( log(Op) );
	nll += log(sigma) * bm1;
	for( i = y1; i <= y2; i++ )
	{
		nll += 0.5 * log(det(Vy(i)));
		nll += (size_count(Op(i))-1) * log(Wy(i));
	}
	nll += 0.5 / sigma2 * ssw;

	RETURN_ARRAYS_DECREMENT();
	return nll;
}

/**
 * Implementing the LN2 version.
**/
dvariable nll_logistic_normal(const dmatrix &O, const dvar_matrix &E, 
                              const double &minp, const double &eps,
                              double &age_tau2, const dvariable &phi)
{
	int i,y1,y2;
	RETURN_ARRAYS_INCREMENT();
	y1 = O.rowmin();
	y2 = O.rowmax();

	// 1) (aggregate || add constant) && compress tails
	dmatrix n_Age  = get_tail_compressed_index(O,minp);
	dmatrix Op     = tail_compress(O,n_Age);
	dvar_matrix Ep = tail_compress(E,n_Age);
	
	
	if( eps )
	{
		// cout<<"adding constant"<<endl;
		add_constant_normalize(Op,eps);
		add_constant_normalize(Ep,eps);
	}
	else
	{
		// cout<<"aggregating cohorts"<<endl;
		aggregate(Op,Ep,minp);
	}

	// 2) Compute relative weights for each year W_y
	dvector Wy = compute_relative_weights(O);
	
	// 3) Compute covariance matrix V_y = K C K'
	dvar3_array Vy = compute_correlation_matrix(n_Age,phi);
	
	// 4) Compute residual differences (w_y)
	dvar_matrix wwy = compute_residual_difference(Op,Ep);

	// 5) Compute weighted sum of squares (wSS)
	dvariable ssw = compute_weighted_sumofsquares(Wy,wwy,Vy);

	// 6) Compute MLE of variance
	double bm1 = size_count(Op) - (y2-y1+1.0);
	dvariable sigma2 = ssw / bm1;
	dvariable sigma  = sqrt(sigma2);
	age_tau2         = value(sigma2);
	
	// 7) Compute nll_logistic_normal
	dvariable nll;
	nll  = 0.5 * log(2.0 * PI) * bm1;
	nll += sum( log(Op) );
	nll += log(sigma) * bm1;
	for( i = y1; i <= y2; i++ )
	{
		nll += 0.5 * log(det(Vy(i)));
		nll += (size_count(Op(i))-1) * log(Wy(i));
	}
	nll += 0.5 / sigma2 * ssw;

	RETURN_ARRAYS_DECREMENT();
	return nll;
}








/**
 * What is done here is that proportions less than minp are pooled with the adjacent
 * younger cohort, and split evenly to each age, and the sample size is halved.
 * This approximates the likelihood of ages 4 and 5 (eg.) as ages 4-5.
**/
void aggregate(dmatrix Op, dvar_matrix Ep, const double &minp)
{
	int i,y1,y2;
	int j,b1,b2;

	y1 = Op.rowmin();
	y2 = Op.rowmax();

	for( i = y1; i <= y2; i++ )
	{
		b1 = Op(i).indexmin();
		b2 = Op(i).indexmax();
		for( j = b1+1; j <= b2; j++ )
		{
			if( Op(i,j) < minp )
			{
				double    tmpO = 0.5*(Op(i,j) + Op(i,j-1));
				dvariable tmpE = 0.5*(Ep(i,j) + Ep(i,j-1));
				Op(i)(j-1,j) = tmpO;
				Ep(i)(j-1,j) = tmpE;
			}
		}
	}
}

dvariable compute_weighted_sumofsquares(const dvector &Wy, 
                                        const dvar_matrix &wwy,
                                        const dvar3_array &V)
{
	int i,y1,y2;
	RETURN_ARRAYS_INCREMENT();
	y1 = wwy.rowmin();
	y2 = wwy.rowmax();

	dvariable ssw = 0;
	for( i = y1; i <= y2; i++ )
	{
		dvar_matrix Vinv = inv(V(i));
		ssw += (wwy(i) * Vinv * wwy(i))/ (Wy(i)*Wy(i));
	}

	RETURN_ARRAYS_DECREMENT();
	return ssw;
}

dvar_matrix compute_residual_difference(const dmatrix &O, const dvar_matrix &E)
{
	int i,y1,y2;
	
	RETURN_ARRAYS_INCREMENT();

	y1 = O.rowmin();
	y2 = O.rowmax();
	ivector b1(y1,y2);
	ivector b2(y1,y2);
	for( i = y1; i <= y2; i++ )
	{
		b1(i) = O(i).indexmin();
		b2(i) = O(i).indexmax();
	}
	
	dvar_matrix ww(y1,y2,b1,b2-1);
	ww.initialize();
	
	for( i = y1; i <= y2; i++ )
	{
		int l = b1(i);
		int u = b2(i);
		
		dvector     t1 = O(i);
		dvar_vector t2 = E(i);
		ww(i) = (log(t1(l,u-1)) - log(t1(u)))
		      - (log(t2(l,u-1)) - log(t2(u)));
	}

	RETURN_ARRAYS_DECREMENT();
	return(ww);
}

/**
 * No Autocorrelation case
**/
dvar3_array compute_correlation_matrix(const dmatrix &n_Age)
{
	int i,y1,y2;
	RETURN_ARRAYS_INCREMENT();
	y1 = n_Age.rowmin();
	y2 = n_Age.rowmax();
	ivector b1(y1,y2);
	ivector b2(y1,y2);
	for( i = y1; i <= y2; i++ )
	{
		b1(i) = min(n_Age(i));
		b2(i) = max(n_Age(i));
	}

	dvar3_array V;
	V.allocate(y1,y2,b1,b2-1,b1,b2-1);
	V.initialize();
	
	// Now compute the correlation matrix C, and set V = K C K'
	// Which is just the 1+identity matrix (i.e. no covariance structure)
	for( i = y1; i <= y2; i++ )
	{
		V(i) = 1 + identity_matrix(b1(i),b2(i)-1);
	}
	RETURN_ARRAYS_DECREMENT();
	return (V);
}
/**
 * AR1 case
**/
dvar3_array compute_correlation_matrix(const dmatrix &n_Age,const dvariable &phi)
{
	int i,y1,y2;
	int j,k;
	RETURN_ARRAYS_INCREMENT();
	y1 = n_Age.rowmin();
	y2 = n_Age.rowmax();
	ivector b1(y1,y2);
	ivector b2(y1,y2);
	for( i = y1; i <= y2; i++ )
	{
		b1(i) = min(n_Age(i));
		b2(i) = max(n_Age(i));
	}

	dvar3_array V;
	V.allocate(y1,y2,b1,b2-1,b1,b2-1);
	V.initialize();
	
	// Now compute the correlation matrix C, and set V = K C K'
	for( i = y1; i <= y2; i++ )
	{
		int l = b1(i);
		int u = b2(i);
		dvar_vector rho(l,u);

		for( j = l, k = 1; j <= u; j++, k++ )
		{
			rho(j) = pow(phi,k);
			// rho(j) = phi / k;
			// rho(j) = mfexp( k * log(phi) );
		}
		
		dvar_matrix C = identity_matrix(l,u);
		for( j = l; j <= u; j++ )
		{
			for( k = l; k <= u; k++ )
			{
				if(j != k) C(j,k) = rho(l-1+abs(j-k));
			}
		}
		
		dmatrix I = identity_matrix(l,u-1);
		dmatrix tK(l,u,l,u-1);
		tK.sub(l,u-1) = I;
		tK(u)         = -1.0;
		dmatrix K = trans(tK);
		
		V(i) = K * C * tK;
		
	}
	RETURN_ARRAYS_DECREMENT();
	return (V);
}



dvector compute_relative_weights(const dmatrix &O)
{
	dvector Wy(O.rowmin(),O.rowmax());
	Wy = rowsum(O);
	double mu = mean(Wy);
	Wy = sqrt(mu / Wy);
	return(Wy);
}



dmatrix get_tail_compressed_index(const dmatrix &O, const double &minp)
{
	int i,y1,y2;
	int j,b1,b2;

	y1 = O.rowmin();
	y2 = O.rowmax();
	b1 = O.colmin();
	b2 = O.colmax();
	dmatrix p(y1,y2,b1,b2);
	ivector n_B1(y1,y2);
	ivector n_B2(y1,y2);
	n_B1.initialize();
	n_B2.initialize();

	for( i = y1; i <= y2; i++ )
	{
		bool blt = true;    //left tail
		n_B1(i) = b1;
		n_B2(i) = b1-1;
		p(i) = O(i)/sum(O(i));
		for( j = b1; j <= b2; j++ )
		{
			if( p(i,j) <= minp && blt ) n_B1(i) ++; else blt=false;
			if( p(i,j) >  minp ) n_B2(i) ++;
		}	
	}
	
	dmatrix n_Age(y1,y2,n_B1,n_B2);
	n_Age.initialize();
	
	for( i = y1; i <= y2; i++ )
	{
		int k = n_B1(i);
		for( j = b1; j <= b2; j++ )
		{
			if( p(i,j) > minp )
			{
				if( k <= n_B2(i) ) n_Age(i,k) = k;
				if( k <  n_B2(i) ) k++;
			}
		}
	}
	
	return (n_Age);
}