#include <admodel.h>
#include "LogisticNormal.h"

/**
	Implementation of the logistic normal negative loglikelihood.

	The overloaded operator () can have the following args:
	()             -no args, anompute the MLE of the variance
	(sig)          -estimate the variance
	(sig,phi1)     -AR1 process
	(sig,phi1,phi2)-AR2 procd cess

	TODO:

	
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

/** Constructor */
logistic_normal::logistic_normal(const dmatrix& _O,const dvar_matrix& _E,
	                			const double _minp,const double _eps)
: minp(_minp),eps(_eps),m_O(_O),m_E(_E)
{
	m_y1 = m_O.rowmin();
	m_y2 = m_O.rowmax();
	m_b1 = m_O.colmin();
	m_b2 = m_O.colmax();

	// 1). Aggregate and compress arrays, add small constant if eps > 0
	if(eps)
	{	
		add_constant_normalize(m_O,eps);
		add_constant_normalize(m_E,eps);
	}
	aggregate_and_compress_arrays();

	
	// Now that we now the dimensions of each array, allocate
	// memory for the covariance matrixes.
	m_V.allocate(m_y1,m_y2,m_nb1,m_nb2-1,m_nb1,m_nb2-1);
	m_V.initialize();


	// Total number of bins minus 1.
	double Y = m_y2-m_y1+1.0;
	m_bm1 = Y*(size_count(m_Op) - Y);
	
	// Relative weights to assign to each year.
	m_Wy = compute_relative_weights(m_O);
	
	// Residuals for use in likelihood calculations
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

	// Get correlation vector rho
	get_rho();

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
 * assuming no autocorrelation and estimates the variance parameter sigma.
**/
dvariable logistic_normal::operator() (const dvariable &sigma2)
{
	m_nll = 0;

	// Get correlation vector rho
	get_rho();
	
	// Construct covariance (m_V)
	compute_correlation_array();

	// Compute weighted sum of squares
	compute_weighted_sumofsquares();

	// estimated variance
	m_sigma2 = sigma2;
	m_sigma = sqrt(m_sigma2);

	// compute negative loglikelihood
	m_nll = negative_log_likelihood();
	return m_nll;
}


/**
 * Returns the negative loglikelihood for the logistic normal distribution
 * assuming AR1 autocorrelation in the composition data.
**/
dvariable logistic_normal::operator() (const dvariable &sigma2,const dvariable &phi)
{
	m_nll = 0;
	cout<<phi<<endl;
	// Get correlation vector rho
	get_rho(phi);

	// Construct covariance (m_V)
	compute_correlation_array();

	// Compute weighted sum of squares
	compute_weighted_sumofsquares();

	// estimated variance
	m_sigma2  = sigma2;// / (1.0-phi);
	m_sigma   = sqrt(m_sigma2);

	// compute negative loglikelihood
	m_nll = negative_log_likelihood();
	return m_nll;
}

/**
 * Returns the negative loglikelihood for the logistic normal distribution
 * assuming AR1 autocorrelation in the composition data.
**/
dvariable logistic_normal::operator() (const dvariable &sigma2,const dvariable &phi,
                                       const dvariable &psi)
{
	m_nll = 0;

	// Get correlation vector rho
	get_rho(phi,psi);
	// cout<<m_rho<<endl;

	// Construct covariance (m_V)
	compute_correlation_array();

	// Compute weighted sum of squares
	compute_weighted_sumofsquares();

	// estimated variance
	m_sigma2  = sigma2;// / (1.0-phi);
	m_sigma   = sqrt(m_sigma2);

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

// void logistic_normal::compute_correlation_array()
// {
// 	int i;
// 	for( i = m_y1; i <= m_y2; i++ )
// 	{
// 		m_V(i) = 1 + identity_matrix(m_nb1(i),m_nb2(i)-1);
// 	}
// }

void logistic_normal::get_rho()
{
	m_rho.allocate(min(m_nb1),max(m_nb2));
	m_rho = 1.0;
}

void logistic_normal::get_rho(const dvariable &phi)
{
	int j,k;
	m_rho.allocate(min(m_nb1),max(m_nb2));
	m_rho.initialize();
	for( j = min(m_nb1), k=1; j <= max(m_nb2); j++, k++ )
	{
		m_rho(j) = pow(phi,k);
	}
}

void logistic_normal::get_rho(const dvariable &phi, const dvariable &psi)
{
	int j;
	int lb = min(m_nb1);
	int ub = max(m_nb2);
	m_rho.allocate(lb,ub);
	m_rho.initialize();
	dvar_vector tmp(lb-1,ub);
	tmp = 1.0;
	dvariable phi2 = -1. + (1. + sfabs(phi)) * psi;
	tmp(lb) = phi /(1.0-phi2);
	for( j = lb+1; j <= ub; j++ )
	{
		tmp(j) = phi*tmp(j-1) + phi2*tmp(j-2);
	}
	for( j = lb; j <= ub; j++ )
	{
		m_rho(j) = tmp(j);
	}
}


/**
 * Compute the vector of correlation matrixes for each year based on the
 * vector of m_rho values.  If there is no correlation then use a simple identity
 * matrix plus 1.  Else, C is the correlation matrix that is used to calculate the
 * covariance matrix m_V(in year i).
**/
void logistic_normal::compute_correlation_array()
{
	int i,j,k;
	

	if ( sum(first_difference(m_rho)) == 0 )
	{
		for( i = m_y1; i <= m_y2; i++ )
		{
			m_V(i) = 1 + identity_matrix(m_nb1(i),m_nb2(i)-1);
		}
		return; 
	}

	for( i = m_y1; i <= m_y2; i++ )
	{
		int l = m_nb1(i);
		int u = m_nb2(i);

		dvar_matrix C = identity_matrix(l,u);
		for( j = l; j <= u; j++ )
		{
			for( k = l; k <= u; k++ )
			{
				if(j != k) C(j,k) = m_rho(l-1+abs(j-k));
			}
		}
		dmatrix I = identity_matrix(l,u-1);
		dmatrix tK(l,u,l,u-1);
		tK.sub(l,u-1) = I;
		tK(u)         = -1;
		dmatrix K = trans(tK);
		m_V(i) = K * C * tK;
	}

}

/**
 * This function aggregates adjacent cohorts if the proportion of each
 * cohort is less than minp. This routine also allocates m_Op and m_Ep.
 *  
**/
void logistic_normal::aggregate_and_compress_arrays()
{
	// Determine the max column index for each year in the array.
	int i,j;
	m_nb1.allocate(m_y1,m_y2);
	m_nb2.allocate(m_y1,m_y2);
	for( i = m_y1; i <= m_y2; i++ )
	{
		m_nb1(i) = m_O(i).indexmin();
		int n = m_nb1(i)-1;
		double sumO = sum(m_O(i));
		for( j = m_nb1(i); j <= m_O(i).indexmax(); j++ )
		{
			double p = m_O(i,j) / sumO;
			if( p > minp ) n++;
		}
		m_nb2(i) = n;
	}

	// Now allocate arrays
	m_Op.allocate(m_y1,m_y2,m_nb1,m_nb2);
	m_Ep.allocate(m_y1,m_y2,m_nb1,m_nb2);
	m_Op.initialize();
	m_Ep.initialize();
	m_nAgeIndex.allocate(m_y1,m_y2,m_nb1,m_nb2);

	// Now aggregate observed and expected proprtions.
	for( i = m_y1; i <= m_y2; i++ )
	{
		dvector     oo = m_O(i) / sum(m_O(i));
		dvar_vector ee = m_E(i) / sum(m_E(i));
		int k = m_nb1(i);
		for( j = m_nb1(i); j <= m_O(i).indexmax(); j++ )
		{
			if( oo(j) <= minp )
			{
				m_Op(i)(k) += oo(j);
				m_Ep(i)(k) += ee(j);
			}
			else
			{
				m_Op(i)(k) += oo(j);
				m_Ep(i)(k) += ee(j);
				if( k <= m_nb2(i) ) m_nAgeIndex(i,k) = k;
				if( k <  m_nb2(i) ) k++;
			}
		}
	}
}

dvector logistic_normal::compute_relative_weights(const dmatrix &O)
{
	dvector Wy(O.rowmin(),O.rowmax());
	Wy = rowsum(O);
	double mu = mean(Wy);
	Wy = sqrt(mu / Wy);
	return(Wy);
}


/**
 * Compute the standardized residuals.
 * res{_by} = [log(O_{by}/\tilde(O)) - log(E_{by}/\tilde(E))] / (W_y G_y^{0.5})
 * Note that this matrix is of the same dimensions as the input matrix m_O.
 * Pseudocode:
 * 1). calculate \tilde(O) & \tilde(E)
**/
void logistic_normal::std_residuals()
{
	m_std_residual.allocate(m_O);
	m_std_residual.initialize();

	int i,j,k;
	for( i = m_y1; i <= m_y2; i++ )
	{
		int l = m_nb1(i);
		int u = m_nb2(i);
		
		// geometric means of each vector.
		double tO = geomean<double>(m_Op(i));
		dvector t1 = m_Op(i)/tO;
		double tE = geomean<double>(value(m_Ep(i)));
		dvector t2 = value(m_Ep(i))/tE;
		
		dmatrix  I = identity_matrix(l,u-1);
		
		dmatrix tF(l,u,l,u-1);
		tF.sub(l,u-1) = I;
		tF(u)         = 1;
		
		dmatrix  Hinv    = inv(I + 1);
		dmatrix FHinv    = tF * Hinv;
		dmatrix     G    = FHinv * value(m_V(i)) * trans(FHinv);
		dvector sd       = sqrt(diagonal(G));
		for( j = m_nb1(i); j <= m_nb2(i); j++ )
		{
			k = m_nAgeIndex(i)(j);
			m_std_residual(i)(k) = (t1(j)-t2(j)) / (sd(j)*m_Wy(i));
		}
	}}


































// dvariable nll_logistic_normal(const dmatrix &O, const dvar_matrix &E, 
//                               const double &minp, const double &eps,
//                               double &age_tau2)
// {
// 	int i,y1,y2;
// 	RETURN_ARRAYS_INCREMENT();
// 	y1 = O.rowmin();
// 	y2 = O.rowmax();

// 	// 1) (aggregate || add constant) && compress tails
// 	dmatrix n_Age  = get_tail_compressed_index(O,minp);
// 	dmatrix Op     = tail_compress(O,n_Age);
// 	dvar_matrix Ep = tail_compress(E,n_Age);
	
	
// 	if( eps )
// 	{
// 		// cout<<"adding constant"<<endl;
// 		add_constant_normalize(Op,eps);
// 		add_constant_normalize(Ep,eps);
// 	}
// 	else
// 	{
// 		// cout<<"aggregating cohorts"<<endl;
// 		aggregate(Op,Ep,minp);
// 	}

// 	// 2) Compute relative weights for each year W_y
// 	dvector Wy = compute_relative_weights(O);
	
// 	// 3) Compute covariance matrix V_y = K C K'
// 	dvar3_array Vy = compute_correlation_matrix(n_Age);
	
// 	// 4) Compute residual differences (w_y)
// 	dvar_matrix wwy = compute_residual_difference(Op,Ep);

// 	// 5) Compute weighted sum of squares (wSS)
// 	dvariable ssw = compute_weighted_sumofsquares(Wy,wwy,Vy);

// 	// 6) Compute MLE of variance
// 	double bm1 = size_count(Op) - (y2-y1+1.0);
// 	dvariable sigma2 = ssw / bm1;
// 	dvariable sigma  = sqrt(sigma2);
// 	age_tau2         = value(sigma2);
	
// 	// 7) Compute nll_logistic_normal
// 	dvariable nll;
// 	nll  = 0.5 * log(2.0 * PI) * bm1;
// 	nll += sum( log(Op) );
// 	nll += log(sigma) * bm1;
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		nll += 0.5 * log(det(Vy(i)));
// 		nll += (size_count(Op(i))-1) * log(Wy(i));
// 	}
// 	nll += 0.5 / sigma2 * ssw;

// 	RETURN_ARRAYS_DECREMENT();
// 	return nll;
// }

// /**
//  * Implementing the LN2 version.
// **/
// dvariable nll_logistic_normal(const dmatrix &O, const dvar_matrix &E, 
//                               const double &minp, const double &eps,
//                               double &age_tau2, const dvariable &phi)
// {
// 	int i,y1,y2;
// 	RETURN_ARRAYS_INCREMENT();
// 	y1 = O.rowmin();
// 	y2 = O.rowmax();

// 	// 1) (aggregate || add constant) && compress tails
// 	dmatrix n_Age  = get_tail_compressed_index(O,minp);
// 	dmatrix Op     = tail_compress(O,n_Age);
// 	dvar_matrix Ep = tail_compress(E,n_Age);
	
	
// 	if( eps )
// 	{
// 		// cout<<"adding constant"<<endl;
// 		add_constant_normalize(Op,eps);
// 		add_constant_normalize(Ep,eps);
// 	}
// 	else
// 	{
// 		// cout<<"aggregating cohorts"<<endl;
// 		aggregate(Op,Ep,minp);
// 	}

// 	// 2) Compute relative weights for each year W_y
// 	dvector Wy = compute_relative_weights(O);
	
// 	// 3) Compute covariance matrix V_y = K C K'
// 	dvar3_array Vy = compute_correlation_matrix(n_Age,phi);
	
// 	// 4) Compute residual differences (w_y)
// 	dvar_matrix wwy = compute_residual_difference(Op,Ep);

// 	// 5) Compute weighted sum of squares (wSS)
// 	dvariable ssw = compute_weighted_sumofsquares(Wy,wwy,Vy);

// 	// 6) Compute MLE of variance
// 	double bm1 = size_count(Op) - (y2-y1+1.0);
// 	dvariable sigma2 = ssw / bm1;
// 	dvariable sigma  = sqrt(sigma2);
// 	age_tau2         = value(sigma2);
	
// 	// 7) Compute nll_logistic_normal
// 	dvariable nll;
// 	nll  = 0.5 * log(2.0 * PI) * bm1;
// 	nll += sum( log(Op) );
// 	nll += log(sigma) * bm1;
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		nll += 0.5 * log(det(Vy(i)));
// 		nll += (size_count(Op(i))-1) * log(Wy(i));
// 	}
// 	nll += 0.5 / sigma2 * ssw;

// 	RETURN_ARRAYS_DECREMENT();
// 	return nll;
// }

// /**
//  * Implementing the LN2 version.
// **/
// dvariable nll_logistic_normal(const dmatrix &O, const dvar_matrix &E, 
//                               const double &minp, const double &eps,
//                               const dvariable &theta, const dvariable &phi)
// {
// 	int i,y1,y2;
// 	RETURN_ARRAYS_INCREMENT();
// 	y1 = O.rowmin();
// 	y2 = O.rowmax();

// 	// 1) (aggregate || add constant) && compress tails
// 	dmatrix n_Age  = get_tail_compressed_index(O,minp);
// 	dmatrix Op     = tail_compress(O,n_Age);
// 	dvar_matrix Ep = tail_compress(E,n_Age);
	
	
// 	if( eps )
// 	{
// 		// cout<<"adding constant"<<endl;
// 		add_constant_normalize(Op,eps);
// 		add_constant_normalize(Ep,eps);
// 	}
// 	else
// 	{
// 		// cout<<"aggregating cohorts"<<endl;
// 		aggregate(Op,Ep,minp);
// 	}

// 	// 2) Compute relative weights for each year W_y
// 	dvector Wy = compute_relative_weights(O);
	
// 	// 3) Compute covariance matrix V_y = K C K'
// 	dvar3_array Vy = compute_correlation_matrix(n_Age,phi);
	
// 	// 4) Compute residual differences (w_y)
// 	dvar_matrix wwy = compute_residual_difference(Op,Ep);

// 	// 5) Compute weighted sum of squares (wSS)
// 	dvariable ssw = compute_weighted_sumofsquares(Wy,wwy,Vy);

// 	// 6) Compute variance
// 	double bm1 = size_count(Op) - (y2-y1+1.0);
// 	dvariable sigma = sqrt(theta) / (1.0-phi);
// 	dvariable sigma2 = square(sigma);
	
	
// 	// 7) Compute nll_logistic_normal
// 	dvariable nll;
// 	nll  = 0.5 * log(2.0 * PI) * bm1;
// 	nll += sum( log(Op) );
// 	nll += log(sigma) * bm1;
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		nll += 0.5 * log(det(Vy(i)));
// 		nll += (size_count(Op(i))-1) * log(Wy(i));
// 	}
// 	nll += 0.5 / sigma2 * ssw;

// 	RETURN_ARRAYS_DECREMENT();
// 	return nll;
// }







// /**
//  * What is done here is that proportions less than minp are pooled with the adjacent
//  * younger cohort, and split evenly to each age, and the sample size is halved.
//  * This approximates the likelihood of ages 4 and 5 (eg.) as ages 4-5.
// **/
// void aggregate(dmatrix Op, dvar_matrix Ep, const double &minp)
// {
// 	int i,y1,y2;
// 	int j,b1,b2;

// 	y1 = Op.rowmin();
// 	y2 = Op.rowmax();

// 	for( i = y1; i <= y2; i++ )
// 	{
// 		b1 = Op(i).indexmin();
// 		b2 = Op(i).indexmax();
// 		for( j = b1+1; j <= b2; j++ )
// 		{
// 			if( Op(i,j) < minp )
// 			{
// 				double    tmpO = 0.5*(Op(i,j) + Op(i,j-1));
// 				dvariable tmpE = 0.5*(Ep(i,j) + Ep(i,j-1));
// 				Op(i)(j-1,j) = tmpO;
// 				Ep(i)(j-1,j) = tmpE;
// 			}
// 		}
// 	}
// }

// dvariable compute_weighted_sumofsquares(const dvector &Wy, 
//                                         const dvar_matrix &wwy,
//                                         const dvar3_array &V)
// {
// 	int i,y1,y2;
// 	RETURN_ARRAYS_INCREMENT();
// 	y1 = wwy.rowmin();
// 	y2 = wwy.rowmax();

// 	dvariable ssw = 0;
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		dvar_matrix Vinv = inv(V(i));
// 		ssw += (wwy(i) * Vinv * wwy(i))/ (Wy(i)*Wy(i));
// 	}

// 	RETURN_ARRAYS_DECREMENT();
// 	return ssw;
// }

// dvar_matrix compute_residual_difference(const dmatrix &O, const dvar_matrix &E)
// {
// 	int i,y1,y2;
	
// 	RETURN_ARRAYS_INCREMENT();

// 	y1 = O.rowmin();
// 	y2 = O.rowmax();
// 	ivector b1(y1,y2);
// 	ivector b2(y1,y2);
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		b1(i) = O(i).indexmin();
// 		b2(i) = O(i).indexmax();
// 	}
	
// 	dvar_matrix ww(y1,y2,b1,b2-1);
// 	ww.initialize();
	
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		int l = b1(i);
// 		int u = b2(i);
		
// 		dvector     t1 = O(i);
// 		dvar_vector t2 = E(i);
// 		ww(i) = (log(t1(l,u-1)) - log(t1(u)))
// 		      - (log(t2(l,u-1)) - log(t2(u)));
// 	}

// 	RETURN_ARRAYS_DECREMENT();
// 	return(ww);
// }

// /**
//  * No Autocorrelation case
// **/
// dvar3_array compute_correlation_matrix(const dmatrix &n_Age)
// {
// 	int i,y1,y2;
// 	RETURN_ARRAYS_INCREMENT();
// 	y1 = n_Age.rowmin();
// 	y2 = n_Age.rowmax();
// 	ivector b1(y1,y2);
// 	ivector b2(y1,y2);
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		b1(i) = min(n_Age(i));
// 		b2(i) = max(n_Age(i));
// 	}

// 	dvar3_array V;
// 	V.allocate(y1,y2,b1,b2-1,b1,b2-1);
// 	V.initialize();
	
// 	// Now compute the correlation matrix C, and set V = K C K'
// 	// Which is just the 1+identity matrix (i.e. no covariance structure)
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		V(i) = 1 + identity_matrix(b1(i),b2(i)-1);
// 	}
// 	RETURN_ARRAYS_DECREMENT();
// 	return (V);
// }
// /**
//  * AR1 case
// **/
// dvar3_array compute_correlation_matrix(const dmatrix &n_Age,const dvariable &phi)
// {
// 	int i,y1,y2;
// 	int j,k;
// 	RETURN_ARRAYS_INCREMENT();
// 	y1 = n_Age.rowmin();
// 	y2 = n_Age.rowmax();
// 	ivector b1(y1,y2);
// 	ivector b2(y1,y2);
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		b1(i) = min(n_Age(i));
// 		b2(i) = max(n_Age(i));
// 	}

// 	dvar3_array V;
// 	V.allocate(y1,y2,b1,b2-1,b1,b2-1);
// 	V.initialize();
	
// 	// Now compute the correlation matrix C, and set V = K C K'
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		int l = b1(i);
// 		int u = b2(i);
// 		dvar_vector rho(l,u);

// 		for( j = l, k = 1; j <= u; j++, k++ )
// 		{
// 			rho(j) = pow(phi,k);
// 			// rho(j) = phi / k;
// 			// rho(j) = mfexp( k * log(phi) );
// 		}
		
// 		dvar_matrix C = identity_matrix(l,u);
// 		for( j = l; j <= u; j++ )
// 		{
// 			for( k = l; k <= u; k++ )
// 			{
// 				if(j != k) C(j,k) = rho(l-1+abs(j-k));
// 			}
// 		}
		
// 		dmatrix I = identity_matrix(l,u-1);
// 		dmatrix tK(l,u,l,u-1);
// 		tK.sub(l,u-1) = I;
// 		tK(u)         = -1.0;
// 		dmatrix K = trans(tK);
		
// 		V(i) = K * C * tK;
		
// 	}
// 	RETURN_ARRAYS_DECREMENT();
// 	return (V);
// }







// dmatrix get_tail_compressed_index(const dmatrix &O, const double &minp)
// {
// 	int i,y1,y2;
// 	int j,b1,b2;

// 	y1 = O.rowmin();
// 	y2 = O.rowmax();
// 	b1 = O.colmin();
// 	b2 = O.colmax();
// 	dmatrix p(y1,y2,b1,b2);
// 	ivector n_B1(y1,y2);
// 	ivector n_B2(y1,y2);
// 	n_B1.initialize();
// 	n_B2.initialize();

// 	for( i = y1; i <= y2; i++ )
// 	{
// 		bool blt = true;    //left tail
// 		n_B1(i) = b1;
// 		n_B2(i) = b1-1;
// 		p(i) = O(i)/sum(O(i));
// 		for( j = b1; j <= b2; j++ )
// 		{
// 			if( p(i,j) <= minp && blt ) n_B1(i) ++; else blt=false;
// 			if( p(i,j) >  minp ) n_B2(i) ++;
// 		}	
// 	}
	
// 	dmatrix n_Age(y1,y2,n_B1,n_B2);
// 	n_Age.initialize();
	
// 	for( i = y1; i <= y2; i++ )
// 	{
// 		int k = n_B1(i);
// 		for( j = b1; j <= b2; j++ )
// 		{
// 			if( p(i,j) > minp )
// 			{
// 				if( k <= n_B2(i) ) n_Age(i,k) = k;
// 				if( k <  n_B2(i) ) k++;
// 			}
// 		}
// 	}
	
// 	return (n_Age);
// }