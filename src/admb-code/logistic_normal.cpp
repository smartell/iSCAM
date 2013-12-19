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


**/


#include <admodel.h>
#include "logistic_normal.h"

logistic_normal::~logistic_normal()
{}

logistic_normal::logistic_normal()
{}

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


dvariable logistic_normal::negative_loglikelihood(const dvariable& tau2)
{

	aggregate_arrays();

	m_nll = 0;

	int i;
	m_sig2=tau2;

	for( i = m_y1; i <= m_y2; i++ )
	{	
		int    nB     = m_nNminp(i);
		double t1     = 0.5 * (nB - 1.0) * log(2.0*PI);
		double t3     = sum( log(m_Oa(i)) );
		dvariable det = pow(nB * m_sig2, 2.*(nB-1.));
		dvariable t5  = 0.5 * log(det);
		double t7     = (nB - 1.0) * log(m_dWy(i));
		dvariable t9  = 0.5 * sum( square(m_w(i)) / (m_sig2 * square(m_dWy(i))) );
		
		m_nll        += t1 + t3 + t5 + t7 + t9;
	}
		
	return(m_nll);
}

dvariable logistic_normal::negative_loglikelihood()
{
	// negative loglikelihood evaluated at the MLE of the variance tau2.
	m_nll = 0;
	cout<<"Integrate over the variance"<<endl;
	int i;
	correlation_matrix();
	m_sig2 = 0;
	dvariable sws = 0;
	for( i = m_y1; i <= m_y2; i++ )
	{
		cout<<m_w(i)<<endl<<endl;
		cout<<m_V(i)<<endl<<endl;
		cout<<m_w(i)*m_V(i)<<endl<<endl;
		cout<<(m_w(i)*m_V(i))*m_w(i)<<endl;
		sws += (m_w(i) * m_V(i)) * m_w(i);		///> equation A11
		exit(1);
	}

	double t1 = 0.5 * log(2.0*PI) * sum( m_nNminp-1 );
	double t2 = sum( log(m_Oa) );

	return(m_nll);
}

dvar_matrix logistic_normal::standardized_residuals()
{
	/*
		Standardized residuals
		nu = [log(O/~O) - log(E/~E)]/(Wy m_sig2^0.5)
		where O & E are the observed and expected proportion vectors
		~O and ~E is the geometric means of each proporition vectors
		Wy is the relative weight each year.
		
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
			m_residual(i)(idx) = w(j);
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

				if(k <=m_nNminp(i)) m_nAgeIndex(i,k) = j;
				if(k < m_nNminp(i)) k++;
			}
		}

		// compressed matrix of log-residuals
		nB = m_nNminp(i);
		m_w(i) = log(m_Oa(i)(1,nB-1)/m_Oa(i,nB))
		        -log(m_Ea(i)(1,nB-1)/m_Ea(i,nB));
		
	}
}


void logistic_normal::correlation_matrix()
{
	// Calculate the covariance matrix for each year, based on ragged arrays m_w;

	int i,j,k;
	m_V.allocate(m_y1,m_y2,1,m_nNminp,1,m_nNminp);
	m_V.initialize();

	for( i = m_y1; i <= m_y2; i++ )
	{
		 dmatrix K = identity_matrix(1,m_nNminp(i));
		 dvar_matrix C(1,m_nNminp(i),1,m_nNminp(i));
		 for( j = 1; j <= m_nNminp(i); j++ )
		 {
		 	 for( k = 1; k <= m_nNminp(i); k++ )
		 	 {
		 	 	C(j,k) = pow(m_rho,abs(double(j)-double(k)));
		 	 }
		 }
		 m_V(i) = K * C * trans(K);
	}	

}













