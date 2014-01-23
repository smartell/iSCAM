#include <admodel.h>
#include "LogisticNormal.h"

/*
	Implementation of the logistic normal negative loglikelihood.
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

/* FIXME Need to rethink this algorithm  TEMPLATE? */
void aggregate(dmatrix M, const double &minp)
{
	int i,y1,y2;
	int j,b1,b2;

	y1 = M.rowmin();
	y2 = M.rowmax();

	dmatrix tmpM = M;
	tmpM.initialize();

	for( i = y1; i <= y2; i++ )
	{
		b1 = M.indexmin(); 
		b2 = M.indexmax();
		int k=b1;
		for( j = b1; j <= b2; j++ )
		{
			if( M(i,j) <= minp )
			{
				tmpM(i)(k) += M(i,j);
			} 
			else
			{
				tmpM(i)(k) += M(i,j);
				if( k < b2 ) k++;
			}
		}
	}
	M = tmpM;
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
		// TODO Need to re-think this one.
		cout<<"aggregating cohorts"<<endl;
		//aggregate(O,minp);
	}
	
	// 2) Compute relative weights for each year W_y
	dvector Wy = compute_relative_weights(O);
	
	// 3) Compute covariance matrix V_y = K C K'
	d3_array Vy = compute_correlation_matrix(n_Age);
	
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

dvariable compute_weighted_sumofsquares(const dvector &Wy, 
                                        const dvar_matrix &wwy,
                                        const d3_array &V)
{
	int i,y1,y2;
	RETURN_ARRAYS_INCREMENT();
	y1 = wwy.rowmin();
	y2 = wwy.rowmax();

	dvariable ssw = 0;
	for( i = y1; i <= y2; i++ )
	{
		dmatrix Vinv = inv(V(i));
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
d3_array compute_correlation_matrix(const dmatrix &n_Age)
{
	int i,y1,y2;
	y1 = n_Age.rowmin();
	y2 = n_Age.rowmax();
	ivector b1(y1,y2);
	ivector b2(y1,y2);
	for( i = y1; i <= y2; i++ )
	{
		b1(i) = min(n_Age(i));
		b2(i) = max(n_Age(i));
	}

	d3_array V;
	V.allocate(y1,y2,b1,b2-1,b1,b2-1);
	V.initialize();
	
	// Now compute the correlation matrix C, and set V = K C K'
	// Which is just the 1+identity matrix (i.e. no covariance structure)
	for( i = y1; i <= y2; i++ )
	{
		V(i) = 1 + identity_matrix(b1(i),b2(i)-1);
	}
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

dmatrix tail_compress(const dmatrix &O,const dmatrix &n_Age)
{
	int i,y1,y2;
	int b1,b2;

	y1 = O.rowmin();
	y2 = O.rowmax();
	b1 = O.colmin();
	b2 = O.colmax();

	dmatrix P;
	P.allocate(n_Age);
	P.initialize();
	
	for( i = y1; i <= y2; i++ )
	{
		ivector ix = n_Age(i);
		dvector px = O(i)/sum(O(i));
		P(i)(min(ix),max(ix)) = px(ix);
		if( min(ix) > b1 )
		{
			P(i)(min(ix)) = sum(px(b1,min(ix)));
		}
		if( max(ix) < b2 )
		{
			P(i)(max(ix)) = sum(px(max(ix),b2));
		}
	}
	return(P);
}

dvar_matrix tail_compress(const dvar_matrix &O,const dmatrix &n_Age)
{
	int i,y1,y2;
	int b1,b2;

	y1 = O.rowmin();
	y2 = O.rowmax();
	b1 = O.colmin();
	b2 = O.colmax();

	dvar_matrix P;
	P.allocate(n_Age);
	P.initialize();
	
	for( i = y1; i <= y2; i++ )
	{
		ivector ix = n_Age(i);
		dvar_vector px = O(i)/sum(O(i));
		P(i)(min(ix),max(ix)) = px(ix);
		if( min(ix) > b1 )
		{
			P(i)(min(ix)) = sum(px(b1,min(ix)));
		}
		if( max(ix) < b2 )
		{
			P(i)(max(ix)) = sum(px(max(ix),b2));
		}
	}
	return(P);
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