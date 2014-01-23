#include <admodel.h>

/*
	Implementation of the logistic normal negative loglikelihood.
*/
dmatrix get_tail_compressed_index(const dmatrix &O, const double &minp);
dmatrix tail_compress(const dmatrix &O,const dmatrix &n_Age);
dvar_matrix tail_compress(const dvar_matrix &O,const dmatrix &n_Age);
dvector compute_relative_weights(const dmatrix &O);
d3_array compute_correlation_matrix(const dmatrix &n_Age);
dvariable nll_logistic_normal(const dmatrix &O, const dvar_matrix &E, 
                              const double &minp, const double &eps);

template <typename T>
void add_constant_normalize(T M, const double &eps)
{
	int i,y1,y2;
	y1 = M.rowmin();
	y2 = M.rowmax();
	for( i = y1; i <= y2; i++ )
	{
		M(i) = M(i) + eps;
		M(i) = M(i) / sum(M(i));
	}
}

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

/* FIXME Need to rethink this algorithm */
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
                              const double &minp, const double &eps)
{
	RETURN_ARRAYS_INCREMENT();
	dvariable nll;

	// 1) (aggregate || add constant) && compress tails
	dmatrix n_Age  = get_tail_compressed_index(O,minp);
	dmatrix Op     = tail_compress(O,n_Age);
	dvar_matrix Ep = tail_compress(E,n_Age);
	if( eps )
	{
		cout<<"adding constant"<<endl;
		add_constant_normalize(O,eps);
		add_constant_normalize(E,eps);
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
	cout<<" Whats up"<<endl;
	d3_array Vy = compute_correlation_matrix(n_Age);

	exit(1);
	RETURN_ARRAYS_DECREMENT();
	return nll;
}

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
	V.allocate(y1,y2,b1,b2,b1,b2);
	V.initialize();
	
	// Now compute the correlation matrix C, and set V = K C K'
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