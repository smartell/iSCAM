#include <admodel.h>
#include "multinomial.h"


dvariable multivariate_t_likelihood(const dmatrix &o, const dvar_matrix &p, 
                                    const dvariable &log_var, const dvariable &log_v,
                                    const dvariable &expon, dvar_matrix &nu)
{
	// kludge to ensure observed and predicted matrixes are the same size
	if(o.colsize()!=p.colsize() || o.rowsize()!=p.rowsize())
	{
		cerr<<"Error in multivariate_t_likelihood, observed and predicted matrixes"
		" are not the same size"<<endl;
		ad_exit(1);
	}
	
	dmatrix O = o;
	dvar_matrix P = p;
	O.colshift(1);
	P.colshift(1);


	int r1 = o.rowmin();
	int r2 = o.rowmax();
	int c1 = o.colmin();
	int c2 = o.colmax();
	dvariable ll=0.0;
	dvariable v   = exp(log_v);
	dvariable var = exp(log_var);
	dvar_matrix COVAR(c1,c2,c1,c2);
	dvar_matrix tmp_nu(r1,r2,c1,c2);
	COVAR.initialize();
	const double pp    = (c2-c1)+1;
	const double lppi2 = 0.5*pp*log(3.1415926535);

	for(int i = r1; i <= r2; i++ )
	{
		for(int j = c1; j <= c2; j++ )
		{
			COVAR(j,j) = 0.001 + pow( P(i,j)*(1.0-P(i,j)),exp(expon) );
		}
		COVAR *=var;

		dvar_vector diff = O(i)/sum(O(i)) - P(i);
		dvariable ln_det = 0.0;
		int sgn = 1;
		dvar_vector e = choleski_solve(COVAR,diff,ln_det,sgn);
		dvariable ln_det_choleski_inv = -ln_det;
		tmp_nu(i) = e/sqrt(v);
		
		ll += -gammln(0.5*(v+pp)) + gammln(0.5*v)
		      + lppi2 + (0.5*pp)*log(v) 
		      + 0.5*(v+pp) * log(1.0+e*e/v) - ln_det_choleski_inv;
	}
	
	tmp_nu.rowshift(nu.rowmin());
	tmp_nu.colshift(nu.colmin());
	nu = tmp_nu;
	return ll;
}

// FUNCTION  dvariable compositional_likelihood(void)
//  v=exp(log_v);
//  var=exp(log_var);
//  dvariable ll=0.0;
//  COVAR.initialize();
//  for (int ii=1;ii<=nyrs;ii++)
//  {

//    for (int i=1;i<=nage;i++)
//    {
//      COVAR(i,i)=.001+pcomp(i)*(1.0-pcomp(i));   // add the .001 so variance never goes to 0.
//    }
//    COVAR*=var;    // scale covariance
//    dvar_vector diff=obs_composition_at_age(ii)-pcomp(ii);
//    dvariable ldet=0.0;
//    int sgn=1;
//    dvar_vector e=choleski_solve(COVAR,diff,ldet,sgn);
//    dvariable ln_det_choleski_inv=-ldet;
//    double pp=nage;
//    const double lppi2=0.5*pp*log(3.1415926535);
//    ll += -gammln(0.5*(v+pp)) + gammln(0.5*v)
//         + lppi2 +(0.5*pp)*log(v)
//         +0.5*(v+pp)*log(1.0+e*e/v)
//         - ln_det_choleski_inv;
//  }
//  return ll;