#include <admodel.h>
#include "multinomial.h"
dvariable mult_likelihood(const dmatrix &o, const dvar_matrix &p, dvar_matrix &nu, 
                          const dvariable &log_vn)
{

	// kludge to ensure observed and predicted matrixes are the same size
	if(o.colsize()!=p.colsize() || o.rowsize()!=p.rowsize())
	{
		cerr<<"Error in multivariate_t_likelihood, observed and predicted matrixes"
		" are not the same size"<<endl;
		ad_exit(1);
	}
	
	dvariable vn = mfexp(log_vn);
	dvariable ff = 0.0;
	int r1 = o.rowmin();
	int r2 = o.rowmax();
	int c1 = o.colmin();
	int c2 = o.colmax();

	for(int i = r1; i <= r2; i++ )
	{
		dvar_vector sobs = vn * o(i)/sum(o(i));  //scale observed numbers by effective sample size.
		ff -= gammln(vn);
		for(int j = c1; j <= c2; j++ )
		{
			if( value(sobs(j)) > 0.0 )
				ff += gammln(sobs(j));
		}
		ff -= sobs * log(TINY + p(i));
		
		dvar_vector o1=o(i)/sum(o(i));
		dvar_vector p1=p(i)/sum(p(i));
		nu(i) = elem_div(o1-p1,sqrt(elem_prod(p1,1.-p1)/vn));


	}
	// exit(1);
	return ff;
}