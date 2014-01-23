#include <admodel.h>
#ifndef __LOGISTIC_NORMAL_H
#define __LOGISTIC_NORMAL_H

dmatrix get_tail_compressed_index(const dmatrix &O, const double &minp);
dmatrix tail_compress(const dmatrix &O,const dmatrix &n_Age);
dvar_matrix tail_compress(const dvar_matrix &O,const dmatrix &n_Age);
dvector compute_relative_weights(const dmatrix &O);
d3_array compute_correlation_matrix(const dmatrix &n_Age);
dvar_matrix compute_residual_difference(const dmatrix &O, const dvar_matrix &E);
dvariable compute_weighted_sumofsquares(const dvector &Wy, const dvar_matrix &wwy,const d3_array &V);
dvariable nll_logistic_normal(const dmatrix &O, const dvar_matrix &E, 
                              const double &minp, const double &eps,
                              double &age_tau2);

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

#endif