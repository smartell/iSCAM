#ifndef __MULTINOMIAL_H
#define __MULTINOMIAL_H

#define TINY     1.e-08

dvariable mult_likelihood(const dmatrix &o, const dvar_matrix &p, dvar_matrix &nu, const dvariable &log_vn);


#endif