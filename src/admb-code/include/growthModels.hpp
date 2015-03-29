#include<admodel.h>

#ifndef GROWTHMODELS_HPP
#define GROWTHMODELS_HPP

template <typename T, typename T1>
inline
T1 vonBertalanffy(const T1 age, const T lmin, const T lmax, const T rho)
{
	T1 len;
	for(int i = age.indexmin(); i <= age.indexmax; i++ )
	{
		len(i) = lmin + (lmax - lmin) 
				* (1.0-pow(rho,i-1.0)) / (1.0-pow(rho,age.indexmax()-1.0));
	}
	return(len);
}

#endif
