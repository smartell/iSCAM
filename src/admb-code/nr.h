// Header file for numerical recipies code.
#include <admodel.h>
#include <dfridr.cpp>

#ifndef _NR_H_
#define _NR_H_

typedef double DP;

namespace NR {
	DP dfridr(DP func(const DP), const DP x, const DP h, DP &err);
}
#endif /* _NR_H_ */
