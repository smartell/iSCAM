// Calculate a penalty on the jaggedness of a selectivity schedule.
// This is done by computing and penalizing second differences.
//
# include <fvar.hpp>
# include <ourlib.hpp>
//
  dvariable pen2diffs(_CONST dvar_vector& selvec, _CONST double& pensd)
{
RETURN_ARRAYS_INCREMENT();
//
int ibot = selvec.indexmin();
int itop = selvec.indexmax();
int i;
dvar_vector diff2vec(ibot,itop);
dvariable penalty = 0.0;
// Need at least three points.
if (itop - ibot > 2) for (i=ibot+2; i<=itop; i++)
   {
   diff2vec(i) = selvec(i) - 2.0*selvec(i-1) + selvec(i-2);
   penalty += sq( diff2vec(i) / pensd );
   }
RETURN_ARRAYS_DECREMENT();
return(penalty);
}
