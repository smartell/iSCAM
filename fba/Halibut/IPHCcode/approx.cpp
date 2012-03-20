// Implement the S-PLUS function approx(xvec, yvec, xout=xout, rule=2).
// This interpolates linearly between the values of xvec and yvec at
// the points xout, and returns yout. If any element of xout lies
// outside the range of xvec, it is set to the corresponding extreme
// value of yvec. xvec is assumed to be sorted in ascending order.
//
# include <fvar.hpp>
//
// First a dvector version.
//
dvector approx(_CONST dvector& xvec, _CONST dvector& yvec,
               _CONST dvector& xout)
{
// Construct return value container.
int xouta = xout.indexmin();
int xoutz = xout.indexmax();
dvector yout(xouta, xoutz);
// Check xvec and yvec for conformity.
int xveca = xvec.indexmin();
int xvecz = xvec.indexmax();
int xvecn = xvec.size();
int yveca = yvec.indexmin();
int yvecz = yvec.indexmax();
int yvecn = yvec.size();
if (xvecn != yvecn)
   {
   cerr << "error in approx: xvec and yvec differ in length" << endl;
   exit(1);
   }
// Check xvec for sortedness.
for (int ix = xveca+1; ix <= xvecz; ix++)
   {
   double xdiff = xvec(ix) - xvec(ix-1);
   if (xdiff < 0.0)
      {
      cerr << "error in approx: xvec decreases" << endl;
      exit(1);
      }
   }
// Interpolate point by point.
for (int iout = xouta; iout <= xoutz; iout++)
   {
   if (xout(iout) <= xvec(xveca))
      {
      yout(iout) = yvec(yveca);
      }
   else if (xout(iout) > xvec(xvecz))
      {
      yout(iout) = yvec(yvecz);
      }
   else for (int ioff = 1; ioff <= xvecn-1; ioff++)
      {
      int ix = xveca + ioff;
      int iy = yveca + ioff;
      if (xout(iout) <= xvec(ix))
         {
         yout(iout) = yvec(iy-1) + (yvec(iy) - yvec(iy-1)) *
                      (xout(iout) - xvec(ix-1))/(xvec(ix) - xvec(ix-1));
         break;
         } // Found slot and stored interpolated y.
      } // End loop on ioff.
   } // End loop on iout.
// Done.
return(yout);
}
//
// A dvar_vector version with the obligatory RETURN_ARRAYS stuff.
dvar_vector approx(_CONST dvar_vector& xvec, _CONST dvar_vector& yvec,
               _CONST dvar_vector& xout)
{
RETURN_ARRAYS_INCREMENT();
// Construct return value container.
int xouta = xout.indexmin();
int xoutz = xout.indexmax();
dvar_vector yout(xouta, xoutz);
// Check xvec and yvec for conformity.
int xveca = xvec.indexmin();
int xvecz = xvec.indexmax();
int xvecn = xvec.size();
int yveca = yvec.indexmin();
int yvecz = yvec.indexmax();
int yvecn = yvec.size();
if (xvecn != yvecn)
   {
   cerr << "error in approx: xvec and yvec differ in length" << endl;
   exit(1);
   }
// Check xvec for sortedness.
for (int ix = xveca+1; ix <= xvecz; ix++)
   {
   dvariable xdiff = xvec(ix) - xvec(ix-1);
   if (xdiff < 0.0)
      {
      cerr << "error in approx: xvec decreases" << endl;
      exit(1);
      }
   }
// Interpolate point by point.
for (int iout = xouta; iout <= xoutz; iout++)
   {
   if (xout(iout) <= xvec(xveca))
      {
      yout(iout) = yvec(yveca);
      }
   else if (xout(iout) > xvec(xvecz))
      {
      yout(iout) = yvec(yvecz);
      }
   else for (int ioff = 1; ioff <= xvecn-1; ioff++)
      {
      int ix = xveca + ioff;
      int iy = yveca + ioff;
      if (xout(iout) <= xvec(ix))
         {
         yout(iout) = yvec(iy-1) + (yvec(iy) - yvec(iy-1)) *
                      (xout(iout) - xvec(ix-1))/(xvec(ix) - xvec(ix-1));
         break;
         } // Found slot and stored interpolated y.
      } // End loop on ioff.
   } // End loop on iout.
RETURN_ARRAYS_DECREMENT();
// Done.
return(yout);
}
//
// A dvariable version with a SAFE flag to bypass list checks.
dvariable approx(_CONST dvar_vector& xvec, _CONST dvar_vector& yvec,
               _CONST dvariable& xout, int SAFE)
{
RETURN_ARRAYS_INCREMENT();
// Construct return value container.
dvariable yout;
// Get indexes of xvec and yvec.
int xveca = xvec.indexmin();
int xvecz = xvec.indexmax();
int xvecn = xvec.size();
int yveca = yvec.indexmin();
int yvecz = yvec.indexmax();
int yvecn = yvec.size();
// Maybe do some safety checks.
if (SAFE != 0)
{
// Check xvec and yvec for conformity.
if (xvecn != yvecn)
   {
   cerr << "error in approx: xvec and yvec differ in length" << endl;
   exit(1);
   }
// Check xvec for sortedness.
for (int ix = xveca+1; ix <= xvecz; ix++)
   {
   dvariable xdiff = xvec(ix) - xvec(ix-1);
   if (xdiff < 0.0)
      {
      cerr << "error in approx: xvec decreases" << endl;
      exit(1);
      }
   }
} // End SAFEty checks.
// Interpolate.
if (xout <= xvec(xveca))
{
yout = yvec(yveca);
}
else if (xout > xvec(xvecz))
{
yout = yvec(yvecz);
}
else for (int ioff = 1; ioff <= xvecn-1; ioff++)
{
int ix = xveca + ioff;
int iy = yveca + ioff;
if (xout <= xvec(ix))
{
yout = yvec(iy-1) + (yvec(iy) - yvec(iy-1)) *
(xout - xvec(ix-1))/(xvec(ix) - xvec(ix-1));
break;
} // Found slot and stored interpolated y.
} // End loop on ioff.
RETURN_ARRAYS_DECREMENT();
// Done.
return(yout);
}
//
// Another dvariable version with a SAFE flag to bypass list checks.
// Here xvec is a dvector.
dvariable approx(_CONST dvector& xvec, _CONST dvar_vector& yvec,
               _CONST dvariable& xout, int SAFE)
{
RETURN_ARRAYS_INCREMENT();
// Construct return value container.
dvariable yout;
// Get indexes of xvec and yvec.
int xveca = xvec.indexmin();
int xvecz = xvec.indexmax();
int xvecn = xvec.size();
int yveca = yvec.indexmin();
int yvecz = yvec.indexmax();
int yvecn = yvec.size();
// Maybe do some safety checks.
if (SAFE != 0)
{
// Check xvec and yvec for conformity.
if (xvecn != yvecn)
   {
   cerr << "error in approx: xvec and yvec differ in length" << endl;
   exit(1);
   }
// Check xvec for sortedness.
for (int ix = xveca+1; ix <= xvecz; ix++)
   {
   dvariable xdiff = xvec(ix) - xvec(ix-1);
   if (xdiff < 0.0)
      {
      cerr << "error in approx: xvec decreases" << endl;
      exit(1);
      }
   }
} // End SAFEty checks.
// Interpolate.
if (xout <= xvec(xveca))
{
yout = yvec(yveca);
}
else if (xout > xvec(xvecz))
{
yout = yvec(yvecz);
}
else for (int ioff = 1; ioff <= xvecn-1; ioff++)
{
int ix = xveca + ioff;
int iy = yveca + ioff;
if (xout <= xvec(ix))
{
yout = yvec(iy-1) + (yvec(iy) - yvec(iy-1)) *
(xout - xvec(ix-1))/(xvec(ix) - xvec(ix-1));
break;
} // Found slot and stored interpolated y.
} // End loop on ioff.
RETURN_ARRAYS_DECREMENT();
// Done.
return(yout);
}
