// Predict the observed surface and/or burn age distn matrix ObsMat
// corresponding to a true age distn matrix TrueMat. Either or both
// kinds of readings may have to be predicted during the range of 
// years in TrueMat.
//
# include <fvar.hpp>
//
void SmearAges(dvar_matrix& TrueMat, dvar_matrix& ObsMat,
               dmatrix& TrueToSurf, dmatrix& TrueToBurn,
               double NA)
//
// In this routine, "b" is true age and is assumed to run from 1 to
// TrueMat.colmax(), e.g. 30 in 2003. Observed age is "a" which runs
// from amin = ObsMat.colmin(), e.g. 6 in 2003, through 
// amax = ObsMat.colmax(), e.g. 25 in 2003. Cells in ObsMat beyond
// the plus age are set to NA.
{
RETURN_ARRAYS_INCREMENT();
//
int FirstYear = TrueMat.rowmin();
int LastYear = TrueMat.rowmax();
//
int bmax = TrueMat.colmax();
int bmin = TrueMat.colmin();
if (bmin != 1)
   {cerr << "SmearAges: TrueMat must start at age 1" << endl;
    exit(1);}
//
int bmin_surf = TrueToSurf.colmin();
int bmax_surf = TrueToSurf.colmax();
if (bmin_surf != bmin || bmax_surf != bmax)
   {cerr << "SmearAges: TrueToSurf does not match TrueMat" << endl;
    exit(1);}
//
int bmin_burn = TrueToBurn.colmin();
int bmax_burn = TrueToBurn.colmax();
if (bmin_burn != bmin || bmax_burn != bmax)
   {cerr << "SmearAges: TrueToBurn does not match TrueMat" << endl;
    exit(1);}
//
int FirstYearObs = ObsMat.rowmin();
int LastYearObs = ObsMat.rowmax();
if (FirstYearObs != FirstYear || LastYearObs != LastYear)
   {cerr << "SmearAges: different year ranges" << endl;
    exit(1);}
//
int amin = ObsMat.colmin();
int amax = ObsMat.colmax();
int aplus_surf = TrueToSurf.rowmax();
int aplus_burn = TrueToBurn.rowmax();
//
// Compute full matrices of predicted surface and burn observations.
dvar_matrix SurfMat(FirstYear,LastYear,1,aplus_surf);
dvar_matrix BurnMat(FirstYear,LastYear,1,aplus_burn);
SurfMat = trans(TrueToSurf * trans(TrueMat));
BurnMat = trans(TrueToBurn * trans(TrueMat));
//
// cout << "TrueMaT" << endl << TrueMat << endl;
// cout << "SurfMaT" << endl << SurfMat << endl;
//
// Fill in ObsMat year by year.
int y;
int a;
for (y=FirstYear; y<=LastYear; y++)
   {
   if (y <= 2001) // Surface.
      {
      for (a=amin; a<=aplus_surf; a++) ObsMat(y,a) = SurfMat(y,a);
      if (amax > aplus_surf)
      for (a=aplus_surf+1; a<=amax; a++) ObsMat(y,a) = NA;
      //
      } else      // Burn.
      {
      for (a=amin; a<=aplus_burn; a++) ObsMat(y,a) = BurnMat(y,a);
      if (amax > aplus_burn)
      for (a=aplus_burn+1; a<=amax; a++) ObsMat(y,a) = NA;
      }
      // Done with this year.
   } // End loop on year.
// cout << "TrueMaT" << endl << TrueMat << endl;
// cout << "SurfMaT" << endl << SurfMat << endl;
// cout << "BurnMaT" << endl << BurnMat << endl;
// cout << "ObsMat" << endl << ObsMat << endl;
    //
RETURN_ARRAYS_DECREMENT();
}
