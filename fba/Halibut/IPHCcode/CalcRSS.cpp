// Calculate a scaled sum of squares from corresponding matrices
// (or vectors) of predictions, observations, and standard
// deviations, leaving out cells where any item is NA (missing).
// The absolute deviation SD = |(PredX - ObsX)/SDX| is raised to
// the power SDevPow so with SDevPow=2 the objective function is
// least squares, with SDevPow=1 it is least absolute deviations
// and so on.
// 
// Begin at column MinJ (element MinJ of a vector)--this to allow for 
// calculating RSS when the predictions are a subset of the ages in the 
// observation and SD matrices. The matrices do not have to have the 
// same starting column index, but they have to have the same ending 
// index, and both must include MinJ. 
//
// End at row MaxI (element MaxI of a vector)--this to allow for
// leaving a year or two out of the RSS calcs at the end.
//
// Return a dvar_vector containing RSS, n (points contributing),
// and mean score ms = RSS/n.
//
# include <fvar.hpp>
//
// First the dvar_matrix version.
dvar_vector CalcRSS(_CONST dvar_matrix& PredX, _CONST dvar_matrix& ObsX,
                    _CONST dvar_matrix& SDX, double SDevPow, int MinJ, 
                    int MaxI, double NA)
{
RETURN_ARRAYS_INCREMENT();
// Determine and check sizes.
int ibot = ObsX.rowmin();
int itop = ObsX.rowmax();
int jbot = ObsX.colmin();
int jtop = ObsX.colmax();
//
if (MinJ < jbot || MinJ > jtop)
   {
   cerr << "error in CalcRSS: start column out of range" << endl;
   exit(1);
   }
if (MaxI < ibot || MaxI > itop)
   {
   cerr << "error in CalcRSS: end row out of range" << endl;
   exit(1);
   }
if (  SDX.rowmin() != ibot ||
      SDX.rowmax() != itop ||
      SDX.colmin() != jbot ||
      SDX.colmax() != jtop  )
   {
   cerr << "error in CalcRSS: dimensions of SD matrix" << endl;
   exit(1);
   }
if (  PredX.rowmin() != ibot ||
      PredX.rowmax() != itop ||
      PredX.colmin() > MinJ  ||
      PredX.colmax() != jtop  )
   {
   cerr << "error in CalcRSS: dimensions of Pred matrix" << endl;
   exit(1);
   }
//
// Run through cells.
int n = 0;
dvariable RSS = 0.0;
dvariable sdev;
for (int i = ibot; i<=MaxI; i++)
   for (int j = MinJ; j<=jtop; j++)
      {
      // Skip if any item is NA.
      if (  PredX(i,j) == NA ||
            ObsX(i,j)  == NA ||
            SDX(i,j)   == NA  ) continue;
      // Bail if SD is zero.
      if (SDX(i,j) == 0.0)
         {
         cerr << "error in CalcRSS: zero SD" << endl;
         exit(1);
         }
      // Add to accumulators.
      sdev = fabs((PredX(i,j) - ObsX(i,j))/SDX(i,j));
      n++;
      RSS += pow(sdev, SDevPow);
      } // End loop on i and j.
// Load return vector.
dvar_vector RSS_vec(1,3);
RSS_vec(1) = RSS;
RSS_vec(2) = n;
RSS_vec(3) = RSS/n;
RETURN_ARRAYS_DECREMENT();
return(RSS_vec);
}
//
// Next the dvar_vector (arg) version.
dvar_vector CalcRSS(_CONST dvar_vector& PredX, _CONST dvar_vector& ObsX,
                    _CONST dvar_vector& SDX, double SDevPow, int MinI, 
                    int MaxI, double NA)
{
RETURN_ARRAYS_INCREMENT();
// Determine and check sizes.
int ibot = ObsX.indexmin();
int itop = ObsX.indexmax();
//
if (MinI < ibot || MinI > itop)
   {
   cerr << "error in CalcRSS: start point out of range" << endl;
   exit(1);
   }
if (MaxI < ibot || MaxI > itop)
   {
   cerr << "error in CalcRSS: end point out of range" << endl;
   exit(1);
   }
if (  SDX.indexmin() != ibot ||
      SDX.indexmax() != itop  )
   {
   cerr << "error in CalcRSS: indices of SD vector" << endl;
   exit(1);
   }
if (  PredX.indexmin() > MinI  ||
      PredX.indexmax() != itop  )
   {
   cerr << "error in CalcRSS: indices of Pred vector" << endl;
   exit(1);
   }
//
// Run through cells.
int n = 0;
dvariable RSS = 0.0;
dvariable sdev;
for (int i = MinI; i<=MaxI; i++)
      {
      // Skip if any item is NA.
      if (  PredX(i) == NA ||
            ObsX(i)  == NA ||
            SDX(i)   == NA  ) continue;
      // Bail if SD is zero.
      if (SDX(i) == 0.0)
         {
         cerr << "error in CalcRSS: zero SD" << endl;
         exit(1);
         }
      // Add to accumulators.
      sdev = fabs((PredX(i) - ObsX(i))/SDX(i));
      n++;
      RSS += pow(sdev, SDevPow);
      } // End loop on i.
// Load return vector.
dvar_vector RSS_vec(1,3);
RSS_vec(1) = RSS;
RSS_vec(2) = n;
RSS_vec(3) = RSS/n;
RETURN_ARRAYS_DECREMENT();
return(RSS_vec);
}
