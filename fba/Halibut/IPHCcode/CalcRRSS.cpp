// CalcRRSS -- calculate a robust residual sum of squares using
//             a variance scaler.
//
// Calculate a scaled sum of squares from corresponding matrices
// (or vectors) of predictions, observations, and standard
// deviations, leaving out cells where any item is NA (missing).
// The externally extimated sampling standard deviation of each 
// observation SDX is scaled up by an externally estimated factor
// Tau that accounts for overdispersion.
// The scaled deviation SDevX = (PredX - ObsX)/(Tau*SDX) may be
// robustified by an asymptotic function RSD = f(SD, ZBend, ZMax)
// which is equal to SD until ZBend and then bends and levels out
// with ZXax as an asymptote. With ZBend = 2.5 and ZMax =3 this
// function closely matches Fournier's robustification device.
// Robustification is performed if the flag Robust is set. With
// the flag unset and Tau=1, this routine calculates the same sum
// of squares as the old CalcRSS. The expectation is that a conventional
// sum of squares will be calculated in the first phase of minimization
// and a robust sum of squares in the second and/or later phases.
//
// Set shape of robustification function (ZMax > ZBend).
# define ZBend 2.5
# define ZMax  3.0
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
dvar_vector CalcRRSS(_CONST dvar_matrix& PredX, _CONST dvar_matrix& ObsX,
                    _CONST dvar_matrix& SDX, double Tau, int Robust, 
                    int MinJ, int MaxI, double NA)
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
double ZDiff = ZMax - ZBend;
for (int i = ibot; i<=MaxI; i++)
   for (int j = MinJ; j<=jtop; j++)
      {
      // Skip if any item is NA.
      if (  PredX(i,j) == NA ||
            ObsX(i,j)  == NA ||
            SDX(i,j)   == NA  ) continue;
      // Bail if SD is zero.
      if (Tau*SDX(i,j) == 0.0)
         {
         cerr << "error in CalcRSS: zero SD" << endl;
         exit(1);
         }
      // Add to accumulators.
      n++;
      sdev = fabs(PredX(i,j) - ObsX(i,j))/(Tau*SDX(i,j));
      if (Robust != 0 && sdev > ZBend) 
         sdev = ZBend + ZDiff * (1 - exp(-(sdev-ZBend)/ZDiff));
      RSS += sdev*sdev;
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
dvar_vector CalcRRSS(_CONST dvar_vector& PredX, _CONST dvar_vector& ObsX,
                    _CONST dvar_vector& SDX, double Tau, int Robust, 
                    int MinI, int MaxI, double NA)
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
double ZDiff = ZMax - ZBend;
for (int i = MinI; i<=MaxI; i++)
      {
      // Skip if any item is NA.
      if (  PredX(i) == NA ||
            ObsX(i)  == NA ||
            SDX(i)   == NA  ) continue;
      // Bail if SD is zero.
      if (Tau*SDX(i) == 0.0)
         {
         cerr << "error in CalcRSS: zero SD" << endl;
         exit(1);
         }
      // Add to accumulators.
      n++;
      sdev = fabs(PredX(i) - ObsX(i))/(Tau*SDX(i));
      if (Robust != 0 && sdev > ZBend) 
         sdev = ZBend + ZDiff * (1 - exp(-(sdev-ZBend)/ZDiff));
      RSS += sdev*sdev;
      } // End loop on i.
// Load return vector.
dvar_vector RSS_vec(1,3);
RSS_vec(1) = RSS;
RSS_vec(2) = n;
RSS_vec(3) = RSS/n;
RETURN_ARRAYS_DECREMENT();
return(RSS_vec);
}
