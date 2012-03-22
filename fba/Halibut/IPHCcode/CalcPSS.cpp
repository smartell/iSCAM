// Calculate a penalty sum of squares from a matrix (or vector)
// of deviations and a SCALAR penalty SD. Return PSS, PSS_n, and
// PSS_rms, just like CalcRSS.
// Also like CalcRSS, allow option to begin at a specified column
// in a matrix (or element in a vector).
//
# include <fvar.hpp>
//
// First the dvar_matrix version.
dvar_vector CalcPSS(_CONST dvar_matrix& DevX, _CONST dvariable& DevSD,
                     int MinJ, double NA)
{
RETURN_ARRAYS_INCREMENT();
// Determine and check sizes.
int ibot = DevX.rowmin();
int itop = DevX.rowmax();
int jbot = DevX.colmin();
int jtop = DevX.colmax();
//
if (MinJ < jbot || MinJ > jtop)
   {
   cerr << "error in CalcPSS: start column out of range" << endl;
   exit(1);
   }
//
// Check DevSD.
if (DevSD == 0.0 || DevSD == NA)
   {
   cerr << "error in CalcPSS: bumDevSD" << endl;
   exit(1);
   }
// Run through cells.
int n = 0;
dvariable PSS = 0.0;
dvariable sdev;
for (int i = ibot; i<=itop; i++)
   for (int j = MinJ; j<=jtop; j++)
      {
      // Skip if any item is NA.
      if (  DevX(i,j) == NA ) continue;
      // Add to accumulators.
      n++;
      sdev = DevX(i,j)/DevSD;
      PSS += sdev * sdev;
      }
// Load return vector.
dvar_vector PSS_vec(1,3);
PSS_vec(1) = PSS;
PSS_vec(2) = n;
PSS_vec(3) = sqrt(PSS/n);
RETURN_ARRAYS_DECREMENT();
return(PSS_vec);
}
//
// Next the dvar_vector (arg) version.
dvar_vector CalcPSS(_CONST dvar_vector& DevX, _CONST dvariable& DevSD,
                    int MinI, double NA)
{
RETURN_ARRAYS_INCREMENT();
// Determine and check sizes.
int ibot = DevX.indexmin();
int itop = DevX.indexmax();
//
if (MinI < ibot || MinI > itop)
   {
   cerr << "error in CalcPSS: start element out of range" << endl;
   exit(1);
   }
//
// Check SD.
if (DevSD == 0.0 || DevSD == NA)
   {
   cerr << "error in CalcPSS: bum DevSD" << endl;
   exit(1);
   }
// Run through cells.
int n = 0;
dvariable PSS = 0.0;
dvariable sdev;
for (int i = MinI; i<=itop; i++)
      {
      // Skip if any item is NA.
      if (  DevX(i) == NA ) continue;
      // Add to accumulators.
      n++;
      sdev = DevX(i)/DevSD;
      PSS += sdev * sdev;
      }
// Load return vector.
dvar_vector PSS_vec(1,3);
PSS_vec(1) = PSS;
PSS_vec(2) = n;
PSS_vec(3) = sqrt(PSS/n);
RETURN_ARRAYS_DECREMENT();
return(PSS_vec);
}
