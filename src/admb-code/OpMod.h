/**
* \file OpMod.h
* \bief This is the operating model code for the population reference model.
* \author Steve Martell
* \date Jun 4, 2013
*
* DETAILED DESCRIPTION OF THE CLASS OperatingModel
*  - There is a single constructor that takes two arguments which are the
*    structs fro the data types (s_iSCAMdata), and the variable types 
*    (s_iSCAMvariables).  Once instantiated, a few of the methods of the 
*    OperatingModel class executed to initialize variables that are used 
*    in the main routines of the OperatingModel class.
* 
*/
#include <admodel.h>
#include "iscam.htp"

// |---------------------------------------------------------------|
// | Structures for constants used in the general operating model. |
// |---------------------------------------------------------------|
#ifndef _SCENARIO_DATA_H
#define _SCENARIO_DATA_H

/** \brief  Structure for operating model constants
  \author  Steven Martell
  \author $LastChangedBy$
  \date `date +%Y-%m-%d`
  \date $LastChangedDate$
  \version $Rev$  \sa

  To be deprecated, now OperatingModel class is derived from model_data
**/
struct s_iSCAMdata
{
  // |------------------|
  // | Model dimensions |
  // |------------------|
  int nStock;               //!< Number of stocks
  int nArea;                //!< Number of management areas
  int nSex;                 //!< Number of distinct sexes in population
  int nSyr;                 //!< Initial year
  int nNyr;                 //!< Terminal year of available data
  int nPyr;                 //!< Terminal year of the projection period
  int nSage;                //!< Initial age class (years)
  int nNage;                //!< Terminal plus group age class (years)
  int nGear;                //!< Number of distinct fleets incl. survey vessels
  int nFleet;               //!< Number of fishing fleets with allocation > 0.

  dvector dAllocation;
  ivector nFleetIndex;

  // Growth & Maturity data
  dvector d_linf;
  dvector d_vonbk;
  dvector d_to;
  dvector d_a;
  dvector d_b;
  dvector d_ah;
  dvector d_gh;

  // Catch Data
  int     nCtNobs;
  dmatrix dCatchData;

  // Survey data
  int      nIt;
  ivector  nItNobs;
  ivector  n_survey_type;
  d3_array* dSurveyData;

  // Composition data
  int      nAgears;
  ivector  nAnobs;
  ivector  nAsage;
  ivector  nAnage;
  d3_array* dA;

  // Empirical weight-at-age
  int      nWtNobs;
  d3_array* d3_wt_avg;
  d3_array* d3_wt_mat;

  // control vector;
  dvector d_iscamCntrl;

  
};
#endif




// |----------------------------------------------------------------------------|
// | Structure for variable parameters that define the operating model dynamics |
// |----------------------------------------------------------------------------|
#ifndef _SCENARIO_PARAMETERS_H
#define _SCENARIO_PARAMETERS_H
struct s_iSCAMvariables 
{
  /* Parameters used in the constructor */
  dvector   d_log_ro;
  dvector   d_steepness;
  dvector   d_log_m;
  dvector   d_log_rbar;
  dvector   d_log_rinit;
  dvector   d_rho;
  dvector   d_varphi;
  dvector   dLog_M_devs;
  dmatrix   dLog_Rbar_devs;
  dmatrix   dLog_Rinit_devs;

  /* Selectivity parameters */
  ivector   nSel_type;
  imatrix   nSel_block;
  d3_array* dSelPars;
  d4_array* d4_log_sel;

  
  /* Mortality rate arrays */
  d3_array* d3_Ft;
  d3_array* d3_Mt;
  d3_array* d3_St;


};
#endif









// |--------------------------------------------------|
// | Class definition for the general operating model |
// |--------------------------------------------------|
#ifndef OPERATINGMODEL_H
#define OPERATINGMODEL_H
class OperatingModel : public model_data
{
private:
  // |------------------|
  // | Model dimensions |
  // |------------------|
  int nStock;        
  int nArea;
  int nSex;
  int nSyr;
  int nNyr;
  int nPyr;
  int nSage;
  int nNage;
  int nGear;
  int nFleet;
  dvector dAge;

 

  // |------------------|
  // | Model variablies |
  // |------------------|
  dvector dRo;
  dvector dSteepness;
  dvector dM;
  dvector dAvgRec;
  dvector dInitRec;
  dvector dRho;
  dvector dVartheta;
  dvector dLog_m_devs;
  dmatrix dLog_rbar_devs;
  dmatrix dLog_init_rec_devs;
  dvector dKappa;
  
  d3_array d3_Ft;
  d3_array d3_Mt;
  

  // Selectivity parameters
  ivector  nSel_type;
  imatrix  nSel_block;
  d3_array d3_selPars;
  d4_array d4_log_sel;
  d5_array d5_logSel;     // deprecate

  // derived variables for SRR 
  dvector m_kappa;
  dvector m_so;
  dvector m_beta;
  dvector m_dSbo;

  d3_array d3_Nt;
  d3_array d3_Zt;
  d3_array d3_St;
  // d3_array d3_wt_avg;
  // d3_array d3_wt_mat;


  // reference points
  dvector m_dFmsy;

  // stock assessment results
  dvector m_est_bo;
  dvector m_est_fmsy;
  dvector m_est_msy;
  dvector m_est_bmsy;
  dvector m_est_sbt;
  dvector m_est_bt;

  int m_nSeed;    /**< Random number seed */
protected:
  dvector m_dTac;
  dmatrix m_dFt;

  // | Catch data for assessment.
  int     m_nCtNobs;
  int     m_nCtNobs_counter;
  dmatrix m_dCatchData;

  // | Relative abundance index.
  ivector m_n_it_nobs;
  d3_array m_d3_survey_data;

  // | Age composition
  ivector m_n_A_nobs;
  d3_array m_d3_A;

  // | Empirical weight-at-age data
  int m_nWtNobs;
  dmatrix m_imp_wt_avg;
  d3_array m_d3_wt_avg;
  d3_array m_d3_wt_mat;

public:
  // OperatingModel(Scenario &cScenario);
  OperatingModel(const s_iSCAMdata&  mse_data, 
                 const s_iSCAMvariables& mse_vars,
                 int argc,
                 char * argv[]);
  

  ~OperatingModel();

  /* data */
  void initializeConstants(const s_iSCAMdata& cS);
  void initializeVariables(const s_iSCAMvariables& cS);
  void calcSelectivities();
  void runScenario(const int& seed);

protected:
  void calcStockRecruitment();
  void conditionReferenceModel();
  void calcReferencePoints();
  void calcTAC();
  void implementFisheries(const int& iyr);
  void calcRelativeAbundance(const int& iyr);
  void updateReferencePopulation(const int& iyr);
  void calcGrowth(const int& iyr);
  void generateStockAssessmentData(const int& iyr);
  void runStockAssessment();

};
#endif