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

  dvector dAllocation;      //!< test        
  ivector nFleetIndex;      //!< test

  // Growth & Maturity data
  dvector d_linf;     //!< test
  dvector d_vonbk;      //!< test
  dvector d_to;     //!< test
  dvector d_a;      //!< test
  dvector d_b;      //!< test
  dvector d_ah;     //!< test
  dvector d_gh;     //!< test

  // Catch Data
  int     nCtNobs;      //!< test
  dmatrix dCatchData;     //!< test

  // Survey data
  int      nIt;     //!< test
  ivector  nItNobs;     //!< test
  ivector  n_survey_type;     //!< test
  d3_array* dSurveyData;      //!< test

  // Composition data
  int      nAgears;     //!< test
  ivector  nAnobs;      //!< test
  ivector  nAsage;      //!< test
  ivector  nAnage;      //!< test
  d3_array* dA;     //!< test

  // Empirical weight-at-age
  int      nWtNobs;     //!< test
  d3_array* d3_wt_avg;      //!< test
  d3_array* d3_wt_mat;      //!< test

  // control vector;
  dvector d_iscamCntrl;     //!< test

  
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
  dvector   d_log_ro;         //!< Equilibrium unfished recruitment.
  dvector   d_steepness;      //!< Steepness of the stock-recruitment relationship.
  dvector   d_log_m;          //!< log of the instantaneous natural mortality rate
  dvector   d_log_rbar;       //!< Annaul average recruitment
  dvector   d_log_rinit;      //!< Initial average recruitment
  dvector   d_rho;            //!< Proportion of variance due to observation error.
  dvector   d_varphi;         //!< Total precision (1/variance)
  dvector   dLog_M_devs;      //!< Annual deviates in natural mortality
  dmatrix   dLog_Rbar_devs;   //!< Annual recruitment deviates
  dmatrix   dLog_Rinit_devs;  //!< Initial recruitment deviates

  /* Selectivity parameters */
  ivector   nSel_type;      //!< Type of selectivity function
  imatrix   nSel_block;     //!< Number of selectivity blocks
  d3_array* dSelPars;       //!< Selectivity parameters
  d4_array* d4_log_sel;     //!< Selectivity coefficients

  
  /* Mortality rate arrays */
  d3_array* d3_Ft;          //!< Annual instantaneous fishing mortality rate
  d3_array* d3_Mt;          //!< Annual instantaneous natural mortality rate
  d3_array* d3_St;          //!< Annual survival rate


};
#endif









// |--------------------------------------------------|
// | Class definition for the general operating model |
// |--------------------------------------------------|
#ifndef OPERATINGMODEL_H
#define OPERATINGMODEL_H
class mse_data
{
private:
  int m_nPyr;

public:
  ~mse_data();
  mse_data(int argc, char * argv[]);

  friend class OperatingModel;

};


class OperatingModel : public model_data, public mse_data
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
  dvector m_dTac;                   //!< Vector of total allowable catch
  dmatrix m_dFt;                    //!< Vector of annual fishing mortality rates

  // | Catch data for assessment.
  int     m_nCtNobs;                //!< total number of catch observations
  int     m_nCtNobs_counter;        //!< counter for number of catch observations
  ivector m_catch_sex_composition;  //!< ivector of sex composition of catch for gear
  ivector m_catch_type;             //!< ivector of catch type for each gear
  dmatrix m_dCatchData;             //!< matrix of catch observations 

  // | Relative abundance index.
  ivector m_n_it_nobs;              //!< total number of survey observations
  ivector m_n_it_counter;           //!< counter for the number of survey observations
  dvector m_survey_q;               //!< vector of survey catchability coeffients
  d3_array m_d3_survey_data;        //!< array of survey data

  // | Age composition
  ivector m_n_A_nobs;               //!< integer vector of number of age composition data
  d3_array m_d3_A;                  //!< Age composition data  

  // | Empirical weight-at-age data
  int m_nWtNobs;                    //!< Number of empirical weight-at-age rows
  dmatrix m_imp_wt_avg;             //!< Average weight-at-age 
  d3_array m_d3_wt_avg;             //!< Average weight-at-age 
  d3_array m_d3_wt_mat;             //!< Average mature weight-at-age 

  // | Selectivity
  d4_array m_d4_log_sel;

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
  void calcSurveyCatchability();
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