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



#ifndef _SCENARIO_DATA_H
#define _SCENARIO_DATA_H
struct s_iSCAMdata
{
  // Model dimensions
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
  ivector  nSurveyType;
  d3_array* dSurveyData;

  // Composition data
  int      nAgears;
  ivector  nAnobs;
  ivector  nAsage;
  ivector  nAnage;
  d3_array* dA;

  // Empirical weight-at-age
  int      nWtNobs;
  d3_array* dWt_avg;
  d3_array* dWt_mat;

  // control vector;
  dvector cntrl;

  
};
#endif



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









#ifndef MODELDATA_H
#define MODELDATA_H
// To be deprecated
class ModelData
{
public:
  // Model dimensions
  int m_nStock;
  int m_nArea;
  int m_nSex;
  int m_nSyr;
  int m_nNyr;
  int m_nPyr;
  int m_nSage;
  int m_nNage;
  int m_nGear;

public:
  ModelData(const int nStock=1, const int nArea=1, const int nSex=1, const int nSyr=1950,
            const int nNyr=1999, const int nPyr=2020, const int nSage=1, 
            const int nNage=10,const int nGear=1)
  :m_nStock(nStock),m_nArea(nArea),m_nSex(nSex),m_nSyr(nSyr),m_nNyr(nNyr),m_nPyr(nPyr),
   m_nSage(nSage),m_nNage(nNage),m_nGear(nGear)
  {};

  
  ~ModelData();

  
  /* getters */
  int     get_nStock()    { return m_nStock;   }
  int     get_nArea()     { return m_nArea;    }
  int     get_nSex()      { return m_nSex;     }
  int     get_nSyr()      { return m_nSyr;     }
  int     get_nNyr()      { return m_nNyr;     }
  int     get_nPyr()      { return m_nPyr;     }
  int     get_nSage()     { return m_nSage;    }
  int     get_nNage()     { return m_nNage;    }
  int     get_nGear()     { return m_nGear;    }

  /* setters */
  void set_nStock(int x)  { m_nStock  = x;   }
  void set_nArea(int x)   { m_nArea   = x;   }
  void set_nSex(int x)    { m_nSex    = x;   }
  void set_nSyr(int x)    { m_nSyr    = x;   }
  void set_nNyr(int x)    { m_nNyr    = x;   }
  void set_nPyr(int x)    { m_nPyr    = x;   }
  void set_nSage(int x)   { m_nSage   = x;   }
  void set_nNage(int x)   { m_nNage   = x;   }
  void set_nGear(int x)   { m_nGear   = x;   }

};
#endif








#ifndef MODELPARAMS_H
#define MODELPARAMS_H
//  To be deprecated
class ModelParams
{
public:
  dvector m_log_ro;
  dvector m_steepness;
  dvector m_log_m;
  dvector m_log_avgrec;
  dvector m_log_initrec;
  dvector m_rho;
  dvector m_vartheta;


  d3_array m_selpars;


public:
  ModelParams(){};
  ModelParams(const dvector log_ro, const dvector steepness,const dvector log_m,
              const dvector log_avgrec, const dvector log_initrec, 
              const dvector rho,const dvector vartheta,
              const d3_array & selpars)
  :m_log_ro(log_ro),m_steepness(steepness),m_log_m(log_m),
   m_log_avgrec(log_avgrec),m_log_initrec(log_initrec), m_rho(rho),m_vartheta(vartheta),
   m_selpars(selpars)
  {};


  ~ModelParams();

  /* getters */
  dvector get_log_ro       (){ return m_log_ro;      }
  dvector get_steepness    (){ return m_steepness;   }
  dvector get_log_m        (){ return m_log_m;       }
  dvector get_log_avgrec   (){ return m_log_avgrec;  }
  dvector get_log_initrec  (){ return m_log_initrec; }
  dvector get_rho          (){ return m_rho;         }
  dvector get_vartheta     (){ return m_vartheta;    }

  /* setters */
  dvector set_log_ro      (dvector t1) { m_log_ro      = t1;}
  dvector set_steepness   (dvector t1) { m_steepness   = t1;}
  dvector set_log_m       (dvector t1) { m_log_m       = t1;}
  dvector set_log_avgrec  (dvector t1) { m_log_avgrec  = t1;}
  dvector set_log_initrec (dvector t1) { m_log_initrec = t1;}
  dvector set_rho         (dvector t1) { m_rho         = t1;}
  dvector set_vartheta    (dvector t1) { m_vartheta    = t1;}

};
#endif









#ifndef SCENARIO_H
#define SCENARIO_H
// To be deprecated
class Scenario: public ModelData, public ModelParams
{
private:
  // ScenarioData       m_ScenarioData;
  // ScenarioParameters m_ScenarioParameters;
  // Scenario();
public:
  // Scenario(const ScenarioData& data, const ScenarioParameters& params);

  Scenario( const int nStock, 
            const int nArea, 
            const int nSex, 
            const int nSyr, 
            const int nNyr, 
            const int nPyr, 
            const int nSage, 
            const int nNage,
            const int nGear,
            const dvector log_ro, 
            const dvector steepness,
            const dvector log_m,
            const dvector log_avgrec,
            const dvector log_initrec,
            const dvector rho,
            const dvector vartheta,
            const d3_array& selpar,
            const ivector& sel_type );


  ~Scenario();
  ivector m_sel_type;
 
};

#endif


#ifndef OPERATINGMODEL_H
#define OPERATINGMODEL_H
class OperatingModel 
{
private:
  // |-----------------|
  // | Model constants |
  // |-----------------|
  
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

  // LINKS TO MANAGE ARRAY INDEXING
  int n_ags;
  int n_ag;
  int n_gs;
  ivector n_area;
  ivector n_group;
  ivector n_sex;
  imatrix pntr_ag;
  imatrix pntr_gs;
  d3_array pntr_ags;

  dvector dAllocation;
  ivector nFleetIndex;


  dvector dAge;
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

  d3_array d3_Ct;

  // Survey data
  int       nIt;
  ivector   nItNobs;
  ivector   nSurveyType;
  d3_array  dSurveyData;

   // Composition data
  int      nAgears;
  ivector  nAnobs;
  ivector  nAsage;
  ivector  nAnage;
  d3_array dA;

   // Empirical weight-at-age
  int      nWtNobs;
  dmatrix  dWt_bar;
  dmatrix  dEt_bar;
  d3_array dWt_avg;
  d3_array dWt_mat;


  // Model variables
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
  dvector nCntrl;
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
  d3_array d3_wt_avg;
  d3_array d3_wt_mat;


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
  dvector m_dFt;
  int     m_nCtNobs;
  dmatrix m_dCatchData;

public:
  // OperatingModel(Scenario &cScenario);
  OperatingModel(const s_iSCAMdata&  mse_data, const s_iSCAMvariables& mse_vars);
  

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
  void updateReferencePopulation(const int& iyr);
  void generateStockAssessmentData(const int& iyr);
  void runStockAssessment();

};
#endif