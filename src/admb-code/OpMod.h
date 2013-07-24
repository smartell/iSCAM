/**
 * \file OpMod.h
 * \author Steve Martell
 * \date Jun 4, 2013
**/

#ifndef _SCENARIO_DATA_H
#define _SCENARIO_DATA_H
class ScenarioData
{
private:
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

  dvector m_dAllocation;

  // Growth & Maturity data
  dvector m_d_linf;
  dvector m_d_vonbk;
  dvector m_d_to;
  dvector m_d_a;
  dvector m_d_b;
  dvector m_d_ah;
  dvector m_d_gh;

  // Catch Data
  int     m_nCtNobs;
  dmatrix m_dCatchData;

  // Survey data
  int      m_nIt;
  ivector  m_nItNobs;
  ivector  m_nSurveyType;
  d3_array m_dSurveyData;

  // Composition data
  int      m_nAgears;
  ivector  m_nAnobs;
  ivector  m_nAsage;
  ivector  m_nAnage;
  d3_array m_dA;

  // Empirical weight-at-age
  int      m_nWtNobs;
  d3_array m_dWt_avg;
  d3_array m_dWt_mat;

public:
  ScenarioData(
      const int&      _nStock,
      const int&      _nArea,
      const int&      _nSex,
      const int&      _nSyr,
      const int&      _nNyr,
      const int&      _nPyr,
      const int&      _nSage,
      const int&      _nNage,
      const int&      _nGear,
      const dvector&  _dAllocation,
      const dvector&  _d_linf,
      const dvector&  _d_vonbk,
      const dvector&  _d_to,
      const dvector&  _d_a,
      const dvector&  _d_b,
      const dvector&  _d_ah,
      const dvector&  _d_gh,
      const int&      _nCtNobs,
      const dmatrix&  _dCatchData,
      const int&      _nIt,
      const ivector&  _nItNobs,
      const ivector&  _nSurveyType,
      const d3_array& _dSurveyData,
      const int&      _nWtNobs,
      const d3_array& _dWt_avg,
      const d3_array& _dWt_mat );
  ~ScenarioData();
  
  /* Friends */
};
#endif




#ifndef _SCENARIO_PARAMETERS_H
#define _SCENARIO_PARAMETERS_H
class ScenarioParameters 
{
private:
  /* Parameters used in the constructor */
  const dvector&  m_dBo;
  const dvector&  m_dSteepness;
  const dvector&  m_dM;
  const dvector&  m_dRbar;
  const dvector&  m_dRinit;
  const double&   m_dRho;
  const double&   m_dVarphi;
  const dvector&  m_dLog_M_devs;
  const dmatrix&  m_dLog_Rbar_devs;
  const dmatrix&  m_dLog_Rinit_devs;
  const d3_array& m_dSelPars;
  const dmatrix&  m_dFt;

public:
  ScenarioParameters(
          const dvector&  _dBo,
          const dvector&  _dSteepness,
          const dvector&  _dM,
          const dvector&  _dRbar,
          const dvector&  _dRinit,
          const double&   _dRho,
          const double&   _dVarphi,
          const dvector&  _dLog_M_devs,
          const dmatrix&  _dLog_Rbar_devs,
          const dmatrix&  _dLog_Rinit_devs,
          const d3_array& _dSelPars,
          const dmatrix&  _dFt
                    );
  ~ScenarioParameters();

  /* Friends */

};
#endif









#ifndef MODELDATA_H
#define MODELDATA_H
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
class Scenario: public ModelData, public ModelParams
{
private:

public:
  Scenario();
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
class OperatingModel //: public Scenario
{
private:
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

  dvector dAge;
  dvector dRo;
  dvector dSteepness;
  dvector dM;
  dvector dAvgRec;
  dvector dInitRec;
  dvector dRho;
  dvector dVartheta;

  dvector dKappa;

  // Selectivity parameters
  ivector  nSel_type;
  d3_array d3_selPars;
  d5_array d5_logSel;


  // Class aggregations
  Scenario m_cScenario;
public:
  OperatingModel(Scenario &cScenario);
  

  ~OperatingModel();

  /* data */
  void initializeVariables(Scenario& cS);
  void calcSelectivities();
};
#endif