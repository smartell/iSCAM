/**
 * \file OpMod.h
 * \author Steve Martell
 * \date Jun 4, 2013
**/

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

public:
  ModelData(const int nStock=1, const int nArea=1, const int nSex=1, const int nSyr=1950,
            const int nNyr=1999, const int nPyr=2020, const int nSage=1, 
            const int nNage=10)
  :m_nStock(nStock),m_nArea(nArea),m_nSex(nSex),m_nSyr(nSyr),m_nNyr(nNyr),m_nPyr(nPyr),
   m_nSage(nSage),m_nNage(nNage)
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

  /* setters */
  void set_nStock(int x)  { m_nStock  = x;   }
  void set_nArea(int x)   { m_nArea   = x;   }
  void set_nSex(int x)    { m_nSex    = x;   }
  void set_nSyr(int x)    { m_nSyr    = x;   }
  void set_nNyr(int x)    { m_nNyr    = x;   }
  void set_nPyr(int x)    { m_nPyr    = x;   }
  void set_nSage(int x)   { m_nSage   = x;   }
  void set_nNage(int x)   { m_nNage   = x;   }

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

public:
  ModelParams(){};
  ModelParams(const dvector log_ro, const dvector steepness,const dvector log_m,
              const dvector log_avgrec, const dvector log_initrec, 
              const dvector rho,const dvector vartheta)
  :m_log_ro(log_ro),m_steepness(steepness),m_log_m(log_m),
   m_log_avgrec(log_avgrec),m_log_initrec(log_initrec), m_rho(rho),m_vartheta(vartheta)
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
            const dvector log_ro, 
            const dvector steepness,
            const dvector log_m,
            const dvector log_avgrec,
            const dvector log_initrec,
            const dvector rho,
            const dvector vartheta);

  ~Scenario();

 
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

  dvector dRo;
  dvector dSteepness;
  dvector dM;
  dvector dAvgRec;
  dvector dInitRec;
  dvector dRho;
  dvector dVartheta;

  // Class aggregations
  Scenario m_cScenario;
public:
  OperatingModel(Scenario &cScenario);
  

  ~OperatingModel();

  /* data */
  void initializeVariables(Scenario& cS);
};
#endif