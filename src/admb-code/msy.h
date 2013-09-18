// Header file for msy.cpp
#include <admodel.h>

#ifndef _MSY_H_
#define _MSY_H_
#define MAXITER 300
#define TOL 1.e-04

#include <admodel.h>
#include <fvar.hpp>

/** \brief  Msy class
	
	The Msy class provides computational support for age-structured models
	in which MSY-based (Maximum Sustainable Yield) reference points are 
	required.  It is specifically designed to deal with cases in which there
	are multiple fishing fleets that differ in selectivity and is capable of
	determining optimal fishing mortality rates for each fleet and optimal 
	allocations for each fleet such that the sum of catches acrosss all fleets
	is maximized.  There is also an option to determine MSY-based reference points
	in cases where there are allocation agreements in place.

	LUCIE'S RULE: the derivative of a sum is the sum of its derivatives.
	Lucie says: there is some nasty calculus in here daddy, are you sure
	you've got it right?
	
	© Copyright 2012 UBC Fisheries Centre - Martell. All Rights Reserved.

	\author Steve Martell
	\author $Martell
	\date 2012-07-29
	\date 2012-08-20
	\version 1.1
	\sa
**/
class Msy
{
private:
	int     m_sage;  /**< youngest age class*/
	int     m_nage;
	int     m_ngear;
	bool    m_FAIL;
	
	double  m_ro;
	double  m_h;
	double  m_M;
	dmatrix m_dM;	// Age/sex specific natural mortality rates
	double  m_rho;  // Fraction of total mortality that occurs before spawning.
	dvector m_wa;
	dmatrix m_dWa;
	dvector m_fa;
	dmatrix m_dFa;	// sex specific fecundity
	dmatrix m_V;	// Selectivity for each gear (rows) at age (col)
	d3_array m_d3_V; // Selectivity for each gear, sex, age.
	dmatrix m_lz;   // survivorship under fished conditions.
	
	double  m_phie;
	double  m_phif;
	dvector m_phiq;
	dvector m_fmsy;
	dvector m_msy;
	double  m_bmsy;
	double  m_rmsy;
	double  m_spr_msy;
	double  m_bo;
	
	dvector m_ye;
	double  m_be;
	double  m_bi;	// spawning biomass at the start of the year.
	double  m_re;
	double  m_spr;
	
	double m_dYe;
	double m_d2Ye;
	
	dvector m_f;	// value of the function to minimize (norm(p))
	dvector m_g;	// gradient
	dvector m_p;	// Newton-Raphson step for iteratively solving for Fmsy
	
	
public:
	
	Msy()		// default constructor
	{
		m_ro = 1;
		m_h  = 0.75;
		m_M  = 0.3;
	}
	Msy(double ro, double h, double m, double rho, dvector wa, dvector fa, dmatrix V);
	Msy(double ro, double h, dmatrix m, double rho, dmatrix wa, dmatrix fa,const d3_array *V);


	~Msy(){}	// destructor
	
	// Getters
	bool     getFail() { return m_FAIL;    }  /**< Flag for convergence */
	
	double     getRo() { return m_ro;      }  /**< Return unfished recruits*/
	double   getPhie() { return m_phie;    }  /**< Return unfished spawning biomass per recruit*/
	dvector  getFmsy() { return m_fmsy;    }  /**< Return vector of fishing mortality rates at MSY*/
	double   getBmsy() { return m_bmsy;    }  /**< Return spawning biomass at MSY*/
	double   getRmsy() { return m_rmsy;    }  /**< Return recruitment at MSY*/
	dvector   getMsy() { return m_msy;     }  /**< Return maximum sustainable yield vector*/
	dvector    getYe() { return m_ye;      }  /**< Return equilibrium yield */
	dvector   getdYe() { return m_f;       }  /**< Return derivatives of catch vector*/
	dvector  getphiq() { return m_phiq;    }  /**< Return vector of per-recruit yields*/
	double   getphif() { return m_phif;    }  /**< Return spawning biomass per recruit*/
	double getSprMsy() { return m_spr_msy; }  /**< Return Spawning Potential Ratio at MSY*/
	double     getBo() { return m_bo;      }  /**< Return unfished spawning biomass*/
	double     getBe() { return m_be;      }  /**< Return equilibrium spawning biomass*/
	double     getBi() { return m_bi;      }  /**< Return start of the year equilirbiurm spawning biomass*/
	double     getRe() { return m_re;      }  /**< Return equilibrium recruits*/
	double    getSpr() { return m_spr;     }  /**< Return Spawning Potential Ratio*/
	dmatrix    getLz() { return m_lz;      }  /**< Return unfished survivorship vector*/
	
	
	// Setters
	void set_ro(double ro)   { m_ro = ro; }
	void  set_m(double m)    { m_M  = m;  }
	void  set_h(double& h)   { m_h  = h;  }
	void set_fa(dvector& fa) { m_fa = fa; }
	void set_wa(dvector& wa) { m_wa = wa; }
	void  set_V(dmatrix& V)  { m_V  = V;  }
	
	// Member functions
	void        calc_phie();
	void        calc_phie(double& _m, dvector& _fa);
	void        calc_phie(const dmatrix& _m, const dmatrix& _fa);
	void          calc_bo(double& _m, dvector& _va);
	void          calc_bo(const dmatrix& _m, const dmatrix& _fa);
	void calc_equilibrium(const dvector& fe);
	void  calcEquilibrium(const dvector& fe);
	void         get_fmsy(dvector& fe);
	void         get_fmsy(dvector& fe, dvector& ak);
	void            print();
	
};

#endif
