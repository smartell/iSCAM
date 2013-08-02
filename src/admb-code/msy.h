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
	
© Copyright 2012 UBC Fisheries Centre - Martell. All Rights Reserved.

	\author Martell UBC Fisheries Centre
	\author $Martell$
	\date 2012-07-29
	\date 2012-08-20
	\version $Rev$	\sa
**/
class Msy
{
private:
	int     m_sage;
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
	bool     getFail() { return m_FAIL;    }
	
	double     getRo() { return m_ro;      }
	double   getPhie() { return m_phie;    }
	dvector  getFmsy() { return m_fmsy;    }
	double   getBmsy() { return m_bmsy;    }
	double   getRmsy() { return m_rmsy;    }
	dvector   getMsy() { return m_msy;     }
	dvector    getYe() { return m_ye;      }
	dvector   getdYe() { return m_f;       }
	dvector  getphiq() { return m_phiq;    }
	double   getphif() { return m_phif;    }
	double getSprMsy() { return m_spr_msy; }
	double     getBo() { return m_bo;      }
	double     getBe() { return m_be;      }
	double     getBi() { return m_bi;      }
	double     getRe() { return m_re;      }
	double    getSpr() { return m_spr;     }
	dmatrix    getLz() { return m_lz;      }
	
	
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
