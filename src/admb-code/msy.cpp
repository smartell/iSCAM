/** \brief  Msy class
	
	The Msy class provides computational support for age-structured models
	in which MSY-based (Maximum Sustainable Yield) reference points are 
	required.
	
© Copyright 2012 UBC Fisheries Centre - Martell. All Rights Reserved.

	\author Martell UBC Fisheries Centre
	\author $LastChangedBy$
	\date 2012-07-29
	\date $LastChangedDate$
	\version $Rev$	\sa
**/

#ifndef _MSY_H_
#define _MSY_H_

#include <admodel.h>
#include <fvar.hpp>

class Msy{
private:
	double m_ro;
	double m_h;
	double m_M;
	dvector m_wa;
	dvector m_fa;
	dmatrix m_V;	// Selectivity for each gear (rows) at age (col)
	
	double m_phie;
	dvector m_delta; // Newton-Raphson step for iteratively solving for Fmsy
	dvector m_fmsy;
	dvector m_msy;
	double m_bmsy;
	double m_rmsy;
	double m_spr_msy;
	double m_bo;
	
	dvector m_ye;
	double m_be;
	double m_re;
	double m_spr;
	
	double m_dYe;
	double m_d2Ye;
	
	
public:
	
	Msy()		// default constructor
	{
		m_ro = 1;
		m_h  = 0.75;
		m_M  = 0.3;
	}
	Msy(double ro, double h, double m, dvector wa, dvector fa, dmatrix V);
	
	~Msy(){}	// destructor
	
	// Getters
	double     getRo() { return m_ro;      }
	double   getPhie() { return m_phie;    }
	dvector  getFmsy() { return m_fmsy;    }
	double   getBmsy() { return m_bmsy;    }
	double   getRmsy() { return m_rmsy;    }
	dvector   getMsy() { return m_msy;     }
	dvector    getYe() { return m_ye;      }
	double getSprMsy() { return m_spr_msy; }
	double     getBo() { return m_bo;      }
	double     getBe() { return m_be;      }
	double     getRe() { return m_re;      }
	double    getSpr() { return m_spr;     }
	
	// Setters
	void set_ro(double ro)   { m_ro = ro; }
	void set_m(double m)     { m_M  = m;  }
	void set_h(double& h)    { m_h  = h;  }
	void set_fa(dvector& fa) { m_fa = fa; }
	void set_wa(dvector& wa) { m_wa = wa; }
	void set_fa(dmatrix& V)  { m_V  = V;  }
	
	// Member functions
	void calc_phie();
	void calc_phie(double& _m, dvector& _fa);
	void calc_bo(double& _m, dvector& _va);
	void calc_equilibrium(dvector& fe);
	
	void get_fmsy(dvector& fe);
	void get_fmsy(dvector& fe, dvector& ak);
};

#endif

// Constructor with arguments (most likely way user will create an Msy object)
Msy::Msy(double ro, double h, double m, dvector wa, dvector fa, dmatrix V)
{
	m_ro = ro;
	m_h  = h;
	m_M  = m;
	m_wa = wa;
	m_fa = fa;
	m_V  = V;
	
	
	calc_phie();
}

// Use Newton-Raphson method to get Fmsy
void Msy::get_fmsy(dvector& fe)
{
	/*
		Iteratively solve the derivative of the catch equation
		to find values of fe that correspond to dYe.df=0
	*/
	int iter = 0;
	m_delta = 1.0;
	do
	{
		calc_equilibrium(fe);
		fe += m_delta;
		iter++;
		//cout<<iter<<"norm(delta) ="<<norm(m_delta)<<endl;
	}
	while ( norm(m_delta)>1.e-12 && iter < 100 );
	m_fmsy= fe;
	m_msy = m_ye;
	m_bmsy = m_be; 
	m_rmsy = m_re;
	m_spr_msy = m_spr;
	
}

// calculate msy value given an allocation for each gear type.
void Msy::get_fmsy(dvector& fe, dvector& ak)
{
	/*
		This algorithm iteratively solves for a vector fe that will
		maximize the sum of catches over all fleets where the catch 
		is allocated (ak, where sum ak=1) to each fleet.
		
		How it works:
		 - Run the model with all fleets having the same Fe
		 - Then calculate fmultiplier (lambda) as the ratio 
		   allocation:yield/sum(yield) and correct Fe's down or up
		   such that the yield proportions matches the allocation.
		 - Run the equilibrium model again to get the derivative 
		   information for the total yield and use Newton-Rhaphson to
		   update the estimate of fbar.
		 - Repeat above steps until derivative for the total yield 
		   approaches zero (1.e-12).
		
	*/
	
	int iter =1;
	ak             = ak/sum(ak);
	double    fbar = mean(fe);
	dvector lambda = ak/mean(ak);
	dvector     fk = elem_prod(fe,lambda);
	
	do
	{
		// Run the model with constant fk = fbar to get
		// lambda multipliers.
		fk     = fbar;
		calc_equilibrium(fk);
		lambda = elem_div(ak,m_ye/sum(m_ye));
		
		// Run the model with fk = lambda*fbar to get
		// derivatives of the total catch equation to update fbar;
		fk     = fbar*lambda;
		calc_equilibrium(fk);
		
		fbar   = fbar - m_dYe/m_d2Ye;
		iter++;
		//cout<<iter<<" fbar "<<fbar<<" dYe "<<fabs(m_dYe)<<endl;
	}
	while ( sqrt(square(m_dYe)) >1.e-12 && iter < 100 );
	fe = fk;
}

// calculate survivorship under fished conditions
void Msy::calc_equilibrium(dvector& fe)
{
	/*
		This is the main engine behind the Msy class.
		Using the private member variables, calculate the equilibrium yield
		for a given fishing mortality rate vector (fe). Also, compute the
		derivative information such that this info can be used to compute 
		the corresponding MSY-based reference points (Fmsy, MSY, Bmsy).
		
	*/
	
	// Indexes for dimensions
	int i,j,k;
	int sage,nage,ngear;
	sage  = m_wa.indexmin();
	nage  = m_wa.indexmax();
	ngear = m_V.rowmax();
	
	// Survivorship for fished conditions.
	dvector lz(sage,nage);
	dmatrix dlz(1,ngear,sage,nage);
	dmatrix d2lz(1,ngear,sage,nage);
	dlz.initialize();
	d2lz.initialize();

	dvector za = m_M + fe * m_V;
	dvector sa = exp(-za);
	dvector oa = 1.-sa;
	dmatrix qa(1,ngear,sage,nage);
	
	for(k=1;k<=ngear;k++)
	{
		qa(k)        = elem_div(elem_prod(elem_prod(m_V(k),m_wa),oa),za);
		dlz(k,sage)  = 0;
		d2lz(k,sage) = 0;
	}
	
	lz(sage) = 1.0;
	for(j=sage+1; j<=nage; j++)
	{
		lz(j)   = lz(j-1)  * sa(j-1);
		if( j==nage )
		{
			lz(j) = lz(j)/oa(j);
		}
		
		for(k=1; k<=ngear; k++)
		{
			dlz(k)(j)  = sa(j-1) * ( dlz(k)(j-1) - lz(j-1)*m_V(k)(j-1) );
			d2lz(k)(j) = sa(j-1) * ( d2lz(k)(j-1)+ lz(j-1)*m_V(k)(j-1)*m_V(k)(j-1) );
			
			if( j==nage ) // + group derivatives
			{
				dlz(k)(j)  = dlz(k)(j)/oa(j) - lz(j-1)*sa(j-1)*m_V(k)(j)*sa(j)/square(oa(j));
				
				double V1  = m_V(k)(j-1);
				double V2  = m_V(k)(j);
				double oa2 = oa(j)*oa(j);
				
				d2lz(k)(j) = d2lz(k)(j)/oa(j) 
							+ 2*lz(j-1)*V1*sa(j-1)*V2*sa(j)/oa2
							+ 2*lz(j-1)*sa(j-1)*V2*V2*sa(j)*sa(j)/(oa(j)*oa2)
							+ lz(j-1)*sa(j-1)*V2*V2*sa(j)/oa2;
			}
		}// gear
		
	}// age
	
	// Incidence functions and associated derivatives
	double      ro = m_ro;
	double   kappa = 4.0*m_h/(1.0-m_h);  // Beverton-Holt model
	double     km1 = kappa-1.0;
	double    phif = lz * m_fa;
	double   phif2 = phif*phif;
	dvector  dphif(1,ngear);
	dvector d2phif(1,ngear);
	dvector   phiq(1,ngear);
	dvector  dphiq(1,ngear);
	dvector d2phiq(1,ngear);
	dvector    dre(1,ngear);
	dvector   d2re(1,ngear);
	
	for(k=1; k<=ngear; k++)
	{
		dphif(k)   = dlz(k)  * m_fa;
		d2phif(k)  = d2lz(k) * m_fa;
		
		// per recruit yield
		phiq(k)    = lz * qa(k);
		
		dvector t1 = elem_div(elem_prod(elem_prod(lz,m_wa),m_V(k)),za);
		dvector t0 = elem_div(oa,za);
		dvector t3 = sa-t0;
		dphiq(k)   = qa(k)*dlz(k) + t1 * t3;
		
		// 2nd derivative for per recruit yield (nasty)
		dvector t2  = 2. * dlz(k);
		dvector V2  = elem_prod(m_V(k),m_V(k));
		dvector t5  = elem_div(elem_prod(m_wa,V2),za);
		dvector t7  = elem_div(m_V(k),za);
		dvector t9  = elem_prod(t5,sa);
		dvector t11 = elem_prod(t5,t0);
		dvector t13 = elem_prod(lz,t5);
		dvector t14 = elem_prod(m_V(k),sa);
		dvector t15 = elem_prod(t7,sa);
		dvector t17 = elem_prod(m_V(k),t0);
		dvector t18 = elem_div(t17,za);
		d2phiq(k)  =  d2lz(k)*qa(k) + t2*t9 - t2*t11 - t13*t14 -2.*t13*t15 + 2.*t13*t18;
		
		
		// 1st & 2nd partial derivatives for recruitment
		dre(k)      = ro/km1*m_phie/(phif2)*dphif(k);
		d2re(k)     = -2.*ro*m_phie*dphif(k)*dphif(k)/(phif2*phif*km1) 
					+ ro*m_phie*d2phif(k)/(phif2*km1);
	}
	
	// Equilibrium calculations
	double    re;
	dvector   ye(1,ngear);
	dvector fstp(1,ngear);
	dvector  dye(1,ngear);
	dmatrix d2ye(1,ngear,1,ngear);
	dmatrix invJ(1,ngear,1,ngear);
	
	re   = ro*(kappa-m_phie/phif)/km1;
	ye   = re*elem_prod(fe,phiq);
	dye  = re*phiq + elem_prod(fe,elem_prod(phiq,dre)) + elem_prod(fe,re*dphiq);
	
	// Caclculate Jacobian matrix (2nd derivative of the catch equation)
	for(j=1; j<=ngear; j++)
	{
		for(k=1; k<=ngear; k++)
		{
			d2ye(j)(k) = fe(j)*phiq(j)*d2re(k) + 2.*fe(j)*dre(k)*dphiq(k) + fe(j)*re*d2phiq(k);
			if(k == j)
			{
				d2ye(j)(k) += 2.*dre(j)*phiq(j)+2.*re*dphiq(j);
			}
		} 
	}
	
	invJ    = -inv(d2ye);
	fstp    = invJ * dye;  //Newton-Raphson step.
	m_delta = fstp;
	m_ye    = ye;
	m_re    = re;
	m_be    = re*phif;
	m_spr   = phif/m_phie;
	m_dYe   = sum(dye);
	m_d2Ye  = sum(diagonal(d2ye));
	
}

// calculate unfished Eggs per recruit
void Msy::calc_phie()
{
	int i,sage,nage;
	sage      = m_fa.indexmin();
	nage      = m_fa.indexmax();
	dvector lx(sage,nage);
	for(i=sage; i<=nage; i++)
	{
		lx(i) = exp( -m_M*(i-sage) );
		if(i==nage) lx(i)/=1.-exp(-m_M);
	}
	m_phie = lx*m_fa;
	m_bo   = m_ro * m_phie;
	//cout<<"This is lx "<<lx<<endl;
	//cout<<"This is phie "<<m_phie<<endl;
}

void Msy::calc_phie(double& _m, dvector& _fa)
{
	int i,sage,nage;
	m_M       = _m;
	m_fa      = _fa;
	sage      = m_fa.indexmin();
	nage      = m_fa.indexmax();
	dvector lx(sage,nage);
	for(i=sage; i<=nage; i++)
	{
		lx(i) = exp( -m_M*(i-sage) );
		if(i==nage) lx(i)/=1.-exp(-m_M);
	}
	m_phie = lx   * m_fa;
	m_bo   = m_ro * m_phie;
}

void Msy::calc_bo(double& _m, dvector& _fa)
{
	calc_phie(_m,_fa);
}