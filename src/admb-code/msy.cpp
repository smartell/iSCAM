#include "msy.h"


// Constructor with arguments (most likely way user will create an Msy object)
Msy::Msy(double ro, double h, double m, double rho, dvector wa, dvector fa, dmatrix V)
{
	m_ro    = ro;
	m_h     = h;
	m_M     = m;
	m_rho   = rho;
	m_wa    = wa;
	m_fa    = fa;
	m_V     = V;
	
	m_FAIL  = false;
	m_sage  = wa.indexmin();
	m_nage  = wa.indexmax();
	m_ngear = V.rowmax();

	m_dM.allocate(1,1,m_sage,m_nage);
	m_dM(1)   = m_M;
	
	m_dWa.allocate(1,1,m_sage,m_nage);
	m_dWa(1)  = wa;
	
	m_dFa.allocate(1,1,m_sage,m_nage);
	m_dFa(1)  = fa;
	
	m_d3_V.allocate(1,1,1,m_ngear,m_sage,m_nage);
	m_d3_V(1) = V;
	
	calc_phie();
}

// Constructor with age/sex-specific mortality fecundity & weight-at-age
Msy::Msy(double ro, double h, dmatrix m, double rho, dmatrix wa, dmatrix fa, const d3_array *V)
{
	m_ro    = ro;
	m_h     = h;
	m_dM    = m;
	m_rho   = rho;
	m_dWa    = wa;
	m_dFa    = fa;

	m_d3_V.allocate(*V);
	m_d3_V  = *V;

	
	m_FAIL  = false;

	calc_phie(m_dM,m_dFa);
}

// Use Newton-Raphson method to get Fmsy
void Msy::get_fmsy(dvector& fe)
{
	/*
		Iteratively solve the derivative of the catch equation
		to find values of fe that correspond to dYe.df=0
		
		Use Newtons method to interatively solve for fe in the range
		of x1 to x2.  If fe is outside the range of x1--x2, then use
		a simple bactrace method method to reduce the step size m_p.
	*/
	int i;
	int iter = 0;
	int n    = size_count(fe);
	double x1, x2;
	dvector fold(1,n);
	x1 = 1.0e-5;
	x2 = 3.0e02;
	m_p      = 1.0;
	// Spawning biomass per recruit for unfished conditions
	// calc_phie(m_M,m_fa);
	// calc_equilibrium(fe);
	calc_phie(m_dM,m_dFa);
	calcEquilibrium(fe);
	do
	{
		iter++;
		fold = m_f;
		// calc_equilibrium(fe);
		calcEquilibrium(fe);
		fe  += m_p;
		
		// check boundary conditions
		for(i = 1; i<=n; i++)
		{
			if( (x1-fe[i])*(fe[i]-x2) < 0.0 ) // backtrack 98% of the newton step.
			{                                 // if outside the boundary conditions.
				fe[i] -= 0.999*m_p[i];       
				cout<<"OOOPS"<<endl;  
			}
		}
		cout<<iter<<" fe "<<fe<<" f "<<m_f<<endl;
		
	}
	while ( norm(m_f) > TOL && iter < MAXITER );
	
	if( iter == MAXITER ) m_FAIL = true;
	
	//if(min(fe)<0) exit(1);
	
	m_fmsy    = fe;
	m_msy     = m_ye;
	m_bmsy    = m_be; 
	m_rmsy    = m_re;
	m_spr_msy = m_spr;
	
}


// Calculate msy value given an allocation for each gear type.
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
		   approaches zero (TOL).
		
	*/
	int i;
	int n = size_count(fe);
	int iter =1;
	ak             = ak/sum(ak);
	double    fbar = mean(fe);
	dvector lambda = ak/mean(ak);
	dvector     fk = elem_prod(fe,lambda);
	dvector p(1,n);
	double x1, x2;
	x1 = 1.0e-5;
	x2 = 3.0e02;
	// Spawning biomass per recruit for unfished conditions
	calc_phie(m_dM,m_dFa);
	
	do
	{
		// Run the model with constant fk = fbar to get
		// lambda multipliers.
		fk     = fbar;
		calcEquilibrium(fk);
		lambda = elem_div(ak,m_ye/sum(m_ye));
		
		// Run the model with fk = lambda*fbar to get
		// derivatives of the total catch equation to update fbar;
		fk     = fbar*lambda;
		calcEquilibrium(fk);
		//cout<<iter<<" fbar "<<fbar<<" dYe "<<fabs(m_dYe)<<" fk "<<fk<<endl;
		
		fbar   = fbar - m_dYe/m_d2Ye;
		iter++;
		
		// Check boundary conditions and reduce step size if necessary.
		if( (x1-fbar)*(fbar-x2)<0.0 )
		{
			cout<<"get_fmsy:: reducing Newton-step"<<endl;
			fbar += 0.98 * m_dYe/m_d2Ye;
			//exit(1);
		}
 		
	}
	while ( sqrt(square(m_dYe)) >TOL && iter < MAXITER );
	fe = fk;

	m_fmsy    = fe;
	m_msy     = m_ye;
	m_bmsy    = m_be; 
	m_rmsy    = m_re;
	m_spr_msy = m_spr;
	
}

// calculate survivorship under fished conditions
void Msy::calc_equilibrium(const dvector& fe)
{
	/* TO BE Deprecated in iSCAM */	
	/*

		This is the main engine behind the Msy class.
		Using the private member variables, calculate the equilibrium yield
		for a given fishing mortality rate vector (fe). Also, compute the
		derivative information such that this info can be used to compute 
		the corresponding MSY-based reference points (Fmsy, MSY, Bmsy).
		
		There are two different survivorship calculations inside this routine:
		lz is the survivorship of the total numbers, and lw is the survivorship
		upto the time of spawning.  The spawning biomass per recruit (phif) depends
		on the survivorship up to the time of spawning.
	*/
	
	// Indexes for dimensions
	int i,j,k;
	int sage,nage,ngear;
	sage  = m_wa.indexmin();
	nage  = m_wa.indexmax();
	ngear = m_V.rowmax();
	
	
	// Survivorship for fished conditions.
	dvector   lz(sage,nage);
	dvector   lw(sage,nage);
	dmatrix  dlz(1,ngear,sage,nage);
	dmatrix d2lz(1,ngear,sage,nage);
	dmatrix  dlw(1,ngear,sage,nage);
	dmatrix d2lw(1,ngear,sage,nage);
	dlz.initialize();
	d2lz.initialize();
	dlw.initialize();
	d2lw.initialize();

	dvector za  = m_M + fe * m_V;
	dvector pza = m_rho*za;
	dvector sa  = exp(-za);
	dvector psa = exp(-pza);
	dvector oa  = 1.-sa;
	dvector poa = 1.-elem_prod(sa,psa);
	dmatrix qa(1,ngear,sage,nage);
	
	for(k=1;k<=ngear;k++)
	{
		qa(k)        = elem_div(elem_prod(elem_prod(m_V(k),m_wa),oa),za);
		dlz(k,sage)  = 0;
		d2lz(k,sage) = 0;
		dlw(k,sage)  = -psa(sage)*m_rho*m_V(k)(sage);
		d2lw(k,sage) =  psa(sage)*square(m_rho)*square(m_V(k)(sage));
	}
	
	lz(sage) = 1.0;
	lw(sage) = 1.0 * psa(sage);
	for(j=sage+1; j<=nage; j++)
	{
		lz(j)   = lz(j-1) * sa(j-1);
		lw(j)   = lz(j)   * psa(j);
		if( j==nage )
		{
			lz(j) = lz(j)/oa(j);
			lw(j) = lz(j-1)*sa(j-1)*psa(j)/oa(j);
		}
		
		for(k=1; k<=ngear; k++)
		{
			// derivatives for survivorship
			dlz(k)(j)  = sa(j-1) * ( dlz(k)(j-1) - lz(j-1)*m_V(k)(j-1) );
			d2lz(k)(j) = sa(j-1) * ( d2lz(k)(j-1)+ lz(j-1)*m_V(k)(j-1)*m_V(k)(j-1) );
			
			// derivatives for spawning survivorship
			dlw(k)(j)  = -lz(j)*m_rho*m_V(k)(j)*psa(j); 
			d2lw(k)(j) =  lz(j)*square(m_rho)*square(m_V(k)(j))*psa(j);
			
			if( j==nage ) // + group derivatives
			{
				dlz(k)(j)  = dlz(k)(j)/oa(j) - lz(j-1)*sa(j-1)*m_V(k)(j)*sa(j)/square(oa(j));
				
				dlw(k)(j)  = -lz(j-1)*sa(j-1)*m_rho*m_V(k)(j)/oa(j)
							- lz(j-1)*psa(j)*m_V(k)(j)*sa(j)/square(oa(j));
				
				double V1  = m_V(k)(j-1);
				double V2  = m_V(k)(j);
				double oa2 = oa(j)*oa(j);
				
				d2lz(k)(j) = d2lz(k)(j)/oa(j) 
							+ 2*lz(j-1)*V1*sa(j-1)*V2*sa(j)/oa2
							+ 2*lz(j-1)*sa(j-1)*V2*V2*sa(j)*sa(j)/(oa(j)*oa2)
							+ lz(j-1)*sa(j-1)*V2*V2*sa(j)/oa2;
				
				d2lw(k)(j) = lz(j-1)*square(m_rho)*square(V2)*psa(j)/oa(j)
							+ 2*lz(j-1)*m_rho*square(V2)*psa(j)*sa(j)/oa2
							+ 2*lz(j-1)*psa(j)*square(V2)*square(sa(j))/(oa(j)*oa2)
							+ lz(j-1)*psa(j)*square(V2)*sa(j)/oa2;
			}
		}// gear
		
	}// age
	
	// Incidence functions and associated derivatives
	double      ro = m_ro;
	double   kappa = 4.0*m_h/(1.0-m_h);  // Beverton-Holt model
	double     km1 = kappa-1.0;
	//double    phif = lz * m_fa;
	double    phif = lw * m_fa;
	double   phif2 = phif*phif;
	dvector  dphif(1,ngear);
	dvector d2phif(1,ngear);
	dvector   phiq(1,ngear);
	dvector  dphiq(1,ngear);
	dvector d2phiq(1,ngear);
	dvector    dre(1,ngear);
	dvector   d2re(1,ngear);
	dvector     t1(sage,nage);
	
	for(k=1; k<=ngear; k++)
	{
		dphif(k)   = dlz(k)  * m_fa;
		d2phif(k)  = d2lz(k) * m_fa;
		//dphif(k)   = dlw(k)  * m_fa;
		//d2phif(k)  = d2lw(k) * m_fa;
		
		// per recruit yield
		phiq(k)    = lz * qa(k);
		if(ngear==1)
		{
			// dphiq = wa*oa*va*dlz/za + lz*wa*va^2*sa/za - lz*wa*va^2*oa/za^2
			t1 = elem_div(elem_prod(elem_prod(lz,m_wa),square(m_V(k))),za);
		}
		else
		{
			// dphiq = wa*oa*va*dlz/za + lz*wa*va*sa/za - lz*wa*va*oa/za^2
			t1 = elem_div(elem_prod(elem_prod(lz,m_wa),m_V(k)),za);
		}
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
		dre(k)      = ro*m_phie*dphif(k)/(phif2*km1);
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
	
	re   = ro*(kappa-m_phie/phif) / km1;
	ye   = re*elem_prod(fe,phiq);
	dye  = re*phiq + elem_prod(fe,elem_prod(phiq,dre)) + elem_prod(fe,re*dphiq);
	
	// Caclculate Jacobian matrix (2nd derivative of the catch equation)
	for(j=1; j<=ngear; j++)
	{
		for(k=1; k<=ngear; k++)
		{
			d2ye(k)(j) = fe(j)*phiq(j)*d2re(k) + 2.*fe(j)*dre(k)*dphiq(k) + fe(j)*re*d2phiq(k);
			if(k == j)
			{
				d2ye(j)(k) += 2.*dre(j)*phiq(j)+2.*re*dphiq(j);
			}
		} 
	}
	// Newton-Raphson step.
	invJ   = -inv(d2ye);
	fstp   = invJ * dye;  
	
	// Set private members
	m_p    = fstp;
	m_ye   = ye;
	m_re   = re;
	m_be   = re*phif;
	m_bi   = re*(lz*m_fa);
	m_spr  = phif/m_phie;
	m_dYe  = sum(dye);
	m_d2Ye = sum(diagonal(d2ye));
	m_g    = diagonal(d2ye);		//Gradient vector
	m_f    = dye;   				//Value of the function to minimize
	
	cout<<"OB1 fe "<<fe<<" ye "<<ye<<" dye "<<dye<<endl;
}

























// calculate survivorship under fished conditions with sex-specific age-schedules
void Msy::calcEquilibrium(const dvector& fe)
{
	/*
		This is the main engine behind the Msy class.
		Using the private member variables, calculate the equilibrium yield
		for a given fishing mortality rate vector (fe). Also, compute the
		derivative information such that this info can be used to compute 
		the corresponding MSY-based reference points (Fmsy, MSY, Bmsy).
		
		There are two different survivorship calculations inside this routine:
		lz is the survivorship of the total numbers, and lw is the survivorship
		upto the time of spawning.  The spawning biomass per recruit (phif) depends
		on the survivorship up to the time of spawning.

		 July 31, 2013.
		Created by Steve Martell.  Copied from calc_equilibrium
		- This routine has been modified from the previous to include sex-spcific 
		  information in the age-schedules. The basic yield equation for each
		  fleet k now looks like this:
		  Ye_k = Re * sum_h Fe_k phiq_{k,h}

		- Add age-specific natural mortality to survivorship calculations.
		
		- The major modification here is to wrap the original loops that relied on
		  selectivity, natural mortality, weight-at-age inside a loop over sex. 
		
		- [] Bug, with ngear > 1 get much lower Re values than if ngear ==1?


	*/
	
	// Indexes for dimensions
	int h,i,j,k;
	int sage,nage,ngear,ngrp;
	sage  = m_dWa.colmin();
	nage  = m_dWa.colmax();
	ngear = m_d3_V(1).rowmax();
	ngrp  = m_dWa.rowmax();
	
	m_lz.allocate(1,ngrp,sage,nage);

	double      ro = m_ro;
	double   kappa = 4.0*m_h/(1.0-m_h);  // Beverton-Holt model
	double     km1 = kappa-1.0;
	double    phif = 0;
	double   phisf = 0;
	double   phif2 = 0;
	
	// Survivorship for fished conditions.
	dvector   lz(sage,nage);
	dvector   lw(sage,nage);
	dvector   za(sage,nage);
	dvector  pza(sage,nage);
	dvector   sa(sage,nage);
	dvector  psa(sage,nage);
	dvector   oa(sage,nage);
	dvector  poa(sage,nage);

	dmatrix lz_m(1,ngrp,sage,nage);
	dmatrix lw_m(1,ngrp,sage,nage);
	dmatrix za_m(1,ngrp,sage,nage);
	dmatrix oa_m(1,ngrp,sage,nage);
	dmatrix sa_m(1,ngrp,sage,nage);

	dmatrix   qa(1,ngear,sage,nage);
	dmatrix  dlz(1,ngear,sage,nage);
	dmatrix d2lz(1,ngear,sage,nage);
	dmatrix  dlw(1,ngear,sage,nage);
	dmatrix d2lw(1,ngear,sage,nage);
	dlz.initialize();
	d2lz.initialize();
	dlw.initialize();
	d2lw.initialize();

	d3_array   qa_m(1,ngrp,1,ngear,sage,nage);
	d3_array  dlz_m(1,ngrp,1,ngear,sage,nage);
	d3_array d2lz_m(1,ngrp,1,ngear,sage,nage);
	d3_array  dlw_m(1,ngrp,1,ngear,sage,nage);
	d3_array d2lw_m(1,ngrp,1,ngear,sage,nage);

	for( h = 1; h <= ngrp; h++ )
	{
		// dvector za  = m_dMt(h) + fe * m_d3_V;
		za = m_dM(h);
		for( k = 1; k <= ngear; k++ )
		{
			za  = za + fe(k) * m_d3_V(h)(k);
		}
		za_m(h) = za;

		pza     = m_rho*za;
		sa      = exp(-za);
		sa_m(h) = sa;
		psa     = exp(-pza);
		oa      = 1.-sa;
		oa_m(h) = oa;
		poa     = 1.-elem_prod(sa,psa);
		for(k=1;k<=ngear;k++)
		{
			qa(k)        = elem_div(elem_prod(elem_prod(m_d3_V(h)(k),m_dWa(h)),oa),za);
			dlz(k,sage)  = 0;
			dlz_m(h,k,sage) = 0;
			d2lz(k,sage) = 0;
			d2lz_m(h,k,sage) = 0;
			dlw(k,sage)  = -psa(sage)*m_rho*m_d3_V(h)(k)(sage);
			dlw_m(h,k,sage) = 0;
			d2lw(k,sage) =  psa(sage)*square(m_rho)*square(m_d3_V(h)(k)(sage));
			dlw_m(h,k,sage) = 0;
			qa_m(h)(k)   = qa(k);
		}
		// cout<<"ngear "<<ngear<<endl;
		lz(sage) = 1.0/ngrp;
		lw(sage) = 1.0/ngrp * psa(sage);
		for(j=sage+1; j<=nage; j++)
		{
			lz(j)   = lz(j-1) * sa(j-1);
			lw(j)   = lz(j)   * psa(j);
			if( j==nage )
			{
				lz(j) = lz(j)/oa(j);
				lw(j) = lz(j-1)*sa(j-1)*psa(j)/oa(j);
			}
			
			for(k=1; k<=ngear; k++)
			{
				// derivatives for survivorship
				dlz(k)(j)  = sa(j-1) * ( dlz(k)(j-1) - lz(j-1)*m_d3_V(h)(k)(j-1) );
				d2lz(k)(j) = sa(j-1) * ( d2lz(k)(j-1)+ lz(j-1)*m_d3_V(h)(k)(j-1)*m_d3_V(h)(k)(j-1) );
				
				// derivatives for spawning survivorship
				dlw(k)(j)  = -lz(j)*m_rho*m_d3_V(h)(k)(j)*psa(j); 
				d2lw(k)(j) =  lz(j)*square(m_rho)*square(m_d3_V(h)(k)(j))*psa(j);
				
				if( j==nage ) // + group derivatives
				{
					dlz(k)(j)  = dlz(k)(j)/oa(j) - lz(j-1)*sa(j-1)*m_d3_V(h)(k)(j)*sa(j)/square(oa(j));
					
					dlw(k)(j)  = -lz(j-1)*sa(j-1)*m_rho*m_d3_V(h)(k)(j)/oa(j)
								- lz(j-1)*psa(j)*m_d3_V(h)(k)(j)*sa(j)/square(oa(j));
					
					double V1  = m_d3_V(h)(k)(j-1);
					double V2  = m_d3_V(h)(k)(j);
					double oa2 = oa(j)*oa(j);
					
					d2lz(k)(j) = d2lz(k)(j)/oa(j) 
								+ 2*lz(j-1)*V1*sa(j-1)*V2*sa(j)/oa2
								+ 2*lz(j-1)*sa(j-1)*V2*V2*sa(j)*sa(j)/(oa(j)*oa2)
								+ lz(j-1)*sa(j-1)*V2*V2*sa(j)/oa2;
					
					d2lw(k)(j) = lz(j-1)*square(m_rho)*square(V2)*psa(j)/oa(j)
								+ 2*lz(j-1)*m_rho*square(V2)*psa(j)*sa(j)/oa2
								+ 2*lz(j-1)*psa(j)*square(V2)*square(sa(j))/(oa(j)*oa2)
								+ lz(j-1)*psa(j)*square(V2)*sa(j)/oa2;
				}
				dlz_m(h,k,j)  =  dlz(k)(j);
				d2lz_m(h,k,j) = d2lz(k)(j);
				dlw_m(h,k,j)  =  dlw(k)(j);
				d2lw_m(h,k,j) = d2lw(k)(j);

			}// gear		
		}// age
		lz_m(h) = lz;
		lw_m(h) = lw;
		phisf += lz * m_dFa(h);
		phif  += lz * m_dFa(h);
		// cout<<h<< " phif "<<phif<<endl;
	}// ngrp
	phif2 = phif*phif;
	m_lz  = lz_m;
	
	// Incidence functions and associated derivatives
	dvector  dphif(1,ngear);   dphif.initialize();
	dvector d2phif(1,ngear);  d2phif.initialize();
	dvector   phiq(1,ngear);    phiq.initialize();
	dvector  dphiq(1,ngear);   dphiq.initialize();
	dvector d2phiq(1,ngear);  d2phiq.initialize();
	dvector    dre(1,ngear);     dre.initialize();
	dvector   d2re(1,ngear);    d2re.initialize();
	dvector     t1(sage,nage);    t1.initialize();
	
	for( h = 1; h <= ngrp; h++ )
	{	
		for(k=1; k<=ngear; k++)
		{
			// cout<<"dlz_m"<<endl;
			// cout<<dlz_m(h)(k) * m_dFa(h)<<endl;
			dphif(k)  += dlz_m(h)(k)  * m_dFa(h);
			d2phif(k) += d2lz_m(h)(k) * m_dFa(h);
			//dphif(k)   += dlw_m(h)(k)  * m_dFa(h);
			//d2phif(k)  += d2lw_m(h)(k) * m_dFa(h);
			
			// per recruit yield
			phiq(k)   +=  lz_m(h) * qa_m(h)(k);
			if(ngear>=1)  // was (if ngear==1), changed during debugging of nfleet>1
			{
				// dphiq = wa*oa*va*dlz/za + lz*wa*va^2*sa/za - lz*wa*va^2*oa/za^2
				t1 = elem_div(elem_prod(elem_prod(lz_m(h),m_dWa(h)),square(m_d3_V(h)(k))),za_m(h));
			}
			else
			{
				// dphiq = wa*oa*va*dlz/za + lz*wa*va*sa/za - lz*wa*va*oa/za^2
				t1 = elem_div(elem_prod(elem_prod(lz_m(h),m_dWa(h)),m_d3_V(h)(k)),za_m(h));
			}
			dvector t0 = elem_div(oa_m(h),za_m(h));
			dvector t3 = sa_m(h)-t0;
			dphiq(k)  += qa_m(h)(k)*dlz_m(h)(k) + t1 * t3;
			
			// 2nd derivative for per recruit yield (nasty)
			dvector t2  = 2. * dlz_m(h)(k);
			dvector V2  = elem_prod(m_d3_V(h)(k),m_d3_V(h)(k));
			dvector t5  = elem_div(elem_prod(m_dWa(h),V2),za_m(h));
			dvector t7  = elem_div(m_d3_V(h)(k),za_m(h));
			dvector t9  = elem_prod(t5,sa_m(h));
			dvector t11 = elem_prod(t5,t0);
			dvector t13 = elem_prod(lz_m(h),t5);
			dvector t14 = elem_prod(m_d3_V(h)(k),sa_m(h));
			dvector t15 = elem_prod(t7,sa_m(h));
			dvector t17 = elem_prod(m_d3_V(h)(k),t0);
			dvector t18 = elem_div(t17,za_m(h));
			d2phiq(k)  += d2lz_m(h)(k)*qa_m(h)(k) 
						+ t2*t9 - t2*t11 - t13*t14 -2.*t13*t15 
						+ 2.*t13*t18;
		}   // gear
	}  // ngrp
			
	for( k = 1; k <= ngear; k++ )
	{
		// 1st & 2nd partial derivatives for recruitment
		dre(k)      = ro*m_phie*dphif(k)/(phif2*km1);
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

	re   = ro*(kappa-m_phie/phif) / km1;
	ye   = re*elem_prod(fe,phiq);
	// dye  = re*phiq + elem_prod(fe,elem_prod(phiq,dre)) + elem_prod(fe,re*dphiq);
	dye  = re*phiq + elem_prod(fe,phiq)*dre + (fe*re)*dphiq;

	cout<<"Re     = "<<re<<endl;
	cout<<"phiq   = "<<phiq<<endl;
	cout<<"fe     = "<<fe<<endl;
	cout<<"dre    = "<<dre<<endl;
	cout<<"dphiq  = "<<dphiq<<endl;
	cout<<"dye    = "<<dye<<endl;
	
	// Caclculate Jacobian matrix (2nd derivative of the catch equation)
	for(j=1; j<=ngear; j++)
	{
		for(k=1; k<=ngear; k++)
		{
			d2ye(k)(j) = fe(j)*phiq(j)*d2re(k) + 2.*fe(j)*dre(k)*dphiq(k) + fe(j)*re*d2phiq(k);
			if(k == j)
			{
				d2ye(j)(k) += 2.*dre(j)*phiq(j)+2.*re*dphiq(j);
			}
		} 
	}
	// Newton-Raphson step.
	invJ   = -inv(d2ye);
	fstp   = invJ * dye;  
	
	// Set private members
	m_p    = fstp;
	m_ye   = ye;
	// cout<<"Re =    "<<re<<endl;
	m_re   = re;
	m_be   = re*phif;
	m_bi   = re*phisf;
	m_phif = phif;
	m_phiq = phiq;
	m_spr  = phif/m_phie;
	m_dYe  = sum(dye);
	m_d2Ye = sum(diagonal(d2ye));
	m_g    = diagonal(d2ye);		//Gradient vector
	m_f    = dye;   				//Value of the function to minimize
	
	// cout<<"mean za = "<<mean(za_m)<<endl;

	// numerical test for dre passed for 2+ fleets.
	// numerical test for dphif failed.
    // numerical test for dphiq passed for 1 fleet, fails on 2+ fleets.


	cout<<"fe "<<setw(8)<<fe<<" re "<<setw(8)<<re<<" ye "<<setw(8)<<ye
	    <<" dye "<<setw(8)<<dye<<endl;
	// cout<<"dlz_m"<<endl<<dlz_m<<endl;

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
		lx(i) = exp( -m_M*(i-sage) -m_rho*m_M );
		if(i==nage) 
			lx(i) /= 1.0 - exp( -m_M );
	}
	m_phie = lx   * m_fa;
	m_bo   = m_ro * m_phie;
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
		lx(i) = exp( -m_M*(i-sage) -m_rho*m_M );
		if(i==nage) 
			lx(i) /= 1.0 - exp( -m_M );
	}
	m_phie = lx   * m_fa;
	m_bo   = m_ro * m_phie;
}

void Msy::calc_phie(const dmatrix& _m, const dmatrix& _fa)
{
	int i,j,sage,nage;
	m_dM    = _m;
	m_dFa   = _fa;
	int r1  = m_dM.rowmin();
	int r2  = m_dM.rowmax();
	sage    = m_dFa.colmin();
	nage    = m_dFa.colmax();
	dvector lx(sage,nage);
	m_phie  = 0;
	for(i=m_dM.rowmin(); i<=m_dM.rowmax(); i++)
	{
		lx.initialize();
		for( j = sage; j <= nage; j++ )
		{
			lx(j) = exp( -m_dM(i,j)*(j-sage) - m_rho*m_dM(i,j) );
			if( j==nage )
				lx(j) /= 1.0 - exp( -m_dM(i,j) );
		}
		
		m_phie += 1./(r2-r1+1)*lx * m_dFa(i);
	}

	// [] TODO: check the discrepency between bo here and bo in calcStockRecruitment.
	m_bo = m_ro * m_phie;	
	
}

void Msy::calc_bo(double& _m, dvector& _fa)
{
	calc_phie(_m,_fa);
}

void Msy::calc_bo(const dmatrix& _m, const dmatrix& _fa)
{
	calc_phie(_m,_fa);
}

void Msy::print()
{
	cout<<"|------------------------------------------|" <<endl;
	cout<<"| Bo   = "<<setw(10)<<m_bo                    <<endl;
	cout<<"| Re   = "<<setw(10)<<m_re                    <<endl;
	cout<<"| Rmsy = "<<setw(10)<<m_rmsy                  <<endl;
	cout<<"| Bmsy = "<<setw(10)<<m_bmsy                  <<endl;
	cout<<"| Fmsy = "<<setw(10)<<m_fmsy                  <<endl;
	cout<<"| SPR  = "<<setw(10)<<m_spr                   <<endl;
	cout<<"| BPR  = "<<setw(10)<<m_phie                  <<endl;
	cout<<"| bpr  = "<<setw(10)<<m_phif                  <<endl;
	cout<<"| MSY  = "<<setw(10)<<m_msy                   <<endl;
	cout<<"| dYe  = "<<setw(10)<<m_f                     <<endl;
	cout<<"| dYes = "<<setw(10)<<sum(m_f)                <<endl;
	cout<<"|------------------------------------------|" <<endl;
}

