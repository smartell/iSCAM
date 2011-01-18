/*
   Routines for solving the Baranov catch equation.
*/                                                 

#include<admodel.h>

	double get_ft(const double& ct, const double& m, const dvector& va, const dvector& ba);
	

/** get_ft
  Solving the baranov catch equation using Newtons method
  Catch is based on weight.
*/	
double get_ft(const double& ct, const double& m, const dvector& va, const dvector& ba)
{
	double ft;
	//initial guess for ft
	ft=ct/(va*(ba*exp(-m/2.)));
	
	for(int i=1;i<=50;i++)
	{
		dvector f = ft*va;
		dvector z = m+f;
		dvector s = exp(-z);
		dvector o = (1.-s);
		
		dvector t1 = elem_div(f,z);
		dvector t2 = elem_prod(t1,o);
		dvector t3 = elem_prod(o,ba);
		//predicted catch
		double pct = t2*ba;
		
		//derivative of catch wrt ft
		double dct = sum(
			elem_div(t3,z) 
			- elem_div(elem_prod(f,t3),square(z))
			+ elem_prod(elem_prod(t1,s),ba));
		
		ft -= (pct-ct)/dct;  //newton step
		if(fabs(pct-ct)<1.e-4) break; //do not use for dvariables
	}
	
	return(ft);
}

/** get_ft
  Solving the baranov catch equation using Newtons method
  -this implementation is for multiple fleets where the catch
   is based on weight.
*/
dvector get_ft(dvector& ct,const double& m, const dmatrix& V,const dvector& ba)
{
	/*  ARGUMENTS:
	   ct is the observed catch for each fleet
	   m is the natural mortality rate
	   va is the vulnerability row fleet, col age va(1,ngear,1,nage)
	   ba is the start of year biomass at age
	*/
	
	int i,a,A;
	double minsurv = 0.05;
	int ng=size_count(ct);	//number of gears
	a=ba.indexmin(); A=ba.indexmax();
	dvector ft(1,ng); ft=0;
	dvector ctmp(1,ng);
	dvector ct_hat(1,ng);	//predicted catch
	dvector dct_hat(1,ng);	//derivative
	dvector zage(a,A);
	dvector sage(a,A);
	dvector ominus(a,A);
	dmatrix F(1,ng,a,A);
	
	
	//initial guess for ft=ct/(0.98 Bt);
	ctmp=ct;
	
	for(i=1;i<=ng;i++)
	{   
		ft(i) = ctmp(i)/(0.98*ba*V(i)*exp(-0.5*m));
		if(1.-ft(i)<minsurv)
		{
			ft(i)=1.-minsurv;
			ctmp=ft(i)*ba*V(i)*exp(-0.5*m);
		}
	}
	ct=ctmp;	//don't do this for the differentiable version.
	
	//now solve baranov catch equation iteratively.
	for(int iter=1; iter<=17; iter++)
	{
		for(i=1;i<=ng;i++)F(i)=ft(i)*V(i);
		zage=m+colsum(F);
		sage=mfexp(-zage);
		ominus=(1.-sage);
		
		for(i=1;i<=ng;i++)
		{   
			dvector t1 = elem_div(F(i),zage);
			dvector t2 = elem_prod(t1,ominus);
			dvector t3 = elem_prod(ominus,ba);
			
			ct_hat(i) = t2*ba;
			
			dct_hat(i) = sum(
						elem_div(t3,zage)
						- elem_div(elem_prod(F(i),t3),square(zage))
						+ elem_prod(elem_prod(t1,sage),ba));
				 
		    ft(i) -= (ct_hat(i)-ctmp(i))/dct_hat(i);
		}
		//cout<<iter<<"\t"<<ft<<"\t"<<ct_hat-ct<<endl;
	}
	//cout<<ft<<"\t\t"<<ct<<"\t\t"<<ctmp<<endl;
	
	return(ft);
}  
