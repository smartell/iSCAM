
#include "baranov.h"

BaranovCatchEquation::BaranovCatchEquation()
{
	// Constructor
}
BaranovCatchEquation::~BaranovCatchEquation()
{
	// Destructor
}



/** \brief Baranov catch equation solution for 1 or more fleets.
	
		The following function solves the Baranov catch equation for multiple fleets using
		a Newton-Raphson alogrithm to find a vector of fishing mortlity rates (ft) that 
		predicts the total catch for each fleet.  The Jacobian matrix is computed and 
		subsequently inverted to compute the Newton step for each fishing rate.
		
		BARANOV CATCH EQUATION:
		\f$ C = \frac{FN(1-exp(-Z))}{Z} \f$

	\author Martell IPHC
	\date 2012-08-05
	\param  ct a vector of observed catches, 1 element for each gear.
	\param  m instananeous natural mortality rate.
	\param  V a matrix of selectivities (row for each gear, col for each age)
	\param  na a vector of numbers-at-age at the start of each year
	
	\return Returns a vector of instantaneous fishing mortality rates.
	\sa
**/
dvector BaranovCatchEquation::getFishingMortality(const dvector &ct, const double &m, const dmatrix &V, const dvector &na)
{
	
	int i,j,its;
	int ngear = V.rowmax()-V.rowmin()+1;

	dvector   ft(1,ngear);
	dvector chat(1,ngear);
	dvector ctmp(1,ngear);
	dvector   fx(1,ngear);
	dmatrix    J(1,ngear,1,ngear);
	dmatrix invJ(1,ngear,1,ngear);
	dvector   ba(na.indexmin(),na.indexmax());
	dmatrix    F(1,ngear,V.colmin(),V.colmax());
	
	
	
	// Initial guess for fishing mortality rates;
	ba        = na; //elem_prod(na,wa);
	double bt = sum(na * exp(-0.5*m));
	ft        = ct / bt;
		
	// Iterative soln for catch equation using Newton-Raphson
	for(its=1; its<=MAXITS; its++)
	{
		for(i=1;i<=ngear;i++) F(i) =ft(i)*V(i);
		
		dvector za = m + colsum(F);
		dvector sa = exp(-za);
		dvector oa = (1.-sa);
		
		for(i=1;i<=ngear;i++)
		{
			for(j=1;j<=ngear;j++)
			{
				if(i==j)
				{
					dvector k1   =  elem_prod(ba,elem_div(V(i),za));
					dvector k2   =  elem_prod(k1,V(i));
					dvector k3   =  elem_div(k2,za);
					double dCdF  = -(k1*oa) - ft(i)*(k2*sa) + ft(i)*(k3*oa);
					J(i)(j)      = dCdF;
					chat(i)      = (ft(i)*k1) * oa;
				}
				else
				{
					dvector t1   = elem_div(elem_prod(ft(i)*ba,V(i)),za);
					dvector t2   = elem_prod(t1,V(j));
					dvector t3   = elem_div(t2,za);
					double dCdF  = -(t2*sa) + (t3*oa);
					J(j)(i)      = dCdF;
				}
			}	
		}
		fx   = ct - chat;
		//The following couts were used to debug the transpose error in the Jacobian.
		
		//cout<<"fx = "<<fx<<endl;
		// cout<<"Jacobian\t"<<"its = "<<its<<"\n"<<J<<endl;
		invJ = -inv(J);
		ft  += fx*invJ;
		
		if( norm(fx) < TOL || max(ft) > MAXF ) break;
	}
	
	for(i=1;i<=ngear;i++) if(ft(i)>MAXF) ft(i) = MAXF;
	
	return (ft);
}



/** \brief Get instantaneous fishing mortality rates based on Baranov catch equation.
	
		The following function solves the Baranov catch equation for multiple fleets using
		a Newton-Raphson alogrithm to find a vector of fishing mortlity rates (ft) that 
		predicts the total catch for each fleet.  The Jacobian matrix is computed and 
		subsequently inverted to compute the Newton step for each fishing rate.
	
	\author Martell IPHC
	\date 2012-08-05
	\param  ct a vector of observed catches based on weight, 1 element for each gear.
	\param  m age-independent instananeous natural mortality rate.
	\param  V a matrix of selectivities (row for each gear, col for each age)
	\param  na a vector of numbers-at-age at the start of each year
	\param  wa a vector of mean weight-at-age.
	\return Returns a vector of instantaneous fishing mortality rates.
	\sa
**/
dvector BaranovCatchEquation::getFishingMortality(const dvector &ct, const double &m, const dmatrix &V, const dvector &na, const dvector &wa)
{
	
	int i,j,its;
	int ngear = V.rowmax()-V.rowmin()+1;

	dvector   ft(1,ngear);
	dvector chat(1,ngear);
	dvector ctmp(1,ngear);
	dvector   fx(1,ngear);
	dmatrix    J(1,ngear,1,ngear);
	dmatrix invJ(1,ngear,1,ngear);
	dvector   ba(na.indexmin(),na.indexmax());
	dmatrix    F(1,ngear,V.colmin(),V.colmax());
	
	
	
	// Initial guess for fishing mortality rates;
	ba        = elem_prod(na,wa);
	double bt = (na * exp(-0.5*m)) * wa;
	ft        = ct / bt;
		
	// Iterative soln for catch equation using Newton-Raphson
	for(its=1; its<=MAXITS; its++)
	{
		for(i=1;i<=ngear;i++) F(i) =ft(i)*V(i);
		
		dvector za = m + colsum(F);
		dvector sa = exp(-za);
		dvector oa = (1.-sa);
		
		for(i=1;i<=ngear;i++)
		{
			for(j=1;j<=ngear;j++)
			{
				if(i==j)
				{
					dvector k1   =  elem_prod(ba,elem_div(V(i),za));
					dvector k2   =  elem_prod(k1,V(i));
					dvector k3   =  elem_div(k2,za);
					double dCdF  = -(k1*oa) - ft(i)*(k2*sa) + ft(i)*(k3*oa);
					J(i)(j)      = dCdF;
					chat(i)      = (ft(i)*k1) * oa;
				}
				else
				{
					dvector t1   = elem_div(elem_prod(ft(i)*ba,V(i)),za);
					dvector t2   = elem_prod(t1,V(j));
					dvector t3   = elem_div(t2,za);
					double dCdF  = -(t2*sa) + (t3*oa);
					J(j)(i)      = dCdF;
				}
			}	
		}
		fx   = ct - chat;
		//The following couts were used to debug the transpose error in the Jacobian.
		
		//cout<<"fx = "<<fx<<endl;
		//cout<<"Jacobian\t"<<"its = "<<its<<"\n"<<J<<endl;
		invJ = -inv(J);
		ft  += fx*invJ;
		
		if( norm(fx) < TOL || max(ft) > MAXF ) break;
	}
	
	for(i=1;i<=ngear;i++) if(ft(i)>MAXF) ft(i) = MAXF;
	
	return (ft);
}


// Vector of fishing mortality rate for catch based on numbers with age-dependent M.
/** \brief Get instantaneous fishing mortality rates based on Baranov catch equation.
	
		The following function solves the Baranov catch equation for multiple fleets using
		a Newton-Raphson alogrithm to find a vector of fishing mortlity rates (ft) that 
		predicts the total catch for each fleet.  The Jacobian matrix is computed and 
		subsequently inverted to compute the Newton step for each fishing rate.
	
	\author Martell IPHC
	\date 2012-08-05
	\param  ct a vector of observed catches based on numbers, 1 element for each gear.
	\param  ma age-dependent instananeous natural mortality rate.
	\param  V a matrix of selectivities (row for each gear, col for each age)
	\param  na a vector of numbers-at-age at the start of each year
	
	\return Returns a vector of instantaneous fishing mortality rates.
	\sa
**/

dvector BaranovCatchEquation::getFishingMortality(const dvector &ct, const dvector &ma, const dmatrix &V, const dvector &na)
{
	
	int i,j,its;
	int ngear = V.rowmax()-V.rowmin()+1;

	dvector   ft(1,ngear);
	dvector chat(1,ngear);
	dvector ctmp(1,ngear);
	dvector   fx(1,ngear);
	dmatrix    J(1,ngear,1,ngear);
	dmatrix invJ(1,ngear,1,ngear);
	dvector   ba(na.indexmin(),na.indexmax());
	dmatrix    F(1,ngear,V.colmin(),V.colmax());
	
	
	
	// Initial guess for fishing mortality rates;
	ba        = na; //elem_prod(na,wa);
	double bt = sum(elem_prod(na , exp(-0.5*ma)));
	ft        = ct / bt;
		
	// Iterative soln for catch equation using Newton-Raphson
	for(its=1; its<=MAXITS; its++)
	{
		for(i=1;i<=ngear;i++) F(i) =ft(i)*V(i);
		
		dvector za = ma + colsum(F);
		dvector sa = exp(-za);
		dvector oa = (1.-sa);
		
		for(i=1;i<=ngear;i++)
		{
			for(j=1;j<=ngear;j++)
			{
				if(i==j)
				{
					dvector k1   =  elem_prod(ba,elem_div(V(i),za));
					dvector k2   =  elem_prod(k1,V(i));
					dvector k3   =  elem_div(k2,za);
					double dCdF  = -(k1*oa) - ft(i)*(k2*sa) + ft(i)*(k3*oa);
					J(i)(j)      = dCdF;
					chat(i)      = (ft(i)*k1) * oa;
				}
				else
				{
					dvector t1   = elem_div(elem_prod(ft(i)*ba,V(i)),za);
					dvector t2   = elem_prod(t1,V(j));
					dvector t3   = elem_div(t2,za);
					double dCdF  = -(t2*sa) + (t3*oa);
					J(j)(i)      = dCdF;
				}
			}	
		}
		fx   = ct - chat;
		//The following couts were used to debug the transpose error in the Jacobian.
		
		//cout<<"fx = "<<fx<<endl;
		//cout<<"Jacobian\t"<<"its = "<<its<<"\n"<<J<<endl;
		invJ = -inv(J);
		ft  += fx*invJ;
		
		if( norm(fx) < TOL || max(ft) > MAXF ) break;
	}
	
	for(i=1;i<=ngear;i++) if(ft(i)>MAXF) ft(i) = MAXF;
	
	return (ft);
}



// Vector of fishing mortality rate for catch based on weight with age-dependent M.
/** \brief Get instantaneous fishing mortality rates based on Baranov catch equation.
	
		The following function solves the Baranov catch equation for multiple fleets using
		a Newton-Raphson alogrithm to find a vector of fishing mortlity rates (ft) that 
		predicts the total catch for each fleet.  The Jacobian matrix is computed and 
		subsequently inverted to compute the Newton step for each fishing rate.
	
	\author Martell IPHC
	\date 2012-08-05
	\param  ct a vector of observed catches based on weight, 1 element for each gear.
	\param  ma age-dependent instananeous natural mortality rate.
	\param  V a matrix of selectivities (row for each gear, col for each age)
	\param  na a vector of numbers-at-age at the start of each year
	\param  wa a vector of mean weight-at-age.
	\return Returns a vector of instantaneous fishing mortality rates.
	\sa
**/
dvector BaranovCatchEquation::getFishingMortality(const dvector &ct, const dvector &ma, const dmatrix &V, const dvector &na, const dvector &wa)
{
	
	int i,j,its;
	int ngear = V.rowmax()-V.rowmin()+1;

	dvector   ft(1,ngear);
	dvector chat(1,ngear);
	dvector ctmp(1,ngear);
	dvector   fx(1,ngear);
	dmatrix    J(1,ngear,1,ngear);
	dmatrix invJ(1,ngear,1,ngear);
	dvector   ba(na.indexmin(),na.indexmax());
	dmatrix    F(1,ngear,V.colmin(),V.colmax());
	
	
	
	// Initial guess for fishing mortality rates;
	ba        = elem_prod(na,wa);
	double bt = elem_prod(na , exp(-0.5*ma)) * wa;
	ft        = ct / bt;
		
	// Iterative soln for catch equation using Newton-Raphson
	for(its=1; its<=MAXITS; its++)
	{
		for(i=1;i<=ngear;i++) F(i) =ft(i)*V(i);
		
		dvector za = ma + colsum(F);
		dvector sa = exp(-za);
		dvector oa = (1.-sa);
		
		for(i=1;i<=ngear;i++)
		{
			for(j=1;j<=ngear;j++)
			{
				if(i==j)
				{
					dvector k1   =  elem_prod(ba,elem_div(V(i),za));
					dvector k2   =  elem_prod(k1,V(i));
					dvector k3   =  elem_div(k2,za);
					double dCdF  = -(k1*oa) - ft(i)*(k2*sa) + ft(i)*(k3*oa);
					J(i)(j)      = dCdF;
					chat(i)      = (ft(i)*k1) * oa;
				}
				else
				{
					dvector t1   = elem_div(elem_prod(ft(i)*ba,V(i)),za);
					dvector t2   = elem_prod(t1,V(j));
					dvector t3   = elem_div(t2,za);
					double dCdF  = -(t2*sa) + (t3*oa);
					J(j)(i)      = dCdF;
				}
			}	
		}
		fx   = ct - chat;
		//The following couts were used to debug the transpose error in the Jacobian.
		
		// cout<<"fx = "<<fx<<"\t ft = "<<max(ft)<<endl;
		// cout<<"Jacobian\t"<<"its = "<<its<<"\n"<<J<<endl;
		invJ = -inv(J);
		ft  += fx*invJ;
		
		if( norm(fx) < TOL || max(ft) > MAXF ) break;
	}
	// | DO NOT DO THE FOLLOWING FOR DIFFERENTIABLE PROBLEMS.
	for(i=1;i<=ngear;i++) if(ft(i)>MAXF) ft(i) = MAXF;
	
	return (ft);
}



/** \brief Instantaneous fishing mortality rate for a single fleet.
	
		Solves the Baranov Catch equation for a single fleet.
	
	\author  Martell, IPHC
	\date Sept 3, 2013.
	\param  ct Total catch in weight
	\param  m  Natural mortality rate
	\param  va age-specific selectivity.
	\param  ba weight-at-age vector.
	\return description of return value
	\sa
**/
double BaranovCatchEquation::get_ft(const double& ct, const double& m, const dvector& va, const dvector& ba)
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
		//if(fabs(pct-ct)<1.e-9) break; //do not use for dvariables
	}
	
	return(ft);
}



/** \brief Get conditional instantaneous fishing mortality rate when sex ratio of catch
           is unknown.
	
		This algorithm adds an additional dimension of sex to the gear-age problem.
		Use this function in cases
	
	\author Steven Martell 
	\date Sept 10, 2013
	\param  ct vector of catches for each gear.
	\param  ma a matrix of instanaenous age-specific natural mortality rates.
	\param  _V a pointer to a d3_array for age-gear-sex-specific selectivity.
	\param  na a matrix of numbers-at-age.
	\return description of return value
	\sa
**/
dvector BaranovCatchEquation::getFishingMortality(const dvector &ct, const dmatrix &ma, const d3_array *_V, const dmatrix &na)
{

	int h,i,j,k,its;
	d3_array V;
	V.allocate(*_V); //V(1,nsex,1,ngear,1,nage)
	V = *_V;

	int ngear = size_count(ct);
	int nsex  = ma.rowmax() - ma.rowmin() + 1;

	dvector   ft(1,ngear);
	dvector chat(1,ngear);
	// dvector ctmp(1,ngear);
	dvector   fx(1,ngear);
	dmatrix    J(1,ngear,1,ngear);
	dmatrix invJ(1,ngear,1,ngear);
	dvector ntmp(1,ngear);
	dmatrix   ba(1,nsex,na.colmin(),na.colmax());
	d3_array   F(1,nsex,1,ngear,na.colmin(),na.colmax());
	m_hCt.allocate(1,nsex,1,ngear);
	m_hCt.initialize();
	
	
	// Initial guess for fishing mortality rates based on Pope's approximation;
	ntmp.initialize();
	for( k = 1; k <= ngear; k++ )
	{
		for( h = 1; h <= nsex; h++ )
		{
			ntmp(k) +=  elem_prod(na(h),V(h)(k)) * ( exp(-0.5*ma(h)) );
		}
	}
	ft = elem_div(ct,ntmp);
	ba = na; //elem_prod(na,wa);
		
	// Iterative soln for catch equation using Newton-Raphson
	for(its=1; its<=MAXITS; its++)
	{
		J.initialize();
		chat.initialize();
		for( h = 1; h <= nsex; h++ )
		{
			 
			for(i=1;i<=ngear;i++) F(h)(i) =ft(i)*V(h)(i);
			
			dvector za = ma(h) + colsum(F(h));
			dvector sa = exp(-za);
			dvector oa = (1.-sa);
			
			for(i=1;i<=ngear;i++)
			{
				for(j=1;j<=ngear;j++)
				{
					if(i==j)
					{
						dvector k1   =  elem_prod(ba(h),elem_div(V(h)(i),za));
						dvector k2   =  elem_prod(k1,V(h)(i));
						dvector k3   =  elem_div(k2,za);
						double dCdF  = -(k1*oa) - ft(i)*(k2*sa) + ft(i)*(k3*oa);
						J(i)(j)     += dCdF;
						chat(i)     += (ft(i)*k1) * oa;
						m_hCt(h)(i)  = (ft(i)*k1) * oa;
					}
					else
					{
						dvector t1   = elem_div(elem_prod(ft(i)*ba(h),V(h)(i)),za);
						dvector t2   = elem_prod(t1,V(h)(j));
						dvector t3   = elem_div(t2,za);
						double dCdF  = -(t2*sa) + (t3*oa);
						J(j)(i)     += dCdF;
					}
				}	
			}
		}  // end of sex
		fx   = ct - chat;
		//The following couts were used to debug the transpose error in the Jacobian.
		
		//cout<<"fx = "<<fx<<endl;
		//cout<<"Jacobian\t"<<"its = "<<its<<"\n"<<J<<endl;
		invJ = -inv(J);
		ft  += fx*invJ;
		
		if( norm(fx) < TOL || max(ft) > MAXF ) break;
	}
	
	for(i=1;i<=ngear;i++) if(ft(i)>MAXF) ft(i) = MAXF;
	
	return (ft);

}


/** \brief Get conditional instantaneous fishing mortality rate when sex ratio of catch
           is unknown.
	
		This algorithm adds an additional dimension of sex to the gear-age problem.
		Use this function in cases
	
	\author Steven Martell 
	\date Sept 10, 2013
	\param  ct vector of catches for each gear.
	\param  ma a matrix of instanaenous age-specific natural mortality rates.
	\param  _V a pointer to a d3_array for age-gear-sex-specific selectivity.
	\param  na a matrix of numbers-at-age.
	\param  wa a matrix of weight-at-age.
	\return description of return value
	\sa
**/
dvector BaranovCatchEquation::getFishingMortality(const dvector &ct, const dmatrix &ma, const d3_array *_V, const dmatrix &na, const dmatrix &wa)
{

	int h,i,j,k,its;
	d3_array V;
	V.allocate(*_V); //V(1,nsex,1,ngear,1,nage)
	V = *_V;

	int ngear = size_count(ct);
	int nsex  = ma.rowmax() - ma.rowmin() + 1;

	dvector   ft(1,ngear);
	dvector chat(1,ngear);
	// dvector ctmp(1,ngear);
	dvector   fx(1,ngear);
	dmatrix    J(1,ngear,1,ngear);
	dmatrix invJ(1,ngear,1,ngear);
	dvector ntmp(1,ngear);
	dmatrix   ba(1,nsex,na.colmin(),na.colmax());
	d3_array   F(1,nsex,1,ngear,na.colmin(),na.colmax());
	m_hCt.allocate(1,nsex,1,ngear);
	m_hCt.initialize();
	
	
	// Initial guess for fishing mortality rates based on Pope's approximation;
	ba = elem_prod(na,wa);
	ntmp.initialize();
	for( k = 1; k <= ngear; k++ )
	{
		for( h = 1; h <= nsex; h++ )
		{
			ntmp(k) +=  elem_prod(na(h),V(h)(k)) * elem_prod(wa(h), exp(-0.5*ma(h)) );
		}
	}
	ft = elem_div(ct,ntmp);
		
	// Iterative soln for catch equation using Newton-Raphson
	for(its=1; its<=MAXITS; its++)
	{
		J.initialize();
		chat.initialize();
		for( h = 1; h <= nsex; h++ )
		{
			 
			for(i=1;i<=ngear;i++) F(h)(i) =ft(i)*V(h)(i);
			
			dvector za = ma(h) + colsum(F(h));
			dvector sa = exp(-za);
			dvector oa = (1.-sa);
			
			for(i=1;i<=ngear;i++)
			{
				for(j=1;j<=ngear;j++)
				{
					if(i==j)
					{
						dvector k1   =  elem_prod(ba(h),elem_div(V(h)(i),za));
						dvector k2   =  elem_prod(k1,V(h)(i));
						dvector k3   =  elem_div(k2,za);
						double dCdF  = -(k1*oa) - ft(i)*(k2*sa) + ft(i)*(k3*oa);
						J(i)(j)     += dCdF;
						chat(i)     += (ft(i)*k1) * oa;
						cout<<"CT = "<<(ft(i)*k1) * oa<<endl;
						m_hCt(h)(i)  = (ft(i)*k1) * oa;
					}
					else
					{
						dvector t1   = elem_div(elem_prod(ft(i)*ba(h),V(h)(i)),za);
						dvector t2   = elem_prod(t1,V(h)(j));
						dvector t3   = elem_div(t2,za);
						double dCdF  = -(t2*sa) + (t3*oa);
						J(j)(i)     += dCdF;
					}
				}	
			}
		}  // end of sex
		fx   = ct - chat;
		//The following couts were used to debug the transpose error in the Jacobian.
		
		cout<<"its = "<<its-1<<" fx = "<<fx<<endl;
		//cout<<"Jacobian\t"<<"its = "<<its<<"\n"<<J<<endl;
		invJ = -inv(J);
		ft  += fx*invJ;
		
		if( norm(fx) < TOL || max(ft) > MAXF ) break;
	}
	
	for(i=1;i<=ngear;i++) if(ft(i)>MAXF) ft(i) = MAXF;
	
	return (ft);

}


/** \brief Get conditional instantaneous fishing mortality rate when sex ratio of catch
           is unknown.
	
		This algorithm adds an additional dimension of sex to the gear-age problem.
		Use this function in cases where catch is based on numbers.
	
	\author Steven Martell 
	\date Sept 10, 2013
	\param  ct vector of catches for each gear.
	\param  ma a matrix of instanaenous age-specific natural mortality rates.
	\param  p_V a pointer to a d3_array for age-gear-sex-specific selectivity.
	\param  na a matrix of numbers-at-age.
	\param  _hCt a matrix of catch by sex(row) by gear(col)
	\return description of return value
	\sa
**/
dvector BaranovCatchEquation::getFishingMortality(const dvector &ct, const dmatrix &ma, const d3_array *p_V, const dmatrix &na, dmatrix &_hCt)
{
	dvector ft = getFishingMortality(ct,ma,p_V,na);
	_hCt = m_hCt;
	return(ft);
}


/** \brief Get conditional instantaneous fishing mortality rate when sex ratio of catch
           is unknown.
	
		This algorithm adds an additional dimension of sex to the gear-age problem.
		Use this function in cases where catch is based on weight.
	
	\author Steven Martell 
	\date Sept 10, 2013
	\param  ct vector of catches for each gear.
	\param  ma a matrix of instanaenous age-specific natural mortality rates.
	\param  p_V a pointer to a d3_array for age-gear-sex-specific selectivity.
	\param  na a matrix of numbers-at-age.
	\param  wa a matrix of weight-at-age.
	\param  _hCt a matrix of catch by sex(row) by gear(col)
	\return description of return value
	\sa
**/
dvector BaranovCatchEquation::getFishingMortality(const dvector &ct, const dmatrix &ma, const d3_array *p_V, const dmatrix &na, const dmatrix &wa, dmatrix &_hCt)
{
	dvector ft = getFishingMortality(ct,ma,p_V,na,wa);
	_hCt = m_hCt;
	return(ft);
}









// /** get_ft
//   Solving the baranov catch equation using Newtons method
//   -this implementation is for multiple fleets where the catch
//    is based on weight.
// */
// dvector get_ft(dvector& ct,const double& m, const dmatrix& V,const dvector& ba)
// {
// 	/*  ARGUMENTS:
// 	   ct is the observed catch for each fleet
// 	   m is the natural mortality rate
// 	   va is the vulnerability row fleet, col age va(1,ngear,1,nage)
// 	   ba is the start of year biomass at age
// 	*/
	
// 	int i,a,A;
// 	double minsurv = 0.05;
// 	int ng=size_count(ct);	//number of gears
// 	a=ba.indexmin(); A=ba.indexmax();
// 	dvector ft(1,ng); ft=0;
// 	dvector ctmp(1,ng);
// 	dvector ct_hat(1,ng);	//predicted catch
// 	dvector dct_hat(1,ng);	//derivative
// 	dvector zage(a,A);
// 	dvector sage(a,A);
// 	dvector ominus(a,A);
// 	dmatrix F(1,ng,a,A);
	
	
// 	//initial guess for ft=ct/(0.98 Bt);
// 	ctmp=ct;
	
// 	for(i=1;i<=ng;i++)
// 	{   
// 		ft(i) = ctmp(i)/(0.98*ba*V(i)*exp(-0.5*m));
// 		if(1.-ft(i)<minsurv)
// 		{
// 			ft(i)=1.-minsurv;
// 			ctmp=ft(i)*ba*V(i)*exp(-0.5*m);
// 		}
// 	}
// 	ct=ctmp;	//don't do this for the differentiable version.
	
// 	//now solve baranov catch equation iteratively.
// 	for(int iter=1; iter<=17; iter++)
// 	{
// 		for(i=1;i<=ng;i++)F(i)=ft(i)*V(i);
// 		zage=m+colsum(F);
// 		sage=mfexp(-zage);
// 		ominus=(1.-sage);
		
// 		for(i=1;i<=ng;i++)
// 		{   
// 			dvector t1 = elem_div(F(i),zage);
// 			dvector t2 = elem_prod(t1,ominus);
// 			dvector t3 = elem_prod(ominus,ba);
			
// 			ct_hat(i) = t2*ba;
			
// 			dct_hat(i) = sum(
// 						elem_div(t3,zage)
// 						- elem_div(elem_prod(F(i),t3),square(zage))
// 						+ elem_prod(elem_prod(t1,sage),ba));
				 
// 		    ft(i) -= (ct_hat(i)-ctmp(i))/dct_hat(i);
// 		}
// 		//cout<<iter<<"\t"<<ft<<"\t"<<ct_hat-ct<<endl;
// 	}
// 	//cout<<ft<<"\t\t"<<ct<<"\t\t"<<ctmp<<endl;
	
// 	return(ft);
// }  

// /** get_ft
//   Solving the baranov catch equation using Newtons method
//   -this implementation is for multiple fleets where the catch
//    is based on weight.
//  in this case, the predicted catch is based on catch-at-numbers * wa
// */
// dvector get_ft(dvector& ct,const double& m, const dmatrix& V,const dvector& na, const dvector& wa)
// {
// 	  ARGUMENTS:
// 	   ct is the observed catch for each fleet
// 	   m is the natural mortality rate
// 	   va is the vulnerability row fleet, col age va(1,ngear,1,nage)
// 	   na is the start of year numbers at age
// 	   wa is the mean weight-at-age

// 	   THis function should be deprecated.
	
// 	cout<<"I'm in the Baranov equation for multiple fleets"<<endl;
// 	int i,a,A;
// 	double minsurv = 0.05;
// 	int ng=size_count(ct);	//number of gears
// 	a=na.indexmin(); A=na.indexmax();
// 	dvector ft(1,ng); ft=0;
// 	dvector ctmp(1,ng);
// 	dvector ct_hat(1,ng);	//predicted catch
// 	dvector dct_hat(1,ng);	//derivative
// 	dvector ba(a,A);		//biomass at age
// 	dvector ca(a,A);		//catch-at-age in numbers
// 	dvector zage(a,A);
// 	dvector sage(a,A);
// 	dvector ominus(a,A);
// 	dmatrix F(1,ng,a,A);
	
	
// 	//initial guess for ft=ct/(0.98 Bt);
// 	ctmp=ct;
// 	ba = elem_prod(na,wa);
// 	for(i=1;i<=ng;i++)
// 	{   
// 		ft(i) = ctmp(i)/(0.98*ba*V(i)*exp(-0.5*m));
// 		if(1.-ft(i)<minsurv)
// 		{
// 			ft(i)=1.-minsurv;
// 			ctmp(i)=ft(i)*ba*V(i)*exp(-0.5*m);
// 		}
// 	}
// 	//ct=ctmp;	//don't do this for the differentiable version.
	
// 	//now solve baranov catch equation iteratively.
// 	for(int iter=1; iter<=17; iter++)
// 	{
// 		for(i=1;i<=ng;i++)F(i)=ft(i)*V(i);
// 		zage=m+colsum(F);
// 		sage=mfexp(-zage);
// 		ominus=(1.-sage);
		
// 		for(i=1;i<=ng;i++)
// 		{   
// 			dvector t1 = elem_div(F(i),zage);
// 			dvector t2 = elem_prod(t1,ominus);
// 			dvector t3 = elem_prod(ominus,na);
// 			ca = elem_prod(t2,na);
			
// 			ct_hat(i) = ca*wa;
			
// 			/*dct_hat(i) = sum(
// 						elem_div(t3,zage)
// 						- elem_div(elem_prod(F(i),t3),square(zage))
// 						+ elem_prod(elem_prod(t1,sage),na));
// 			*/
// 			dvector t4 = ft(i)*square(V(i));
			
// 			dvector t5 = elem_div(elem_prod(elem_prod(V(i),ominus),na),zage)
// 				- elem_div(elem_prod(elem_prod(t4,ominus),na),square(zage))
// 				+ elem_div(elem_prod(elem_prod(t4,sage),na),zage);
// 			dct_hat(i) = t5*wa;
				 
// 		    ft(i) -= (ct_hat(i)-ctmp(i))/dct_hat(i);
// 		}
// 		//cout<<iter<<"\t"<<ft<<"\t"<<ct_hat-ctmp<<endl;
// 		//SJDM, this algorithm does converge niceley for multiple fleets
// 	}
// 	//cout<<ft<<"\t\t"<<ct<<"\t\t"<<ctmp<<endl;
// 	//cout<<ct<<endl;
// 	return(ft);
// }  

