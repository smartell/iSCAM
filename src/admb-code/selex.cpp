/** \brief  class Selex
	
	Selectivity function used in age and size structured models.  This is meant to be
	a single repository for various selectivity functions used in fisheries assessment
	models.
	
&copy; Copyright 2012 UBC Fisheries Centre - . All Rights Reserved.

	\author  Steven Martell
	\author $Martell$
	\date 2012-09-04
	\date $LastChangedDate$
	\version $Rev$
	\sa
**/


#include <admodel.h>
#include <contrib.h>
#include "Selex.h"


Selex::~Selex(){}

Selex::Selex()
{}


Selex::Selex(const dvector& _x)
: m_x(_x)
{
	cout<<"In the dvar_matrix constructor:"<<endl;
	
	
};

// void Selex::fill_selex_array(dvar_matrix& log_sel)
// {
// 	/*
// 		Fill in the selectivity matrix based on selType and selPars defined during
// 		instantiation of the class.
// 	*/


// 	int i,j;
// 	int r1 = log_sel.rowmin();
// 	int r2 = log_sel.rowmax();
// 	int c1 = log_sel.colmin();
// 	int c2 = log_sel.colmax();
// 	for( i = r1; i <= r2; i++ )
// 	{
// 		switch (m_selType)
// 		{
// 			default:
// 				log_sel(i) = 0;
// 				break;

// 			case 1:    // logistic
// 				log_sel(i) = log( logistic( m_dvar_selPar(1)) );
// 				break;

// 			case 2:
// 				for( j = c1; j < c2; j++ )
// 				{
// 					log_sel(i,j) = m_dvar_selPar(1,j-c1+1);
// 				}
// 				log_sel(i)(c2) = log_sel(i)(c2-1);
// 				break;
// 		}
// 	}

// 	cout<<m_selType<<endl;
// 	cout<<r1<<"\t"<<r2<<"\t"<<c1<<"\t"<<c2<<endl;
// 	cout<<m_dvar_selPar<<endl;
// }

// // ------------------------------------------------------------------------------------ //
// class Selex{
// private:
// 	dvariable m_mu;
// 	dvariable m_sd;
		
// public:
// 	~Selex() {}                                   // destructor
		
// 	Selex( const dvariable& mu = 0., const dvariable& sd = 1.0) // default constructor
// 	{
// 		m_mu = mu;
// 		m_sd = sd;
// 	}
	
// 	// Logistic function
// 	dvector     logistic( const dvector& x,const double& mu, const double& sd );
// 	dvar_vector logistic( const dvector& x,const dvariable& mu, const dvariable& sd );
	
	
// 	// Exponential logistic
// 	dvector     eplogis(const dvector& x, const double& x1, const double& x2, const double& gamma);
// 	dvar_vector eplogis(const dvar_vector& x, const dvariable& x1, const dvariable& x2, const dvariable& gamma);
	
	
// 	// Linear Interoplation
// 	dvector     linapprox(const dvector& x, const dvector& y, const dvector& xout);
// 	dvar_vector linapprox(const dvector& x, const dvar_vector& y, const dvector& xout);
	
// 	dvariable GetMu()  { return m_mu; }
// 	dvariable GetSd()  { return m_sd; }
// };
// ------------------------------------------------------------------------------------ //
// Logistic function                                                                    //
// ------------------------------------------------------------------------------------ //
dvar_vector Selex::logistic( const dvector& x,const dvariable& mu, const dvariable& sd )
{
	return 1./(1.+mfexp(-(x-mu)/sd) );
}


dvector Selex::logistic( const dvector& x,const double& mu, const double& sd )
{
	return 1./(1.+mfexp(-(x-mu)/sd) );
}
 
dvar_vector Selex::logistic(const dvar_vector& log_selpar)
{
	RETURN_ARRAYS_INCREMENT();
	dvariable mu = mfexp(log_selpar(1));
	dvariable sd = mfexp(log_selpar(2));
	RETURN_ARRAYS_DECREMENT();
	return 1./(1.+mfexp(-(m_x-mu)/sd) );
}

void Selex::logistic( const dmatrix& theta, const ivector& blk, dmatrix& log_sel)
{
	/* 
		Logistic selectivity for MSE operating model.
		Args: theta     (matrix of selectivity parameters)
		      blk       (time blocks, start year.)
		      log_sel   (matrix of selectivity coefficients)

	*/	

	int i,j;
	int cntr = 1;
	int r1 = log_sel.rowmin();
	int r2 = log_sel.rowmax();
	int c1 = log_sel.colmin();
	int c2 = log_sel.colmax();
	double p1,p2;
	double tiny = 1.e-30;

	log_sel.initialize();
	j = 0;
	for(i = r1; i<= r2; i++)
	{
		if( i == blk(cntr))
		{
			j++;
			if(cntr < blk.indexmax()) cntr++;
		}
		p1 = mfexp( theta(j,1) );
		p2 = mfexp( theta(j,2) );
		log_sel(i) = log( plogis(m_x, p1, p2) + tiny );
		log_sel(i)-= log( mean(mfexp(log_sel(i))) );
	}
}

void Selex::logistic( const dmatrix& theta, const ivector& blk, const dmatrix& len, dmatrix& log_sel)
{
	/*
		Logistic selectivity based on mean length of individuals at age.
		Args: theta   (matrix of selectivity parameters)
		      blk     (time blocks, start year)
		      log_sel (matrix of selectivity coefficients)
	*/
    cout<<"Selex::logistic HELP"<<endl;
    if(len.colmin() != log_sel.colmin() || len.colmax() != log_sel.colmax())
    {
    	cerr<<"Length & log_sel matrix dimensions do not match"<<endl;
    	ad_exit(1);
    }
	int i,j;
	int cntr = 1;
	int r1 = log_sel.rowmin();
	int r2 = log_sel.rowmax();
	int c1 = log_sel.colmin();
	int c2 = log_sel.colmax();
	double p1,p2;
	double tiny = 1.e-30;

	log_sel.initialize();
	j = 0;
	for(i = r1; i<= r2; i++)
	{
		cout<<len(i)<<endl;
		if( i == blk(cntr))
		{
			j++;
			if(cntr < blk.indexmax()) cntr++;
		}
		p1 = mfexp( theta(j,1) );
		p2 = mfexp( theta(j,2) );
		log_sel(i) = log( plogis(len(i), p1, p2) + tiny );
		log_sel(i)-= log( mean(mfexp(log_sel(i))) );
	}
}

// |---------------------------------------------------------------------------------|
// | Selectivity Coefficients
// |---------------------------------------------------------------------------------|
// |
void Selex::selcoeff( const dmatrix& theta, const ivector& blk, dmatrix& log_sel)
{
	/*
		Selectivity coefficients for each age class, where it is assumed that the
		last two age classes have the same selectivity coefficient.
	*/
	int i,j,k;
	int cntr = 1;
	int r1 = log_sel.rowmin();
	int r2 = log_sel.rowmax();
	int c1 = log_sel.colmin();
	int c2 = log_sel.colmax();

	log_sel.initialize();
	k = 0;
	for( i = r1; i <= r2; i++ )
	{
		if( i == blk(cntr))
		{
			k++;
			if(cntr < blk.indexmax()) cntr++;
		}
		for( j = c1; j < c2; j++ )
		{
			log_sel(i,j) = theta(k)(j-c1+1);
		}
		log_sel(i,c2) = log_sel(i,c2-1);
		log_sel(i)   -= log( mean(mfexp(log_sel(i))) );
	}

}

// ------------------------------------------------------------------------------------ //
// Exponential Logistic                                                                 //
// ------------------------------------------------------------------------------------ //
// dvar_vector Selex::eplogis( const dvector& x, const dvariable& x1, 
//                             const dvariable& x2, const dvariable& gamma )
// {
// 	//exponential logistic based on Grant Thompson (1994) Paper, CJFAS.
// 	/*
// 	A modified version of the exponential logistic presented in Thompson's 1994 paper in CJFAS
	
// 	Here the arguemnts x1 and x2 represnt the inflection points for the ascending 
// 	and desciending limb of the selectivity function.  Gamma is descending limb 
// 	parameter where gamma=0 is logistic, and gamma <1.0 dome=shaped and gamma==1 is undefined.
	
// 	*/
// 	RETURN_ARRAYS_INCREMENT();
// 	dvariable k1,k2,alpha,beta;
// 	dvariable t1 = 2.-4.*gamma+2.*gamma*gamma;
// 	dvariable t3 = 1.+2.*gamma-2*gamma*gamma;
// 	dvariable t5 = sqrt(1. + 4.*gamma - 4.*gamma*gamma);
	
// 	k1    = log(t1/(t3+t5));
// 	k2    = log(t1/(t3-t5));
// 	beta  = (k1*x2-x1*k2)/(k1-k2);
// 	alpha = k2/(x2-beta);
// 	dvar_vector sx = (1./(1.-gamma))*pow((1.-gamma)/gamma,gamma)*elem_div(exp(alpha*gamma*(beta-x)),1.+exp(alpha*(beta-x)));
	
// 	RETURN_ARRAYS_DECREMENT();
// 	return sx;
// }

// dvector Selex::eplogis( const dvector& x, const double& x1, 
//                         const double& x2, const double& gamma )
// {
// 	//exponential logistic based on Grant Thompson (1994) Paper, CJFAS.
// 	/*
// 	A modified version of the exponential logistic presented in Thompson's 1994 paper in CJFAS
	
// 	Here the arguemnts x1 and x2 represnt the inflection points for the ascending 
// 	and desciending limb of the selectivity function.  Gamma is descending limb 
// 	parameter where gamma=0 is logistic, and gamma <1.0 dome=shaped and gamma==1 is undefined.
	
// 	*/
// 	double k1,k2,alpha,beta;
// 	double t1 = 2.-4.*gamma+2.*gamma*gamma;
// 	double t3 = 1.+2.*gamma-2*gamma*gamma;
// 	double t5 = sqrt(1. + 4.*gamma - 4.*gamma*gamma);
	
// 	k1    = log(t1/(t3+t5));
// 	k2    = log(t1/(t3-t5));
// 	beta  = (k1*x2-x1*k2)/(k1-k2);
// 	alpha = k2/(x2-beta);
// 	dvector sx = (1./(1.-gamma))*pow((1.-gamma)/gamma,gamma)*elem_div(exp(alpha*gamma*(beta-x)),1.+exp(alpha*(beta-x)));
	
// 	return sx;
// }


// ------------------------------------------------------------------------------------ //
// Linear Interpolation using approx function from R libraries                          //
// ------------------------------------------------------------------------------------ //
static double approx1(const double& v, const dvector& x, const dvector& y)
{
    /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
    int i, j, ij;

    //if(!n) return R_NaN;

    i = x.indexmin();
    j = x.indexmax() - 1;

    /* handle out-of-domain points */
    if(v < x[i]) return min(y);
	if(v > x[j]) return max(y);

    /* find the correct interval by bisection */
    while(i < j - 1) 
	{ /* x[i] <= v <= x[j] */
		ij = (i + j)/2; /* i+1 <= ij <= j-1 */
		if(v < x[ij]) j = ij;
		else i = ij;
		/* still i < j */
    }
    /* provably have i == j-1 */

    /* interpolation */
    if(v == x[j]) return y[j];
    if(v == x[i]) return y[i];
    /* impossible: if(x[j] == x[i]) return y[i]; */

	/* linear */
	return y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i]));
}/* approx1() */

static dvariable approx1(const double& v, const dvector& x, const dvar_vector& y)
{
    /* Approximate  y(v),  given (x,y)[i], i = 0,..,n-1 */
    int i, j, ij;

    //if(!n) return R_NaN;

    i = x.indexmin();
    j = x.indexmax() - 1;

    /* handle out-of-domain points */
    if(v < x[i]) return min(y);
	if(v > x[j]) return max(y);

    /* find the correct interval by bisection */
    while(i < j - 1) 
	{ /* x[i] <= v <= x[j] */
		ij = (i + j)/2; /* i+1 <= ij <= j-1 */
		if(v < x[ij]) j = ij;
		else i = ij;
		/* still i < j */
    }
    /* provably have i == j-1 */

    /* interpolation */
    if(v == x[j]) return y[j];
    if(v == x[i]) return y[i];
    /* impossible: if(x[j] == x[i]) return y[i]; */

	/* linear */
	return y[i] + (y[j] - y[i]) * ((v - x[i])/(x[j] - x[i]));
}/* approx1() */


dvector Selex::linapprox(const dvector& x, const dvector& y, const dvector& xout)
{
	// Piece-wise linear approximation for n points in xout between min(x) and max(y):
	int i,j,k;
	int n = xout.indexmax() - xout.indexmin() + 1;
	double v;
	dvector yout(xout.indexmin(),xout.indexmax());
	
	for(k = xout.indexmin(); k <= xout.indexmax(); k++)
	{
		v       = xout(k);
		yout(k) = approx1(v,x,y);
	}
	return yout;
}

dvar_vector Selex::linapprox(const dvector& x, const dvar_vector& y, const dvector& xout)
{
	// Piece-wise linear approximation for n points in xout between min(x) and max(y):
	int i,j,k;
	int n = xout.indexmax() - xout.indexmin() + 1;
	double v;
	RETURN_ARRAYS_INCREMENT();
	dvar_vector yout(xout.indexmin(),xout.indexmax());
	for(k = xout.indexmin(); k <= xout.indexmax(); k++)
	{
		v       = xout(k);
		yout(k) = approx1(v,x,y);
	}
	RETURN_ARRAYS_DECREMENT();
	return yout;
}
