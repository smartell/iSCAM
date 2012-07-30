// refpoints.cpp
// Author: Steven Martell
// An algorithm for calculating MSY-based reference points.


#include <admodel.h>

typedef struct {
	double ro;       // Unfished recruits
	double kappa;    // Recruitment compensation of the stock-recruitment function
	double m;        // Instantaneous natural mortality rate
	dvector fe;      // Fishing mortality rates
	dvector wa;      // Average weight-at-age
	dvector fa;      // Average fecundity-at-age
	dmatrix V;       // Selectivity for each gear (rows) at age (col)
	} params;

void calcEquilibrium(params theta)
{
	cout<<"in calcEquilibrium"<<endl;
	cout<<theta.m<<endl;
	// Calculate the equilibrium values to determine MSY-based reference points.
	int i,j,k;
	int sage,nage,ngear;
	sage  = theta.wa.indexmin();
	nage  = theta.wa.indexmax();
	ngear = theta.V.rowmax();
	// Survivorship for unfished conditions.
	dvector lx(sage,nage);
	lx.initialize();
	lx(sage) = 1.0;
	
	for(i=sage+1; i<=nage; i++)
	{
		lx(i) = lx(i-1)*exp(-theta.m);
		if(i==nage) lx(i) /= (1.-exp(-theta.m));
	}
	
	double phie = lx * theta.fa;
	
	// Survivorship for fished conditions.
	dvector lz(sage,nage);
	dmatrix dlz(1,ngear,sage,nage);
	dmatrix d2lz(1,ngear,sage,nage);
	dmatrix V = theta.V;
	dlz.initialize();
	d2lz.initialize();
	dvector fe = theta.fe;
	dvector za = theta.m + fe * V;
	dvector sa = exp(-za);
	dvector oa = 1.-sa;
	dmatrix qa = theta.V;
	
	for(k=1;k<=ngear;k++)
	{
		qa(k)        = elem_div(elem_prod(elem_prod(theta.V(k),theta.wa),oa),za);
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
			dlz(k)(j)  = sa(j-1) * ( dlz(k)(j-1) - lz(j-1)*V(k)(j-1) );
			d2lz(k)(j) = sa(j-1) * ( d2lz(k)(j-1)+ lz(j-1)*V(k)(j-1)*V(k)(j-1) );
			
			if( j==nage ) // + group derivatives
			{
				dlz(k)(j)  = dlz(k)(j)/oa(j) - lz(j-1)*sa(j-1)*V(k)(j)*sa(j)/square(oa(j));
				
				double V1  = V(k)(j-1);
				double V2  = V(k)(j);
				double oa2 = oa(j)*oa(j);
				
				d2lz(k)(j) = d2lz(k)(j)/oa(j) + 2*lz(j-1)*V1*V1*sa(j-1)/oa(j)
							+ 2*lz(j-1)*V1*sa(j-1)*V2*sa(j)/oa2
							+ 2*lz(j-1)*sa(j-1)*V2*V2*sa(j)/(oa[i]*oa2)
							+ lz(j-1)*sa(j-1)*V2*V2*sa(j)/oa2;
			}
		}
		
	}
	
	// Incidence functions and associated derivatives
	double ro   = theta.ro;
	double km1  = theta.kappa-1;
	double    phif = lz * theta.fa;
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
		dphif(k)   = dlz(k)  * theta.fa;
		d2phif(k)  = d2lz(k) * theta.fa;
		
		phiq(k)    = lz * qa(k);
		
		dvector t1 = elem_div(elem_prod(elem_prod(lz,theta.wa),V(k)),za);
		dvector t3 = elem_div(sa-oa,za);
		dphiq(k)   = qa(k)*dlz(k) + t1 * t3;
		
		dvector t7  = elem_div(elem_prod(theta.wa,V(k)),za);
		dvector t9  = elem_prod(t7,oa);
		dvector t11 = elem_prod(elem_prod(t7,V(k)),sa);
		dvector t13 = elem_div(elem_prod(t9,V(k)),za);
		dvector t15 = elem_prod(t11,V(k));
		dvector t17 = elem_div(t15,za);
		dvector t19 = elem_div(t17,za);
		d2phiq(k)   = d2lz(k)*t9 + 2.*dlz(k)*t11 
					- 2.*dlz(k)*t13 - lz*t15 
					- 2.*lz*t17 + 2.*lz*t19;
		
		dre(k)      = ro/km1*phie/(phif2)*dphif(k);
		d2re(k)     = -2.*ro*phie*d2phif(k)/(phif2*phif*km1) 
					+ ro*phie*d2phif(k)/(phif2*km1);
	}
	
	
	// Equilibrium calculations
	double re;
	dvector   ye(1,ngear);
	dvector  dye(1,ngear);
	dmatrix d2ye(1,ngear,1,ngear);
	
	re   = ro*(theta.kappa-phie/phif)/km1;
	ye   = re*elem_prod(fe,phiq);
	dye  = re*phiq + elem_prod(fe,elem_prod(phiq,dre)) + elem_prod(fe,re*dphiq);
	
	// Caclculate Jacobian matrix
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
	cout<<"Re\n"<<re<<endl;
	cout<<"Ye\n"<<ye<<endl;
	cout<<"dYe\n"<<dye<<endl;
	cout<<"Jacobian\n"<<d2ye<<endl;
	cout<<"Inv(J)\n"<<-inv(d2ye)<<endl;
}
