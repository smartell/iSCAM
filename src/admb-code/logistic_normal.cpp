#include <admodel.h>
#include "logistic_normal.h"

logistic_normal::~logistic_normal()
{}

logistic_normal::logistic_normal()
{}

logistic_normal::logistic_normal(const dmatrix& _O, const dvar_matrix& _E)
: m_O(_O), m_E(_E)
{
	/*
	O - observed numbers-at-age or proportions-at-age.
	E - expected numbers-at-age or porportions-at-age.

	Matrix rows correspond to years, cols age or length bins.
	*/
	m_b1 = m_O.colmin();
	m_b2 = m_O.colmax();
	m_B  = m_b2 - m_b1 + 1;

	m_y1 = m_O.rowmin();
	m_y2 = m_O.rowmax();
	m_Y  = m_y2 - m_y1 + 1;

	m_dWy.allocate(m_y1,m_y2);
	m_dWy.initialize();

	m_Op.allocate(m_y1,m_y2,m_b1,m_b2);
	m_Ep.allocate(m_y1,m_y2,m_b1,m_b2);
	m_Ox.allocate(m_y1,m_y2,m_b1,m_b2);
	m_Ex.allocate(m_y1,m_y2,m_b1,m_b2);

	int i;
	for( i = m_y1; i <= m_y2; i++ )
	{
		// Ensure proportions-at-age/size sume to 1.
		m_Op(i) = m_O(i) / sum( m_O(i) );
		m_Ep(i) = m_E(i) / sum( m_E(i) );

		// Logistic transformation
		m_Ox(i) = exp(m_Op(i)) / sum(exp(m_Op(i)));
		m_Ex(i) = exp(m_Ep(i)) / sum(exp(m_Ep(i)));
	}
	
	// Calculate mean weighting parameters for each year.
	dvector dNy = rowsum(m_O);	
	double dN   = mean(dNy);
	m_dWy       = sqrt( dN / dNy );

	
	// Minimum proportion to pool into adjacent cohort.
	m_dMinimumProportion = 0;
	aggregate_arrays();
}


dvariable logistic_normal::negative_loglikelihood()
{
	double t1 = 0.5 * m_Y * (m_B -1.0) * log(2.0*PI);
	double t3 = sum( log(m_Op) );

	m_nll = 0;

		cout<<"Here I am in the likelihood"<<endl;
		cout<<t1<<endl;
		cout<<t3<<endl;
		
	return(0);
}


void logistic_normal::aggregate_arrays()
{
	/*
		- Aggregate minimum proportions in to adjacent cohorts (tail compression)
		- This routine populates m_Oa and m_Ea as ragged objects to be used for
		  likelihood calculations.
	*/
	m_nNminp.allocate(m_y1,m_y2);
	

	int k;
	for(int i = m_y1; i <= m_y2; i++ )
	{
		int n = 0;  // number of non-zero observations in each year
		for(int j = m_b1; j <= m_b2; j++ )
		{
			if( m_Op(i,j) > m_dMinimumProportion )
			{
				n ++;
			}
		}
		m_nNminp(i) = n;
	}

	// Ragged arrays with tail compression and zeros omitted.
	m_nAgeIndex.allocate(m_y1,m_y2,1,m_nNminp);
	m_Oa.allocate(m_y1,m_y2,1,m_nNminp);
	m_Ea.allocate(m_y1,m_y2,1,m_nNminp);
	m_nAgeIndex.initialize();
	m_Oa.initialize();
	m_Ea.initialize();

	for(int i = m_y1; i <= m_y2; i++ )
	{
		dvector     oo = m_Op(i);
		dvar_vector pp = m_Ep(i);
		k = 1;
		for(int j = m_b1; j <= m_b2; j++ )
		{
			if( oo(j) <= m_dMinimumProportion )
			{
				m_Oa(i,k) += oo(j);
				m_Ea(i,k) += pp(j);
			}
			else
			{
				m_Oa(i,k) += oo(j);
				m_Ea(i,k) += pp(j);

				if(k <=m_nNminp(i)) m_nAgeIndex(i,k) = j;
				if(k < m_nNminp(i)) k++;
			}
		}
	}
}














