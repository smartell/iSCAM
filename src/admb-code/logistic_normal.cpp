#include <admodel.h>
#include "logistic_normal.h"

logistic_normal::~logistic_normal()
{}

logistic_normal::logistic_normal()
{}

logistic_normal::logistic_normal(const dmatrix& _O, const dvar_matrix _E)
: m_O(_O), m_E(_E)
{
	/*
	O - observed numbers-at-age or proportions-at-age.
	E - expected numbers-at-age or porportions-at-age.

	Matrix rows correspond to years, cols age or length bins.
	*/
	m_b1 = m_O.colmin();
	m_b2 = m_O.colmax();
	m_y1 = m_O.rowmin();
	m_y2 = m_O.rowmax();

	int i;
	// Ensure proportions-at-age/size sume to 1.
	for( i = m_y1; i <= m_y2; i++ )
	{
		m_Op(i) = m_O(i) / sum(m_O(i));
		m_Ep(i) = m_E(i) / sum(m_E(i));
	}

	// Minimum proportion to pool into adjacent cohort.
	m_dMinimumProportion = 0;
}


dvariable logistic_normal::negative_loglikelihood()
{
	m_nll = 0;
	return(0);
}