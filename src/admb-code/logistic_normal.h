/**
	This is a class for implementing the logistic normal negative
	loglikelihoood for composition data.
*/

#include <admodel.h>

#ifndef LOGISTIC_NORMAL_H
#define LOGISTIC_NORMAL_H

class logistic_normal
{
private:
	int         m_b1;
	int         m_b2;
	int         m_B;

	int         m_y1;
	int         m_y2;
	int         m_Y;

	double m_dMinimumProportion;

	dvariable   m_nll;

	dvector     m_dWy;		// relative weight for each year.

	dmatrix     m_O;
	dmatrix     m_Op;

	dvar_matrix m_E;
	dvar_matrix m_Ep;

public:
	~logistic_normal();
	logistic_normal();
	logistic_normal(const dmatrix& _O, const dvar_matrix _E);

	/* data */
	dvariable negative_loglikelihood();


	/* setters */
	void set_MinimumProporiton(double &p) {m_dMinimumProportion = p;}
};


#endif