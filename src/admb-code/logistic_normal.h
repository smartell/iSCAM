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
	dvariable   m_sig2;

	ivector     m_nNminp;
	dvector     m_dWy;		// relative weight for each year.

	imatrix     m_nAgeIndex;// Index for aggregated residuals
	dmatrix     m_O;		// Raw Data
	dmatrix     m_Op;		// Proportions
	dmatrix     m_Ox;		// Logistic transform.
	dmatrix     m_Oa;		// Aggregated matrix for tail compression and zeros

	dvar_matrix m_E;
	dvar_matrix m_Ep;
	dvar_matrix m_Ex;
	dvar_matrix m_Ea;		// Aggregated matrix for tail compression and zeros

	dvar_matrix m_residual;

public:
	~logistic_normal();
	logistic_normal();
	logistic_normal(const dmatrix& _O, const dvar_matrix& _E);

	/* data */
	dvariable negative_loglikelihood(const dvariable& tau2);
	dvar_matrix standardized_residuals();
	void aggregate_arrays();


	/* setters */
	void set_MinimumProporiton(double &p) {m_dMinimumProportion = p;}
};


#endif