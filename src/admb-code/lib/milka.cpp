// milka.cpp

#include <admodel.h>
#include "milka.h"

using namespace mse;

OperatingModel::~OperatingModel(){}

/**
 * @brief Operating Model constructor
 * @details Default constructor for the operating model.
 * 
 * @param _md Model Data Struct
 * @param _mv Model Variables Struct
 * 
 */
OperatingModel::OperatingModel(const ModelData &_md, const ModelVariables &_mv)
:md(_md), mv(_mv)
{
	cout<<"In the constructor"<<endl;
	cout<<mv.log_ro<<endl;
}

void OperatingModel::initParameters()
{
	m_dRo        = exp(mv.log_ro);
	m_dSteepness = mv.steepness;
	m_dM         = exp(mv.m);
	m_dRho       = mv.rho;
	m_dVarphi    = sqrt(1.0/mv.varphi);
	m_dSigma     = sqrt(m_dRho) * m_dVarphi;
	m_dTau       = sqrt(1.0-m_dRho)*m_dVarphi;

	for(int ih = 1; ih <= md.n_ag; ih++ )
	{
		m_dRbar  = exp(mv.log_rbar(ih));
		m_dRinit = exp(mv.log_rinit(ih));
	}

	switch(int(md.d_iscamCntrl(2)))
	{
		case 1:
			//Beverton-Holt model
			m_dKappa = elem_div(4.*m_dSteepness,(1.-m_dSteepness));
			break;
		case 2:
			//Ricker model
			m_dKappa = pow((5.*m_dSteepness),1.25);
		break;
	}

}