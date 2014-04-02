// milka.cpp

#include <admodel.h>
#include "milka.h"
#include "contrib.h"

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
	int pyr = 2030;

	// Initialize Catch array
	m_nCtNobs = md.nCtNobs + (pyr-md.nNyr+1)*md.n_ags*md.nFleet;
	m_dCatchData.allocate(1,m_nCtNobs,1,7);
	m_dCatchData.initialize();
	m_dCatchData.sub(1,md.nCtNobs) = md.dCatchData;

	// Initialize Log_selectivity array (4d array)
	d4_logSel.allocate(1,md.nGear,1,md.n_ags,md.nSyr,pyr,md.nSage,md.nNage);
	d4_logSel.initialize();

	m_log_sel_par.allocate(*mv.d3_log_sel_par);
	m_log_sel_par = *mv.d3_log_sel_par;

	cout<<m_log_sel_par(1)<<endl;

	initParameters();
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



/**
 * @brief Routine for calculating selectivities.
 * @details Selectivity models.
 * 
 */
void OperatingModel::calcSelectivities(const ivector& isel_type)
{


	int ig,i,j,k,byr,bpar,kgear;
	double tiny=1.e-10;
	double p1,p2,p3;
	dvector age_dev=md.age;
	dmatrix t1;
	dmatrix   tmp(md.nSyr,md.nNyr-1,md.nSage,md.nNage);
	dmatrix  tmp2(md.nSyr,md.nNyr,md.nSage,md.nNage);
	dmatrix ttmp2(md.nSage,md.nNage,md.nSyr,md.nNyr);

	for( kgear = 1; kgear <= md.nGear; kgear++ )
	{
		k = kgear;
		if(md.i_sel_phz(k) < 0)
		{
			k = abs(md.i_sel_phz(kgear));
		}

		for( ig = 1; ig <= md.n_ags; ig++ )
		{
			tmp.initialize(); tmp2.initialize();
			dvector iy(1,md.n_yr_nodes(k));
			dvector ia(1,md.n_age_nodes(k));
			byr  = 1;
			bpar = 0; 

			switch(isel_type(k))
			{
				case 1: //logistic selectivity (2 parameters)
					for(i=md.nSyr; i<=md.nNyr; i++)
					{
						if( i == md.n_sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < md.nSelBlocks(k) ) byr++;
						}
						p1 = mfexp(m_log_sel_par(k,bpar,1));
						p2 = mfexp(m_log_sel_par(k,bpar,2));
						d4_logSel(kgear)(ig)(i) = log( plogis(md.age,p1,p2)+tiny );
					}
					break;

				case 6:	// fixed logistic selectivity
					p1 = mfexp(m_log_sel_par(k,1,1));
					p2 = mfexp(m_log_sel_par(k,1,2));
					for(i=md.nSyr; i<=md.nNyr; i++)
					{
						d4_logSel(kgear)(ig)(i) = log( plogis(md.age,p1,p2) );
					}
					break;

				case 2:	// age-specific selectivity coefficients
					for(i=md.nSyr; i<=md.nNyr; i++)
					{
						if( i == md.n_sel_blocks(k,byr) )
						{
							bpar ++;
							if( byr < md.nSelBlocks(k) ) byr++;
						}
						for(j=md.nSage;j<=md.nNage-1;j++)
						{
							d4_logSel(k)(ig)(i)(j)   = m_log_sel_par(k)(bpar)(j-md.nSage+1);
						}
						d4_logSel(kgear)(ig)(i,md.nNage) = d4_logSel(k)(ig)(i,md.nNage-1);
					}
					break;

				case 3:	// cubic spline 
					for(i=md.nSyr; i<md.nNyr; i++)
					{
						if( i==md.n_sel_blocks(k,byr) )
						{
							bpar ++;	
							d4_logSel(k)(ig)(i)=cubic_spline( m_log_sel_par(k)(bpar),md.age );
							if( byr < md.nSelBlocks(k) ) byr++;
						}
						d4_logSel(kgear)(ig)(i+1) = d4_logSel(k)(ig)(i);
					}
					break;

				case 4:	// time-varying cubic spline every year				
					for(i=md.nSyr; i<=md.nNyr; i++)
					{
						d4_logSel(kgear)(ig)(i) = cubic_spline(m_log_sel_par(k)(i-md.nSyr+1),md.age);
					}
					break;


			} // end of switch
		}


	}
}


dvector cubic_spline(const dvector& spline_coffs, const dvector& la)
  {
	/*interplolation for selectivity coefficeients*/
	int nodes=size_count(spline_coffs);
	dvector ia(1,nodes);
	ia.fill_seqadd(0,1./(nodes-1));
	dvector fa = (la-min(la))/(max(la)-min(la));
	vcubic_spline_function ffa(ia,spline_coffs);

	return(value(ffa(fa)));

  }










