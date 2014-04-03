// milka.cpp
/**
 * @Milka Source code for Operating Model
 * @author Steven Martell & Catarina Wor
 * @details The default constuctort uses model_data as the base
 * class for the OperatingModel class.  
 * 
 * The OperatingModel class has a major function called 
 * runScenario:
 * 
 * runScenario:
 * 		|- readMSEcontrols
 * 		|- initParameters
 * 			|- surveyQ
 * 			|- stock-recruitment parameters
 * 		|- conditionReferenceModel
 * 		|- setRandomVariables
 * 		|- | getReferencePointsAndStockStatus
 * 		   | calculateTAC
 * 		   | allocateTAC
 * 		   | implementFisheries
 * 		   		|- calcSelectivity
 * 		   		|- calcRetentionDiscards
 * 		   		|- calcTotalMortality
 * 		   | updateReferenceModel
 * 		   | writeDataFile
 * 		   | runStockAssessment
 * 		|- |
 * 		|- writeSimulationVariables
 * 		|- calculatePerformanceMetrics
 * 		
 */


#include <admodel.h>
#include "milka.h"
#include <contrib.h>

// Destructor
OperatingModel::~OperatingModel(){}

// Constructor
OperatingModel::OperatingModel(ModelVariables _mv,int argc,char * argv[])
:model_data(argc,argv), mv(_mv)
{
	cout<<"Inheritance version using model_data as base class"<<endl;
	cout<<"Ngroup "<<ngroup<<endl;
	cout<<"Catch Data\n"<<dCatchData<<endl;
	cout<<"d3 Survey Data\n"<<d3_survey_data<<endl;
	cout<<"eof "<<eof<<endl;

}


void OperatingModel::runScenario()
{
	readMSEcontrols();

	initParameters();

	conditionReferenceModel();

	setRandomVariables();

	for(int i = nyr+1; i <= nyr+10; i++ )
	{
		 
		getReferencePointsAndStockStatus();

		calculateTAC();

		allocateTAC();

		implementFisheries();



		updateReferenceModel();

		writeDataFile();

		runStockAssessment();
	}

}

/**
 * @brief Read control file for Management Strategy Evaluation.
 * @details [long description]
 */
void OperatingModel::readMSEcontrols()
{
	if(verbose) cout<<"MSE Control file\n"<<ProjControlFile<<endl;

	cifstream ifs(ProjControlFile);
	ifs>>m_nPyr;

	cout<<m_nPyr<<endl;

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

	for(int ih = 1; ih <= n_ag; ih++ )
	{
		m_dRbar  = exp(mv.log_rbar(ih));
		m_dRinit = exp(mv.log_rinit(ih));
	}

	switch(int(d_iscamCntrl(2)))
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

void OperatingModel::conditionReferenceModel()
{

}

void OperatingModel::setRandomVariables()
{

}


void OperatingModel::getReferencePointsAndStockStatus()
{
	// read iscam.mse file to get this information.
}


void OperatingModel::calculateTAC()
{

}

void OperatingModel::allocateTAC()
{

}

void OperatingModel::implementFisheries()
{

}

void OperatingModel::updateReferenceModel()
{

}

void OperatingModel::writeDataFile()
{

}

void OperatingModel::runStockAssessment()
{

}
