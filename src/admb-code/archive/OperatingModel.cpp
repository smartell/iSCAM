/**
 * @ingroup Milka
 * @author Steve Martell & Catarina Wor
 * @details This is the main operating model to be used in conjuction with iSCAM.  In the
 * FINALS_SECTION of iSCAM, with the command line option '-mse <random number seed>'. The
 * routine called runMSE in iSCAM, first instantiates the omData class, then the omVariables
 * class which are aggregated into the OperatingModel class.
 * 
 * <br> The user will set the variables for each respective class then call the 
 * OperatingModel::runScenario routine which has the following structure:
 * 
 * PSUEDOCODE:<br>
 * <br> - read operating model controls from an mse File.
 * <br> - initialize variables and arrays.
 * <br> - initialize stock-recruitment model.
 * <br> - condition reference population on historical assessment (omVariables).
 * <br> | -- calculate reference points
 * <br> | -- calculate harevest control rule & tac
 * <br> | -- apportion tac across nAreas
 * <br> | -- implement harvest on reference population (ct < bt)
 * <br> | -- update reference population.
 * <br> | -- generate & update data for annual assessment model.
 * <br> | -- run stock assessment procedures.
 * <br> | -- repeat for n simulation years.
 * <br> - write simulation variables to output file.
 * <br> - calculate and write performance statistics to an output file.
 * <br><br>
 * 
 * 
 */


#include <admodel.h>
#include "milka.h"




using namespace mse;
/**
 * @brief Constructor for operating model
 * @details Default constructor for operating model based on aggregation with the omData
 * class and the omVariables class objects used to initialize the model.
 * 
 * @param md Model Data
 * @param mv Model Variables
 */
OperatingModel::OperatingModel(const omData& md, const omVariables& mv, const int& seed)
:m_data(md), m_vars(mv), m_nSeed(seed)
{
	cout<<"In the constructor"<<endl;
	//cout<<"nArea\n"<<m_nArea<<endl;
	exit(1);	
}



OperatingModel::~OperatingModel(){}

void OperatingModel::runScenario(const int &seed)
{
	
}



