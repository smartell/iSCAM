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
	cout<<mv.log_ro<<endl;
}

