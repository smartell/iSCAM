/**
 * \file OpMod.h
 * \author Steve Martell
 * \date Jun 4, 2013
**/

#ifndef MODELDATA_H
#define MODELDATA_H
class ModelData
{
public:
  ModelData();
  ~ModelData();

  /* data */
};
#endif

#ifndef MODELPARAMS_H
#define MODELPARAMS_H
class ModelParams
{
public:
  ModelParams();
  ~ModelParams();

  /* data */
};
#endif

#ifndef SCENARIO_H
#define SCENARIO_H
class Scenario
{
public:
  Scenario();
  ~Scenario();

  /* data */
};

#endif


#ifndef OPERATINGMODEL_H
#define OPERATINGMODEL_H
class OperatingModel : public Scenario
{
private:
  // Model dimensions
  int nStock;
  int nArea;
  int nSex;
  int nSyr;
  int nNyr;
  int nPyr;
  int nSage;
  int nNage;

  // Class aggregations
  Scenario m_cScenario;
public:
  OperatingModel(Scenario &cScenario);
  ~OperatingModel();

  /* data */
};
#endif