/**
 * \file OpMod.h
 * \author Steve Martell
 * \date Jun 4, 2013
**/



#include <admodel.h>

#ifndef OM_MODEL_DATA_H
#define OM_MODEL_DATA_H

class OM_model_data
{
public:
  int narea;
  int ngroup;
  int nsex;
  int syr;
  int nyr;
  int sage;
  int nage;
  int ngear;
  // dvector age;
  // int n_ags;
  // int n_ag;
  // ivector i_area;
  // ivector i_group;
  // ivector i_sex;
  // imatrix pntr_ag;
  // d3_array pntr_ags;
  // int nfleet;
  dvector allocation;
  // ivector fsh_flag;
  // ivector ifleet;
  dvector linf;
  dvector vonbk;
  dvector to;
  dvector a;
  dvector b;
  dvector ah;
  dvector gh;
  // dmatrix la;
  // dmatrix wa;
  // dmatrix ma;
  int n_ct_obs;
  dmatrix catch_data;
  // d3_array catch_array;
  // int ft_count;
  int nit;
  ivector nit_nobs;
  ivector survey_type;
  d3_array survey_data;
  // dmatrix it_wt;
  int na_gears;
  ivector na_nobs;
  ivector a_sage;
  ivector a_nage;
  d3_array A;
  // d3_array A_obs;
  int n_wt_nobs;
  d3_array wt_avg;
  // dmatrix wt_bar;
  // d3_array wt_avg;
  // d3_array wt_dev;
  // d3_array wt_mat;
  // data_int eof;

public:
	OM_model_data();
	OM_model_data(
	            const int nArea,
	   			const int nGroup,
	   			const int nSex,
	   			const int nSyr,
	   			const int nNyr,
	   			const int nSage,
	   			const int nNage,
	   			const int nGear,
	   			dvector allocation,
	   			dvector linf,
	   			dvector vonbk,
	   			dvector to,
	   			dvector a,
	   			dvector b,
	   			dvector ah,
	   			dvector gh,
	   			const int n_ct_obs,
	   			dmatrix& catch_data,
	   			const int nit,
	   			ivector nit_nobs,
	   			ivector survey_type,
	   			d3_array survey_data,
	   			const int na_gears,
	   			ivector na_nobs,
	   			ivector a_sage,
	   			ivector a_nage,
	   			d3_array A,
	   			int n_wt_nobs,
	   			d3_array wt_avg
	   			);

	~OM_model_data();

	/* data */
};

#endif

#ifndef OM_MODEL_PARAMETERS_H
#define OM_MODEL_PARAMETERS_H
class OM_model_parameters : public OM_model_data
{
public:
	OM_model_parameters();
	~OM_model_parameters();

	/* data */
};
#endif


