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
  dmatrix theta;
  // param_init_bounded_matrix_vector sel_par;
  // param_init_bounded_vector log_ft_pars;
  // param_init_bounded_matrix init_log_rec_devs;
  // param_init_bounded_matrix log_rec_devs;
  // param_init_bounded_vector log_m_nodes;
  // param_number prior_function_value;
  // param_number likelihood_function_value;
  // objective_function_value objfun;
  // param_number m_bar;
  // param_number rho;
  // param_number varphi;
  // param_number sig;
  // param_number tau;
  // param_vector ro;
  // param_vector bo;
  // param_vector sbo;
  // param_vector kappa;
  // param_vector steepness;
  // param_vector so;
  // param_vector beta;
  // param_vector m;
  // param_vector log_avgrec;
  // param_vector log_recinit;
  // param_vector q;
  // param_vector ct;
  // param_vector eta;
  // param_vector log_m_devs;
  // param_matrix ft;
  // param_matrix log_rt;
  // param_matrix nlvec;
  // param_matrix epsilon;
  // param_matrix it_hat;
  // param_matrix qt;
  // param_matrix sbt;
  // param_matrix rt;
  // param_matrix delta;
  // param_3array F;
  // param_3array M;
  // param_3array Z;
  // param_3array S;
  // param_3array N;
  // param_3array A_hat;
  // param_3array A_nu;
  // param_4array log_sel;
  // param_stddev_vector sd_depletion;
public:
  OM_model_parameters();
	OM_model_parameters(
                      dmatrix theta
                      );
	~OM_model_parameters();

	/* data */
};
#endif


