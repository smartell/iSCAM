#include <admodel.h>
#include "OpMod.h"
// |---------------------------------------------------------------------------------|
// | OM_model_data
// |---------------------------------------------------------------------------------|
// |
OM_model_data::~OM_model_data(){}


OM_model_data::OM_model_data(
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
	   						d3_array wt_avg)
:	narea(nArea),
	ngroup(nGroup),
	nsex(nSex),
	syr(nSyr),
	nyr(nNyr),
	sage(nSage),
	nage(nNage),
	ngear(nGear),
	allocation(allocation),
	linf(linf),
	vonbk(vonbk),
	to(to),
	a(a),
	b(b),
	ah(ah),
	gh(gh),
	n_ct_obs(n_ct_obs),
	catch_data(catch_data),
	nit(nit),
	nit_nobs(nit_nobs),
	survey_type(survey_type),
	survey_data(survey_data),
	na_gears(na_gears),
	na_nobs(na_nobs),
	a_sage(a_sage),
	a_nage(a_nage),
	A(A),
	n_wt_nobs(n_wt_nobs),
	wt_avg(wt_avg)
{
	cout<<"Read in OM model data"<<endl;
	cout<<A(1)<<endl;
}


// |---------------------------------------------------------------------------------|
// | OM_model_parameters
// |---------------------------------------------------------------------------------|
// |

OM_model_parameters::~OM_model_parameters()
{

}

OM_model_parameters::OM_model_parameters(
                                         dmatrix theta
                                         )
:theta(theta)
{
	cout<<"Read in OM model parameters"<<endl;
}
