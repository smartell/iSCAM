---
layout: post
title:  "IPHC Branches"
date:   2014-03-23 06:50:11
categories: iscam overview
---

#IPHC branches
The master branch on the github repository is dated and contains the original source code that was developed for a single-sex, single stock model.  Branches that begin with the prefix IPHC (short for International Pacific Halibut Commission) contain code for a much more complex population models, observation models, and statistical criterion.  The IPHC branch contains the most stable release of this new model.

##IPHC-developer
New developments in iSCAM at this time are all occurring in the IPHC developer branch.  New features that do not exists in the master branch include:

- dimensions for management or regulatory areas.
- dimensions for the number of stocks or functional groups.
- dimensions for the number of sexes (1 sex, or 2 sex models).

An excerpt from the data template:		
-
	## ------------------------------------------------------------------------- ##
	## MODEL DIMENSIONS                                                          ##
	## ------------------------------------------------------------------------- ##
	1        # -number of areas            (narea)
	1        # -number of groups or stocks (ngroup)
	2        # -number of sexes            (nsex)
	1968     # -first year of data         (syr)
	1979     # -last year of data          (nyr)
	3        # -age of youngest age class  (sage)
	9        # -age of plus group          (nage)
	1        # -number of gears            (ngear)
	##



###1 Spatial considerations

The transition to a spatially explicit model is still under development. Many of the data structures and nearly all of the existing code has been ported over to accommodate the increased dimensionality.  However, movement among spatial has not been fully implemented and is currently work in progress.  More details will come on the movement kernels; but, in short the kernels will can be age/size-dependent or independent and is parameterized using N + 1 parameters for each class, where N is the number of spatial areas.  These parameters can be informed based on relative changes in spatial CPUE data as well as tag-recapture information.

###2 Number of groups or stocks
It is possible to jointly fit the data from two different species that are harvested in the same fishery. A classic example is the trawl fisheries for hake of the west coast of South Africa which harvest two different species of hake.  In such a case, each species would be assigned a unique group number.  The input data file for example, now requires that each catch observation be assigned a unique year, gear, area, group, sex, type, and value (see below).  Different combinations of gear, area, group, can be used to identify unique fisheries in certain areas for certain species by specific gear types.  Similar input information is also required for relative abundance indices, and composition information, that would allow the model to be fitted to data for a given species in a given area, etc.  

	## ------------------------------------------------------------------------- ##
	## TIME SERIES CATCH DATA                                                    ##
	## ------------------------------------------------------------------------- ##
	## Observed catch from all gears, areas, and sex                             ##
	## sex: 1=female, 2=male, 0=asexual
	## Type of catch: legend                                                     ##
	##               1 = catch in weight                                         ##
	##               2 = catch in numbers                                        ##
	##               3 = catch in spawn (roe)                                    ##
	## 
	## n_ct_obs
	   12
	## Year gear area group sex type value
	   1968 1    1    1     0   1    1864
	   1969 1    1    1     0   1    1943
	   1970 1    1    1     0   1    1656
	   1971 1    1    1     0   1    2903
	   1972 1    1    1     0   1    2697
	   1973 1    1    1     0   1    2203
	   1974 1    1    1     0   1    1266
	   1975 1    1    1     0   1    1245
	   1976 1    1    1     0   1    2110
	   1977 1    1    1     0   1    2670
	   1978 1    1    1     0   1    2879
	   1979 1    1    1     0   1    2916
	##

One motivation for jointly fitting data from multiple species, is often the density of one species can affect the catchability of another species, and in this integrated framework we can jointly model changes in commercial catchability associated with changes in stock density (e.g., Meaghan Darcy, PhD Thesis, Buzz Holling disc equation).

###3 Sex
iSCAM now accommodates both single sex and two sex models.  In the case of single sex models, all the data structures must have sex = 0 to denote that the data represents both male and female, or just male only, or female only.  In the case of two sex models, the data structures can be both sexed or unsexed.  The convention is to always use sex=1 for the female index, sex=2 for the male index, and sex=0 to denote sex combined data.  It is possible to fit a model with no sex specific data, however, the growth curves and maturity-at-age information must be specified for each sex.

##Logistic Normal Likelihood
Two likelihood options for composition data were available in the previous versions of iSCAM. These were the Multivariate Logistic likelihood based on papers by Jon Schnute and Laura Richards, and the widely used Multinomial likelihood.  There have been a number of recent papers in the literature that have questioned the statistical assumptions behind the multinomial distribution for age and size composition data, and noted that this distribution does not allow for within year correlations, and estimating the variance scaling (or effective sample size) can be subjective or conditional on structural assumptions and process error assumptions.  The multivariate logistic likelihood addresses some of these issues, namely the data are self weighting based on the MLE estimate of the variance.  Recently, Chris Francis (New Zealand) has done a lot of work in the area of data weighting and correlations in composition data, and has come up with an alternative likelihood based on the Logistic normal distribution. 

In iSCAM, three alternative forms of the logistic normal are available:

* No autocorrelation and use the conditional maximum likelihood estimate of the variance
* Estimate the variance and assume no autocorrelation
* Estimate the variance and assume a lag-1 AR for the correlation structure.
* A four option with estimated variance and a AR2 process for the correlations.




The following header file highlights the public functions available for a new class object that I recently developed to implement Chris Fracis's ideas for the Logisitc Normal likelihood.  It is currently undergoing testing by me and Dave Fournier in his Multifan CL software.  We are both quite happy with the performance so far.

{% highlight cpp %}
class logistic_normal
{
private:
	int m_y1;
	int m_y2;
	int m_b1;
	int m_b2;
	double minp;
	double eps;
	double m_bm1;
	ivector m_nb1;
	ivector m_nb2;
	dvector m_Wy;

	dvariable m_wss;	///> weighted sum of squares  
	dvariable m_nll;
	dvariable m_sigma;
	dvariable m_sigma2;

	imatrix m_nAgeIndex;

	dmatrix m_O;
	dmatrix m_Op;
	dmatrix m_Oz;
	dmatrix m_nAidx;
	dmatrix m_std_residual;

	dvar_vector m_rho;

	dvar_matrix m_E;
	dvar_matrix m_Ep;
	dvar_matrix m_Ez;
	dvar_matrix m_ww;
	dvar3_array m_V;

	logistic_normal();
	dvariable negative_log_likelihood();
	dvector compute_relative_weights(const dmatrix &O);
	
	void get_rho();
	void get_rho(const dvariable &phi);
	void get_rho(const dvariable &phi, const dvariable &psi);

	void std_residuals();

	void compute_correlation_array();
	void compute_likelihood_residuals();
	void compute_weighted_sumofsquares();
	void aggregate_and_compress_arrays();

	friend class logistic_student_t;

public:
	// Constructor
	logistic_normal(const dmatrix& _O,const dvar_matrix& _E,
	                const double _minp,const double _eps=0);

	// Four alternative methods for calculating the nll.
	dvariable operator() ();
	dvariable operator() (const dvariable &sigma2);
	dvariable operator() (const dvariable &sigma2 ,const dvariable &phi);
	dvariable operator() (const dvariable &sigma2 ,const dvariable &phi,
	                      const dvariable &psi);


	// Return the estimated (or mle) of the variance
	double      get_sigma () { return value(m_sigma ); }
	double      get_sigma2() { return value(m_sigma2); }
	dmatrix     get_standardized_residuals() { std_residuals(); return m_std_residual; }
	

};

{% endhighlight %}