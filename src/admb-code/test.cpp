	#include <contrib.h>
	#include <msy.cpp>
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <test.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
 A = 20;
 ngear = 3;
  age.allocate(1,A);
  ah.allocate(1,ngear);
  gh.allocate(1,ngear);
  wa.allocate(1,A);
  fa.allocate(1,A);
  V.allocate(1,ngear,1,A);
		age.fill_seqadd(1,1);
		ah = 6.;
		gh = 0.5 * ah;
		wa = pow((1.-exp(-0.2*age)),3);
		fa = wa;
		fa(1,4) = 0;
		int k;
		// Selectivities for each gear.
		// Chang the ah gh values here to not how they impact MSY calcs.
		for(k=1;k<=ngear;k++)
		{
			V(k) = log(plogis(age,ah(k)+double(k),gh(k)));  // progressive increase in age-50% harvest
		}
		for(k=1;k<=ngear;k++)
		{
			V(k) = exp(V(k) - log(mean(exp(V(k)))) );
		}
		cout<<"fa\n"<<fa<<endl;
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  dum.allocate("dum");
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::userfunction(void)
{
  f =0.0;
}

void model_parameters::final_calcs()
{
	double ro    = 500.0;
	double h     = 0.75;
	double m     = 0.3;
	double rho   = 0.00;
	int i;
	dvector fe(1,ngear);
	Msy cMSY(ro,h,m,rho,wa,fa,V);
	/* Lets have a look at the unfished spawning biomass per recruit */
	double phie = cMSY.getPhie();
	cout<<"\tUnfished spawning biomass per recruit = "<<phie<<endl;
	/* Now specify an initial guess for Fmsy and get Fmsy from Msy class */
	dvector fmsy(1,ngear);
	fmsy = m*h/0.8;  // initial guess for fmsy
	cMSY.get_fmsy(fmsy);
	//cMSY.get_fmsy_safe(fmsy,1e-3,5.0);
	/* Now print out other associated reference points */
	cout<<"\t---------------------------------------"<<endl;
	cout<<setprecision(5)<<endl;
	cout<<"\t M = 0.3 and h = 0.75"<<endl;
	cout<<"\tUnfished spawning biomass             = "<<cMSY.getBo()<<endl;
	cout<<"\tF_MSY values for each fishing fleet   ="<<fmsy<<endl;
	cout<<"\tMSY values for each fishing fleet     ="<<cMSY.getMsy()<<endl;
	cout<<"\tSpawning biomass at MSY (Bmsy)        = "<<cMSY.getBmsy()<<endl;
	cout<<"\tRecruitment at MSY (Rmsy)             = "<<cMSY.getRmsy()<<endl;
	cout<<"\tSpawning potential ratio at MSY       = "<<cMSY.getSprMsy()<<endl;
	exit(1);
	/* Now change M and steepness and recalculate reference poinots */
	h = 0.80;
	m = 0.20;
	cMSY.set_h(h);
	cMSY.set_m(m);
	cMSY.calc_bo(m,fa);  //necessary to update phie due to change in m.
	fmsy = 0.6*m/ngear;  // initial guess for fmsy
	cMSY.get_fmsy(fmsy);
	/* Now print out other associated reference points */
	cout<<"\t---------------------------------------"<<endl;
	cout<<setprecision(4)<<endl;
	cout<<"\t M = 0.2 and h = 0.80"<<endl;
	cout<<"\tUnfished spawning biomass             = "<<cMSY.getBo()<<endl;
	cout<<"\tF_MSY values for each fishing fleet   ="<<fmsy<<endl;
	cout<<"\tMSY values for each fishing fleet     ="<<cMSY.getMsy()<<endl;
	cout<<"\tSpawning biomass at MSY (Bmsy)        = "<<cMSY.getBmsy()<<endl;
	cout<<"\tRecruitment at MSY (Rmsy)             = "<<cMSY.getRmsy()<<endl;
	cout<<"\tSpawning potential ratio at MSY       = "<<cMSY.getSprMsy()<<endl;
	/* Now run the model with a fixed allocation among gears */
	dvector allocation(1,ngear);
	allocation.fill_seqadd(1,1);
	allocation/=sum(allocation);
	cout<<"\t---------------------------------------"<<endl;
	cMSY.get_fmsy(fmsy,allocation);
	cout<<setprecision(4)<<endl;
	cout<<"\t M = 0.2 and h = 0.80, & allocations  "<<endl;
	cout<<"\tAllocaton for each gear               ="<<allocation<<endl;
	cout<<"\tUnfished spawning biomass             = "<<cMSY.getBo()<<endl;
	cout<<"\tFishing rate for each fishing fleet   ="<<fmsy<<endl;
	cout<<"\tYield values for each fishing fleet   ="<<cMSY.getYe()<<endl;
	cout<<"\tSpawning biomass at given allocation  = "<<cMSY.getBe()<<endl;
	cout<<"\tRecruitment at given allocaiton       = "<<cMSY.getRe()<<endl;
	cout<<"\tSpawning potential ratio              = "<<cMSY.getSpr()<<endl;
	cout<<"\t---------------------------------------"<<endl;
	cout<<"fe\t ye \t dye \t re"<<endl;
	for(i=1;i<=100;i++)
	{
		fe = double(i-1.)/100;
		cMSY.calc_equilibrium(fe);
		cout<<fe<<"\t"<<cMSY.getYe()<<cMSY.getdYe()<<"\t"<<cMSY.getRe()<<endl;
	}
	cout<<"\t---------------------------------------"<<endl;
	ro     =  3598.1;
	h      =  0.62303;
	m      =  0.3;
	/* The data below are from the Flak Lake Model (catage demo in ADMB)*/
	dvector waa  =   ("{1.0121, 1.3858, 1.725, 2.0174, 2.2609, 2.4589, 2.6171}");
	dvector faa  =   ("{0.072924, 0.66935, 1.5841, 2.0026, 2.2595, 2.4588, 2.6171}");
	dvector tV   =   ("{0.0041512, 0.084884, 0.60249, 1.8459, 1.7443, 1.3592, 1.3592}");
	dmatrix VV(1,ngear,1,7);
	VV(1) = tV;
	cout<<"Flak Lake Trout"<<endl;
	Msy cFLK(ro,h,m,rho,waa,faa,VV);
	rho = 0.5;
	Msy cFLK5(ro,h,m,rho,waa,faa,VV);
	fe = 0.1;
	cout<<"fe\t ye1 \t ye2\t dye1 \t dye2\t re"<<endl;
	for(i=1;i<=100;i++)
	{
		fe = double(i-1.)/100;
		cFLK.calc_equilibrium(fe);;
		cFLK5.calc_equilibrium(fe);
		cout<<fe<<"\t"<<cFLK.getYe()<<cFLK5.getYe()<<cFLK.getdYe()<<cFLK5.getdYe()<<"\t"<<cFLK.getRe()<<endl;
	}
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
}

void model_parameters::preliminary_calculations(void){
  admaster_slave_variable_interface(*this);
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
