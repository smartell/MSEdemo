	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	#undef COUT
	#define COUT(object) cout<<fixed<<#object "\n"<<object<<endl;
	#include <iostream>
	#include <iomanip>
	using namespace std;
	#include <admodel.h>
	#include <time.h>
	#include <statsLib.h>
	#include "MSYReferencePoints.h"
	// #include "Scenario.h"
	#include "OperatingModel.h"
	// #include "harvestControlRule.h"
	
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
	
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <om.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
		int opt, on;
		do_mse = 0;
		rseed  = 999;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-mse",opt))>-1)
		{
			do_mse = 1;
			rseed=atoi(ad_comm::argv[on+1]);
			COUT(rseed);
		}
  agek.allocate("agek");
  syr.allocate("syr");
  nyr.allocate("nyr");
  iyr.allocate(syr,nyr,"iyr");
  ct.allocate(syr,nyr,"ct");
  it.allocate(syr,nyr,"it");
  n_hcr.allocate("n_hcr");
  n_pyr.allocate("n_pyr");
  sEstimator.allocate("sEstimator");
}

void model_parameters::initializationfunction(void)
{
  log_bo.set_initial_value(8.0);
  h.set_initial_value(0.9);
  s.set_initial_value(0.9);
  log_sigma.set_initial_value(4.65);
  log_tau.set_initial_value(4.65);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_bo.allocate(0,10,1,"log_bo");
  h.allocate(0.2,1.0,1,"h");
  s.allocate(0.0,1.0,1,"s");
  log_sigma.allocate(2,"log_sigma");
  log_tau.allocate(3,"log_tau");
  wt.allocate(syr,nyr,-15,15,3,"wt");
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  bo.allocate("bo");
  #ifndef NO_AD_INITIALIZE
  bo.initialize();
  #endif
  ro.allocate("ro");
  #ifndef NO_AD_INITIALIZE
  ro.initialize();
  #endif
  sig.allocate("sig");
  #ifndef NO_AD_INITIALIZE
  sig.initialize();
  #endif
  tau.allocate("tau");
  #ifndef NO_AD_INITIALIZE
  tau.initialize();
  #endif
  a.allocate("a");
  #ifndef NO_AD_INITIALIZE
  a.initialize();
  #endif
  b.allocate("b");
  #ifndef NO_AD_INITIALIZE
  b.initialize();
  #endif
  reck.allocate("reck");
  #ifndef NO_AD_INITIALIZE
  reck.initialize();
  #endif
  q.allocate("q");
  #ifndef NO_AD_INITIALIZE
  q.initialize();
  #endif
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  bt.allocate(syr,nyr+1,"bt");
  #ifndef NO_AD_INITIALIZE
    bt.initialize();
  #endif
  rt.allocate(syr,nyr,"rt");
  #ifndef NO_AD_INITIALIZE
    rt.initialize();
  #endif
  ft.allocate(syr,nyr,"ft");
  #ifndef NO_AD_INITIALIZE
    ft.initialize();
  #endif
  epsilon.allocate(syr,nyr,"epsilon");
  #ifndef NO_AD_INITIALIZE
    epsilon.initialize();
  #endif
  nll.allocate(1,8,"nll");
  #ifndef NO_AD_INITIALIZE
    nll.initialize();
  #endif
  sd_dep.allocate("sd_dep");
}

void model_parameters::userfunction(void)
{
  f =0.0;
	initialize_model();
	population_dynamics();
	observation_model();
	calc_objective_function();
}

void model_parameters::initialize_model(void)
{
	bo               = mfexp(log_bo);
	ro               = bo*(1.-s);
	reck             = 4.*h/(1.-h);
	a                = reck*ro/bo;
	b                = (reck-1.0)/bo;
	rt(syr,syr+agek) = ro * exp(wt(syr,syr+agek));
	bt(syr)          = bo;
	sig              = sqrt(1.0/mfexp(log_sigma));
	tau              = sqrt(1.0/mfexp(log_tau));
}

void model_parameters::population_dynamics(void)
{
	int i;
	dvariable btmp;
	fpen.initialize();
	for(i=syr;i<=nyr;i++)
	{
		ft(i) = -log((-ct(i)+bt(i))/bt(i));
		if(i-syr > agek)
		{
			rt(i) = a*bt(i-agek)/(1.+b*bt(i-agek)) * exp(wt(i));	
		}
		btmp    = s*bt(i) + rt(i) - ct(i);
		bt(i+1) = posfun(btmp,0.1,fpen);
	}
	sd_dep = bt(nyr)/bo;
}

void model_parameters::observation_model(void)
{
	dvar_vector zt = log(it) - log(bt(syr,nyr));
	q              = exp(mean(zt));
	epsilon        = zt - mean(zt);
}

void model_parameters::calc_objective_function(void)
{
	nll.initialize();
	dvariable isig2 = mfexp(log_sigma);
	dvariable itau2 = mfexp(log_tau);
	// negative loglikelihoods and priors stored in nll vector
	nll(1) = dnorm(epsilon,sig);
	nll(2) = dbeta(s,30.01,10.01);
	nll(3) = dnorm(log_bo,log(3000),1.0);
	nll(4) = dbeta((h-0.2)/0.8,1.01,1.01);
	nll(5) = dgamma(isig2,1.01,1.01);
	if(active(log_tau))
	{
		nll(6) = dgamma(itau2,1.01,1.01);
		nll(7) = dnorm(wt,tau);
	}
	else
	{
		nll(7) = dnorm(wt,1.0);
	}
	// objective function + penalty
	if(fpen>0 && !mc_phase()) cout<<"Fpen = "<<fpen<<endl;
	f = sum(nll) + 100000.*fpen;
}

void model_parameters::calcReferencePoints(void)
{
	// Numerically checking my calculus.
	int i,j;
	dvector fe(1,100);
	double be,ce;
	for( i = 1; i <= 1000; i++ )
	{
		be = value(bo);
		fe(i) = double((i-1.)/(100.-1.))*1.2;
		for( j = 1; j <= 100; j++ )
		{
			ce = be * (1.-exp(-fe(i)));
			be = value(s*be + a*be/(1.+b*be)) - ce;
		}
		cout<<setprecision(5)<<fe(i)<<"\t"<<ce<<"\t"<<be<<endl;
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
	msy_reference_points cMSY(value(reck),value(s),value(bo));
	REPORT(bo);
	REPORT(h);
	REPORT(s);
	REPORT(sig);
	REPORT(tau);
	double fmsy = cMSY.get_fmsy();
	double bmsy = cMSY.get_bmsy();
	double msy  = cMSY.get_msy();
	REPORT(fmsy);
	REPORT(bmsy);
	REPORT(msy);
	REPORT(ft);
	REPORT(wt);
	REPORT(bt);
	REPORT(ct);
	REPORT(epsilon);
	// print mle estimates of key parameters for MSE
	ofstream ofs("mse.par");
	ofs<<bo<<"\n"<<reck<<"\n"<<s<<"\n"<<bt(nyr+1)<<endl;		
}

void model_parameters::final_calcs()
{
	if(do_mse)
	{
		cout<<"Running MSE"<<endl;
		run_mse();
	}
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;
	// Check calculus: A-O-K
	// calcReferencePoints();
}

void model_parameters::run_mse(void)
{
	// This is the entire management strategy evaluaiton routine.  
	// So far I use 3 class objects to this via OOP.
	// 1) The scenario class: -parameters & data for operating model
	// 2) The harvestControlRule class: Use FORTY_TEN, FIXED_HARVEST_RATE, FIXED_ESCAPMENT
	// 3) The Operating model class: call .runMSEScenario to run the simulation.
	// 4) The OPertating model calls msyrefPoints.h to calculate Fmsy etc.
	// Scenario class
	Scenario cScenario1(agek,n_pyr,rseed,value(bo),value(h),value(s),
	                    value(q),value(sig),value(tau),value(ft),
	                    value(wt),it,ct);
	// Harvest control rule
	// int e_hcr = HarvestControlRule::FORTY_TEN;
	// int e_hcr = HarvestControlRule::FIXED_ESCAPEMENT;
	// int e_hcr = HarvestControlRule::FIXED_ESCAPEMENT_CAP;
	// int e_hcr = HarvestControlRule::FIXED_HARVEST_RATE;
	// int e_hcr = HarvestControlRule::CONDITIONAL_CONSTANT_CATCH;
	int e_hcr = n_hcr;
	HarvestControlRule c_hcr(e_hcr);
	// Estimator class (allow user defined estimator)
	// Operating model class
	OperatingModel cOM(cScenario1,c_hcr);
	cOM.runMSEscenario(cScenario1);
	ofstream ofs("OM.rep",ios::app);
	ofs<<"t_bo\n"  << cOM.get_bo()   <<endl;
	ofs<<"t_bmsy\n"<< cOM.get_bmsy() <<endl;
	ofs<<"t_fmsy\n"<< cOM.get_fmsy() <<endl;
	ofs<<"t_msy\n" << cOM.get_msy()  <<endl;
	ofs<<"t_bt\n"  << cOM.get_bt()   <<endl;
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
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 
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
