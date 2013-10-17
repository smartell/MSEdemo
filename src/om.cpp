	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	#undef COUT
	#define COUT(object) cout<<fixed<<#object "\n"<<object<<endl;
	#include <iostream>
	#include <iomanip>
	using namespace std;
	#include <admodel.h>
	#include <time.h>
	//#include <statsLib.h>
	#include "LRGS.h"
	#include "MSYReferencePoints.h"
	// #include "Scenario.h"
	#include "OperatingModel.h"
	// #include "harvestControlRule.h"
	
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	// sLRGSparameters sPars;
	sLRGSdata data;
	
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
		if( ( on=option_match(ad_comm::argc,ad_comm::argv,"-mse",opt) ) >-1 )
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
  nScenario.allocate("nScenario");
  n_hcr.allocate("n_hcr");
  n_pyr.allocate("n_pyr");
  n_flg_perfect_information.allocate("n_flg_perfect_information");
  iuu_rate.allocate("iuu_rate");
  min_tac.allocate("min_tac");
  sEstimator.allocate("sEstimator");
 data.syr  = syr;
 data.nyr  = nyr;
 data.agek = agek;
 data.ct   = ct;
 data.it   = it;
 cout<<data.it<<endl;
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
	sLRGSparameters sPars;
	sPars.log_bo = log_bo;
	sPars.h  = h;
	sPars.s = s;
	sPars.log_sigma = log_sigma;
	sPars.log_tau  = log_tau;
	sPars.wt = wt;
	bo               = mfexp(log_bo);
	sig              = sqrt(1.0/mfexp(log_sigma));
	tau              = sqrt(1.0/mfexp(log_tau));
	reck             = 4.*h/(1.-h);
	// Testing out a new class called LRGS for doing all of the 
	// model calculations.
	// LRGS cLRGSmodel(syr,nyr,agek,bo,h,s,sig,tau,ct,it,wt);
	// Test cTest;
	LRGS cLRGSmodel(data,sPars);
	cLRGSmodel.initialize_model();
	cLRGSmodel.population_dynamics();
	cLRGSmodel.observation_model();
	epsilon = cLRGSmodel.get_epsilon();
	sd_dep  = cLRGSmodel.get_depletion();
	bt      = cLRGSmodel.get_bt();	
	ft      = cLRGSmodel.get_ft();
	q       = cLRGSmodel.get_q();
	// cout<<q<<endl;
	// initialize_model();
	// population_dynamics();
	// observation_model();
	calc_objective_function();
}

void model_parameters::calc_Reference_Points()
{
	int i;
	int j;
	dvector fe(1,100);
	double be,ce;
	for( i = 1; i <= 100; i++ )
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

void model_parameters::initialize_model()
{
	rt.initialize();
	bt.initialize();
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

void model_parameters::population_dynamics()
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

void model_parameters::observation_model()
{
	int i;
	dvar_vector zt = log(it) - log(bt(syr,nyr));
	q              = exp(mean(zt));
	epsilon        = zt - mean(zt);
}

void model_parameters::run_mse()
{
	int j;
	Scenario cScenario1(agek,nScenario,n_pyr,n_flg_perfect_information,rseed,value(bo),
	                    value(h),value(s),iuu_rate,
	                    value(q),value(sig),value(tau),value(ft),
	                    value(wt),it,ct,min_tac);
	int e_hcr = n_hcr;
	HarvestControlRule c_hcr(e_hcr);
	EstimatorClass cEstimator(sEstimator);
	OperatingModel cOM(cScenario1,cEstimator,c_hcr);
	cOM.runMSEscenario(cScenario1);
	ofstream ofs("OM.rep",ios::app);
	ofs<<"t_bo\n"      << cOM.get_bo()       <<endl;
	ofs<<"t_est_bo\n"  << cOM.get_est_bo()   <<endl;
	ofs<<"t_bmsy\n"    << cOM.get_bmsy()     <<endl;
	ofs<<"t_fmsy\n"    << cOM.get_fmsy()     <<endl;
	ofs<<"t_msy\n"     << cOM.get_msy()      <<endl;
	ofs<<"t_bt\n"      << cOM.get_bt()       <<endl;
	ofs<<"t_aav\n"     << cOM.get_aav()      <<endl;
	ofs.close();
	// |--------------------------------------------|
	// | NOW RUN THE MODEL WITH PERFECT INFORMATION |
	// |--------------------------------------------|
	Scenario cScenarioP(agek,nScenario,n_pyr,0,rseed,value(bo),
	                    value(h),value(s),iuu_rate,
	                    value(q),value(sig),value(tau),value(ft),
	                    value(wt),it,ct);
	OperatingModel cOMP(cScenarioP,cEstimator,c_hcr);
	cOMP.runMSEscenario(cScenarioP);
	ofstream ofs1("OM.rep",ios::app);
	ofs1<<"p_bo\n"  << cOMP.get_est_bo()   <<endl;
	ofs1<<"p_bmsy\n"<< cOMP.get_bmsy()     <<endl;
	ofs1<<"p_fmsy\n"<< cOMP.get_fmsy()     <<endl;
	ofs1<<"p_msy\n" << cOMP.get_msy()      <<endl;
	ofs1<<"p_bt\n"  << cOMP.get_bt()       <<endl;
	ofs1<<"p_ct\n"  << cOMP.get_hat_ct()   <<endl;
	ofs1<<"p_aav\n" << cOMP.get_aav()      <<endl;
	ofs1.close();
}

void model_parameters::calc_objective_function()
{
	nll.initialize();
	dvariable isig2 = mfexp(log_sigma);
	dvariable itau2 = mfexp(log_tau);
	// No comments
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
	if(fpen>0 && !mc_phase()) cout<<"Fpen = "<<fpen<<endl;
	f = sum(nll) + 100000.*fpen;
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
	REPORT(q);
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
