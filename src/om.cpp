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
	#include "lrgsOM.h"
	#include "MSYReferencePoints.h"
	// #include "Scenario.h"
	// #include "OperatingModel.h"
	// #include "harvestControlRule.h"
	
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	//sLRGSparameters sPars;  //causes segmenation fault.
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
  ngear.allocate("ngear");
  nEpochs.allocate("nEpochs");
  nIt_nobs.allocate(1,nEpochs,"nIt_nobs");
  catch_data.allocate(syr,nyr,0,ngear,"catch_data");
  it_data.allocate(1,nEpochs,1,nIt_nobs,1,4,"it_data");
  eof.allocate("eof");
 cout<<eof<<endl;
 if(eof != -999){cout<<"Error reading data"<<endl; exit(1);}
  iyr.allocate(syr,nyr);
  ct.allocate(syr,nyr,1,ngear);
  it_yr.allocate(1,nEpochs,1,nIt_nobs);
  epoch.allocate(1,nEpochs,1,nIt_nobs);
  it.allocate(1,nEpochs,1,nIt_nobs);
  cv.allocate(1,nEpochs,1,nIt_nobs);
		iyr   = ivector(column(catch_data,0));
		ct    = trans(trans(catch_data).sub(1,ngear));
		for(int i = 1; i <= nEpochs; i++ )
		{
			it_yr(i) = ivector(column(it_data(i),1));
			epoch(i) = ivector(column(it_data(i),2));
			it(i)    = column(it_data(i),3);
			cv(i)    = column(it_data(i),4);		
		}
 ad_comm::change_datafile_name("om.ctl");
  npar=8;
  d_PC.allocate(1,8,1,7,"d_PC");
 data.syr      = syr;
 data.nyr      = nyr;
 data.agek     = agek;
 data.ngear    = ngear;
 data.nIt_nobs = nIt_nobs;
 data.ct       = ct;
 data.it       = it;
 data.rseed    = rseed;
 data.it_yr    = it_yr;
 data.nEpochs  = nEpochs;
 data.epoch    = epoch;
 data.cv       = cv;
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
 ir=1; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
  log_bo.allocate(dlb,dub,iphz,"log_bo");
 log_bo = dval;
 ir=2; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
  log_b1.allocate(dlb,dub,iphz,"log_b1");
 log_b1 = dval;
 ir=3; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
  h.allocate(dlb,dub,iphz,"h");
 h = dval;
 ir=4; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
  s.allocate(dlb,dub,iphz,"s");
 s = dval;
 ir=5; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
  gamma.allocate(dlb,dub,iphz,"gamma");
 gamma = dval;
 ir=6; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
  log_sigma.allocate(dlb,dub,iphz,"log_sigma");
 log_sigma = dval;
 ir=7; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
  log_tau.allocate(dlb,dub,iphz,"log_tau");
 log_tau = dval;
 ir=8; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
  wt.allocate(syr,nyr,dlb,dub,iphz,"wt");
 wt = dval;
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  bo.allocate("bo");
  #ifndef NO_AD_INITIALIZE
  bo.initialize();
  #endif
  b1.allocate("b1");
  #ifndef NO_AD_INITIALIZE
  b1.initialize();
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
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  q.allocate(1,nEpochs,"q");
  #ifndef NO_AD_INITIALIZE
    q.initialize();
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
  delta.allocate(syr,nyr,"delta");
  #ifndef NO_AD_INITIALIZE
    delta.initialize();
  #endif
  epsilon.allocate(1,nEpochs,1,nIt_nobs,"epsilon");
  #ifndef NO_AD_INITIALIZE
    epsilon.initialize();
  #endif
  negloglike.allocate(1,2,1,nEpochs,"negloglike");
  #ifndef NO_AD_INITIALIZE
    negloglike.initialize();
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
	sPars.log_b1 = log_b1;
	sPars.h  = h;
	sPars.s = s;
	sPars.log_sigma = log_sigma;
	sPars.log_tau  = log_tau;
	sPars.wt = wt;
	sPars.gamma = &gamma;
	bo               = mfexp(log_bo);
	b1               = mfexp(log_b1);
	sig              = sqrt(1.0/mfexp(log_sigma));
	tau              = sqrt(1.0/mfexp(log_tau));
	reck             = 4.*h/(1.-h);
	// Testing out a new class called LRGS for doing all of the 
	// model calculations.
	// LRGS cLRGSmodel(syr,nyr,agek,bo,h,s,sig,tau,ct,it,wt);
	// Test cTest;
	fpen = 0;
	LRGS cLRGSmodel(data,sPars);
	cLRGSmodel.initialize_model();
	cLRGSmodel.population_dynamics();
	cLRGSmodel.observation_model();
	cLRGSmodel.calc_negative_loglikelihoods();
	epsilon = cLRGSmodel.get_epsilon();
	sd_dep  = cLRGSmodel.get_depletion();
	bt      = cLRGSmodel.get_bt();	
	ft      = cLRGSmodel.get_ft();
	q       = cLRGSmodel.get_q();
	negloglike = cLRGSmodel.get_nll();
	delta   = cLRGSmodel.get_delta();
	fpen    = cLRGSmodel.get_fpen();
	//cout<<"Fpen " <<fpen<<endl;
	calc_objective_function();
}

void model_parameters::calc_objective_function()
{
	nll.initialize();
	dvariable isig2 = mfexp(log_sigma);
	dvariable itau2 = mfexp(log_tau);
	// No comments
	for(int i = 1; i <= nEpochs; i++ )
	{
		//nll(1) += dnorm(epsilon(i),sig);
	}
	//nll(2) = dbeta(s,30.01,10.01);
	//The following is based on E(x) = exp(-0.15), Sig2 = (.15*CV)^2, where assumed CV=0.02
	nll(2) = dbeta(s,347.3694,56.21626);
	nll(3) = dnorm(log_bo,log(500),5.0);
	nll(3)+= dnorm(log_b1,log(500),5.0);
	nll(4) = dbeta((h-0.2)/0.8,1.01,1.01);
	nll(5) = dgamma(isig2,1.01,1.01);
	if(active(log_tau))
	{
		nll(6) = dgamma(itau2,1.01,1.01);
		//nll(7) = dnorm(wt,tau);
	}
	else
	{
		//nll(7) = dnorm(wt,1.0);
	}
	if(fpen>0 && !mc_phase()) cout<<"Fpen = "<<fpen<<endl;
	f = sum(nll) + 100000.*fpen + sum(negloglike);
}

void model_parameters::mse2()
{
	//lrgsOM OM;
	//lrgsOM OM("S1.scn");
	sLRGSparameters sPars;
	sPars.log_bo    = log_bo;
	sPars.log_b1    = log_b1;
	sPars.h         = h;
	sPars.s         = s;
	sPars.log_sigma = log_sigma;
	sPars.log_tau   = log_tau;
	sPars.wt        = wt;
	sPars.gamma     = &gamma;
	lrgsOM cOM(data,sPars,"Scenario.scn","ManagementProcedure.mp");
	ofstream ofs("OM.rep",ios::app);
	ofs<<"t_bo\n"      << cOM.get_bo()       <<endl;
	ofs<<"t_est_bo\n"  << cOM.get_est_bo()   <<endl;
	ofs<<"t_bmsy\n"    << cOM.get_bmsy()     <<endl;
	ofs<<"t_fmsy\n"    << cOM.get_fmsy()     <<endl;
	ofs<<"t_msy\n"     << cOM.get_msy()      <<endl;
	ofs<<"t_bt\n"      << cOM.get_bt()       <<endl;
	ofs<<"t_aav\n"     << cOM.get_aav()      <<endl;
	ofs.close();
	cout<<"Running the new Operating Model Class"<<endl;
}

void model_parameters::run_mse()
{
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
	REPORT(b1);
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
	REPORT(delta);
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
		//run_mse();
		mse2();
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
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
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
