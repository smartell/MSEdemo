/// @file OM.tpl
/// @author Steve Martell 




///
/// \def REPORT(object)
/// \brief Prints name and value of \a object on ADMB report %ofstream file.
///

///
/// \def COUT(object)
/// \brief Screen dump using cout<<"object\n"<<object<<endl;
///


/// \remarks
///  ******************************************************************
///  | OM, short for Operating model.
///  |
///  | Created by Martell on 2013-05-14.
///  | Copyright (c) 2013. All rights reserved.
///  | Comments: This operation model is conditioned on the LRGS model
///  |           based on chapter 10 in the Ecological Detective.
///  ******************************************************************

DATA_SECTION
	
	// |---------------------------------------------------------------------------------|
	// | COMMAND LINE OPTIONS
	// |---------------------------------------------------------------------------------|
	// |
	int on;   		/// flag for mse */
	int do_mse;
	int rseed;
	
	LOC_CALCS
		int opt, on;
		do_mse = 0;
		rseed  = 999;
		if( ( on=option_match(ad_comm::argc,ad_comm::argv,"-mse",opt) ) >-1 )
		{
			do_mse = 1;
			rseed=atoi(ad_comm::argv[on+1]);
			COUT(rseed);
		}
	END_CALCS

	init_int agek;
	init_int syr;
	init_int nyr;
	init_ivector iyr(syr,nyr);
	init_vector ct(syr,nyr);
	init_vector it(syr,nyr);

	// |---------------------------------------------------------------------------------|
	// | MANAGEMENT STRATEGY EVALUATION COMMANDS
	// |---------------------------------------------------------------------------------|
	// |
	// Scenario 1 = stationarity
	// Scenario 2 = pdo
	init_int nScenario;
	// Harvest control rule
	init_int n_hcr;
	init_int n_pyr;
	init_adstring sEstimator;
	// !! COUT(sEstimator);
	// !! exit(1);
	// !! sLRGSdata data;
	!! data.syr  = syr;
	!! data.nyr  = nyr;
	!! data.agek = agek;
	!! data.ct   = ct;
	!! data.it   = it;
	!! cout<<data.it<<endl;
	
INITIALIZATION_SECTION
	log_bo     8.0;
	h          0.9;
	s          0.9;
	log_sigma  4.65;
	log_tau    4.65;
	

PARAMETER_SECTION
	init_bounded_number log_bo(0,10,1);
	init_bounded_number h(0.2,1.0,1);
	init_bounded_number s(0.0,1.0,1);
	init_number log_sigma(2);
	init_number log_tau(3);
	init_bounded_dev_vector wt(syr,nyr,-15,15,3);
	

	objective_function_value f;

	number bo;	
	number ro;	
	number sig;
	number tau;	
	number a;	
	number b;	
	number reck;
	number q;	
	number fpen;

	vector bt(syr,nyr+1);	
	vector rt(syr,nyr);
	vector ft(syr,nyr);	
	vector epsilon(syr,nyr);
	
	vector nll(1,8);		

	sdreport_number sd_dep;
PROCEDURE_SECTION

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
	cout<<q<<endl;
	

	// initialize_model();
	// population_dynamics();
	// observation_model();
	calc_objective_function();
	


///
/// @brief Reference points check
/// @author Steve Martell
/// @remarks Checking my calculus.
///
FUNCTION void calc_Reference_Points()
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

///
/// @brief A simple function
/// @author Steve Martell
/// @remarks Yet to be documented.
///
FUNCTION void initialize_model()
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





///
/// @brief Population dynamics for lagged recruit growth survival model
/// @author Steve Martell
/// @remarks
/// \f$ b_{t+1} = s*b_{t}+r_{t} - c_{t} \f$	
/// Taken from the Ecological Detective
///
FUNCTION void population_dynamics()
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

///
/// @brief Observation model
/// @author Steve Martell
/// @remarks
/// Assumes observation errors are log-normal and based on the 
/// conditional maximum likelihood estimate of q.
///
FUNCTION void observation_model()
	int i;
	dvar_vector zt = log(it) - log(bt(syr,nyr));
	q              = exp(mean(zt));
	epsilon        = zt - mean(zt);

///
/// @brief WTF
/// @author Steve Martell
/// @remarks Runs Management Strategy Evaluation based on class OperatingModel
///
FUNCTION void run_mse()
	int j;
	
	Scenario cScenario1(agek,nScenario,n_pyr,rseed,value(bo),value(h),value(s),
	                    value(q),value(sig),value(tau),value(ft),
	                    value(wt),it,ct);

	int e_hcr = n_hcr;
	HarvestControlRule c_hcr(e_hcr);
	EstimatorClass cEstimator(sEstimator);


	OperatingModel cOM(cScenario1,cEstimator,c_hcr);
	cOM.runMSEscenario(cScenario1);

	ofstream ofs("OM.rep",ios::app);
	ofs<<"t_bo\n"  << cOM.get_bo()   <<endl;
	ofs<<"t_bmsy\n"<< cOM.get_bmsy() <<endl;
	ofs<<"t_fmsy\n"<< cOM.get_fmsy() <<endl;
	ofs<<"t_msy\n" << cOM.get_msy()  <<endl;
	ofs<<"t_bt\n"  << cOM.get_bt()   <<endl;

	ofs.close();



///
/// @brief Objective Function	
/// @author Steve Martell
/// @remarks Based on the negative log likelihood
///
FUNCTION void calc_objective_function()
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







///
/// @brief Management Strategy Evaluation
/// @author Steve Martell
/// @remarks Based on OperatingModel class
///
FUNCTION void wtf()
	int i;
	Scenario cScenario1(agek,nScenario,n_pyr,rseed,value(bo),value(h),value(s),
	                    value(q),value(sig),value(tau),value(ft),
	                    value(wt),it,ct);
	exit(1);
	int e_hcr = n_hcr;
	HarvestControlRule c_hcr(e_hcr);

	EstimatorClass cEstimator(sEstimator);


	OperatingModel cOM(cScenario1,cEstimator,c_hcr);
	cOM.runMSEscenario(cScenario1);


	ofstream ofs("OM.rep",ios::app);
	ofs<<"t_bo\n"  << cOM.get_bo()   <<endl;
	ofs<<"t_bmsy\n"<< cOM.get_bmsy() <<endl;
	ofs<<"t_fmsy\n"<< cOM.get_fmsy() <<endl;
	ofs<<"t_msy\n" << cOM.get_msy()  <<endl;
	ofs<<"t_bt\n"  << cOM.get_bt()   <<endl;

	ofs.close();



REPORT_SECTION
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


TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 

GLOBALS_SECTION
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

	
FINAL_SECTION
	
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

