//  ******************************************************************
//  | OM, short for Operating model.
//  |
//  | Created by Martell on 2013-05-14.
//  | Copyright (c) 2013. All rights reserved.
//  | Comments: This operation model is conditioned on the LRGS model
//  |           based on chapter 10 in the Ecological Detective.
//  ******************************************************************


DATA_SECTION
	
	// |---------------------------------------------------------------------------------|
	// | COMMAND LINE OPTIONS
	// |---------------------------------------------------------------------------------|
	// |
	int on;
	int do_mse;
	!! if ( (on=option_match(argc,argv,"-mse"))>-1) do_mse=1;else do_mse=0;

	init_int agek;
	init_int syr;
	init_int nyr;
	init_ivector iyr(syr,nyr);
	init_vector ct(syr,nyr);
	init_vector it(syr,nyr);

INITIALIZATION_SECTION
	log_bo     8.0;
	h          0.9;
	s          0.9;
	log_sigma  4.65;
	log_tau    4.65;
	log_ft    -2.30;

PARAMETER_SECTION
	init_number log_bo;
	init_bounded_number h(0.2,1.0,2);
	init_bounded_number s(0.0,1.0);
	init_number log_sigma(3);
	init_number log_tau(3);
	// init_number mean_log_ft;
	init_bounded_vector log_ft(syr,nyr,-30,10,1);
	init_bounded_dev_vector wt(syr,nyr,-15,15,2);
	// random_effects_vector wt(syr,nyr,2);

	objective_function_value f;

	number bo;	
	number ro;	
	number sig;
	number tau;	
	number a;	
	number b;	
	number reck;
	number q;	

	vector bt(syr,nyr+1);	
	vector rt(syr,nyr);
	vector ft(syr,nyr);	
	vector hat_ct(syr,nyr);	
	vector epsilon(syr,nyr);
	vector nu(syr,nyr);
	vector nll(1,8);		

	sdreport_number sd_dep;

PROCEDURE_SECTION
	initialize_model();
	population_dynamics();
	observation_model();
	calc_objective_function();

FUNCTION initialize_model
	bo      = mfexp(log_bo);
	ro      = bo*(1.-s);
	reck    = 4*h/(1.-h);
	a       = reck*ro/bo;
	b       = (reck-1.0)/bo;
	rt      = ro * exp(wt);
	bt(syr) = bo;
	sig     = sqrt(1.0/mfexp(log_sigma));
	tau     = sqrt(1.0/mfexp(log_tau));
	ft      = mfexp( log_ft );

FUNCTION population_dynamics
	int i;
	dvariable m = -log(s);
	dvariable z;
	for(i=syr;i<=nyr;i++)
	{
		if(i-syr > agek)
		{
			rt(i) = a*bt(i-agek)/(1.+b*bt(i-agek)) * exp(wt(i));	
		}
		z = ft(i) + m;
		hat_ct(i) = bt(i) *ft(i)/z *(1.-mfexp(-z));
		bt(i+1)   = s*bt(i) + rt(i) - hat_ct(i);
	}
	sd_dep = bt(nyr)/bo;

FUNCTION observation_model
	double tiny = 1.0;
	dvar_vector zt = log(it) - log(bt(syr,nyr));
	q = exp(mean(zt));
	epsilon = zt - mean(zt);

	nu = log(ct+tiny)-log(hat_ct+tiny);

FUNCTION calc_objective_function
	nll.initialize();
	nll(1) = dnorm(epsilon,sig);
	// nll(2) = dbeta(s,15,3.321);
	// nll(3) = dbeta((h-0.2)/0.8,5.0,1.666);
	nll(2) = dbeta(s,1.01,1.01);
	nll(3) = dbeta((h-0.2)/0.8,1.01,1.01);
	nll(4) = dnorm(nu,0.01);
	// Phased penalty on mean fishing mortality rate.
	if(!last_phase())
	{
		nll(5) = dnorm(mean(log_ft),log(0.2),0.15);
	}
	else if(last_phase())
	{
		nll(5) = dnorm(mean(log_ft),log(0.2),2.00);
	}
	nll(6) = dnorm(wt,tau);
	dvariable isig2 = mfexp(log_sigma);
	dvariable itau2 = mfexp(log_tau);
	nll(7) = dgamma(isig2,1.01,1.01);
	nll(8) = dgamma(itau2,1.01,1.01);
	f = sum(nll);

REPORT_SECTION
	REPORT(bo);
	REPORT(h);
	REPORT(s);
	REPORT(sig);
	REPORT(tau);
	REPORT(ft);
	REPORT(wt);
	REPORT(bt);
	REPORT(nu);
	REPORT(epsilon);

TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 

GLOBALS_SECTION
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
	#include "OperatingModel.h"
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;


	
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

FUNCTION run_mse
	// Scenario class
	int pyr = 30;
	Scenario cScenario1(agek,pyr,value(bo),value(h),value(s),
	                    value(sig),value(tau),value(ft),
	                    value(wt));
	
	// Operating model class
	operatingModel cOM(cScenario1);
	cOM.populationModel(cScenario1);





