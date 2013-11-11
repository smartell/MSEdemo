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
	init_int ngear;
	init_int nEpochs;
	init_ivector nIt_nobs(1,nEpochs);
	init_matrix catch_data(syr,nyr,0,ngear);
	init_3darray   it_data(1,nEpochs,1,nIt_nobs,1,4);
	init_int eof;
	!! cout<<eof<<endl;
	!! if(eof != -999){cout<<"Error reading data"<<endl; exit(1);}

	ivector iyr(syr,nyr);
	matrix   ct(syr,nyr,1,ngear);
	imatrix it_yr(1,nEpochs,1,nIt_nobs);
	imatrix epoch(1,nEpochs,1,nIt_nobs);
	matrix     it(1,nEpochs,1,nIt_nobs);
	matrix     cv(1,nEpochs,1,nIt_nobs);


	LOC_CALCS
		iyr   = ivector(column(catch_data,0));
		ct    = trans(trans(catch_data).sub(1,ngear));
		for(int i = 1; i <= nEpochs; i++ )
		{
			it_yr(i) = ivector(column(it_data(i),1));
			epoch(i) = ivector(column(it_data(i),2));
			it(i)    = column(it_data(i),3);
			cv(i)    = column(it_data(i),4);		
		}
	END_CALCS

	// |---------------------------------------------------------------------------------|
	// | Read control file for turning on and off parameters & setting priors            |
	// |---------------------------------------------------------------------------------|

	!! ad_comm::change_datafile_name("om.ctl");
	int npar;
	!!  npar=8;
	init_matrix   d_PC(1,8,1,7);
	

	int ir;
	int iphz;

	number dval;
	number dlb;
	number dub;


	// |---------------------------------------------------------------------------------|
	// | MANAGEMENT STRATEGY EVALUATION COMMANDS
	// |---------------------------------------------------------------------------------|
	// |
	// Scenario 1 = stationarity
	// Scenario 2 = pdo
	//init_int nScenario;
	//// Harvest control rule
	//init_int n_hcr;
	//init_int n_pyr;
	//init_int n_flg_perfect_information;
	//init_number iuu_rate;
	//init_number min_tac;
	//init_adstring sEstimator;
	// !! COUT(sEstimator);
	// !! exit(1);
	// !! sLRGSdata data;
	!! data.syr            = syr;
	!! data.nyr            = nyr;
	!! data.agek           = agek;
	!! data.ngear          = ngear;
	!! data.nIt_nobs       = nIt_nobs;
	!! data.ct             = ct;
	!! data.it             = it;
	!! data.rseed          = rseed;
	!! data.it_yr          = it_yr;
	!! data.nEpochs        = nEpochs;
	!! data.epoch          = epoch;
	!! data.cv             = cv;
	!! data.prior_controls = d_PC;
	!! cout<<data.it<<endl;
	
INITIALIZATION_SECTION
	log_bo     8.0;
	h          0.9;
	s          0.9;
	log_sigma  4.65;
	log_tau    4.65;
	

PARAMETER_SECTION
	!! ir=1; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
	init_bounded_number log_bo(dlb,dub,iphz);
	!! log_bo = dval;

	!! ir=2; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
	init_bounded_number log_b1(dlb,dub,iphz);
	!! log_b1 = dval;

	!! ir=3; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
	init_bounded_number h(dlb,dub,iphz);
	!! h = dval;

	!! ir=4; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
	init_bounded_number s(dlb,dub,iphz);
	!! s = dval;

	!! ir=5; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
	init_bounded_number gamma(dlb,dub,iphz);
	!! gamma = dval;

	!! ir=6; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
	init_bounded_number log_sigma(dlb,dub,iphz);
	!! log_sigma = dval;

	!! ir=7; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
	init_bounded_number log_tau(dlb,dub,iphz);
	!! log_tau = dval;

	!! ir=8; dval=d_PC(ir,1); dlb=d_PC(ir,2); dub=d_PC(ir,3); iphz=int(d_PC(ir,4));
	init_bounded_dev_vector wt(syr,nyr,dlb,dub,iphz);
	!! wt = dval;
	

	objective_function_value f;

	number bo;
	number b1;	
	number ro;	
	number sig;
	number tau;	
	number a;	
	number b;	
	number reck;
	number fpen;

	vector q(1,nEpochs);	
	vector bt(syr,nyr+1);	
	vector rt(syr,nyr);
	vector ft(syr,nyr);	
	vector delta(syr,nyr);
	matrix epsilon(1,nEpochs,1,nIt_nobs);
	matrix negloglike(1,2,1,nEpochs);

	
	vector nll(1,8);
	vector prior(1,8);

	sdreport_number sd_dep;
PROCEDURE_SECTION
	
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
	cLRGSmodel.observation_model_q_random_walk();
	cLRGSmodel.calc_negative_loglikelihoods();
	cLRGSmodel.calc_prior_densities();

	epsilon    = cLRGSmodel.get_epsilon();
	sd_dep     = cLRGSmodel.get_depletion();
	bt         = cLRGSmodel.get_bt();	
	ft         = cLRGSmodel.get_ft();
	q          = cLRGSmodel.get_q();
	negloglike = cLRGSmodel.get_nll();
	prior      = cLRGSmodel.get_prior_pdf();
	delta      = cLRGSmodel.get_delta();
	fpen       = cLRGSmodel.get_fpen();
	//cout<<"Fpen " <<fpen<<endl;
	
	calc_objective_function();
	
///
/// @brief Objective Function	
/// @author Steve Martell
/// @remarks Based on the negative log likelihood
///
FUNCTION void calc_objective_function()
	nll.initialize();
	//dvariable isig2 = mfexp(log_sigma);
	//dvariable itau2 = mfexp(log_tau);

	// No comments
	for(int i = 1; i <= nEpochs; i++ )
	{
		//nll(1) += dnorm(epsilon(i),sig);
	}
	//nll(2) = dbeta(s,30.01,10.01);
	//The following is based on E(x) = exp(-0.15), Sig2 = (.15*CV)^2, where assumed CV=0.02
	//nll(2) = dbeta(s,347.3694,56.21626);
	//nll(3) = dnorm(log_bo,log(500),5.0);
	//nll(3)+= dnorm(log_b1,log(500),5.0);
	//nll(4) = dbeta((h-0.2)/0.8,1.01,1.01);
	//nll(5) = dgamma(isig2,1.01,1.01);
	if(active(log_tau))
	{
		//nll(6) = dgamma(itau2,1.01,1.01);
		//nll(7) = dnorm(wt,tau);
	}
	else
	{
		//nll(7) = dnorm(wt,1.0);
	}


	
	if(fpen>0 && !mc_phase()) cout<<"Fpen = "<<fpen<<endl;
	f = 100000.*fpen + sum(negloglike) + sum(prior);


FUNCTION void mse2()
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


///
/// @brief WTF
/// @author Steve Martell
/// @remarks Runs Management Strategy Evaluation based on class OperatingModel
///
FUNCTION void run_mse()
//	int j;
//	
//	Scenario cScenario1(agek,nScenario,n_pyr,n_flg_perfect_information,rseed,value(bo),
//	                    value(h),value(s),iuu_rate,
//	                    value(q),value(sig),value(tau),value(ft),
//	                    value(wt),it,ct,min_tac);
//
//	int e_hcr = n_hcr;
//	HarvestControlRule c_hcr(e_hcr);
//	EstimatorClass cEstimator(sEstimator);
//
//
//	OperatingModel cOM(cScenario1,cEstimator,c_hcr);
//	cOM.runMSEscenario(cScenario1);
//
//	ofstream ofs("OM.rep",ios::app);
//	ofs<<"t_bo\n"      << cOM.get_bo()       <<endl;
//	ofs<<"t_est_bo\n"  << cOM.get_est_bo()   <<endl;
//	ofs<<"t_bmsy\n"    << cOM.get_bmsy()     <<endl;
//	ofs<<"t_fmsy\n"    << cOM.get_fmsy()     <<endl;
//	ofs<<"t_msy\n"     << cOM.get_msy()      <<endl;
//	ofs<<"t_bt\n"      << cOM.get_bt()       <<endl;
//	ofs<<"t_aav\n"     << cOM.get_aav()      <<endl;
//
//	ofs.close();
//
//	// |--------------------------------------------|
//	// | NOW RUN THE MODEL WITH PERFECT INFORMATION |
//	// |--------------------------------------------|
//
//	Scenario cScenarioP(agek,nScenario,n_pyr,0,rseed,value(bo),
//	                    value(h),value(s),iuu_rate,
//	                    value(q),value(sig),value(tau),value(ft),
//	                    value(wt),it,ct);
//	
//
//	OperatingModel cOMP(cScenarioP,cEstimator,c_hcr);
//	cOMP.runMSEscenario(cScenarioP);
//
//	ofstream ofs1("OM.rep",ios::app);
//	ofs1<<"p_bo\n"  << cOMP.get_est_bo()   <<endl;
//	ofs1<<"p_bmsy\n"<< cOMP.get_bmsy()     <<endl;
//	ofs1<<"p_fmsy\n"<< cOMP.get_fmsy()     <<endl;
//	ofs1<<"p_msy\n" << cOMP.get_msy()      <<endl;
//	ofs1<<"p_bt\n"  << cOMP.get_bt()       <<endl;
//	ofs1<<"p_ct\n"  << cOMP.get_hat_ct()   <<endl;
//	ofs1<<"p_aav\n" << cOMP.get_aav()      <<endl;
//
//	ofs1.close();





REPORT_SECTION
	msy_reference_points cMSY(value(reck),value(s),value(bo));

	REPORT(bo);
	REPORT(b1);
	REPORT(h);
	REPORT(s);
	REPORT(q);
	REPORT(sig);
	REPORT(tau);
	REPORT(prior);
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
	REPORT(it_data);

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

	
FINAL_SECTION
	
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

