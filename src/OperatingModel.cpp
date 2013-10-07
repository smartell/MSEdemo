/**
 * \file OperatingModel.cpp
 * \author Steve Martell
**/


#include <admodel.h>
#include "OperatingModel.h"


// #include "MSYReferencePoints.h"
// #include "Scenario.h"

/** \brief Run the Operating model scenario
 *	
 *		This function runs the operating model conditional on the scenario class that is 
 *		passed as an arguement.
 *	
 *	\author  Steve Martell
 *	\date May 29, 2013
 *	\param  a Scenario class object.
 *	\return void
 *	\sa
**/
void OperatingModel::runMSEscenario(const Scenario &cScenario)
{
	// This routine reconstructs the population dynamics based on the Scenario class
	int i;
	double ro, reck, a, b;
	// [ ] TO DO, clean this up and declare this as private member variables.
	int   agek = m_agek;
	double bo  = m_bo;
	double h   = m_h;
	double s   = m_s;
	double q   = m_q;
	double sig = m_sig;
	double tau = m_tau;
	dvector ft = m_ft;
	dvector wt = m_wt;
	dvector it = m_it;

	ro   = bo*(1.-s);
	reck = 4*h/(1.-h);
	a    = reck*ro/bo;
	b    = (reck-1.0)/bo;

	// |---------------------------------------------------------------------------------|
	// | TRUE REFERENCE POINTS
	// |---------------------------------------------------------------------------------|
	// |
	msy_reference_points cRefPoints(reck,s,bo);
	m_fmsy = cRefPoints.get_fmsy();
	m_bmsy = cRefPoints.get_bmsy();
	m_msy  = cRefPoints.get_msy();

	dvector bt(m_syr,m_nyr+m_pyr+1);
	dvector rt(m_syr,m_nyr+m_pyr);
	dvector hat_ct(m_syr,m_nyr+m_pyr);
	dvector hat_it(m_syr,m_nyr+m_pyr);

	// |---------------------------------------------------------------------------------|
	// | CONDITION REFERENCE MODEL BASED ON SCENARIO INFORMATION
	// |---------------------------------------------------------------------------------|
	// |
	bt.initialize();
	bt(m_syr) = bo;
	rt(m_syr,m_syr+agek) = ro * exp(wt(m_syr,m_syr+agek));
	hat_it(m_syr,m_nyr)  = it;
	for(i = m_syr; i <= m_nyr; i++)
	{

		if(i-m_syr > agek)
		{
			rt(i) = a*bt(i-agek)/(1.+b*bt(i-agek)) * exp(wt(i));	
		}
		// hat_ct(i) = bt(i) * (1.-mfexp(-ft(i)));
		hat_ct(i) = bt(i) * ft(i);
		bt(i+1)   = s*bt(i) + rt(i) - hat_ct(i);
	}
	hat_ct(m_syr,m_nyr) = m_ct;
	// | END OF CONDITIONING PERIOD

	// |---------------------------------------------------------------------------------|
	// | PROJECTION PERIOD
	// |---------------------------------------------------------------------------------|
	// | -1. Set Random numbers based on random number seed.
	// | -2. Calculate TAC based on harvest control rule and biomass estimate.
	// | -3. Implement harvest on reference population.
	// | -4. Update reference population.
	// | -5. Update observation models & write data files.
	// | -6. Conduct stock assessment & update Bo, reck and s parameters.
	// | -7. Repeat steps 2-7 for pyr's 
	// | -8. Compute performance measures.
	// | -9. Added sin curve to mimic non-stationarity
	// |
	int pyr1 = m_nyr+1;
	int pyr2 = m_nyr+m_pyr;
	dvector pyr(pyr1,pyr2);
	pyr.fill_seqadd(pyr1,1);
	int seed = m_rng;
	
	random_number_generator rng(seed);
	dvector rt_dev(pyr1,pyr2);
	dvector it_dev(pyr1,pyr2);
	dvector pdo_dev(pyr1,pyr2);

	// | Recruitment Scenario
	double lambda;
	if( m_nScenario==1 )
	{
		lambda = 0;
	}
	else if( m_nScenario==2 )
	{
		lambda = 0.4;
	}
	pdo_dev = sin((pyr-double(pyr1))/double(pyr2-pyr1)*4.0*PI);


	rt_dev.fill_randn(rng);
	it_dev.fill_randn(rng);

	rt_dev = rt_dev*tau - 0.5*tau*tau;
	it_dev = it_dev*sig - 0.5*sig*sig;

	double fmsy,bmsy,msy, tac, frate;
	double est_bo   = bo;
	double est_reck = reck;
	double est_s    = s;
	double est_bt   = bt(pyr1);
	

	for( i = pyr1; i <= pyr2; i++ )
	{
		// -2. Calculate TAC based on HCR & Biomass estimate.
		msy_reference_points cRefPoints(est_reck,est_s,est_bo);
		fmsy = cRefPoints.get_fmsy();
		msy  = cRefPoints.get_msy();
		bmsy = cRefPoints.get_bmsy();
		m_bmsy = bmsy;
		m_msy  = msy;
		m_fmsy = fmsy;
		tac  = m_cHCR.getTac(est_bt,fmsy,msy,bmsy,est_bo);

		// -3. Implement harvest on reference population, add implentation errors
		//     Watch out here if tac > bt(i), then need to set frate to some arbitrary max

		if( tac < bt(i) )
		{
			frate = -log((-tac+bt(i))/bt(i));
		}
		else
		{
			frate = 0.8;
		}
		hat_ct(i) = bt(i) * (1.-mfexp(-frate));


		// -4. Update reference population
		if(i-m_syr > agek)
		{
			rt(i) = a*bt(i-agek)/(1.+b*bt(i-agek)) * exp(rt_dev(i) + lambda*pdo_dev(i));	
		}
		bt(i+1) = s*bt(i) + rt(i) - hat_ct(i);

		// -5. Update observation models and write data files.
		// cout<<"Q = "<<q<<" "<<m_q<<endl;
		hat_it(i) = q*bt(i)*exp(it_dev(i));
		write_data_file(i,hat_ct(m_syr,i),hat_it(m_syr,i));

		// -6. Conduct assessment and update parameters
		// system("./OM -ind MSE.dat -nox -est > NUL");
		// If using perfect information, then don't run assessment.
		if(m_flg_perfect_information)
		{
			m_cEstimator.runEstimator();

			ifstream ifs("mse.par");
			ifs>>est_bo;
			ifs>>est_reck;
			ifs>>est_s;
			ifs>>est_bt;
			m_est_bo = est_bo;
		}
		else if( !m_flg_perfect_information )
		{
			cout<<"|-- PERFECT INFORMATION --|"<<endl;
			est_bo = m_bo;
			est_reck = m_reck;
			est_s = m_s;
			est_bt = bt(i+1);
			m_est_bo = est_bo;
		}


		// -   Screen dump so you can watch the progress.
		cout<<setprecision(3);
		cout<<"|---------------------------|"<<endl;
		cout<<"| - Year      "<<i<<endl;
		cout<<"|---------------------------|"<<endl;
		cout<<"| - est Bo    "<<est_bo<<endl;
		cout<<"| - est Bt    "<<est_bt<<endl;
		cout<<"| - bmsy      "<<bmsy<<endl;
		cout<<"| - bt(i)     "<<bt(i)<<endl;
		cout<<"| - fmsy      "<<fmsy<<endl;
		cout<<"| - frate     "<<frate<<endl;
		cout<<"| - msy       "<<msy<<endl;
		cout<<"| - tac       "<<tac<<endl;
		cout<<"| - hat_ct(i) "<<hat_ct(i)<<endl;
		cout<<"|---------------------------|"<<endl;
		
	} // 7. repeat steps 2-6
	m_bt = bt(m_syr,m_nyr+m_pyr);
	m_hat_ct = hat_ct;

	// 8. Calculate performance measures
	// AAV = |Ct - Ct-1| / Ct
	dvector AAV(m_syr,pyr2);
	AAV.initialize();
	for( i = pyr1+1; i <= pyr2; i++ )
	{
		 AAV(i)  = sum(fabs(first_difference(hat_ct(pyr1,i))));
		 AAV(i) /= sum(hat_ct(pyr1,i));
		 // cout<<AAV(i)<<endl;
	}
	m_AAV = AAV;
	
}

void OperatingModel::write_data_file(const int &nyr, const dvector &ct,const dvector& it)
{
	ofstream ofs("MSE.dat");
	ofs<<m_agek<<endl;
	ofs<<m_syr<<endl;
	ofs<<nyr<<endl;
	ivector iyr(m_syr,nyr);
	iyr.fill_seqadd(m_syr,1);
	ofs<<iyr<<endl;
	ofs<<ct<<endl;
	ofs<<it<<endl;

}

