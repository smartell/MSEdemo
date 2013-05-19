#ifndef OPERATING_MODEL_H
#define OPERATING_MODEL_H

#include <admodel.h>
#include "scenario.h"
#include "msyrefPoints.h"
#include "harvestControlRule.h"
/**
 * This class represents an operating model to be used in Management Strategy Evaluation.
 * 
 * @file OperatingModel.h
 * @author Steven Martell <stevem@iphc.int>
 * 
 */
class OperatingModel
{
private:
	int m_syr;
	int m_nyr;
	int m_pyr;
	int m_rng;

	int     m_agek;
	double  m_bo;
	double  m_h;
	double  m_s;
	double  m_q;
	double  m_sig;
	double  m_tau;
	dvector m_ft;
	dvector m_wt;
	dvector m_it;
	dvector m_ct;
	dvector m_bt;


	Scenario m_cScenario;
	HarvestControlRule m_HCR;

	OperatingModel();

public:
	OperatingModel(const Scenario &cScenario, const HarvestControlRule &cHCR)
	: m_cScenario(cScenario), m_HCR(cHCR)
	{
		m_syr = cScenario.m_ft.indexmin();
		m_nyr = cScenario.m_ft.indexmax();
		m_pyr = cScenario.m_pyr;
		m_rng = cScenario.m_rseed;

		m_agek = cScenario.m_agek;
		m_bo   = cScenario.m_bo;
		m_h    = cScenario.m_h;
		m_s    = cScenario.m_s;
		m_q    = cScenario.m_q;
		m_sig  = cScenario.m_sig;
		m_tau  = cScenario.m_tau;
		m_ft   = cScenario.m_ft;
		m_wt   = cScenario.m_wt;
		m_it   = cScenario.m_it;
		m_ct   = cScenario.m_ct;
	}

	// destructor:
	~OperatingModel() {}

	// getters:
	dvector get_bt() { return m_bt; }

	// member functions
	void runMSEscenario(const Scenario &cScenario);
	void write_data_file(const int &nyr, const dvector &ct,const dvector& it);
};

#endif


// Put this code in a cpp file eventually.
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
		hat_ct(i) = bt(i) * (1.-mfexp(-ft(i)));
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
	// |
	int pyr1 = m_nyr+1;
	int pyr2 = m_nyr+m_pyr;
	int seed = m_rng;
	
	random_number_generator rng(seed);
	dvector rt_dev(pyr1,pyr2);
	dvector it_dev(pyr1,pyr2);

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
		tac = m_HCR.getTac(est_bt,fmsy,msy,bmsy,est_bo);

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
			rt(i) = a*bt(i-agek)/(1.+b*bt(i-agek)) * exp(rt_dev(i));	
		}
		bt(i+1) = s*bt(i) + rt(i) - hat_ct(i);

		// -5. Update observation models and write data files.
		hat_it(i) = q*bt(i)*exp(it_dev(i));
		write_data_file(i,hat_ct(m_syr,i),hat_it(m_syr,i));

		// -6. Conduct assessment and update parameters
		system("./OM -ind MSE.dat -nox -est > NUL");
		ifstream ifs("mse.par");
		ifs>>est_bo;
		ifs>>est_reck;
		ifs>>est_s;
		ifs>>est_bt;


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
	m_bt = bt;

	// 8. Calculate performance measures
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

