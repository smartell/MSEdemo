#ifndef OPERATINGMODEL_H
#define OPERATINGMODEL_H

#include <admodel.h>
#include "scenario.h"

class operatingModel
{
private:
	int m_syr;
	int m_nyr;
	int m_pyr;

	dvector m_bt;

	Scenario m_cScenario;

	operatingModel();

public:
	operatingModel(const Scenario &cScenario)
	: m_cScenario(cScenario)
	{
		m_syr = cScenario.m_ft.indexmin();
		m_nyr = cScenario.m_ft.indexmax();
		m_pyr = cScenario.m_pyr;
	}

	~operatingModel() {}

	// member functions
	void populationModel(const Scenario &cScenario);
};

#endif


// Put this code in a cpp file eventually.
void operatingModel::populationModel(const Scenario &cScenario)
{
	// This routine reconstructs the population dynamics based on the Scenario class
	int i;
	double ro, reck, a, b;
	int   agek = cScenario.m_agek;
	double bo  = cScenario.m_bo;
	double h   = cScenario.m_h;
	double s   = cScenario.m_s;
	double sig = cScenario.m_sig;
	double tau = cScenario.m_tau;
	dvector ft = cScenario.m_ft;
	dvector wt = cScenario.m_wt;

	ro   = bo*(1.-s);
	reck = 4*h/(1.-h);
	a    = reck*ro/bo;
	b    = (reck-1.0)/bo;

	dvector bt(m_syr,m_nyr+m_pyr);
	dvector rt(m_syr,m_nyr+m_pyr);
	dvector hat_ct(m_syr,m_nyr+m_pyr);

	// |---------------------------------------------------------------------------------|
	// | CONDITION REFERENCE MODEL BASED ON SCENARIO INFORMATION
	// |---------------------------------------------------------------------------------|
	// |
	bt.initialize();
	bt(m_syr) = bo;
	rt(m_syr,m_syr+agek) = ro * exp(wt(m_syr,m_syr+agek));
	for(i = m_syr; i <= m_nyr; i++)
	{

		if(i-m_syr > agek)
		{
			rt(i) = a*bt(i-agek)/(1.+b*bt(i-agek)) * exp(wt(i));	
		}
		hat_ct(i) = bt(i) * (1.-mfexp(-ft(i)));
		bt(i+1)   = s*bt(i) + rt(i) - hat_ct(i);
	}

	// | END OF CONDITIONING PERIOD

	// |---------------------------------------------------------------------------------|
	// | PROJECTION PERIOD
	// |---------------------------------------------------------------------------------|
	// | -1. Set Random numbers based on random number seed.
	// | -2. Calculate TAC based on harvest control rule and biomass estimate.
	// | -3. Implement harvest on reference population.
	// | -4. Update reference population.
	// | -5. Update observation models & write data files.
	// | -6. Conduct stock assessment.
	// | -7. Repeat steps 2-7 for pyr's 
	// | -8. Compute performance measures.
	// |
	int pyr1 = m_nyr+1;
	int pyr2 = m_nyr+m_pyr;
	int seed = -12345;
	
	random_number_generator rng(seed);
	dvector rt_dev(pyr1,pyr2);
	dvector it_dev(pyr1,pyr2);

	rt_dev.fill_randn(rng);
	it_dev.fill_randn(rng);

	rt_dev = rt_dev*tau - 0.5*tau*tau;
	it_dev = it_dev*sig - 0.5*sig*sig;

	for( i = pyr1; i <= pyr2; i++ )
	{
		// -2. Calculate TAC based on HCR & Biomass estimate.
	}
	cout<<bt<<endl;
	cout<<bo<<endl;
}


