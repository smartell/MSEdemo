#ifndef OPERATINGMODEL_H
#define OPERATINGMODEL_H

#include <admodel.h>
#include "scenario.h"

class operatingModel
{
private:
	int m_syr;
	int m_nyr;

	
	Scenario m_cScenario;

	operatingModel();

public:
	operatingModel(const Scenario &cScenario)
	: m_cScenario(cScenario)
	{
		m_syr = cScenario.m_ft.indexmin();
		m_nyr = cScenario.m_ft.indexmax();
	}

	~operatingModel() {}

	// member functions
	void populationModel(const Scenario &cScenario);
};

#endif


void operatingModel::populationModel(const Scenario &cScenario)
{
	// This routine reconstructs the population dynamics based on the Scenario class
	int i;
	double ro, reck, a, b, z, m;
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
	m    = -log(s);

	dvector bt(m_syr,m_nyr+1);
	dvector rt(m_syr,m_nyr);
	dvector hat_ct(m_syr,m_nyr);
	bt(m_syr) = bo;
	rt = ro;
	for(i = m_syr; i <= m_nyr; i++)
	{
		cout<<i<<" "<<agek<<endl;
		if(i-m_syr > agek)
		{
			rt(i) = a*bt(i-agek)/(1.+b*bt(i-agek)) * exp(wt(i));	
		}
		z = ft(i) + m;
		hat_ct(i) = bt(i) *ft(i)/z *(1.-mfexp(-z));
		bt(i+1)   = s*bt(i) + rt(i) - hat_ct(i);
	}

	cout<<bt<<endl;
	cout<<bo<<endl;
}


