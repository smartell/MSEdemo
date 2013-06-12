#include <admodel.h>

#ifndef _LRGS_H
#define _LRGS_H
/**
 * \breif Lagged Recruitment Growth Survival Model
 * \author Steve Martell
**/

class LRGS
{
private:
	int m_syr;
	int m_nyr;
	int m_agek;
	dvariable   m_bo;
	dvariable   m_h;
	dvariable   m_ro;
	dvariable   m_s;
	dvariable   m_reck;
	dvariable   m_a;
	dvariable   m_b;
	dvariable   m_sig;
	dvariable   m_tau;
	dvariable   m_fpen;
	dvariable   m_q;
	dvector     m_ct;
	dvector     m_it;
	dvar_vector m_wt;
	dvar_vector m_rt;
	dvar_vector m_bt;
	dvar_vector m_ft;
	dvar_vector m_epsilon;
	

public:
	LRGS(const int& syr,
	     const int& nyr,
	     const int& agek,
	     dvariable& bo,
	     dvariable& h, 
	     dvariable& s, 
	     dvariable& sig, 
	     dvariable& tau,
	     dvector& ct,
	     dvector& it,
	     dvar_vector& wt);
	~LRGS(){};


	void initialize_model();
	void population_dynamics();
	void observation_model();

	dvar_vector get_epsilon() {return m_epsilon;}
	dvariable   get_depletion() {return m_bt(m_nyr)/m_bo;}
	/* data */
};

#endif