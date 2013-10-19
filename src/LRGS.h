#include <admodel.h>
#include "om.htp"
#ifndef _LRGS_H
#define _LRGS_H
/**
 * \breif Lagged Recruitment Growth Survival Model
 * \author Steve Martell
**/

struct sLRGSdata
{
	int syr;
	int nyr;
	int agek;
	int rseed;
	dvector ct;
	dvector it;
	// friend class LRGS;
};

struct sLRGSparameters
{
	dvariable log_bo;
	dvariable h;
	dvariable s;
	dvariable log_sigma;
	dvariable log_tau;
	dvar_vector wt;
};

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
	
	// sLRGSdata       m_data;
	// sLRGSparameters m_pars;
public:
	// friend class model_parameters;
	// friend class model_data;
	~LRGS(){};
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

	LRGS(sLRGSdata& data,sLRGSparameters& pars);

	void initialize_model();
	void population_dynamics();
	void observation_model();  

	dvar_vector get_epsilon()   {return m_epsilon;       }
	dvar_vector get_bt()        {return m_bt;            }
	dvar_vector get_ft()        {return m_ft;            }
	dvariable   get_depletion() {return m_bt(m_nyr)/m_bo;}
	dvariable   get_q()         {return m_q;             }
	/* data */
};

#endif

// class Test: private model_parameters
// {

// 	~Test(){};
// 	Test()
// 	{
// 		cout<<"IN the test class"<<endl;
// 		cout<<syr<<endl;
// 	}


// };