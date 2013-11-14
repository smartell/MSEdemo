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
	int ngear;
	int nEpochs;
	ivector nIt_nobs;
	dmatrix ct;
	dmatrix it;
	imatrix it_yr;
	imatrix epoch;
	dmatrix cv;
	dmatrix prior_controls;
	// friend class LRGS;
};

struct sLRGSparameters
{
	dvariable log_bo;
	dvariable log_b1;
	dvariable h;
	dvariable s;
	dvariable *gamma;
	dvariable log_sigma;
	dvariable log_tau;
	dvar_vector wt;
	dvar_vector log_beta;
};

class LRGS 
{
private:
	int m_syr;
	int m_nyr;
	int m_agek;
	int m_ngear;
	int m_nEpochs;
	ivector m_nIt_nobs;
	dvariable   m_bo;
	dvariable   m_b1;
	dvariable   m_h;
	dvariable   m_ro;
	dvariable   m_s;
	dvariable   m_gamma;
	dvariable   m_reck;
	dvariable   m_a;
	dvariable   m_b;
	dvariable   m_sig;
	dvariable   m_tau;
	dvariable   m_fpen;
	dmatrix     m_ct;
	dmatrix     m_it;
	imatrix     m_it_yr;
	imatrix     m_epoch;
	dmatrix     m_cv;
	dmatrix     m_prior_controls;
	dvar_vector m_wt;
	dvar_vector m_q;
	dvar_vector m_rt;
	dvar_vector m_bt;
	dvar_vector m_ft;
	dvar_vector m_delta;
	dvar_vector m_theta;    // vector of parameters in the control file with priors.
	dvar_vector m_prior_pdf;
	dvar_vector m_log_beta;
	dvar_matrix m_epsilon;
	dvar_matrix m_nll;

	
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
	     dmatrix& ct,
	     dmatrix& it,
	     imatrix& it_yr,
	     imatrix& epoch,
	     dmatrix& cv,
	     dvar_vector& wt);

	LRGS(sLRGSdata& data,sLRGSparameters& pars);

	void initialize_model();
	void population_dynamics();
	void observation_model();
	void observation_model_q_random_walk();
	void calc_negative_loglikelihoods();
	void calc_prior_densities();

	dvar_matrix get_epsilon()   {return m_epsilon;       }
	dvar_matrix get_nll()       {return m_nll;           }
	dvar_vector get_prior_pdf() {return m_prior_pdf;     }
	dvar_vector get_bt()        {return m_bt;            }
	dvar_vector get_ft()        {return m_ft;            }
	dvariable   get_depletion() {return m_bt(m_nyr)/m_bo;}
	dvariable   get_fpen()      {return m_fpen;          }
	dvar_vector get_q()         {return m_q;             }
	dvar_vector get_delta()     {return m_delta;         }
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