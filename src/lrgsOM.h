#ifndef HCR_H
#define HCR_H

class HarvestControlRule
{
private:

public:
	~HarvestControlRule();
	HarvestControlRule(const double &limit, const double &threshold,
                                       const double &target);

	double operator( )(const double& bt);
};

#endif



#ifndef _LRGS_OPERATING_MODEL_H
#define _LRGS_OPERATING_MODEL_H

#include <admodel.h>
#include "LRGS.h"



class lrgsOM
{
private:
	int m_syr;
	int m_nyr;
	int m_agek;
	int m_rseed;


	double m_hr;
	double m_log_bo;
	double m_h;
	double m_s;
	double m_log_sigma;
	double m_log_tau;

	dvector m_wt;
	dvector m_ct;
	dvector m_it;

	int m_pyr;
	double m_bo;
	double m_ro;
	double m_reck;
	double m_sig;
	double m_tau;
	double m_so;
	double m_beta;
	double m_q;

	dvector m_bt;
	dvector m_rt;
	dvector m_chat;
	dvector m_ft;
	dvector m_ihat;

	// MSY-based variables
	double m_fmsy;
	double m_bmsy;
	double m_msy;

	

	// Depeltion based referece variables
	double m_limit;
	double m_threshold;
	double m_target;
	double m_ftarget; // Should correspond to fishing rate at depletion target.
	double m_fmax;    // Should correspond to fishing rate at threshold.
	double m_catch_floor;
	
	

	lrgsOM();
public:
	~lrgsOM();
	lrgsOM(const adstring s_file);
	lrgsOM(const sLRGSdata& data,const sLRGSparameters& pars,
	       const adstring s_Scfile,const adstring s_Mpfile);

	

	void readScenarioInput(const adstring s_file);
	void readProcedureInput(const adstring s_file);
	void conditionOperatingModel();
	void runOperatingModel();
	void calcReferencePoints();
};


#endif

