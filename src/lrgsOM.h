#ifndef HCR_H
#define HCR_H

class HarvestControlRule
{
private:
	double m_limit;
	double m_threshold;
	double m_target;
	double m_ftarget;
	double m_fmax;
	double m_catch_floor;

public:
	~HarvestControlRule();
	HarvestControlRule(const double &limit, const double &threshold,
	                   const double &target,const double &ftarget);
	HarvestControlRule(const double &limit, const double &threshold,
                       const double &target, const double &ftarget,
                       const double &fmax, const double &catch_floor);

	double operator( )(const double& bt,const double &bo);
	double operator( )(const double& bt,const double &bo, double& f_rate);
};

#endif



#ifndef _LRGS_OPERATING_MODEL_H
#define _LRGS_OPERATING_MODEL_H

#include <admodel.h>
#include "LRGS.h"
#include "EstimatorClass.h"



class lrgsOM
{
private:
	int m_syr;
	int m_nyr;
	int m_agek;
	int m_rseed;


	double m_log_bo;
	double m_h;
	double m_s;
	double m_log_sigma;
	double m_log_tau;

	dvector m_wt;
	dvector m_ct;
	dvector m_it;

	int m_pyr;
	double m_phi;
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
	dvector m_aav;

	// MSY-based variables
	double m_fmsy;
	double m_bmsy;
	double m_msy;

	// model name
	adstring m_sAssessmentModel;
	double m_est_bo;
	double m_est_reck;
	double m_est_s;
	double m_est_bt;
	double m_est_fmsy;
	double m_est_bmsy;
	double m_est_msy;

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

	
	double get_bo()     { return m_bo;     }
	double get_est_bo() { return m_est_bo; }
	double get_bmsy()   { return m_bmsy;   }
	double get_fmsy()   { return m_fmsy;   }
	double get_msy()    { return m_msy;    }
	dvector get_bt()     { return m_bt;     }
	dvector get_aav()    { return m_aav;    }

	void readScenarioInput(const adstring s_file);
	void readProcedureInput(const adstring s_file);
	void conditionOperatingModel();
	void runOperatingModel();
	void calcReferencePoints();
	void calcReferencePoints(const double &bo, const double & reck,
                                 const double &s, double &fmsy, double &bmsy, double &msy);
	void write_data_file(const int &nyr, const dvector &ct,const dvector& it);
	void read_parameter_estimates(const adstring &sParFile);
	void print_mse(const int &i);
};


#endif

