/**
 * \file OperatingModel.h
 * \author Steve Martell
**/

#ifndef OPERATING_MODEL_H
#define OPERATING_MODEL_H
/**
 * This class represents an operating model to be used in Management Strategy Evaluation.
 * 
 * @file OperatingModel.h
 * @author Steven Martell <stevem@iphc.int>
 * 
 */

#include <admodel.h>
#include "Scenario.h"
#include "MSYReferencePoints.h"
#include "HarvestControlRule.h"
#include "EstimatorClass.h"

class OperatingModel 
{
private:
	int m_syr;
	int m_nyr;
	int m_pyr;
	int m_rng;
	int m_nScenario;
	int m_flg_perfect_information;

	int     m_agek;
	double  m_bo;
	double  m_est_bo;
	double  m_h;
	double  m_reck;
	double  m_s;
	double  m_q;
	double  m_sig;
	double  m_tau;
	double  m_fmsy;
	double  m_bmsy;
	double  m_msy;
	double  m_iuu_rate;
	double  m_mintac;
	dvector m_ft;
	dvector m_wt;
	dvector m_it;
	dvector m_ct;
	dvector m_bt;
	dvector m_hat_ct;
	dvector m_AAV;


	Scenario m_cScenario;
	HarvestControlRule m_cHCR;
	EstimatorClass	m_cEstimator;

	OperatingModel();

public:
	// OperatingModel(Scenario &cScenario, EstimatorClass &cEstimator, const HarvestControlRule &cHCR)
	// : m_cScenario(cScenario),m_cEstimator(cEstimator),m_cHCR(cHCR)
	// {
	// 	dvector ft = cScenario.get_ft();
	// 	m_syr = ft.indexmin();
	// 	m_nyr = ft.indexmax();
	// 	m_pyr = cScenario.get_pyr();
	// 	m_rng = cScenario.get_rseed();

	// 	m_agek      = cScenario.get_agek();
	// 	m_nScenario = cScenario.get_nScenario();
	// 	m_bo        = cScenario.get_bo();
	// 	m_h         = cScenario.get_h();
	// 	m_s         = cScenario.get_s();
	// 	m_q         = cScenario.get_q();
	// 	m_sig       = cScenario.get_sig();
	// 	m_tau       = cScenario.get_tau();
	// 	m_ft        = cScenario.get_ft();
	// 	m_wt        = cScenario.get_wt();
	// 	m_it        = cScenario.get_it();
	// 	m_ct        = cScenario.get_ct();
	// }

	OperatingModel(Scenario &cScenario, EstimatorClass &cEstimator, const HarvestControlRule &cHCR)
	: m_cScenario(cScenario),m_cEstimator(cEstimator),m_cHCR(cHCR)
	{
		dvector ft = cScenario.get_ft();
		m_syr = ft.indexmin();
		m_nyr = ft.indexmax();
		m_pyr = cScenario.get_pyr();
		m_rng = cScenario.get_rseed();

		m_agek      = cScenario.get_agek();
		m_nScenario = cScenario.get_nScenario();
		m_flg_perfect_information = cScenario.get_nInformation();
		m_bo        = cScenario.get_bo();
		m_h         = cScenario.get_h();
		m_reck      = (4.0*m_h)/(1.-m_h);
		m_s         = cScenario.get_s();
		m_q         = cScenario.get_q();
		m_sig       = cScenario.get_sig();
		m_tau       = cScenario.get_tau();
		m_ft        = cScenario.get_ft();
		m_wt        = cScenario.get_wt();
		m_it        = cScenario.get_it();
		m_ct        = cScenario.get_ct();
		m_iuu_rate  = cScenario.get_iuu();
		m_mintac    = cScenario.get_mintac();
	}

	// destructor:
	~OperatingModel() {}

	
	// getters:
	double  get_bmsy()   { return m_bmsy;   }
	double  get_fmsy()   { return m_fmsy;   }
	double  get_msy()    { return m_msy;    }
	double  get_bo()     { return m_bo;     }
	double  get_est_bo() { return m_est_bo; }
	dvector get_bt()     { return m_bt;     }
	dvector get_hat_ct() { return m_hat_ct; }
	dvector get_aav()    { return m_AAV;    }



	// member functions
	void runMSEscenario(const Scenario &cScenario);
	void write_data_file(const int &nyr, const dvector &ct,const dvector& it);
};

#endif


