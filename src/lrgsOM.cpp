#include <admodel.h>
#include "lrgsOM.h"


lrgsOM::~lrgsOM()
{}

lrgsOM::lrgsOM()
{
	cout<<"Default constructor"<<endl;
	// readScenarioInput("S1.scn");
}

lrgsOM::lrgsOM(const adstring s_file)
{
	readScenarioInput(s_file);
}

lrgsOM::lrgsOM(const sLRGSdata& data,const sLRGSparameters& pars,
               const adstring s_Scfile,const adstring s_Mpfile)
: m_syr(data.syr),
  m_nyr(data.nyr),
  m_agek(data.agek),
  m_rseed(data.rseed),
  m_ct(data.ct),
  m_it(data.it),
  m_log_bo(value(pars.log_bo)),
  m_h(value(pars.h)),
  m_s(value(pars.s)),
  m_log_sigma(value(pars.log_sigma)),
  m_log_tau(value(pars.log_tau)),
  m_wt(value(pars.wt))
{
	readScenarioInput(s_Scfile);

	readProcedureInput(s_Mpfile);

	conditionOperatingModel();

	runOperatingModel();

	cout<<"Start year       "<<m_syr<<endl;
	cout<<"Projection yeear "<<m_pyr<<endl;
	cout<<"Age k            "<<m_agek<<endl;
	cout<<"Threshold dep.   "<<m_threshold<<endl;
	// cout<<m_wt<<endl;
}


void lrgsOM::readScenarioInput(const adstring s_file)
{
	cout<<s_file<<endl;
	cifstream ifs(s_file);
	ifs >> m_pyr;
}

void lrgsOM::readProcedureInput(const adstring s_file)
{
	cifstream ifs(s_file);
	ifs >> m_hr;
	ifs >> m_limit;
	ifs >> m_threshold;
	ifs >> m_target;
	ifs >> m_ftarget; // Should correspond to fishing rate at depletion target.
	ifs >> m_fmax;    // Should correspond to fishing rate at threshold.
	ifs >> m_catch_floor;
}

void lrgsOM::conditionOperatingModel()
{
	int i;
	// Population variables.
	m_bo   = exp(m_log_bo);
	m_ro   = m_bo * (1.0-m_s);
	m_reck = 4.0*m_h/(1.0-m_h);
	m_sig  = exp(m_log_sigma);
	m_tau  = exp(m_log_tau);
	m_so   = m_reck * (1.0 - m_s);
	m_beta = (m_reck - 1.0)/m_bo;

	// Population dynamics
	m_chat.allocate(m_syr,m_pyr);   m_chat.initialize();
	m_ihat.allocate(m_syr,m_pyr);   m_ihat.initialize();
	m_bt.allocate(m_syr,m_pyr+1);   m_bt.initialize();
	m_rt.allocate(m_syr,m_pyr);     m_rt.initialize();
	m_ft.allocate(m_syr,m_pyr);     m_ft.initialize();

	m_bt(m_syr) = m_bo;
	m_rt(m_syr,m_syr+m_agek) = m_ro * exp(m_wt(m_syr,m_syr+m_agek));
	for( i = m_syr; i <= m_nyr; i++ )
	{
		if( i-m_syr > m_agek )
		{
			m_rt(i)  = m_so * m_bt(i-m_agek) / (1.0 + m_beta*m_bt(i-m_agek)); 
			m_rt(i) *= exp(m_wt(i));
		}
		m_bt(i+1) = m_s*m_bt(i) + m_rt(i) - m_ct(i);
	}
	// Verified biomass with assessment model.  SJDM

	m_chat(m_syr,m_nyr) = m_ct(m_syr,m_nyr);
	m_ihat(m_syr,m_nyr) = m_it(m_syr,m_nyr);
	m_ft(m_syr,m_nyr)   = elem_div(m_ct(m_syr,m_nyr),m_bt(m_syr,m_nyr));
	
	// Catchability coefficient
	m_q = exp(mean( log(m_it) - log(m_bt(m_syr,m_nyr)) ));
	
	cout<<"q\n"<<m_q<<endl;
}

void lrgsOM::runOperatingModel()
{
	// |------------------------------------------|
	// | Run the operating model into the future. |
	// |------------------------------------------|
	// | - Set Random numbers
	// | - Calcaulate Reference points
	// | - Calculate TAC based on harvest control rule.
	// | - Implement fisheries
	// | - Update population dynamics
	// | - Catch and effort statistics for reporting
	// | - Write data file
	// | - Conduct stock assessment.

	int i;
	// | - Set Random numbers
	random_number_generator rng(m_rseed);
	dvector rt_dev(m_nyr+1,m_pyr);
	dvector it_dev(m_nyr+1,m_pyr);

	rt_dev.fill_randn(rng);
	it_dev.fill_randn(rng);


	rt_dev = rt_dev*m_tau - 0.5*m_tau*m_tau;
	it_dev = it_dev*m_sig - 0.5*m_sig*m_sig;
	cout<<"Biomass limit = "<<m_limit<<endl;

	HarvestControlRule cTAC(m_limit,m_threshold,m_target);

	for( i = m_nyr; i <= m_pyr; i++ )
	{
		// | - Calculate reference points
		calcReferencePoints();

		// | - Calculate TAC based on harvest control rule parameters.
		m_chat(i) = cTAC(m_bt(i));	
	}
	cout<<m_chat<<endl;
	
	
}

void lrgsOM::calcReferencePoints()
{
	m_bmsy = (m_bo*(-1.+sqrt(m_reck))/(m_reck-1.)); 
	m_msy  = (m_bmsy * (-1. + m_s) * (-1. + m_reck)*(-m_bo + m_bmsy)) 
			 / (m_bo + (-1. + m_reck)*m_bmsy);
	
	m_fmsy = m_msy / m_bmsy;
	//cout<<"MSY\n"<<m_msy<<endl;
}



// HARVEST CONTROL RULE FUNCTIONS
HarvestControlRule::~HarvestControlRule()
{}

HarvestControlRule::HarvestControlRule(const double &limit, const double &threshold,
                                       const double &target)
{
	cout<<"Limit"<<limit<<endl;
}

double HarvestControlRule::operator( )(const double& bt)
{
	return(4);
}
