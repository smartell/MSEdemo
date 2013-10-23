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
  m_ngear(data.ngear),
  m_nEpochs(data.nEpochs),
  m_nIt_nobs(data.nIt_nobs),
  m_rseed(data.rseed),
  m_ct(data.ct),
  m_it(data.it),
  m_it_yr(data.it_yr),
  m_cv(data.cv),
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
	ifs >> m_phi;
	ifs >> m_sAssessmentModel;
	ifs >> m_limit;
	ifs >> m_threshold;
	ifs >> m_target;
	ifs >> m_ftarget; // Should correspond to fishing rate at depletion target.
	ifs >> m_fmax;    // Should correspond to fishing rate at threshold.
	ifs >> m_catch_floor;
	cout<<"Assessment Model"<<m_sAssessmentModel<<endl;


	// Constructor, where if ftarget==0, the calculate the corresponding
	// value of ftarget that results in btarget.
	if( m_ftarget == 0 )
	{
		double t1 = (m_reck-1.0);
		double t2 = (m_bo-m_target)*(-1.0+m_s)*t1;
		double t3 = m_bo+t1*m_target;
		m_ftarget = t2 / t3;
		cout<<"Limit =" <<m_limit<<endl;
		cout<<"The new Ftarget is = "<<m_ftarget<<endl;
	}
}

void lrgsOM::conditionOperatingModel()
{
	int i,j;
	// Population variables.
	m_bo   = exp(m_log_bo);
	m_ro   = m_bo * (1.0-m_s);
	m_reck = 4.0*m_h/(1.0-m_h);
	m_sig  = sqrt(1.0/mfexp(m_log_sigma)); 
	m_tau  = sqrt(1.0/mfexp(m_log_tau));
	m_so   = m_reck * (1.0 - m_s);
	m_beta = (m_reck - 1.0)/m_bo;
	// Population dynamics
	m_chat.allocate(m_syr,m_pyr,1,m_ngear);   m_chat.initialize();
	ivector extra = m_nIt_nobs+(m_pyr-m_nyr);
	m_ihat.allocate(1,m_nEpochs,1,extra);   m_ihat.initialize();
	m_cvhat.allocate(1,m_nEpochs,1,extra);   m_cvhat.initialize();
	m_it_yr_hat.allocate(1,m_nEpochs,1,extra); m_it_yr_hat.initialize();
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
		m_bt(i+1) = m_s*m_bt(i) + m_rt(i) - sum(m_ct(i));

		m_ft(i) = sum(m_ct(i)) / m_bt(i);
	}
	// Verified biomass with assessment model.  SJDM

	for( i = m_syr; i <= m_nyr; i++ )
	{
		for( j = 1; j <= m_ngear; j++ )
		{
			m_chat(i,j) = m_ct(i,j);
		}		
	}
	//m_ft(m_syr,m_nyr)   = elem_div(m_ct(m_syr,m_nyr),m_bt(m_syr,m_nyr));
	
	// Catchability coefficient
	m_q.allocate(1,m_nEpochs);
	for( i = 1; i <= m_nEpochs; i++ )
	{
		m_ihat(i)(1,m_nIt_nobs(i)) = m_it(i)(1,m_nIt_nobs(i));
		m_cvhat(i)(1,m_nIt_nobs(i)) = m_cv(i)(1,m_nIt_nobs(i));
		m_it_yr_hat(i)(1,m_nIt_nobs(i)) = m_it_yr(i)(1,m_nIt_nobs(i));

		ivector iyr = m_it_yr(i);
		dvector zt  = log(m_it(i)) - log(m_bt(iyr).shift(1));
		m_q(i)      = exp(mean(zt));
	}
	cout<<"Malloc"<<endl;
	// m_q = exp(mean( log(m_it) - log(m_bt(m_syr,m_nyr)) ));
	
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

	int i,j,k;
	// | - Set Random numbers
	random_number_generator rng(m_rseed);
	dvector rt_dev(m_nyr+1,m_pyr);
	dmatrix it_dev(1,m_nEpochs,m_nyr+1,m_pyr);
	dvector et_dev(m_nyr+1,m_pyr);

	rt_dev.fill_randn(rng);
	it_dev.fill_randn(rng);
	et_dev.fill_randn(rng);


	rt_dev = rt_dev*m_tau - 0.5*m_tau*m_tau;
	it_dev = it_dev*m_sig - 0.5*m_sig*m_sig;
	et_dev = et_dev*m_phi - 0.5*m_phi*m_phi;
	cout<<"Biomass limit = "<<m_limit<<endl;

	read_parameter_estimates("mse.par");

	calcReferencePoints();
	calcReferencePoints(m_est_bo,m_est_reck,m_est_s,m_est_fmsy,m_est_bmsy,m_est_msy);
	EstimatorClass cAssessmentModel(m_sAssessmentModel);
	HarvestControlRule cTAC(m_limit,m_threshold,m_target,m_ftarget);
	// HarvestControlRule cTAC(m_limit,m_threshold,m_target,m_fmsy);

	for( i = m_nyr+1; i <= m_pyr; i++ )
	{
		// | - Calculate reference points
		calcReferencePoints(m_est_bo,m_est_reck,m_est_s,m_est_fmsy,m_est_bmsy,m_est_msy);

		// | - Calculate TAC based on harvest control rule parameters.
		// TODO: replace with estimated biomass and bo from assessment.
		m_chat(i)  = cTAC(m_est_bt,m_est_bo,m_est_fmsy);

		// | - Implement fishery
		m_chat(i) *= exp( et_dev(i) );
		for( j = 1; j <= m_ngear; j++ )
		{
			if( m_chat(i,j) >= m_bt(i) )
			{
				m_chat(i,j) = 0.8/m_ngear * m_bt(i);
			}
		}


		// | - Update population dynamics
		m_rt(i)   = m_so * m_bt(i-m_agek) / (1.0 + m_beta*m_bt(i-m_agek));
		m_rt(i)   = m_rt(i) * exp(rt_dev(i));
		m_bt(i+1) = m_s*m_bt(i) + m_rt(i) - sum(m_chat(i));

		// | - Catch and effort statistics
		for( j = 1; j <= m_nEpochs; j++ )
		{

			k = m_nIt_nobs(j) + (i-m_nyr);
			m_ihat(j,k) = m_q(j) * m_bt(i) * exp(it_dev(j,i));
			m_it_yr_hat(j,k) = i;
			m_cvhat(j,k) = m_sig / m_ihat(j,k);
		}

		// | - Write new data file
		write_data_file(i,m_chat.sub(m_syr,i),m_ihat);

		// | - Run Assessment Model
		cAssessmentModel.runEstimator();

		// | - Get estimated parameters
		read_parameter_estimates("mse.par");
		
		// | - Screen dump
		if( !(i % 5) ) print_mse(i);
	}
	

	// | - Calculate average annual variation in catch.
	m_aav.allocate(m_syr,m_pyr);
	m_aav.initialize();
	for( i = m_syr+5; i <= m_pyr; i++ )
	{
		m_aav(i)  = sum(fabs(first_difference(rowsum(m_chat.sub(i-5,i)))));
		m_aav(i) /= sum(m_chat.sub(i-5,i));
	}
	
}

void lrgsOM::print_mse(const int &i)
{
	// -   Screen dump so you can watch the progress.
	cout<<setprecision(3)<<setw(11);
	cout<<"|---------------------------------|"                      <<endl;
	cout<<"| - Year    "<<i<<"                  |"                   <<endl;
	cout<<"|---------------------------------|"                      <<endl;
	cout<<"|           True  "<<"\t" <<"Estimated |"                 <<endl;
	cout<<"| - Bo      "<<m_bo<<"\t" <<m_est_bo                      <<endl;
	cout<<"| - Biomass "<<m_bt(i)<<"\t"<<m_est_bt                    <<endl;
	cout<<"| - F       "<<m_chat(i)/m_bt(i)<<"\t"<<m_chat(i)/m_est_bt<<endl;
	cout<<"| - Fmsy    "<<m_fmsy<<"\t"<<m_est_fmsy                   <<endl;
	cout<<"| - Bmsy    "<<m_bmsy<<"\t" <<m_est_bmsy                  <<endl;
	cout<<"| - msy     "<<m_msy<<"\t"<<m_est_msy                     <<endl;
	cout<<"| - Catch   "<<m_chat(i)                                  <<endl;
	cout<<"|---------------------------------|\n"                    <<endl;
}

void lrgsOM::calcReferencePoints()
{
	m_bmsy = (m_bo*(-1.+sqrt(m_reck))/(m_reck-1.)); 
	m_msy  = (m_bmsy * (-1. + m_s) * (-1. + m_reck)*(-m_bo + m_bmsy)) 
			 / (m_bo + (-1. + m_reck)*m_bmsy);
	
	m_fmsy = m_msy / m_bmsy;
	//cout<<"MSY\n"<<m_msy<<endl;
}


void lrgsOM::calcReferencePoints(const double &bo, const double & reck,const double &s,
                                 double &fmsy, double &bmsy, double &msy)
{
	bmsy = (bo*(-1.+sqrt(reck))/(reck-1.)); 
	msy  = (bmsy * (-1. + s) * (-1. + reck)*(-bo + bmsy)) 
			 / (bo + (-1. + reck)*bmsy);
	
	fmsy = msy / bmsy;
	//cout<<"MSY\n"<<m_msy<<endl;
}

void lrgsOM::write_data_file(const int &nyr, const dmatrix &ct,const dmatrix& it)
{
	int extra = (nyr-m_nyr);

	ofstream ofs("MSE.dat");
	ofs<<m_agek<<endl;
	ofs<<m_syr<<endl;
	ofs<<nyr<<endl;
	ofs<<m_ngear<<endl;
	ofs<<m_nEpochs<<endl;
	ofs<<m_nIt_nobs+extra<<endl;

	ofs<<"# Catch Data"<<endl;
	int i,j;
	ivector iyr(m_syr,nyr);
	iyr.fill_seqadd(m_syr,1);
	for( i = m_syr; i <= nyr; i++ )
	{
		ofs<<iyr(i)<<"\t"<<ct(i)<<endl;
	}


	ofs<<"#CPUE data"<<endl;
	ofs<<"#Year Epoch It CV"<<endl;
	for( j = 1; j <= m_nEpochs; j++ )
	{
		for( i = 1; i <= m_nIt_nobs(j) + extra; i++ )
		{
			ofs<<m_it_yr_hat(j,i)<<"\t"<<j<<"\t"<<m_ihat(j,i)<<" \t"<<m_cvhat(j,i)<<endl;;
		}
		
	}
	// ofs<<iyr<<endl;
	// ofs<<ct<<endl;
	// ofs<<it<<endl;

	ofs<<"#eof -999\n"<<-999<<endl;
}

void lrgsOM::read_parameter_estimates(const adstring &sParFile)
{
	ifstream ifs(sParFile);
	ifs >> m_est_bo;
	ifs >> m_est_reck;
	ifs >> m_est_s;
	ifs >> m_est_bt;
}

// HARVEST CONTROL RULE FUNCTIONS
HarvestControlRule::~HarvestControlRule()
{}

HarvestControlRule::HarvestControlRule(const double &limit, const double &threshold,
                                       const double &target, const double &ftarget)
: m_limit(limit), m_threshold(threshold), m_target(target), m_ftarget(ftarget)
{

}

HarvestControlRule::HarvestControlRule(const double &limit, const double &threshold,
                   const double &target, const double &ftarget,
                   const double &fmax, const double &catch_floor)
: m_limit(limit), m_threshold(threshold), m_target(target),
  m_ftarget(ftarget), m_fmax(fmax), m_catch_floor(catch_floor)
{
	
}


double HarvestControlRule::operator( )(const double& bt,const double &bo)
{
	// User defined harvest control rule.
	double tac = 0;
	double f_rate=0;
	double bstatus = bt/bo;
	if( bstatus <= m_limit )
	{
		f_rate = 0;
	}

	if( bstatus >= m_threshold )
	{
		f_rate = m_ftarget;
	}

	if( bstatus > m_limit && bstatus <= m_threshold )
	{
		f_rate = m_ftarget * ( bstatus - m_limit ) 
				            /( m_threshold - m_limit );
	}

	
	tac = f_rate * bt;
	return(tac);
}

double HarvestControlRule::operator( )(const double& bt,const double &bo, 
                                       double& f_rate)
{
	// User defined harvest control rule.
	double tac = 0;
	//double f_rate=0;
	double bstatus = bt/bo;
	if( bstatus <= m_limit )
	{
		f_rate = 0;
	}

	if( bstatus >= m_threshold )
	{
		f_rate = m_ftarget;
	}

	if( bstatus > m_limit && bstatus <= m_threshold )
	{
		f_rate = m_ftarget * ( bstatus - m_limit ) 
				            /( m_threshold - m_limit );
	}

	
	tac = f_rate * bt;
	return(tac);
}
























