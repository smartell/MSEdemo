#include <admodel.h>
#include <contrib.h>
#include <LRGS.h>


// constructor
LRGS::LRGS(const int& syr,
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
           dvar_vector& wt)
:	m_syr(syr),
	m_nyr(nyr),
	m_agek(agek),
	m_bo(bo),
	m_h(h),
	m_s(s),
	m_sig(sig),
	m_tau(tau),
	m_ct(ct),
	m_it(it),
	m_it_yr(it_yr),
	m_epoch(epoch),
	m_cv(cv),
	m_wt(wt)
{
	cout<<"IN CONSTRUCTOR"<<endl;
}

LRGS::LRGS(sLRGSdata& data,sLRGSparameters& pars)
{
	// cout<<"The other constructor"<<endl;
	m_syr      = data.syr;
	m_nyr      = data.nyr;
	m_agek     = data.agek;
	m_ct       = data.ct;
	m_it       = data.it;
	m_it_yr    = data.it_yr;
	m_nIt_nobs = data.nIt_nobs;
	m_epoch    = data.epoch;
	m_cv       = data.cv;
	m_ngear    = data.ngear;
	m_nEpochs  = data.nEpochs;

	m_bo    = mfexp(pars.log_bo);
	m_b1    = mfexp(pars.log_b1);
	m_h     = pars.h;
	m_s     = pars.s;
	m_gamma = *pars.gamma;
	m_sig   = sqrt(1.0/mfexp(pars.log_sigma));
	m_tau   = sqrt(1.0/mfexp(pars.log_tau));
	m_wt    = pars.wt;
}


void LRGS::initialize_model()
{
	m_rt.allocate(m_syr,m_nyr);
	m_ft.allocate(m_syr,m_nyr);
	m_bt.allocate(m_syr,m_nyr+1);
	m_rt.initialize();
	m_ft.initialize();
	m_bt.initialize();
	// bo               = mfexp(log_bo);
	m_ro               = m_bo*(1.-m_s);
	m_reck             = 4.*m_h/(1.-m_h);
	m_a                = m_reck*m_ro/m_bo;
	m_b                = (m_reck-1.0)/m_bo;
	m_rt(m_syr,m_syr+m_agek) = m_ro * exp(m_wt(m_syr,m_syr+m_agek));
	m_bt(m_syr)          = m_b1;
	// m_sig              = sqrt(1.0/mfexp(log_sigma));
	// m_tau              = sqrt(1.0/mfexp(log_tau));

	 // cout<<"OK to HERE"<<endl;


}


void LRGS::population_dynamics()
{
	// cout<<"In population dynamics"<<endl;
	int i;
	dvariable btmp;
	m_fpen.initialize();

	for(i=m_syr;i<=m_nyr;i++)
	{
		// m_ft(i) = -log((-m_ct(i)+m_bt(i))/m_bt(i));
		// m_ft(i) = m_ct(i) / m_bt(i);
		m_ft(i) = sum(m_ct(i)) / m_bt(i);
		if(i-m_syr > m_agek)
		{
			m_rt(i) = m_a*m_bt(i-m_agek)/(1.+m_b*m_bt(i-m_agek)) * exp(m_wt(i));	
		}
		
		btmp    = m_s*m_bt(i) + m_rt(i) - sum(m_ct(i));
		m_bt(i+1) = posfun(btmp,0.01,m_fpen);
		// m_ft(i) = -log( (m_bt(i)*(1.-m_s) + m_bt(i+1) - m_rt(i))/m_bt(i) );
	}
		// cout<<" fpen = "<<m_fpen <<endl;
	// sd_dep = m_bt(nyr)/m_bo;
}

void LRGS::observation_model()
{
	// SM CHanges Oct 23, to accomodate new data structures & multiple surveys.
	int i;
	m_epsilon.allocate(1,m_nEpochs,1,m_nIt_nobs);
	m_epsilon.initialize();

	m_q.allocate(1,m_nEpochs);
	m_q.initialize();

	for( i = 1; i <= m_nEpochs; i++ )
	{
		 ivector iyr    = m_it_yr(i);
		 dvar_vector zt = log(m_it(i)) - log(m_bt(iyr).shift(1));
		 m_q(i)         = exp(mean(zt));
		 m_epsilon(i)   = zt - mean(zt);
	}

	// dvar_vector zt = log(m_it) - log(m_bt(m_syr,m_nyr));
	// m_q            = exp(mean(zt));
	// m_epsilon      = zt - mean(zt);

}

void LRGS::calc_negative_loglikelihoods()
{
	int i;
	m_nll.allocate(1,2,1,m_nEpochs); m_nll.initialize();
	
	// Likelihood component for the residuals.
	for( i = 1; i <= m_nEpochs; i++ )
	{
		dvar_vector std = m_sig + sqrt(log( 1.0 + square(m_cv(i)) ));
		m_nll(1)(i)     = dnorm(m_epsilon(i),std);
	}
	

	// Likelihood component for the process error components.
	// If gamma is not active (i.e. no autocorrelation)
	m_delta.allocate(m_syr,m_nyr);  m_delta.initialize();
	if( (m_gamma==0.0) || (m_gamma==1.0) )
	{
		m_delta = m_wt;
		m_nll(2,1) = dnorm(m_delta,m_tau);
	}
	else
	{
		m_delta(m_syr) = square(m_wt(m_syr));
		for( i = m_syr+1; i <= m_nyr; i++ )
		{
			m_delta(i) = log(m_rt(i)) - (1.0-m_gamma)*(log(m_rt(i))-m_wt(i)) 
			           - m_gamma*(log(m_rt(i-1)));
		}
		m_nll(2,1) = dnorm(m_delta,m_tau);
	}


	// cout<<m_nll<<endl;
}


