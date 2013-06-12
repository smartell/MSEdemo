#include <admodel.h>
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
           dvector& ct,
           dvector& it,
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
	m_wt(wt)
{
	// cout<<"IN CONSTRUCTOR"<<endl;
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
	m_bt(m_syr)          = m_bo;
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
		m_ft(i) = -log((-m_ct(i)+m_bt(i))/m_bt(i));
		if(i-m_syr > m_agek)
		{
			m_rt(i) = m_a*m_bt(i-m_agek)/(1.+m_b*m_bt(i-m_agek)) * exp(m_wt(i));	
		}
		
		btmp    = m_s*m_bt(i) + m_rt(i) - m_ct(i);
		m_bt(i+1) = posfun(btmp,0.1,m_fpen);
	}
	// sd_dep = m_bt(nyr)/m_bo;
}

void LRGS::observation_model()
{

	int i;
	dvar_vector zt = log(m_it) - log(m_bt(m_syr,m_nyr));
	m_q            = exp(mean(zt));
	m_epsilon      = zt - mean(zt);
}
