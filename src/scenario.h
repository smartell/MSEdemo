#ifndef SCENARIO_H
#define SCENARIO_H

class Scenario
{
private:
	double  m_bo;
	double  m_h;
	double  m_s;
	double  m_sig;
	double  m_tau;
	dvector m_ft;

	Scenario();
public:
	// default constructor
	Scenario(const double& bo,const double& h,const double& s,
	         const double& sig,const double tau, const dvector& ft)
		:m_bo(bo), m_h(h), m_s(s), m_sig(sig), m_tau(tau), m_ft(ft)
	{}

	// getters
	double getBo()   { return m_bo;   }
	double geth()    { return m_h;    }
	double gets()    { return m_s;    }
	double getsig()  { return m_sig;  }
	double gettau()  { return m_tau;  }
	dvector getft()  { return m_ft;   }

	// Make the operatingModel class a friend so it can access private members of scenario
	friend class operatingModel;
};

#endif