#ifndef SCENARIO_H
#define SCENARIO_H

class Scenario
{
private:
	int     m_agek;
	double  m_bo;
	double  m_h;
	double  m_s;
	double  m_sig;
	double  m_tau;
	dvector m_ft;
	dvector m_wt;

	Scenario();
public:
	// default constructor
	Scenario(const int& agek, const double& bo,const double& h,const double& s,
	         const double& sig,const double tau,const dvector& ft,const dvector &wt)
		:m_agek(agek),m_bo(bo), m_h(h), m_s(s), m_sig(sig), m_tau(tau), m_ft(ft), m_wt(wt)
	{}

	~Scenario() {}
	// getters
	int     getagek() { return m_agek;  }
	double  getBo()   { return m_bo;   }
	double  geth()    { return m_h;    }
	double  gets()    { return m_s;    }
	double  getsig()  { return m_sig;  }
	double  gettau()  { return m_tau;  }
	dvector getft()   { return m_ft;   }
	dvector getwt()   { return m_wt;   }

	// setters
	void  setBo(double v1)   { m_bo = v1;   }
	void  seth(double v1)    { m_h = v1;    }
	void  sets(double v1)    { m_s = v1;    }
	void  setsig(double v1)  { m_sig = v1;  }
	void  settau(double v1)  { m_tau = v1;  }
	void  setft(dvector v1)  { m_ft = v1;   }
	void  setwt(dvector v1)  { m_wt = v1;   }

	// Make the operatingModel class a friend so it can access private members of scenario
	friend class operatingModel;
};

#endif