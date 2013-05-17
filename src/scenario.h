#ifndef SCENARIO_H
#define SCENARIO_H

class Scenario
{
private:
	int     m_agek;
	int     m_pyr;  // number of projection years
	double  m_bo;
	double  m_h;
	double  m_s;
	double  m_q;
	double  m_sig;
	double  m_tau;
	dvector m_ft;
	dvector m_wt;
	dvector m_it;
	dvector m_ct;

	Scenario();
public:
	// default constructor
	Scenario(const int& agek, const int& pyr, const double& bo,const double& h,
	         const double& s, const double& q, const double& sig,const double tau,
	         const dvector& ft, const dvector &wt, const dvector &it,const dvector &ct)
		:m_agek(agek),m_pyr(pyr),m_bo(bo), m_h(h), m_s(s), m_q(q), m_sig(sig), m_tau(tau), 
		m_ft(ft), m_wt(wt), m_it(it), m_ct(ct)
	{}

	~Scenario() {}
	// getters
	int     getagek() { return m_agek; }
	double  getBo()   { return m_bo;   }
	double  geth()    { return m_h;    }
	double  gets()    { return m_s;    }
	double  getq()    { return m_q;    }
	double  getsig()  { return m_sig;  }
	double  gettau()  { return m_tau;  }
	dvector getft()   { return m_ft;   }
	dvector getwt()   { return m_wt;   }
	dvector getit()   { return m_it;   }
	dvector getct()   { return m_ct;   }

	// setters
	void  setpyr(int v1)     { m_pyr = v1;  }
	void  setBo(double v1)   { m_bo = v1;   }
	void  seth(double v1)    { m_h = v1;    }
	void  sets(double v1)    { m_s = v1;    }
	void  setq(double v1)    { m_q = v1;    }
	void  setsig(double v1)  { m_sig = v1;  }
	void  settau(double v1)  { m_tau = v1;  }
	void  setft(dvector v1)  { m_ft = v1;   }
	void  setwt(dvector v1)  { m_wt = v1;   }
	void  setit(dvector v1)  { m_it = v1;   }
	void  setct(dvector v1)  { m_ct = v1;   }

	// Make the operatingModel class a friend so it can access private members of scenario
	friend class operatingModel;
};

#endif