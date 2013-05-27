// This header file uses an example of inheritance.
// The Scenario class is derived from the population class and inherits its member
// variables (bo,h,s).

#ifndef POPULATION_H
#define POPULATION_H
class Population
{
protected:      // use protected so derived classes can access these variables.
	double m_bo;
	double m_h;
	double m_s;

public:
	Population(const double bo=1.0, const double h=0.75, const double s=0.85)
	: m_bo(bo),m_h(h),m_s(s)
	{}

	~Population(){}
};

#endif


#ifndef SCENARIO_H
#define SCENARIO_H

class Scenario: public Population 
{
private:
	int     m_agek;
	int     m_pyr;  // number of projection years
	int     m_rseed;// Random number seed.
	// double  m_bo;
	// double  m_h;
	// double  m_s;
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
	Scenario(const int& agek, const int& pyr, const int& rseed, double& bo,const double& h,
	         const double& s, const double& q, const double& sig,const double tau,
	         const dvector& ft, const dvector &wt, const dvector &it,const dvector &ct)
		:m_agek(agek),m_pyr(pyr),m_rseed(rseed),Population(bo,h,s), m_q(q),
		 m_sig(sig), m_tau(tau),m_ft(ft), m_wt(wt), m_it(it), m_ct(ct)
	{}

	~Scenario();
	// getters
	int     get_agek() { return m_agek; }
	int     get_pyr()  { return m_pyr;  }
	int     get_rseed(){ return m_rseed;}
	double  get_bo()   { return m_bo;   }
	double  get_h()    { return m_h;    }
	double  get_s()    { return m_s;    }
	double  get_q()    { return m_q;    }
	double  get_sig()  { return m_sig;  }
	double  get_tau()  { return m_tau;  }
	dvector get_ft()   { return m_ft;   }
	dvector get_wt()   { return m_wt;   }
	dvector get_it()   { return m_it;   }
	dvector get_ct()   { return m_ct;   }

	// setters
	void  set_pyr(int v1)     { m_pyr = v1;  }
	void  set_bo(double v1)   { m_bo = v1;   }
	void  set_h(double v1)    { m_h = v1;    }
	void  set_s(double v1)    { m_s = v1;    }
	void  set_q(double v1)    { m_q = v1;    }
	void  set_sig(double v1)  { m_sig = v1;  }
	void  set_tau(double v1)  { m_tau = v1;  }
	void  set_ft(dvector v1)  { m_ft = v1;   }
	void  set_wt(dvector v1)  { m_wt = v1;   }
	void  set_it(dvector v1)  { m_it = v1;   }
	void  set_ct(dvector v1)  { m_ct = v1;   }

	// Make the operatingModel class a friend so it can access private members of scenario
	// friend class OperatingModel;
};

#endif