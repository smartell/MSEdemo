#include <admodel.h>
#undef MAXITER
#define MAXITER  17

#ifndef MSY_REFERENCE_POINTS_H
#define MSY_REFERENCE_POINTS_H


class msy_reference_points
{
private:
	double m_k;
	double m_s;
	double m_bo;

	double m_fmsy;
	double m_msy;
	double m_bmsy;

public:
	~msy_reference_points(){}

	// default constructor
	msy_reference_points(const double &k, const double &s, const double &bo)
		:m_k(k), m_s(s), m_bo(bo)
		{
			calc_rp();
		}

	// getters
	double get_fmsy()  { return m_fmsy; }
	double get_bmsy()  { return m_bmsy; }
	double get_msy()   { return m_msy;  }
	// member functions
	void calc_rp();

	// Make the operatingModel class a friend so it can access private members of scenario
	friend class operatingModel;

};

#endif

void msy_reference_points::calc_rp()
{
	m_bmsy = (m_bo*(-1.+sqrt(m_k))/(m_k-1.)); 
	m_msy  = (m_bmsy*(m_s-1.)*(m_k-1.)*(m_bmsy-m_bo)) / (m_bo + (m_k-1.)*m_bmsy);
	// m_msy  = (m_bmsy*(1.-m_s) + m_k*(1.-m_s)*m_bmsy) / (1.+m_bmsy*(m_k-1.)/m_bo);

	if( m_msy<m_bmsy )
	{
		m_fmsy = -log((-m_msy+m_bmsy)/m_bmsy);
	}
	else
	{
		m_fmsy = m_msy/m_bmsy;
	}

	// cout<<"Bmsy = "<<m_bmsy<<endl;
	// cout<<"Bmsy/Bo = "<<m_bmsy/m_bo<<endl;
	// cout<<m_msy<<endl;
	// cout<<m_fmsy<<endl;	
}