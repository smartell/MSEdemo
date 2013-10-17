/**
 * \file Scenario.h
 * \author Steve Martell
**/

#ifndef POPULATION_H
#define POPULATION_H
/**
\brief Population class
\author Steve Martell
\remarks Contains member variables for unfished biomass, steepness and survival.
*/
class Population
{
private:
	double m_bo; //!< Unfished biomass
	double m_h;  //!< Steepness of Beverton Holt Model
	double m_s;  //!< Survival growth coefficient

public:
	/** Constructor for population class */
	Population(const double bo=1.0, const double h=0.75, const double s=0.85)
	: m_bo(bo),m_h(h),m_s(s)
	{}

	~Population(){}
	friend class Scenario;
};

#endif


#ifndef SCENARIO_H
#define SCENARIO_H
/** \brief  Scenario class
	
	class object for a cadidate management procedure
	
&copy; Copyright 2013 UBC Fisheries Centre - . All Rights Reserved.

	\author  Steven Martell
	\author $LastChangedBy$
	\date 2013-05-27
	\date $LastChangedDate$
	\version $Rev$
	\sa
**/

class Scenario: public Population
{
private:
	int     m_agek;
	int     m_nScenario;  // Stationary or non-stationary stock-recruitment relationship
	int     m_pyr;  // number of projection years
	int     m_flg_perfect_information;  /// flag for perfect information.
	int     m_rseed;// Random number seed.
	// double  m_bo;
	// double  m_h;
	// double  m_s;
	double  m_q;
	double  m_sig;
	double  m_tau;
	double  m_iuu_rate;
	double  m_mintac;
	dvector m_ft;
	dvector m_wt;
	dvector m_it;
	dvector m_ct;


	Scenario(); //!< default constructor
public:
	/**
	Constructor for the scenario class.
	*/
	Scenario(const int& agek, const int& _nScenario, const int& pyr, const int& rseed, double& bo,const double& h,
	         const double& s,const double& iuu_rate, const double& q, const double& sig,const double tau,
	         const dvector& ft, const dvector &wt, const dvector &it,const dvector &ct)
		:m_agek(agek),m_nScenario(_nScenario),m_pyr(pyr),m_rseed(rseed),Population(bo,h,s), m_q(q),
		 m_sig(sig), m_tau(tau),m_ft(ft), m_wt(wt), m_it(it), m_ct(ct),m_iuu_rate(iuu_rate)
	{}

	Scenario(const int& agek, const int& _nScenario, const int& pyr, const int& _flg_perfect_information,
	         const int& rseed, double& bo,const double& h,
	         const double& s,const double& iuu_rate,const double& q, const double& sig,const double tau,
	         const dvector& ft, const dvector &wt, const dvector &it,const dvector &ct)
		:m_agek(agek),m_nScenario(_nScenario),m_pyr(pyr),m_rseed(rseed),Population(bo,h,s), m_q(q),
		 m_sig(sig), m_tau(tau),m_ft(ft), m_wt(wt), m_it(it), m_ct(ct),
		 m_flg_perfect_information(_flg_perfect_information),m_iuu_rate(iuu_rate)
	{}

	Scenario(const int& agek, const int& _nScenario, const int& pyr, const int& _flg_perfect_information,
	         const int& rseed, double& bo,const double& h,
	         const double& s,const double& iuu_rate,const double& q, const double& sig,const double tau,
	         const dvector& ft, const dvector &wt, const dvector &it,const dvector &ct,
	         const double &mintac)
		:m_agek(agek),m_nScenario(_nScenario),m_pyr(pyr),m_rseed(rseed),Population(bo,h,s), m_q(q),
		 m_sig(sig), m_tau(tau),m_ft(ft), m_wt(wt), m_it(it), m_ct(ct),
		 m_flg_perfect_information(_flg_perfect_information),m_iuu_rate(iuu_rate),
		 m_mintac(mintac)
	{}

	~Scenario();
	// getters
	int     get_agek() { return m_agek; } //!< get
	int     get_nScenario(){ return m_nScenario; }
	int     get_nInformation() {return m_flg_perfect_information; }
	int     get_pyr()   { return m_pyr;  } //!< get
	int     get_rseed() { return m_rseed;} //!< get
	double  get_bo()    { return m_bo;   } //!< get
	double  get_h()     { return m_h;    } //!< get
	double  get_s()     { return m_s;    } //!< get
	double  get_q()     { return m_q;    } //!< get
	double  get_sig()   { return m_sig;  } //!< get
	double  get_tau()   { return m_tau;  } //!< get
	double  get_iuu()   { return m_iuu_rate; }
	double  get_mintac(){ return m_mintac;   }
	dvector get_ft()    { return m_ft;   } //!< get
	dvector get_wt()    { return m_wt;   } //!< get
	dvector get_it()    { return m_it;   } //!< get
	dvector get_ct()    { return m_ct;   } //!< get

	// setters
	void  set_pyr(int v1)     { m_pyr = v1;  } //!< set projections years
	void  set_bo(double v1)   { m_bo = v1;   } //!< set unfished biomass
	void  set_h(double v1)    { m_h = v1;    } //!< set steepness
	void  set_s(double v1)    { m_s = v1;    } //!< set survival and growth
	void  set_q(double v1)    { m_q = v1;    } //!< set catchability
	void  set_sig(double v1)  { m_sig = v1;  } //!< set observation error
	void  set_tau(double v1)  { m_tau = v1;  } //!< set process error
	void  set_ft(dvector v1)  { m_ft = v1;   } //!< set fishing mortality
	void  set_wt(dvector v1)  { m_wt = v1;   } //!< set recruitment deviations
	void  set_it(dvector v1)  { m_it = v1;   } //!< set relative abundance
	void  set_ct(dvector v1)  { m_ct = v1;   } //!< set catch

	// Make the operatingModel class a friend so it can access private members of scenario
	// friend class OperatingModel;
};

#endif