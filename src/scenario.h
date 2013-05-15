#ifndef SCENARIO_H
#define SCENARIO_H

class Scenario
{
private:
	double m_bo;
	double m_h;

public:
	// default constructor
	Scenario(const double& bo,const double& h)
		:m_bo(bo), m_h(h)
	{}

	// Make the operatingModel class a friend so it can access private members of scenario
	friend class operatingModel;
};

#endif