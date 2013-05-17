#include <admodel.h>
#ifndef HARVESTCONTROLRULE_H
#define HARVESTCONTROLRULE_H

class harvestControlRule
{
private:
	int m_enum;
	// double (*m_pRuleType); // a pointer to appropriate function to get the tac.
public:

	enum enumHCR 
	{
		FORTY_TEN,
		FIXED_ESCAPEMENT,
		FIXED_HARVEST_RATE,
		FAO_PA_COMPLIANT
	};

	~harvestControlRule(){}
	harvestControlRule(int &eHCR)
	:m_enum(eHCR)
	{
		cout<<"m_enum = "<<m_enum<<endl;	
	}

	// Prototypes
	// m_HCR.getTac(bt(i),fmsy,msy,bmsy,bo);
	double getTac(const double &bt, const double &fmsy, const double &msy,
	              const double &bmsy, const double &bo);
	double FortyTen(const double &bt, const double &bo, const double &fmsy);
	double FixedHarvestRate();

	// Make the operatingModel class a friend so it can access private members of scenario
	friend class operatingModel;
};


#endif

double harvestControlRule::getTac(const double &bt, const double &fmsy, const double &msy,
	              				  const double &bmsy, const double &bo)
{
	double tac = 0;
	switch (m_enum)
	{
		case FORTY_TEN:
			tac = FortyTen(bt, bo, fmsy);
			cout<<"Forty Ten tac = "<<tac<<endl;
			break;
		case FIXED_ESCAPEMENT:
			cout<<"FIXED_ESCAPEMENT"<<endl;
			break;
		case FIXED_HARVEST_RATE:
			cout<<"FIXED_HARVEST_RATE"<<endl;
			break;
		case FAO_PA_COMPLIANT:
			cout<<"FAO_PA_COMPLIANT"<<endl;
			break;
	}
	return tac;
}

double harvestControlRule::FortyTen(const double &bt, const double &bo, const double &fmsy)
{
	double dt = bt/bo;
	double ft = fmsy;
	if( dt <= 0.1 )
	{
		ft = 0;	
	} 
	else if( dt > 0.1 && dt <= 0.4 )
	{
		ft = fmsy * (dt-0.1)/(0.4-0.3);
	}
	else if( dt > 0.4 )
	{
		ft = fmsy;
	}

	return bt*(1.-exp(-ft));
}

double harvestControlRule::FixedHarvestRate()
{
	cout<<"Hi I'm donig fixed harvest rates"<<endl;
	return 0;
}