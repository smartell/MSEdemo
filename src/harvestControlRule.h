#include <admodel.h>
#ifndef HARVEST_CONTROL_RULE_H
#define HARVEST_CONTROL_RULE_H

/*! A harvest control rule class */
class HarvestControlRule
{
private:
	int m_enum;
	// double (*m_pRuleType); // a pointer to appropriate function to get the tac.
public:
	/** A enum type to select the appropriate harvest control rule.
	 *
	 */
	enum enumHCR 
	{
		FORTY_TEN,           /**< 40:10 harvest contrl rule */
		FIXED_ESCAPEMENT,
		FIXED_ESCAPEMENT_CAP,
		FIXED_HARVEST_RATE,
		CONDITIONAL_CONSTANT_CATCH,
		FAO_PA_COMPLIANT
	};

	~HarvestControlRule(){}
	HarvestControlRule(int &eHCR)
	:m_enum(eHCR)
	{
		cout<<"m_enum = "<<m_enum<<endl;	
	}

	// Prototypes
	double getTac(const double &bt, const double &fmsy, const double &msy,
	              const double &bmsy, const double &bo);
	double FortyTen(const double &bt, const double &bo, const double &fmsy);
	double FixedHarvestRate(const double &bt, const double &fmsy);
	double FixedEscapement(const double &bt, const double &bmsy);
	double FixedEscapementCap(const double &bt, const double &bmsy, const double &msy);
	double ConditionalConstantCatch(const double& bt, const double& bmsy, const double& msy, const double& fmsy);
};


#endif

