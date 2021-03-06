#include <admodel.h>
#ifndef HARVEST_CONTROL_RULE_H
#define HARVEST_CONTROL_RULE_H

/*! \brief A harvest control rule class 
	\author Steven Martell
	\remarks Contains an enum type that defines the various harvest control rules 
	that can be used
*/
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
		FORTY_TEN,           		/**< 40:10 harvest contrl rule               */
		FIXED_ESCAPEMENT,			/**< Fixed escapement rule                   */
		FIXED_ESCAPEMENT_CAP,		/**< Fixed escapement rule with a cap at MSY */
		FIXED_HARVEST_RATE,			/**< Fixed harvest rate rule based on Fmys   */
		CONDITIONAL_CONSTANT_CATCH,	/**< IPHC's Conditional constant catch rule  */
		THIRTY_TWENTY,				/**< 30:20 harvest control rule              */
		FIXED_HR_DELTA,             /**< Fixed catch with 15% maximum changed    */
		FLOOR_THIRTY_TWENTY,        /**< 30:20 rule with a minimum catch floor.  */
		FAO_PA_COMPLIANT   			/**< Not implemented yet                     */
	};

	~HarvestControlRule(){}
	HarvestControlRule(int &eHCR)
	:m_enum(eHCR)
	{
		// cout<<"m_enum = "<<m_enum<<endl;	
	}

	// Prototypes
	double getTac(const double &bt, const double &fmsy, const double &msy,
	              const double &bmsy, const double &bo,
	              const double &delta, const double &ptac, const double &mintac);
	double FortyTen(const double &bt, const double &bo, const double &fmsy);
	double ThirtyTwenty(const double &bt, const double &bo, const double &fmsy);
	double FixedHarvestRate(const double &bt, const double &fmsy);
	double FixedEscapement(const double &bt, const double &bmsy);
	double FixedEscapementCap(const double &bt, const double &bmsy, const double &msy);
	double ConditionalConstantCatch(const double& bt, const double& bmsy, const double& msy, const double& fmsy);
	double FixedHarvestRateDelta(const double& bt,const double& fmsy,const double &delta, const double &ptac);
	double FloorThirtyTwenty(const double &bt, const double &bo, const double &fmsy, const double &mintac);
};


#endif

