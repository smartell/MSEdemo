#include <admodel.h>
#include "HarvestControlRule.h"

/**
 * \brief Return the TAC given estimates of current biomass and reference points.
 * \author Steve Martell
 * \remarks Uses a swithc statement based on the type of harvest control rule to call the
 *          appropriate function to calculate the TAC.
**/
double HarvestControlRule::getTac(const double &bt, const double &fmsy, const double &msy,
	              				  const double &bmsy, const double &bo, 
	              				  const double &delta, const double &ptac,
	              				  const double &mintac)
{
	/** \param bt preseason biomass forecast                    */
	/** \param fmsy estimate of Fmsy                            */
	/** \param msy  estimate of MSY                             */
	/** \param bmsy estimate of Bmsy                            */
	/** \param bo estimate of unfished spawning biomass         */
	double tac = 0;
	switch (m_enum)
	{
		case FORTY_TEN:
			tac = FortyTen(bt, bo, fmsy);
			cout<<"Forty Ten tac = "<<tac<<endl;
			break;
		case FIXED_ESCAPEMENT:
			tac = FixedEscapement(bt,bmsy);
			cout<<"Fixed Escapement tac = "<<tac<<endl;
			break;
		case FIXED_ESCAPEMENT_CAP:
			tac = FixedEscapementCap(bt,bmsy,msy);
			cout<<"Fixed Escapement cap tac = "<<tac<<endl;
			break;
		case FIXED_HARVEST_RATE:
			tac = FixedHarvestRate(bt,fmsy);
			cout<<"Fixed harvest rate tac = "<<tac<<endl;
			break;
		case CONDITIONAL_CONSTANT_CATCH:
			tac = ConditionalConstantCatch(bt,bmsy,msy,fmsy);
			cout<<"Conditional constant catch = "<<tac<<endl;
			break;
		case THIRTY_TWENTY:
			tac = ThirtyTwenty(bt,bo,fmsy);
			cout<<"Thirty Twenty tac = "<<tac<<endl;
			break;
		case FIXED_HR_DELTA:
			tac = FixedHarvestRateDelta(bt,fmsy,delta,ptac);
			cout<<"Fixed Harvest Rate Delta tac = "<<tac<<endl;
			break;
		case FAO_PA_COMPLIANT:
			cout<<"FAO_PA_COMPLIANT"<<endl;
			break;
		case FLOOR_THIRTY_TWENTY:
			cout<<"FLOOR_THIRTY_TWENTY"<<endl;
			tac = FloorThirtyTwenty(bt,bo,fmsy,mintac);
			break;

	}

	// Put in SUFD adjustment here
	// 50% decrease, 33% increase

	// Evaluation framework, need to be able
	// to tune the MP to achieve the best 
	// tradeoff in the performance measures.


	// Floor 
	return tac;
}


/**
	\brief Implement the 40:10 harvest control rule.
	\author Steve Martell
	\param <bt> biomass available at time t
	\param <bo> theoretical unfished biomass
	\param <fmsy> Fishing mortality rate that generates Maximum Sustainable Yield.

	This function implements the 40:10 harvest control rule that is commonly used by the
	Pacific Fisheries Management Council.
*/
double HarvestControlRule::FortyTen(const double &bt, const double &bo, const double &fmsy)
{
	double dt = bt/bo;
	double ft = fmsy;
	if( dt <= 0.1 )
	{
		ft = 0;	
	} 
	else if( dt > 0.1 && dt <= 0.4 )
	{
		ft = fmsy * (dt-0.1)/(0.4-0.1);
	}
	else if( dt > 0.4 )
	{
		ft = fmsy;
	}

	return bt*(1.-exp(-ft));
}

/**
 * \brief 30:20 harvest control rule
 * \author Steve Martell
 * \remarks 
 * \param bt preseson available biomass forecast.
 */
double HarvestControlRule::ThirtyTwenty(const double &bt, const double &bo, const double &fmsy)
{
	double dt = bt/bo;
	double ft = fmsy;
	if( dt <= 0.2 )
	{
		ft = 0;	
	} 
	else if( dt > 0.2 && dt <= 0.3 )
	{
		ft = fmsy * (dt-0.2)/(0.3-0.2);
	}
	else if( dt > 0.3 )
	{
		ft = fmsy;
	}

	return bt*(1.-exp(-ft));
}

double HarvestControlRule::FloorThirtyTwenty(const double &bt, const double &bo, const double &fmsy, const double &mintac)
{
	
	double dt = bt/bo;
	double ft = fmsy;
	if( dt <= 0.2 )
	{
		ft = 0;	
	} 
	else if( dt > 0.2 && dt <= 0.3 )
	{
		ft = fmsy * (dt-0.2)/(0.3-0.2);
	}
	else if( dt > 0.3 )
	{
		ft = fmsy;
	}
	double tac = bt*(1.-exp(-ft));
	tac <= mintac? tac=mintac: tac=tac;
	return tac;
}


double HarvestControlRule::FixedHarvestRate(const double &bt, const double &fmsy)
{
	double tac = bt * (1. - exp(-fmsy));
	return tac;
}

double HarvestControlRule::FixedHarvestRateDelta(const double &bt, const double &fmsy,
                                                 const double &delta, const double &ptac)
{
	double tac = bt * (1. - exp(-fmsy));
	double dc = (tac-ptac)/ptac;
	cout<<"tac  = "<<tac<<endl;
	cout<<"ptac = "<<ptac<<endl;
	cout<<"dc   = "<<dc<<endl;
	if( fabs(dc) > delta )
	{
		int v=(int)dc;  v > 0 ? v=1 : v=-1;
		tac = (v*delta) * ptac + ptac;
	cout<<"Steve'O has done some stupid stuff "<<v*delta<<"\t"<<tac<<endl;
	}
	return tac;
}

double HarvestControlRule::FixedEscapement(const double &bt, const double &bmsy)
{
	double tac;
	if(bt > bmsy)
	{
		tac = (-bmsy + bt);
	}
	else
	{
		tac = 0;
	}
	return tac;
}

double HarvestControlRule::FixedEscapementCap(const double &bt, const double &bmsy, const double &msy)
{
	double tac;
	if(bt > bmsy)
	{
		tac = min((-bmsy + bt),msy);
	}
	else
	{
		tac = 0;
	}
	return tac;
}

/**
	\brief Implement the Conditional Constant Catch harvest control rule proposed by the IPHC.
	\author Steve Martell
	\param bt available biomass
	\param bmsy Biomass at MSY
	\param msy  Maximum sustainable yield
	\param fmsy Fishing mortality rate that achieves MSY.

 * For this rule, fish at Fmsy if bt>0.8Bmsy
 * and reduce the TAC to msy of the tac > MSY.
 * 
 If the biomass is less than 0.8Bmsy and greater
 then 0.4Bmsy, then set fishing mortality rate as
 a liner function of depletion. 

 If the biomass is less than 0.4Bmsy, the  set the tac=0
 */
double HarvestControlRule::ConditionalConstantCatch(const double& bt, const double& bmsy, const double& msy, const double& fmsy)
{
	double tac;
	double f = 0;
	if( bt >= 0.8*bmsy)
	{
		f = fmsy;
	}
	else if( bt < 0.8*bmsy && bt >= 0.4*bmsy )
	{
		f = (bt-0.4*bmsy)/(0.8-0.4)*fmsy;
	}
	tac = f * bt;
	if(tac > msy)
	{
		tac = msy;
	}
	return tac;
}
