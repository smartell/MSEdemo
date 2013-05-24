#include <admodel.h>
#include "HarvestControlRule.h"
double HarvestControlRule::getTac(const double &bt, const double &fmsy, const double &msy,
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
		case FAO_PA_COMPLIANT:
			cout<<"FAO_PA_COMPLIANT"<<endl;
			break;
	}
	return tac;
}

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
		ft = fmsy * (dt-0.1)/(0.4-0.3);
	}
	else if( dt > 0.4 )
	{
		ft = fmsy;
	}

	return bt*(1.-exp(-ft));
}

double HarvestControlRule::FixedHarvestRate(const double &bt, const double &fmsy)
{
	double tac = bt * (1. - exp(-fmsy));
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
