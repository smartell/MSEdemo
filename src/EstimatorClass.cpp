#include <admodel.h>
#include "EstimatorClass.h"

/// Destructor
EstimatorClass::~EstimatorClass()
{}

/// Constructor
EstimatorClass::EstimatorClass(adstring model)
: m_model(model)
{
	// cout<<"I'm in the Estimator class constructor"<<endl;
	// cout<<m_model<<endl;
	// system(m_model);
}


/**
\brief Function that runs the user defined estimator.
\author Steve Martell
*/
void EstimatorClass::runEstimator()
{
	adstring arg;
	arg = m_model+" -ind MSE.dat -nox -est > NUL";
	system(arg);
	// cout<<"Finished running the estimator"<<endl;
}