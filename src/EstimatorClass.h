#include <admodel.h>

#ifndef ESTIMATOR_H
#define ESTIMATOR_H
/**
\brief Calls the user defined stock assessment model.
\author Steve Martell
*/
class EstimatorClass
{
public:
	EstimatorClass(adstring model);
	~EstimatorClass();

	/* data */
	adstring m_model;
	void runEstimator();

};

#endif