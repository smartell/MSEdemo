#ifndef OPERATINGMODEL_H
#define OPERATINGMODEL_H

#include <admodel.h>
#include "scenario.h"

class operatingModel
{
private:
	Scenario m_cScenario;

	operatingModel();

public:
	operatingModel(const Scenario &cScenario)
	: m_cScenario(cScenario)
	{}
};

#endif
