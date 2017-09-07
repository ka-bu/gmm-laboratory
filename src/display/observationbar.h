#include "../base.h"

#ifndef OBSERVATIONBAR_H
#define OBSERVATIONBAR_H

#include <vector>

#include <AntTweakBar.h>

#include "../settings/observationsettings.h"

class ObservationBar
{
public:
	ObservationBar(ObservationSettings& settings);

	virtual ~ObservationBar()
	{
	}

	static void TW_CALL buttonCallback(void* clientData);

	void setSolutionRange(int min, int max);
	void setComponentRange(int min, int max);
	void setProjectionRange(int min, int max);
	
	TwBar* bar;
};

#endif