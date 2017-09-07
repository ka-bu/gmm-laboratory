#include "../base.h"

#ifndef DISPLAYBAR_H
#define DISPLAYBAR_H

#include "../settings/displaysettings.h"

#include <vector>
#include <AntTweakBar.h>


class DisplayBar
{
public:
	DisplayBar(DisplaySettings& settings);

	virtual ~DisplayBar()
	{
	}
	
	TwBar* bar;
};

#endif
