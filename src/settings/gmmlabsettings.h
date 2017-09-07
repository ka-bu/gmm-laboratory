#ifndef GMMLABSETTINGS_H
#define GMMLABSETTINGS_H

#include "../base.h"

#include <string>

struct GMMLabSettings
{
	std::string mode = "DEFAULT";
	bool computeCosts = true;
	std::string outputFile = "output";
};

#endif
