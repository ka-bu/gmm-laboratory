#ifndef COMMONSETTINGS_H
#define COMMONSETTINGS_H

#include "../base.h"

#include <string>

struct CommonSettings
{
	bool verbose = false;
	
	// display
	bool ndisplay = false;

	// input
	std::string outputDir;
	std::string datafile, initfile, algofile; 
	
	// store/display intermediate solutions
	bool intermediate = false;
	
	// output
	bool csv = false;
	bool dataset = false;
	bool wgt = false;
	bool dist = false;
	bool cls = false;
	bool gmm = false;
	bool time = false;
	bool difftruth = false;
	bool printStatistics = false;
	
	bool reinit = false; // reinit second algo with last solution of first algo in every diff step
	
	// cost computation
	std::string costmeasure = "GAUSSIAN";
};

#endif
