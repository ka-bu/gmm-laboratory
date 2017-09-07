#include "base.h"

#include "gmmlab.h"
#include "display/gmmdisplay.h"
#include "settings/settings.h"
#include "settings/configparser.h"
#include "data/data_description.h"
#include "algorithm/description/algo_description.h"
#include "algorithm/description/init_description.h"


const static std::string GMMLAB_MODE_DEFAULT = "DEFAULT";
const static std::string GMMLAB_MODE_LOOKAT = "LOOKAT";



int main(int argc, char* argv[])
{  

#ifndef NDEBUG
	std::cout << std::endl
		<< " /-----------------------------------------\\" << std::endl
		<< "(  DEBUG mode is on: asserts are processed  )" << std::endl
		<< " \\-----------------------------------------/" << std::endl;
#endif


		
	configparser::parseGMMLabConfiguration(argc, argv);
	
	GMMLab gmmlab(true, gmmlabSettings().computeCosts, gmmlabSettings().outputFile);
	
	clock_t t;

	if(gmmlabSettings().mode == GMMLAB_MODE_DEFAULT)
	{
		std::vector<DataDescription*> dataDescriptors = DataDescription::createDataDescriptions(commonSettings().datafile);
		std::vector<InitDescription*> initDescriptors = InitDescription::createInitDescriptions(commonSettings().initfile);
		std::vector<AlgoDescription*> algoDescriptors = AlgoDescription::createAlgoDescriptions(commonSettings().algofile);
		
		if(dataDescriptors.size()!=1 || initDescriptors.size()!=1 || algoDescriptors.size()!=1)
		{
			std::cout << "configuration contains #datasets = " << dataDescriptors.size() << ", #inits = " << initDescriptors.size() << ", #algos = " << algoDescriptors.size() << std::endl;
			gmmlab_throw("Please define exactly one dataset, initialization and algorithm for DEFAULT mode.")
		}
		

		t = clock();
		std::cout << "precomputing..." << std::endl;
		gmmlab.runDefaultTests(dataDescriptors, initDescriptors, algoDescriptors);
		std::cout << "time with overhead: " << secondsSince(t) << " seconds." << std::endl;
	}
	else if(gmmlabSettings().mode == GMMLAB_MODE_LOOKAT)
	{
		
		std::vector<DataDescription*> dataDescriptors = DataDescription::createDataDescriptions(commonSettings().datafile);
		if(dataDescriptors.size()!=1)
			gmmlab_throw("Please define exactly one dataset in LOOKAT mode.")
			
		gmmlab.runDataGenerationOnly(dataDescriptors);

	}
	else
		gmmlab_throw("Unknown gmmlab mode");
	


	if (!commonSettings().ndisplay)
	{
		GMMDisplay display(gmmlab.truth, gmmlab.input, gmmlab.sampleSets, gmmlab.solutions, gmmlab.tags, gmmlab.runtimes, gmmlab.costs);
		display.mainloop();
	}


	return 0;
}
