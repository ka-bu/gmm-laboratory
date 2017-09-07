
#include "base.h"
#include "gmmlab.h"

#include "algorithm/description/init_description.h"
#include "algorithm/description/algo_description.h"
#include "data/data_description.h"

#include "settings/settings.h"
#include "settings/configparser.h"
#include "ioutil/data_ioutil.h"

#include <boost/filesystem.hpp>


const std::string TESTLAB_MODE_DEFAULT = "DEFAULT";
const std::string TESTLAB_MODE_DATA = "DATA";


int main(int argc, char* argv[])
{

#ifndef NDEBUG
	std::cout << std::endl
		<< " /-----------------------------------------\\" << std::endl
		<< "(  DEBUG mode is on: asserts are processed  )" << std::endl
		<< " \\-----------------------------------------/" << std::endl;
#endif
	
	
	
	
	configparser::parseTestLabConfiguration(argc, argv);
	
	GMMLab gmmlab(false, false);
	
	if(testlabSettings().mode == TESTLAB_MODE_DEFAULT)
	{
		std::cout << "--> Running Default Tests:" << std::endl << std::endl;
		
		std::vector<DataDescription*> dataDescriptors = DataDescription::createDataDescriptions(commonSettings().datafile);
		std::vector<InitDescription*> initDescriptors = InitDescription::createInitDescriptions(commonSettings().initfile);
		std::vector<AlgoDescription*> algoDescriptors = AlgoDescription::createAlgoDescriptions(commonSettings().algofile);
		
		gmmlab.runDefaultTests(dataDescriptors, initDescriptors, algoDescriptors);
	}
	else if(testlabSettings().mode == TESTLAB_MODE_DATA)
	{
		std::cout << "--> Generating Data Sets:" << std::endl << std::endl;
		commonSettings().verbose = true;
		DataIOUtil::runDataGeneration(commonSettings().datafile, commonSettings().outputDir);
		
	}

	
	
	return 0;
}
