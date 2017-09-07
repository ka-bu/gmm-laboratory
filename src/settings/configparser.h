#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

#include "../base.h"

#include <boost/program_options.hpp>

namespace configparser
{
	void parseCommmonConfiguration(boost::program_options::options_description& config);
	void finalizeCommonConfiguration();
	void parseGMMLabConfiguration(int, char*[]);
	void parseTestLabConfiguration(int, char*[]);
}


#endif