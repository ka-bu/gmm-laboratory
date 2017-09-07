#ifndef GMM_IOUTIL_H
#define GMM_IOUTIL_H

#include "../base.h"
#include "ioutil.h"
#include "../algorithm/gmm/utils/gmmutil.h"

#include <vector>

/**
 * @brief 
 */
namespace GMMIOUtil
{

	void appendToGMM(std::string const path, Parameters const& desc, std::string const tag);
	void loadFromGMM(std::string const path, std::vector<Parameters>& gmms, std::vector<std::string>& tags);
	//void loadGMMs(std::string const pattern, std::vector<GMMDesc>& gmms, std::vector<std::string>& tags);
	
	std::string const GMM_FILE_EXTENSION = ".gmm";
};

#endif