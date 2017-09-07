#ifndef GONZALEZ2GMM_ID_H
#define GONZALEZ2GMM_ID_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"

namespace initializer
{
		
	/**
	 * runs Gonzalez algorithm and calls kmeansutil::meansToGMM.
	 */
	Parameters gonzalezToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen);
}

/**
 * @brief Descriptor for Initialization
 */
class Gonzalez2GMMID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	Gonzalez2GMMID(uint32_t = 1);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
};

#endif
