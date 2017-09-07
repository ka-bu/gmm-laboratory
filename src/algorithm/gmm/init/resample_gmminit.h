#ifndef RESAMPLE_GMMINIT_H
#define RESAMPLE_GMMINIT_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"

namespace initializer
{
	/**
	 * first computes 2*d adaptive spherical initializations (independent from each other), 
	 * then computes an adaptive spherical initialization with respect to the means of these
	 * initializations.
	 */
	Parameters resampleGMMInit(commonutil::DataSet const& input, std::string initmethod, idx_type k, std::mt19937& gen);	
}

/**
 * @brief Descriptor for Initialization
 */
class ResampleGMMInitID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	ResampleGMMInitID(std::string initmethod, uint32_t seed = 1);
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	std::string initmethod;
};

#endif
