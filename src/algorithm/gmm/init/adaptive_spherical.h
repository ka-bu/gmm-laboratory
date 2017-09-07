#ifndef ADAPTIVE_SPHERICAL_H
#define ADAPTIVE_SPHERICAL_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"

namespace initializer
{
	/**
	 * adaptive initialization similar to Dasgupta&Schulman.
	 */
	Parameters adaptiveSphericalGMM(commonutil::DataSet const& dataset, idx_type k, std::mt19937& gen);
}

/**
 * @brief Descriptor for Initialization inspired by Dasgupta & Schulman
 */
class AdaptiveSphericalID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	AdaptiveSphericalID(uint32_t = 1);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const override;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
};

#endif
