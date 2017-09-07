#ifndef UNIFORM_SPHERICAL_H
#define UNIFORM_SPHERICAL_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"


namespace initializer
{
	/**
	 * initialization according to Dasgupta&Schulman but without pruning.
	 */
	Parameters uniformSphericalGMM(commonutil::DataSet const& dataset, idx_type k, std::mt19937& gen);
}

/**
 * @brief Descriptor for Initialization by Dasgupta & Schulman
 */
class UniformSphericalID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	UniformSphericalID(uint32_t = 1);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
};

#endif
