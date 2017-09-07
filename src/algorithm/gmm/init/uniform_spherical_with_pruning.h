#ifndef UNIFORM_SPHERCIAL_WITH_PRUNING_H
#define UNIFORM_SPHERCIAL_WITH_PRUNING_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"

namespace initializer
{
	/**
	 * Gonzales-like pruning according to Dasgupta&Schulman.
	 */
	Parameters prune(Parameters& desc, idx_type k, std::mt19937& gen);
	
	/**
	 * initialization according to Dasgupta&Schulman.
	 * 
	 * @param oversamplingFactor determines how many components are sampled before pruning, i.e., oversamplingFactor*k*ln(oversamplingFactor*k)
	 */
	Parameters uniformSphericalGMMwithPruning(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, unsigned int oversamplingFactor = 2);
}


/**
 * @brief Descriptor for Initialization by Dasgupta & Schulman with pruning!
 */
class UniformSphericalWithPruningID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	UniformSphericalWithPruningID(unsigned int oversamplingFactor, uint32_t seed = 1);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	unsigned int oversamplingFactor;
};

#endif
