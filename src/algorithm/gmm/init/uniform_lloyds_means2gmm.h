#ifndef UNIFORM_LLOYDS_MEANS2GMM_H
#define UNIFORM_LLOYDS_MEANS2GMM_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"


namespace initializer
{
	/**
	 * uses uniformMeans, then runs substeps many steps of the classical Lloyds algorithm and finally uses means2GMM.
	 */
	Parameters uniformLloydsMeans2GMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, unsigned int substeps);

}

/**
 * @brief Descriptor for Initialization
 */
class UniformLloydsMeans2GMMID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	UniformLloydsMeans2GMMID(unsigned int substeps, uint32_t seed = 1);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	unsigned int substeps;
};

#endif
