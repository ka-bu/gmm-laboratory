#ifndef ADAPTIVE_LLOYDS_MEANS2GMM_H
#define ADAPTIVE_LLOYDS_MEANS2GMM_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>

namespace initializer
{
	/**
	 * uses adaptiveMeans, then runs substeps many steps of the classical Lloyds algorithm and finally uses means2GMM.
	 */
	Parameters adaptiveLloydsMeans2GMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, unsigned int substeps);
}


/**
 * @brief Descriptor for Initialization
 */
class AdaptiveLloydsMeans2GMMID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	AdaptiveLloydsMeans2GMMID(unsigned int substeps, uint32_t seed = 1);
	virtual Parameters compute(commonutil::DataSet const& input, unsigned int k) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	unsigned int substeps;
};

#endif
