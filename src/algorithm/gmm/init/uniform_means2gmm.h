#ifndef UNIFORM_MEANS2GMM_H
#define UNIFORM_MEANS2GMM_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"

namespace initializer
{
	/**
	 * computes means according to initutil::uniformMeans, then uses these means to
	 * compute a GMM by using kmeans::meansToGMM.
	 */
	Parameters uniformMeansToGMM(commonutil::DataSet const& input, idx_type k,  std::mt19937& gen);
}

/**
 * @brief Descriptor for Initialization
 */
class UniformMeans2GMMID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	UniformMeans2GMMID(uint32_t = 1);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
};

#endif
