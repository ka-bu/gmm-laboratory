#ifndef ALTERNATING_ADAPTIVEMEANS_MEANS2GMM_H
#define ALTERNATING_ADAPTIVEMEANS_MEANS2GMM_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"

namespace initializer
{
	/**
	 * alternately performs two steps: Chooses a point wrt gmmutil::adaptiveDensities and
	 * calls kmeansutil::meansToGMM wrt the sampled point and the means of the current gmm.
	 */
	Parameters alternatingAdaptiveMeanAndMeanToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type factor, bool spherical);
}

/**
 * @brief Descriptor for Initialization
 */
class AlternatingAdaptiveMeansAndMeans2GMMID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	AlternatingAdaptiveMeansAndMeans2GMMID(fp_type factor, bool spherical, uint32_t seed);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	fp_type factor;
	bool spherical;
};

#endif
