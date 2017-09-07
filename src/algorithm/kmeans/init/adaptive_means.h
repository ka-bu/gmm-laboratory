#ifndef ADAPTIVE_MEANS_H
#define ADAPTIVE_MEANS_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>

namespace initializer
{
	/**
	 * chooses k means according to kmeans++.
	 */
	Matrix adaptiveMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen);
	
	/**
	 * computes means according to initutil::adaptiveMeans, then uses these means
	 * as centers of a GMM by using kmeansutil::wrapMeans.
	 */
	Parameters wrappedAdaptiveMeans(commonutil::DataSet const& dataset, idx_type k, std::mt19937& gen, bool computeWeightAndCovar = 1);
}

/**
 * @brief Descriptor for Initialization
 */
class AdaptiveMeansID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	AdaptiveMeansID(uint32_t = 1);
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
};

#endif
