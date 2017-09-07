#ifndef ADAPTIVE_MEANS2GMM_H
#define ADAPTIVE_MEANS2GMM_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>


namespace initializer
{
	/**
	 * computes means according to initutil::adaptiveMeans, then uses these means to
	 * compute a GMM by using kmeansutil::meansToGMM.
	 */
	Parameters adaptiveMeansToGMM(commonutil::DataSet const& dataset, idx_type k, std::mt19937& gen);
}

/**
 * @brief Descriptor for Initialization
 */
class AdaptiveMeans2GMMID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	AdaptiveMeans2GMMID(uint32_t = 1);
	virtual Parameters compute(commonutil::DataSet const& input, unsigned int k) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
};

#endif
