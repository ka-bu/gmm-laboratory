#ifndef RESAMPLE_ADAPTIVE_MEANS_H
#define RESAMPLE_ADAPTIVE_MEANS_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../../gmm/utils/gmmutil.h"

namespace initializer
{
	/**
	 * computes 2*d adaptive initializations according to kmeans++ (independent of each other),
	 * then samples k centers from these 2*d*k sampled centers according to kmeans++
	 */
	Matrix resampleAdaptiveMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen);
}

/**
 * @brief Descriptor for Initialization
 */
class ResampleAdaptiveMeansID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	ResampleAdaptiveMeansID(uint32_t = 1);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
};

#endif
