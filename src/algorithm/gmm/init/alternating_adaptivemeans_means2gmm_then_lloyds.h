#ifndef ALTERNATING_ADAPTIVE_MEANS2GMM_THEN_LLOYDS_H
#define ALTERNATING_ADAPTIVE_MEANS2GMM_THEN_LLOYDS_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"
#include "../utils/gmmutil.h"

namespace initializer
{
	/**
	 * 
	 */
	Parameters alternatingAdaptiveMeanAndMeanToGMMThenLloyds(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type factor, bool spherical, unsigned int substeps);
}	


/**
 * @brief Descriptor for Initialization
 */
class AlternatingAdaptiveMeans2GMMThenLloydsID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	AlternatingAdaptiveMeans2GMMThenLloydsID(fp_type factor, bool spherical, unsigned int substeps, uint32_t seed = 1);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	unsigned int substeps;
	fp_type factor;
	bool spherical;
};

#endif
