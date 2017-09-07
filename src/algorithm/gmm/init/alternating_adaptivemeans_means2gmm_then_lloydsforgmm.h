#ifndef ALTERNATING_ADAPTIVE_MEANS2GMM_THEN_LLOYDSFORGMM_H
#define ALTERNATING_ADAPTIVE_MEANS2GMM_THEN_LLOYDSFORGMM_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"
#include "../utils/gmmutil.h"

namespace initializer
{
	/**
	 * 
	 */
	Parameters alternatingAdaptiveMeanAndMeanToGMMThenLloydsForGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type factor, bool spherical, unsigned int substeps, bool sphericalLloyds);
	
}

/**
 * @brief Descriptor for Initialization
 */
class AlternatingAdaptiveMeans2GMMThenLloydsForGMMID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	AlternatingAdaptiveMeans2GMMThenLloydsForGMMID(fp_type factor, bool spherical, unsigned int substeps, bool sphericalLloyds, uint32_t seed = 1);	
	virtual Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	unsigned int substeps;
	fp_type factor;
	bool spherical;
	bool sphericalLloyds;
};

#endif
