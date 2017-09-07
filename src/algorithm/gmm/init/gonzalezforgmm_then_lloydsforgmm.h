#ifndef GONZALEZFORGMM_THEN_LLOYDSFORGMM_H
#define GONZALEZFORGMM_THEN_LLOYDSFORGMM_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"
#include "../utils/gmmutil.h"


namespace initializer
{
	/**
	 * 
	 */
	Parameters gonzalezForGMMThenLloydsForGMM(const commonutil::DataSet& input, idx_type k, std::mt19937& gen, bool use2GMM, bool spherical, fp_type sampleFactor, unsigned int substeps, bool sphericalLloyds);

}

/**
 * @brief Descriptor for Initialization
 */
class GonzalezForGMMThenLloydsForGMMID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	GonzalezForGMMThenLloydsForGMMID(bool use2GMM, bool spherical, fp_type sampleFactor, unsigned int substeps, bool sphericalLloyds, uint32_t seed = 1);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	unsigned int substeps;
	bool use2GMM;
	fp_type sampleFactor;
	bool sphericalLloyds;
	bool spherical;
};

#endif
