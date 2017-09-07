#ifndef GONZALEZ_FORGMM_H
#define GONZALEZ_FORGMM_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"

namespace initializer
{
	/**
	 * an version of the Gonzalez algorithm that draws means wrt maximum likelihood and uses means2gmm to estimate the covariances.
	 */
	Parameters gonzalezForGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool use2GMM, bool spherical, fp_type sampleFactor);	
}


/**
 * @brief Descriptor for Initialization
 */
class GonzalezForGMMID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	GonzalezForGMMID(uint32_t seed, bool useGMM, bool spherical, fp_type sampleSizeFactor);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	bool use2GMM;
	fp_type sampleSizeFactor;
	bool spherical;
};

#endif
