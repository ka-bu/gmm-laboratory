#ifndef GONZALEZ_KWEDLO_H
#define GONZALEZ_KWEDLO_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"

namespace initializer
{
	/**
	 * an version of the Gonzalez algorithm that draws means wrt maximum likelihood and uses random covariances.
	 */
	Parameters gonzalezKwedlo(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type sampleFactor);
}

/**
 * @brief Descriptor for Initialization
 */
class GonzalezKwedloID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	GonzalezKwedloID(uint32_t seed, fp_type sampleSizeFactor);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
	fp_type sampleSizeFactor;
};

#endif
