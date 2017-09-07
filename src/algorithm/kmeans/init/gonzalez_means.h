#ifndef GONZALEZ_MEANS_H
#define GONZALEZ_MEANS_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>


namespace initializer
{
	/**
	 * implements Gonazalez algorithm.
	 */
	Matrix gonzalez(commonutil::DataSet const& input, idx_type k, std::mt19937& gen);
	

	
	/**
	 * runs Gonzalez algorithm and wraps the Means using kmeansutil::wrapMeans.
	 */
	Parameters wrappedGonzalez(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar);
}

/**
 * @brief Descriptor for Initialization
 */
class GonzalezMeansID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	GonzalezMeansID(uint32_t = 1);
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
};

#endif