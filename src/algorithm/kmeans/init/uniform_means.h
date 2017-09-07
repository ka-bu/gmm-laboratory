#ifndef UNIFORM_MEANS_H
#define UNIFORM_MEANS_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>

namespace initializer
{
	/**
	 * chooses k means uniformly at random from the given input set.
	 */
	Matrix uniformMeans(commonutil::DataSet const& input, idx_type k,  std::mt19937& gen,  std::vector<idx_type> const& indices_of_points_to_be_used = std::vector<idx_type>());
	
	/**
	 * computes means according to initutil::uniformMeans, then uses these means
	 * as centers of a GMM by using kmeans::wrapMeans. 
	 */
	Parameters wrappedUniformMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, bool computeWeightAndCovar = 1);
}

/**
 * @brief Descriptor for Initialization
 */
class UniformMeansID : public InitDescription
{
public:
	static const std::string CLASSTAG;
	UniformMeansID(uint32_t = 1);	
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring);

private:
	uint32_t seed;
};

#endif
