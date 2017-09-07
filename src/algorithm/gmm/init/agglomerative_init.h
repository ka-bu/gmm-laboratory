#ifndef AGGLOMERATIVE_INIT_H
#define AGGLOMERATIVE_INIT_H

#include "../../../base/parameters.h"
#include "../../description/init_description.h"

#include <vector>
#include "../utils/gmmutil.h"


namespace initializer
{
	/**
	 * computes average linkage clustering. 
	 * TODO: only works with unweighted points
	 */
	Matrix agglomerativeMeans(commonutil::DataSet const&, idx_type, bool = false);
	
	/**
	 * runs agglomerativeMeans and calls kmeansuitl::meansToGMM wrt the input set and the returned means.
	 */
	Parameters sampleAgglomerativeMeansToGMM(const commonutil::DataSet& input, const idx_type k, const fp_type sampleSizeFactor, const bool precompute, std::mt19937& gen);	
}

/**
 * @brief Descriptor for Initialization
 */
class AgglomerativeInitID : public InitDescription
{
public:
	const static std::string CLASSTAG;
	AgglomerativeInitID(fp_type sampleSizeFactor, uint32_t = 1);
	Parameters compute(commonutil::DataSet const&, unsigned int) const;
	static void parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring, bool first = false);
	
private:
	uint32_t seed;
	fp_type sampleSizeFactor;
};

#endif
