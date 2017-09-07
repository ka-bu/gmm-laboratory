#include "adaptive_lloyds_means2gmm.h"

#include "../utils/gmmutil.h"
#include "../../kmeans/lloyds.h"
#include "../../kmeans/init/adaptive_means.h"
#include "../../../base.h"


Parameters initializer::adaptiveLloydsMeans2GMM(const commonutil::DataSet& input, idx_type k, std::mt19937& gen, unsigned int substeps)
{
	Matrix adaptiveMeans = initializer::adaptiveMeans(input, k, gen);
	Lloyds lloyds(input, false, gen);
	lloyds.init(adaptiveMeans);
	lloyds.run(substeps);
	return(gmmutil::meansToGMM(input, lloyds.getGMMDesc().means, false, false));
}

const std::string AdaptiveLloydsMeans2GMMID::CLASSTAG = "AdaptiveLloydsMeans2GMM";

AdaptiveLloydsMeans2GMMID::AdaptiveLloydsMeans2GMMID(unsigned int substeps, uint32_t s) : seed(s), substeps(substeps)
{
	std::stringstream sstream;
	sstream << "AdaptiveLloydsMeans2GMM" << "_s" << substeps << "_i" << s ;
	this->nametag = sstream.str();
}

Parameters AdaptiveLloydsMeans2GMMID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::adaptiveLloydsMeans2GMM(input, k, gen, substeps);
}

void AdaptiveLloydsMeans2GMMID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> initSubstepsList = readIntList(iss);
	std::vector<unsigned int> sList = readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		for(idx_type sst = 0; sst < initSubstepsList.size(); ++sst)
			descriptions.push_back(new AdaptiveLloydsMeans2GMMID(initSubstepsList[sst], sList[s]));
}
