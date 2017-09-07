#include "uniform_lloyds_means2gmm.h"

#include "../../kmeans/lloyds.h"
#include "../../kmeans/init/uniform_means.h"
#include "../../../base.h"

const std::string UniformLloydsMeans2GMMID::CLASSTAG = "UniformLloydsMeans2GMM";

UniformLloydsMeans2GMMID::UniformLloydsMeans2GMMID(unsigned int substeps, uint32_t s) : seed(s), substeps(substeps)
{
	std::stringstream sstream;
	sstream << "UniformLloydsMeans2GMM" << "_s" << substeps << "_i" << s ;
	this->nametag = sstream.str();
}

void UniformLloydsMeans2GMMID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> initSubstepsList = readIntList(iss);
	std::vector<unsigned int> sList = readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		for(idx_type sst = 0; sst < initSubstepsList.size(); ++sst)
			descriptions.push_back(new UniformLloydsMeans2GMMID(initSubstepsList[sst], sList[s]));
}


Parameters UniformLloydsMeans2GMMID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::uniformLloydsMeans2GMM(input, k, gen, substeps);
}


Parameters initializer::uniformLloydsMeans2GMM(const commonutil::DataSet& input, idx_type k, std::mt19937& gen, unsigned int substeps)
{
	Matrix uniformMeans = initializer::uniformMeans(input, k, gen);
	Lloyds lloyds(input, false, gen);
	lloyds.init(uniformMeans);
	lloyds.run(substeps);
	return(gmmutil::meansToGMM(input, lloyds.getGMMDesc().means, false, false));
}
