#include "gonzalez_lloyds_means2gmm.h"

#include "../../kmeans/lloyds.h"
#include "../../kmeans/init/gonzalez_means.h"
#include "../../../base.h"

const std::string GonzalezLloydsMeans2GMMID::CLASSTAG = "GonzalezLloydsMeans2GMM";

GonzalezLloydsMeans2GMMID::GonzalezLloydsMeans2GMMID(unsigned int substeps, uint32_t s) : seed(s), substeps(substeps)
{
	std::stringstream sstream;
	sstream << "GonzalezLloydsMeans2GMM" << "_s" << substeps << "_i" << s ;
	this->nametag = sstream.str();
}

void GonzalezLloydsMeans2GMMID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> initSubstepsList = readIntList(iss);
	std::vector<unsigned int> sList = readIntList(iss);
	for (idx_type s=0; s<sList.size(); ++s)
		for(idx_type sst = 0; sst < initSubstepsList.size(); ++sst)
			descriptions.push_back(new GonzalezLloydsMeans2GMMID(initSubstepsList[sst], sList[s]));
}


Parameters GonzalezLloydsMeans2GMMID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::gonzalezLloydsMeans2GMM(input, k, gen, substeps);
}

Parameters initializer::gonzalezLloydsMeans2GMM(const commonutil::DataSet& input, idx_type k, std::mt19937& gen, unsigned int substeps)
{
	Matrix gonzalezMeans = initializer::gonzalez(input, k, gen);
	Lloyds lloyds(input, false, gen);
	lloyds.init(gonzalezMeans);
	lloyds.run(substeps);
	return(gmmutil::meansToGMM(input, lloyds.getGMMDesc().means, false, false));
}



