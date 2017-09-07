#include "adaptive_means2gmm.h"

#include "../utils/gmmutil.h"
#include "../../kmeans/init/adaptive_means.h"
#include "../../../base.h"

Parameters initializer::adaptiveMeansToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	return gmmutil::meansToGMM(input, initializer::adaptiveMeans(input,k,gen), false, false);
}

const std::string AdaptiveMeans2GMMID::CLASSTAG = "AdaptiveMeans2GMM";

AdaptiveMeans2GMMID::AdaptiveMeans2GMMID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "AdaptiveMeans2GMM" << "_i" << s;
	this->nametag = sstream.str();
}

Parameters AdaptiveMeans2GMMID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::adaptiveMeansToGMM(input, k, gen);
}

void AdaptiveMeans2GMMID::parseDescriptionString(std::vector< InitDescription* > & descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		descriptions.push_back(new AdaptiveMeans2GMMID(sList[s]));
}

