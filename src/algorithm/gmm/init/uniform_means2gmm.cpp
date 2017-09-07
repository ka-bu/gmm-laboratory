#include "uniform_means2gmm.h"

#include "../../../base.h"
#include "../../kmeans/init/uniform_means.h"


const std::string UniformMeans2GMMID::CLASSTAG = "UniformMeans2GMM";

UniformMeans2GMMID::UniformMeans2GMMID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "UniformMeans2GMM" << "_i" << s;
	this->nametag = sstream.str();
}

void UniformMeans2GMMID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		descriptions.push_back(new UniformMeans2GMMID(sList[s]));
}


Parameters UniformMeans2GMMID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::uniformMeansToGMM(input, k, gen);
}


Parameters initializer::uniformMeansToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	return gmmutil::meansToGMM(input, initializer::uniformMeans(input,k,gen), false, false);
}