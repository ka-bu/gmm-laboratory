#include "gonzalez2gmm.h"

#include "../../kmeans/init/gonzalez_means.h"
#include "../../../base.h"

const std::string Gonzalez2GMMID::CLASSTAG = "Gonzalez2GMM";

Gonzalez2GMMID::Gonzalez2GMMID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "Gonzalez2GMM" << "_i" << s;
	this->nametag = sstream.str();
}

Parameters Gonzalez2GMMID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::gonzalezToGMM(input, k, gen);
}

void Gonzalez2GMMID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	for (idx_type s=0; s<sList.size(); ++s)
		descriptions.push_back(new Gonzalez2GMMID(sList[s]));
}


Parameters initializer::gonzalezToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{	
	return gmmutil::meansToGMM(input, initializer::gonzalez(input, k, gen), false, false);
}
