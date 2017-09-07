#include "resample_gmminit.h"

#include <boost/algorithm/string.hpp>

#include "adaptive_spherical.h"


#include "../../../base.h"
#include "../../../settings/configparser.h"

const std::string ResampleGMMInitID::CLASSTAG = "ResampleGMMInit";

ResampleGMMInitID::ResampleGMMInitID(std::string initmethod, uint32_t s) : seed(s), initmethod(initmethod)
{
	std::stringstream sstream;
	boost::to_lower(initmethod);
	sstream << "ResampleGMMInit" << "_" << initmethod << "_i" << s;
	this->nametag = sstream.str();
}

void ResampleGMMInitID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::string initSubInitMethod = readString(iss);
	std::vector<unsigned int> sList = readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		descriptions.push_back(new ResampleGMMInitID(initSubInitMethod, sList[s]));
}


Parameters ResampleGMMInitID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::resampleGMMInit(input, this->initmethod, k, gen);
}

Parameters initializer::resampleGMMInit(commonutil::DataSet const& input, std::string initmethod, idx_type k, std::mt19937& gen)
{	
	idx_type d = input.points.rows();
	idx_type initializations = 2*d;
	
	initializer::check(k, d, input.points.cols());
	
	commonutil::DataSet sample_dataset;
	sample_dataset.points= Matrix::Zero(d,initializations*k);	
	for(idx_type i=0; i<initializations; ++i)
		if(initmethod ==  AdaptiveSphericalID::CLASSTAG)
				sample_dataset.points.block(0, i*k, d, k) =  initializer::adaptiveSphericalGMM(input, k, gen).means;
		else
			gmmlab_throw("Unsupported initmethod.");
	
	sample_dataset.weights = Vector::Ones(sample_dataset.points.cols());
	
	return initializer::adaptiveSphericalGMM(sample_dataset, k, gen);
}
