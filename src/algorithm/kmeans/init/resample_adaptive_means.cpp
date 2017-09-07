#include "resample_adaptive_means.h"

#include "adaptive_means.h"
#include "../../../settings/settings.h"
#include "../../../base.h"
#include "../utils/kmeansutil.h"

const std::string ResampleAdaptiveMeansID::CLASSTAG = "ResampleAdaptiveMeans";

ResampleAdaptiveMeansID::ResampleAdaptiveMeansID(uint32_t s) : seed(s)
{
	std::stringstream sstream;
	sstream << "ResampleAdaptiveMeans" << "_i" << s;
	this->nametag = sstream.str();
}

void ResampleAdaptiveMeansID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		descriptions.push_back(new ResampleAdaptiveMeansID(sList[s]));
}


Parameters ResampleAdaptiveMeansID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return kmeansutil::wrapMeans(input, initializer::resampleAdaptiveMeans(input, k, gen), !commonSettings().ndisplay);
}


Matrix initializer::resampleAdaptiveMeans(commonutil::DataSet const& input, idx_type k, std::mt19937& gen)
{
	idx_type d = input.points.rows();
	idx_type initializations = 2*d;
	
	initializer::check(k, d, input.points.cols());
	
	commonutil::DataSet sample_dataset;
	sample_dataset.points= Matrix::Zero(d,initializations*k);
	for(idx_type i=0; i<initializations; ++i)
	  sample_dataset.points.block(0, i*k, d, k) =  initializer::adaptiveMeans(input, k, gen);
	
	sample_dataset.weights = Vector::Ones(sample_dataset.points.cols());
	
	return initializer::adaptiveMeans(sample_dataset, k, gen);
}
