#include "samplingalgorithm.h"

#include "../../base/linalgutil.h"


SamplingAlgorithm::SamplingAlgorithm(commonutil::DataSet const& ds, bool v, uint32_t s, fp_type r)
	: RandomAlgorithm(ds,v,s), ratio(r)
{
	this->samples.points.resize(this->input->points.rows(),0);
	this->samples.weights.resize(0);
}

commonutil::DataSet const& SamplingAlgorithm::getSamples() const
{
	return this->samples;
}

SamplingAlgorithm* SamplingAlgorithm::toSamplingAlgorithm(Algorithm* a)
{
	return dynamic_cast<SamplingAlgorithm*>(a);
}

