#include "alternating_adaptivemeans_means2gmm.h"

#include "../../../gmmlab.h"

#include "../../../base.h"

const std::string AlternatingAdaptiveMeansAndMeans2GMMID::CLASSTAG = "Alt_AdaptMean+Means2GMM";

AlternatingAdaptiveMeansAndMeans2GMMID::AlternatingAdaptiveMeansAndMeans2GMMID(fp_type factor, bool spherical, uint32_t s) : seed(s), factor(factor), spherical(spherical)
{
	std::stringstream sstream;
	sstream << CLASSTAG;
	if(spherical)
		sstream << "_spherical";
	sstream << "_f" << factor << "_i" << s ;
	this->nametag = sstream.str();
}

Parameters AlternatingAdaptiveMeansAndMeans2GMMID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::alternatingAdaptiveMeanAndMeanToGMM(input, k, gen, factor, spherical);
}

void AlternatingAdaptiveMeansAndMeans2GMMID::parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<fp_type> initNonUniformFactorList = readFpTypeList(iss);
	bool initSphericalFlag = (readInt(iss)==1);
	std::vector<unsigned int> sList = readIntList(iss);
	for(idx_type f = 0; f<initNonUniformFactorList.size(); ++f)
		for(idx_type s=0; s<sList.size(); ++s)
			descriptions.push_back(new AlternatingAdaptiveMeansAndMeans2GMMID(initNonUniformFactorList[f], initSphericalFlag, sList[s]));
}

Parameters initializer::alternatingAdaptiveMeanAndMeanToGMM(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, fp_type factor, bool spherical)
{  
	idx_type d = input.points.rows();
	idx_type n = input.points.cols();
	
	initializer::check(k, d, n);

	Parameters desc;
	Vector sample;
	
	for (idx_type i=0; i<k; ++i)
	{
		if (i==0)
			// draw first point uniformly
			sample = input.points.col(commonutil::randomIndex(input.weights, gen)); 
			// ... yeah, it's useless, but now it stays here because otherwise the seeds are screwed up and we get different results
		else
		{
// clock_t t;
			// draw next point w.r.t. current mixture
			Vector densities = gmmutil::adaptiveDensities(input, desc, factor);
			sample = input.points.col(commonutil::randomIndex(densities, gen)); 
			
// std::cout << "adaptiveDensities computation in time = " << secondsSince(t) << std::endl;
		}
		
		Matrix tmpMeans = Matrix::Zero(d,i+1);
		if(i>0)
			tmpMeans.block(0,0,d,i) = desc.means;
		tmpMeans.col(i) = sample;
		
		desc = gmmutil::meansToGMM(input, tmpMeans, spherical, false);
		
	}

	return desc;
}
