#include "alternating_adaptivemeans_means2gmm_then_lloyds.h"

#include "alternating_adaptivemeans_means2gmm.h"
#include "../../kmeans/lloyds.h"
#include "../../../base.h"

const std::string AlternatingAdaptiveMeans2GMMThenLloydsID::CLASSTAG = "Alt_AdaptMean+Means2GMM_Lloyds";

AlternatingAdaptiveMeans2GMMThenLloydsID::AlternatingAdaptiveMeans2GMMThenLloydsID(fp_type factor, bool spherical, unsigned int substeps, uint32_t s) 
	: seed(s), factor(factor), substeps(substeps), spherical(spherical)
{
	std::stringstream sstream;
	sstream << "Alt_AdaptMean+Means2GMM_Lloyds";
	if(spherical)
		sstream << "_spherical";
	sstream << "_f" << factor << "_s" << substeps << "_i" << s ;
	this->nametag = sstream.str();
}

Parameters AlternatingAdaptiveMeans2GMMThenLloydsID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::alternatingAdaptiveMeanAndMeanToGMMThenLloyds(input, k, gen, factor, spherical, substeps);
}

void AlternatingAdaptiveMeans2GMMThenLloydsID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<fp_type> initNonUniformFactorList = readFpTypeList(iss);
	bool initSphericalFlag = (readInt(iss)==1);
	std::vector<unsigned int> initSubstepsList = readIntList(iss);
	std::vector<unsigned int> sList = Description::readIntList(iss);
	
	for(idx_type s=0; s<sList.size(); ++s)
		for(idx_type sst = 0; sst < initSubstepsList.size(); ++sst)
			for(idx_type f = 0; f<initNonUniformFactorList.size(); ++f)
				descriptions.push_back(new AlternatingAdaptiveMeans2GMMThenLloydsID(initNonUniformFactorList[f], initSphericalFlag, initSubstepsList[sst], sList[s]));
}

Parameters initializer::alternatingAdaptiveMeanAndMeanToGMMThenLloyds(const commonutil::DataSet& input, idx_type k, std::mt19937& gen, fp_type factor,bool spherical, unsigned int substeps)
{
	Parameters adaptive = initializer::alternatingAdaptiveMeanAndMeanToGMM(input, k, gen, factor, spherical);
	Lloyds lloyds(input, false, gen, -1.);
	lloyds.init(adaptive.means);
	lloyds.run(substeps);
	return(gmmutil::meansToGMM(input, lloyds.getGMMDesc().means, false, false));
}
