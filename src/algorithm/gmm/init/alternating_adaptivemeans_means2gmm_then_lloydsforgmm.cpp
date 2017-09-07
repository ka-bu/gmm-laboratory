#include "alternating_adaptivemeans_means2gmm_then_lloydsforgmm.h"

#include "alternating_adaptivemeans_means2gmm.h"
#include "../lloydsforgmm.h"
#include "../../../base.h"

const std::string AlternatingAdaptiveMeans2GMMThenLloydsForGMMID::CLASSTAG = "Alt_AdaptMean+Means2GMM_LLoydsForGMM";

AlternatingAdaptiveMeans2GMMThenLloydsForGMMID::AlternatingAdaptiveMeans2GMMThenLloydsForGMMID(fp_type factor, bool spherical, unsigned int substeps, bool sphericalLloyds, uint32_t s) 
	: seed(s), factor(factor), substeps(substeps), spherical(spherical), sphericalLloyds(sphericalLloyds)
{
	std::stringstream sstream;
	sstream << "Alt_AdaptMean+Means2GMM_";
	if(sphericalLloyds)
		sstream << "Sph";
	sstream << "LloydsForGMM";
	if(spherical)
		sstream << "_spherical";
	sstream << "_f" << factor << "_s" << substeps << "_i" << s ; // beware: factor is not used for lloyds algo!
	this->nametag = sstream.str();
}

Parameters AlternatingAdaptiveMeans2GMMThenLloydsForGMMID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::alternatingAdaptiveMeanAndMeanToGMMThenLloydsForGMM(input, k, gen, factor, spherical, substeps, sphericalLloyds);
}

void AlternatingAdaptiveMeans2GMMThenLloydsForGMMID::parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<fp_type> initNonUniformFactorList = readFpTypeList(iss);
	bool initSphericalFlag = (readInt(iss)==1);
	std::vector<unsigned int> initSubstepsList = readIntList(iss);
	bool initSphericalLloydsFlag = (readInt(iss)==1);
	std::vector<unsigned int> sList = readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		for(idx_type sst = 0; sst < initSubstepsList.size(); ++sst)
			for(idx_type f = 0; f<initNonUniformFactorList.size(); ++f)
				descriptions.push_back(new AlternatingAdaptiveMeans2GMMThenLloydsForGMMID(initNonUniformFactorList[f], initSphericalFlag, initSubstepsList[sst], initSphericalLloydsFlag, sList[s]));

}

Parameters initializer::alternatingAdaptiveMeanAndMeanToGMMThenLloydsForGMM(const commonutil::DataSet& input, idx_type k, std::mt19937& gen, fp_type factor, bool spherical, unsigned int substeps, bool sphericalLloyds)
{	
	Parameters adaptive = initializer::alternatingAdaptiveMeanAndMeanToGMM(input, k, gen, factor, spherical);
	LloydsForGMM lloyds(input, false, gen, -1., sphericalLloyds);
	lloyds.init(adaptive);
	lloyds.run(substeps);
	return(lloyds.getGMMDesc());
}