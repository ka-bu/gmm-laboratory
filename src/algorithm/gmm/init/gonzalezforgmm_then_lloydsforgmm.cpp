#include "gonzalezforgmm_then_lloydsforgmm.h"

#include "gonzalez_forgmm.h"
#include "../lloydsforgmm.h"
#include "../../../base.h"

const std::string GonzalezForGMMThenLloydsForGMMID::CLASSTAG = "GonzalezForGMM_LloydsForGMM";

GonzalezForGMMThenLloydsForGMMID::GonzalezForGMMThenLloydsForGMMID(bool use2GMM, bool spherical, fp_type sampleFactor,unsigned int substeps, bool sphericalLloyds, uint32_t s) 
	: seed(s), use2GMM(use2GMM), spherical(spherical), sampleFactor(sampleFactor), substeps(substeps), sphericalLloyds(sphericalLloyds)
{
	std::stringstream sstream;
	sstream << "GonzalezForGMM_";
	if(sphericalLloyds)
		sstream << "Sph";
	sstream << "LloydsForGMM";
	if(spherical)
		sstream << "_spherical";
	sstream << "_2mm" << use2GMM << "_f" << sampleFactor << "_s" << substeps << "_i" << s ;
	this->nametag = sstream.str();
}

void GonzalezForGMMThenLloydsForGMMID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> initUse2GMMList = readIntList(iss);
	bool initSphericalFlag = (readInt(iss)==1);
	std::vector<fp_type> initSampleSizeFactorList = readFpTypeList(iss);
	std::vector<unsigned int>initSubstepsList = readIntList(iss);
	bool initSphericalLloydsFlag = (readInt(iss)==1);
	std::vector<unsigned int>sList = readIntList(iss);
			
	for(idx_type s = 0; s<sList.size(); ++s)
		for(idx_type u = 0; u<initUse2GMMList.size(); ++u)
			for(idx_type sf = 0; sf<initSampleSizeFactorList.size(); ++sf)
				for(idx_type sst =0; sst<initSubstepsList.size(); ++sst)
					descriptions.push_back(new GonzalezForGMMThenLloydsForGMMID(initUse2GMMList[u]>0, initSphericalFlag, initSampleSizeFactorList[sf], initSubstepsList[sst], initSphericalLloydsFlag, sList[s]));
			
}


Parameters GonzalezForGMMThenLloydsForGMMID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::gonzalezForGMMThenLloydsForGMM(input, k, gen, use2GMM, spherical, sampleFactor, substeps, sphericalLloyds);
}


Parameters initializer::gonzalezForGMMThenLloydsForGMM(const commonutil::DataSet& input, idx_type k, std::mt19937& gen, bool use2GMM, bool spherical, fp_type sampleFactor, unsigned int substeps, bool sphericalLloyds)
{
	Parameters gonzalez = initializer::gonzalezForGMM(input, k, gen, use2GMM, spherical, sampleFactor);
	LloydsForGMM lloyds(input, false, gen, -1., sphericalLloyds);
	lloyds.init(gonzalez);
	lloyds.run(substeps);
	return(lloyds.getGMMDesc());
}