#include "restart_init.h"


#include "alternating_adaptivemeans_means2gmm.h"
#include "../emforgmm.h"
#include "../../../base.h"
#include "../../../settings/configparser.h"
#include <boost/algorithm/string.hpp>

const std::string RestartInitID::CLASSTAG = "RestartInit";

RestartInitID::RestartInitID(std::string initmethod, unsigned int subruns, unsigned int substeps, uint32_t s, fp_type factor) : seed(s), subruns(subruns), substeps(substeps), factor(factor), initmethod(initmethod)
{
	std::stringstream sstream;
	if(this->initmethod == AlternatingAdaptiveMeansAndMeans2GMMID::CLASSTAG)
		sstream << "RestartInit_" << initmethod << "_f" << factor << "_subruns" << subruns << "_substeps" << substeps << "_i" << s;
	else
		sstream << "RestartInit_" << initmethod << "_subruns" << subruns << "_substeps" << substeps << "_i" << s;
	this->nametag = sstream.str();
}

void RestartInitID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	std::istringstream iss(descriptionstring);
	std::string initSubInitMethod = readString(iss);
	std::vector<unsigned int> initSubrunsList = readIntList(iss);
	std::vector<unsigned int> initSubstepsList = readIntList(iss);
	std::vector<fp_type> initNonUniformFactorList;
	if(initSubInitMethod == AlternatingAdaptiveMeansAndMeans2GMMID::CLASSTAG)
		initNonUniformFactorList = readFpTypeList(iss);
	std::vector<unsigned int> sList = readIntList(iss);
	
	for(idx_type s=0; s<sList.size(); ++s)
		for(idx_type sr = 0; sr<initSubrunsList.size(); ++sr)
			for(idx_type sst = 0; sst<initSubstepsList.size(); ++sst)
			{
				if(initSubInitMethod == AlternatingAdaptiveMeansAndMeans2GMMID::CLASSTAG)
					for(idx_type f = 0; f< initNonUniformFactorList.size(); ++f)
						descriptions.push_back(new RestartInitID(initSubInitMethod, initSubrunsList[sr], initSubstepsList[sst], sList[s], initNonUniformFactorList[f]));
				else
					descriptions.push_back(new RestartInitID(initSubInitMethod, initSubrunsList[sr], initSubstepsList[sst], sList[s]));
			}
}


Parameters RestartInitID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	std::mt19937 gen(seed);
	return initializer::restartInit(input, k, gen, this->initmethod, subruns, substeps, factor);
}

Parameters initializer::restartInit(commonutil::DataSet const& input, idx_type k, std::mt19937& gen, std::string initmethod, idx_type subruns, idx_type substeps, fp_type factor)
{
	gmmlab_throw("To be implemented!");
		
	initializer::check(k, input.points.rows(), input.points.cols());
	
	if(subruns == 0)
		subruns = k;
	
	Algorithm* em = new EMforGMM(input, 0, gen);
	Parameters tmp;
	Parameters bestGMMDesc;
	fp_type minCosts;
// 	for(idx_type t=0; t < subruns; ++t)
// 	{
// 		switch(initmethod){
// 			case UNIFORM_MEANS2GMM:
// 				tmp = initutil::uniformMeansToGMM(input, k, gen);
// 				break;
// 			case ADAPTIVE_MEANS2GMM:
// 				tmp = initutil::adaptiveMeansToGMM(input, k, gen);
// 				break;
// 			case GONZALEZ2GMM:
// 				tmp = initutil::gonzalezToGMM(input, k, gen);
// 				break;
// 			case UNIFORM_SPHERICAL:
// 				tmp = initutil::uniformSphericalGMM(input, k, gen);
// 				break;
// 			case ADAPTIVE_SPHERICAL:
// 				tmp = initutil::adaptiveSphericalGMM(input, k, gen);
// 				break;
// 			case EXPERIMENTAL_ADAPTIVE:
// 				tmp = initutil::experimentalAdaptiveGMM(input, k, gen);
// 				break;
// 			case ALTERNATELY_ADAPTIVEMEANS_MEANS2GMM:
// 				tmp = initutil::alternatingAdaptiveMeanAndMeanToGMM(input, k, gen, factor, false);
// 				break;
// 			default:
// 				gmmlab_throw("Unsupported initmethod.");
// 		}		
// 		em->init(tmp);
// 		em->run(substeps);
// 		tmp = em->getGMMDesc();
// 		fp_type tmpCosts = gmmutil::nll(input, tmp);
// 		if(t==0 || (tmpCosts < minCosts))
// 		{
// 			minCosts = tmpCosts;
// 			bestGMMDesc = tmp;
// 		}
// 	}
	delete em;
	
	//std::cout << "restartUniformMeans2GMMEM = \n" << bestGMMDesc << std::endl;
	
	return bestGMMDesc;
}
