#include "alternating_adaptivemeans_means2gmm_initalgo.h"

#include "../init/alternating_adaptivemeans_means2gmm.h"
#include "../lloydsforgmm.h"
#include "../../../base.h"
#include "../../../settings/settings.h"

const std::string AlternatingAdaptiveMeans2GMMAD::CLASSTAG = "Alt_AdaptMean+Means2GMM";

AlternatingAdaptiveMeans2GMMAD::AlternatingAdaptiveMeans2GMMAD(unsigned int numSteps, bool gonzalez_mode, fp_type nonUniformFactor, uint32_t s) 
	:  gonzalez_mode(gonzalez_mode), nonUniformFactor(nonUniformFactor), seed(s)
{
	this->steps = numSteps;
	std::stringstream sstream;
	sstream << "Alt_AdaptMean+Means2GMM_" << "_f" << nonUniformFactor << "_i" << s ;
	this->nametag = sstream.str();
}

Algorithm* AlternatingAdaptiveMeans2GMMAD::create(commonutil::DataSet const& input) const
{
	std::mt19937 gen(seed);
	return new AlternatingAdaptiveMeans2GMMInitAlgo(input, commonSettings().verbose, gonzalez_mode, nonUniformFactor);
}

void AlternatingAdaptiveMeans2GMMAD::parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first)
{
	std::istringstream iss(descriptionstring);
	unsigned int steps = readInt(iss);
	bool gonzalez_mode = (readInt(iss)==1);
	std::vector<fp_type> initNonUniformFactorList = readFpTypeList(iss);
	std::vector<unsigned int> sList = readIntList(iss);
	for(idx_type s=0; s<sList.size(); ++s)
		for(idx_type f=0; f<initNonUniformFactorList.size(); ++f)
			descriptions.push_back(new AlternatingAdaptiveMeans2GMMAD(steps, gonzalez_mode, initNonUniformFactorList[f], sList[s]));

}

// =================================================================================================================================

AlternatingAdaptiveMeans2GMMInitAlgo::AlternatingAdaptiveMeans2GMMInitAlgo(commonutil::DataSet const& input, bool verbose, bool gonzalez_mode, fp_type nonUniformFactor, uint32_t seed)
 : RandomAlgorithm(input,verbose,seed), nonUniformFactor(nonUniformFactor) 
{
	
}

void AlternatingAdaptiveMeans2GMMInitAlgo::run(unsigned int numSteps)
{
	idx_type d = this->input->points.rows();
	idx_type n = this->input->points.cols();
	
	
	for(idx_type step=0; step<numSteps; ++step)
	{
		idx_type k = this->desc.components();
		Vector sample; 
		
		if (k==0)
			// draw first point uniformly
			sample = this->input->points.col(commonutil::randomIndex(this->input->weights, gen)); 
			// ... yeah, it's useless, but now it stays here because otherwise the seeds are screwed up and we get different results
		else
		{
	// clock_t t;
			// draw next point w.r.t. current mixture
			Vector densities = gmmutil::adaptiveDensities(*this->input, desc, nonUniformFactor);
			if(!gonzalez_mode)
				sample = this->input->points.col(commonutil::randomIndex(densities, gen));
			else
			{
				idx_type index;
				densities.maxCoeff(&index);
				sample = this->input->points.col(index); 
			}
				
	// std::cout << "adaptiveDensities computation in time = " << secondsSince(t) << std::endl;
		}
			
			
		
		Matrix tmpMeans = Matrix::Zero(d,k+1);
		if(k>0)
			tmpMeans.block(0,0,d,k) = desc.means;
		tmpMeans.col(k) = sample;
			
		this->desc = gmmutil::meansToGMM(*this->input, tmpMeans, true, false);
			
	}
}

