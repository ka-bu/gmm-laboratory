#include "restartlloyds.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <ctime>

#include "lloyds.h"

#include "init/adaptive_means.h"
#include "../../settings/settings.h"
#include "../../base/parameters.h"
#include "utils/kmeansutil.h"

const std::string RestartLloydsAD::CLASSTAG = "RestartLloyds";

RestartLloydsAD::RestartLloydsAD(unsigned int steps, uint32_t seed, fp_type lloydsRatio, unsigned int substeps) : 
						         seed(seed), subSteps(substeps), lloydsRatio(lloydsRatio)
{
	this->steps = steps;
	std::stringstream sstream;	
	sstream << "RestartLloyds" << "_s" << steps << "_a" << seed << "_llr" << lloydsRatio << "_st" << subSteps;
	this->nametag = sstream.str();
}

Algorithm* RestartLloydsAD::create(const commonutil::DataSet& input) const
{
	return new RestartLloyds(input, commonSettings().verbose, this->subSteps, this->seed, this->lloydsRatio);
}

void RestartLloydsAD::parseDescriptionString(std::vector< AlgoDescription* >& descriptions, const std::string descriptionstring, const bool first)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> stList = readIntList(iss);
	fp_type llr = readFpType(iss);
	std::vector<unsigned int> substList = readIntList(iss);
	std::vector<unsigned int> sdList = readIntList(iss);
	
	if (first)	
		for (std::size_t st=0; st<stList.size(); ++st)
			for (std::size_t sd=0; sd<sdList.size(); ++sd)
				for (std::size_t subst=0; subst<substList.size(); ++subst)
					descriptions.push_back(new RestartLloydsAD(stList[st], sdList[sd], llr, substList[subst]));
	else
		gmmlab_throw("next unsupported");
}


// =================================================================================================================================


RestartLloyds::RestartLloyds(commonutil::DataSet const& ds, bool v,  unsigned int subst, uint32_t s, fp_type lloydsRatio) 
				:  RandomAlgorithm(ds, v, s), subalgoSteps(subst), lloydsRatio(lloydsRatio)
{
}

void RestartLloyds::init(Parameters const& desc)
{
	Algorithm::init(desc);
	if(this->subalgoSteps == 0) 
	  this->subalgoSteps = 2*this->desc.means.cols();
	this->minCost = kmeansutil::kmeanscost(*this->input,this->desc.means);
	this->justInitialized = true;
}

void RestartLloyds::run(unsigned int numSteps)
{
  
    //std::cout << "substeps = " << this->subalgoSteps << std::endl;
  	
    clock_t start, end;
    start = clock();

    idx_type k = this->desc.means.cols();

    if (k==0)
        gmmlab_throw("No or empty initial solution given.");

    idx_type n = this->input->points.cols();
    idx_type d = this->input->points.rows();

    if (n==0 || d==0)
        gmmlab_throw("Input is empty.");    
    
    //std::uniform_int_distribution<uint32_t> uid(0, 100000);
    for (unsigned int step=0; step<numSteps; ++step)
    {
     
        // Usually, we initialize our algorithms with seeds 1,2,3,... . If RestartLLoyds would just increase its seed
        // by one and reuse it for Lloyds, then the different calls of RestartLLoyds would mainly compute the same,
        // i.e. most of their calls of Lloyds would use the same seeds.
        // Thus, here we increase the seed by 1000.
        
		//uint32_t newSeed = uid(this->gen);	
		this->seed += 1000;
		
		Algorithm* lloyds = new Lloyds(*this->input, commonSettings().verbose, this->seed, this->lloydsRatio);
		
		
		if(!this->justInitialized)
			// compute a new initial solution
			lloyds->init(initializer::wrappedAdaptiveMeans(*this->input, k, this->gen, 1));
		else
		{
			// use initial solution computed by restartlloyds::init
			this->justInitialized = false; 
			lloyds->init(this->desc);
		}
		
		lloyds->run(this->subalgoSteps);
		
		this->change = true;
		
		fp_type curCost = kmeansutil::kmeanscost(*this->input,lloyds->getGMMDesc().means);
		if(curCost < this->minCost)
		{
			this->minCost = curCost;
			this->desc = lloyds->getGMMDesc();	
		}
		
		delete lloyds;
    }
        
    end = clock();
    this->runtime += double(end-start)/CLOCKS_PER_SEC;
}



RestartLloyds* RestartLloyds::toRestartLloyds(Algorithm* a)
{
    return dynamic_cast<RestartLloyds*>(a);
}

