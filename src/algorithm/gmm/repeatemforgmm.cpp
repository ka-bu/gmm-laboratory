#include "repeatemforgmm.h"

#include <algorithm>
#include <iostream>
#include <cmath>

#include "init/adaptive_means2gmm.h"
#include "../../base/linalgutil.h"
#include "../../settings/settings.h"



const std::string RepeatEMforGMMAD::CLASSTAG = "RepeatEM";


RepeatEMforGMMAD::RepeatEMforGMMAD(unsigned int st, uint32_t sd, bool spherical, unsigned int substeps) : seed(sd), spherical(spherical), substeps(substeps)
{
	this->steps = st;
	std::stringstream sstream;
	sstream << "RepeatEM" << "_s" << st << (spherical?"_sph":"") << "_subst" << substeps << "_a" << sd;
	this->nametag = sstream.str();
}


void RepeatEMforGMMAD::parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first)
{
	std::istringstream iss(descriptionstring);
	std::vector<unsigned int> stList = Description::readIntList(iss);
	bool spherical = (Description::readInt(iss) == 1);
	std::vector<unsigned int> substList = Description::readIntList(iss);
	std::vector<unsigned int> sdList = Description::readIntList(iss);
	
	if (first)
	{
		for (std::size_t st=0; st<stList.size(); ++st)
			for(std::size_t subst=0; subst<substList.size(); ++subst)
				for (std::size_t sd=0; sd<sdList.size(); ++sd)
					descriptions.push_back(new RepeatEMforGMMAD(stList.at(st), sdList.at(sd), spherical, substList.at(subst)));
	}
	else
	{
		if (stList.size()!=1 || sdList.size()!=1 || substList.size() !=1)
			gmmlab_throw("multiple values forbidden for next:...!");
		descriptions.push_back(new RepeatEMforGMMAD(stList.at(0), sdList.at(0), spherical, substList.at(0)));
		
	}
}

Algorithm* RepeatEMforGMMAD::create(commonutil::DataSet const& input) const
{
	return new RepeatEMforGMM(input, commonSettings().verbose, this->seed, this->substeps, this->spherical);
}

// =================================================================================================================================


RepeatEMforGMM::RepeatEMforGMM(commonutil::DataSet const& ds, bool v, uint32_t s, unsigned int substeps, bool spherical)
	: RandomAlgorithm(ds,v,s), spherical(spherical), substeps(substeps)
{
	emforgmm = new EMforGMM(ds, v, this->gen, spherical);
}

RepeatEMforGMM::RepeatEMforGMM(commonutil::DataSet const& ds, bool v, std::mt19937& gen, unsigned int substeps, bool spherical)
	: RandomAlgorithm(ds,v,0), spherical(spherical), substeps(substeps)
{
	this->gen = gen;
	emforgmm = new EMforGMM(ds, v, this->gen, spherical);
}

void RepeatEMforGMM::setIndicesOfInputPointsToBeUsed(std::vector< idx_type > indices)
{
	this->emforgmm->setIndicesOfInputPointsToBeUsed(indices);
}



void RepeatEMforGMM::run(unsigned int numSteps)
{
	for(int st=0; st<numSteps; ++st)
	{
		this->gen.seed(++this->seed);
		Parameters init = initializer::adaptiveMeansToGMM(*this->input, this->desc.components(), this->gen);
		this->emforgmm->init(init);
		this->emforgmm->run(substeps);
		this->desc = this->emforgmm->getGMMDesc();
		fp_type cost = gmmutil::gmmNLL(this->input->points, this->desc).sum();
		if(best_desc.empty() || cost < best_cost)
		{
			best_desc = this->desc;
			best_cost = cost;
		}
		else if(cost > best_cost)
		{
			this->desc = best_desc;
		}
	}
}


