#include "cost_outil.h"

CostOUtil::CostOUtil(): IOUtil(".csv")
{

}


void CostOUtil::header()
{
	this->filestream << "round" << SEPARATOR << "runtime" << SEPARATOR << "cost";
}

void CostOUtil::store(std::size_t round, double runtime, fp_type costs)
{
	this->filestream << round << SEPARATOR << runtime << SEPARATOR << costs << std::endl;
}

void CostOUtil::storeMarker(std::size_t mainrnd, double maintime, std::size_t minrnd, double mintime)
{
	if(!is_open())
		gmmlab_throw("File not open.");
	this->filestream << mainrnd << SEPARATOR << maintime << std::endl;
	this->filestream << minrnd << SEPARATOR << mintime << std::endl;
}

void CostOUtil::storeMarkerVector(std::vector<std::size_t> endRounds, std::vector<double> endTimes, std::size_t minrnd, double mintime)
{
	if(!is_open())
		gmmlab_throw("File not open.");
	if (endRounds.size()!=endTimes.size())
		gmmlab_throw("round and time vectors have different length!");
		
	for (std::size_t i=0; i<endRounds.size(); ++i)
		this->filestream << endRounds[i] << SEPARATOR << endTimes[i] << std::endl;
	this->filestream << minrnd << SEPARATOR << mintime << std::endl;
}