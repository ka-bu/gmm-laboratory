#include "algo_description.h"

#include <iostream>
#include <fstream>
#include <boost/algorithm/string/erase.hpp>

#include "../base/algorithm.h"
#include "../gmm/emforgmm.h"
#include "../gmm/repeatemforgmm.h"
#include "../gmm/lloydsforgmm.h"
#include "../base/emptyalgorithm.h"
#include "../kmeans/lloyds.h"
#include "../kmeans/restartlloyds.h"
#include "../gmm/init_algo/alternating_adaptivemeans_means2gmm_initalgo.h"


AlgoDescription* AlgoDescription::getNext() const
{
	return this->next;
}

void AlgoDescription::setNext(AlgoDescription* next)
{
	this->next = next;
}

unsigned int AlgoDescription::getSteps() const
{
	return this->steps;
}

void AlgoDescription::parseDescriptionString(std::vector< AlgoDescription* > descriptions, std::string descriptionstring, bool first)
{
	std::cout << "Don't know how to parse " << descriptionstring << " yet." << std::endl;
	gmmlab_throw("Not implemented yet.");
}


std::vector< AlgoDescription* > AlgoDescription::createAlgoDescriptions(std::string filename)
{
	std::vector<AlgoDescription*> descriptions;
	std::cout << "Creating algo descriptors from file: " << filename << std::endl;
	std::ifstream filestream;
	filestream.open(filename.c_str());

	std::string line;
	std::vector<AlgoDescription*> last;
	
	bool empty = true;
	while(getline(filestream, line))
	{
		bool first = false;
		std::istringstream iss(line);
		std::string entry;
		
		entry = Description::readUpperString(iss);
		boost::erase_all(entry, " ");
		if (entry=="FIRST")
		{
			first = true;
			last.clear();
		}
		else if (entry=="NEXT")
		{
			if (last.empty())
				gmmlab_throw("next:... before first:... forbidden!");
		}
		else if (entry.empty())
		{
			last.clear();
			continue; // skip comments
		}
		
		empty = false;

		std::string algo_id = Description::readString(iss);
		std::string algo_description;
		getline(iss, algo_description);
		
		// ------------------------------------------------------------------------------------------
		
		/** create description by id string **/
		
		// TODO: register the other algorithms
		
		std::vector<AlgoDescription*> ad;
		if(algo_id == EMforGMMAD::CLASSTAG)
			EMforGMMAD::parseDescriptionString(ad, algo_description, first);
		else if(algo_id == RepeatEMforGMMAD::CLASSTAG)
			RepeatEMforGMMAD::parseDescriptionString(ad, algo_description, first);
		else if(algo_id == LloydsForGMMAD::CLASSTAG)
			LloydsForGMMAD::parseDescriptionString(ad, algo_description, first);
		else if(algo_id == EmptyAlgorithmAD::CLASSTAG)
			EmptyAlgorithmAD::parseDescriptionString(ad, algo_description, first);
		else if(algo_id == LloydsAD::CLASSTAG)
			LloydsAD::parseDescriptionString(ad, algo_description, first);
		else if(algo_id == RestartLloydsAD::CLASSTAG)
			RestartLloydsAD::parseDescriptionString(ad, algo_description, first);
		else if(algo_id == AlternatingAdaptiveMeans2GMMAD::CLASSTAG)
			AlternatingAdaptiveMeans2GMMAD::parseDescriptionString(ad, algo_description, first);
		else
		{
			std::cout << "AlgoDescription doesn't know the id " << algo_id << " and cannot create a description of " << algo_id << ":" << algo_description << std::endl;
			gmmlab_throw("Unknown algo_id.");
		}
		
		
		// ------------------------------------------------------------------------------------------
		
		
		if(first)
		{
			descriptions.reserve(descriptions.size()+ad.size());
			std::move(ad.begin(), ad.end(), std::inserter(descriptions,descriptions.end()));
			last.reserve(last.size()+ad.size());
			std::move(ad.begin(), ad.end(), std::inserter(last, last.end()));
		}
		else
		{
			if(ad.size()!=1)
				gmmlab_throw("There was more than one algorithm given as 'next'.");
			for (std::size_t i=0; i<last.size(); ++i)
			{
				last.at(i)->setNext(ad.back());
				last.at(i) = ad.back();
			}
		}
		
		ad.clear();
		
	}
	
	if(empty)
		std::cout << "Warning: algo description file is empty." << std::endl;
	
	return(descriptions);
	
}
