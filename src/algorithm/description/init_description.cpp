#include "init_description.h"

#include <fstream>

#include "../kmeans/init/uniform_means.h"
#include "../kmeans/init/adaptive_means.h"
#include "../gmm/init/empty.h"
#include "../gmm/init/adaptive_means2gmm.h"
#include "../gmm/init/adaptive_spherical.h"
#include "../gmm/init/agglomerative_init.h"
#include "../gmm/init/alternating_adaptivemeans_means2gmm.h"
#include "../gmm/init/alternating_adaptivemeans_means2gmm_then_lloyds.h"
#include "../gmm/init/alternating_adaptivemeans_means2gmm_then_lloydsforgmm.h"
#include "../gmm/init/gonzalez2gmm.h"
#include "../gmm/init/gonzalez_forgmm.h"
#include "../gmm/init/gonzalez_kwedlo.h"
#include "../gmm/init/gonzalez_lloyds_means2gmm.h"
#include "../gmm/init/gonzalezforgmm_then_lloydsforgmm.h"
#include "../gmm/init/resample_gmminit.h"
#include "../gmm/init/restart_init.h"
#include "../gmm/init/uniform_lloyds_means2gmm.h"
#include "../gmm/init/uniform_spherical.h"
#include "../gmm/init/uniform_spherical_with_pruning.h"
#include "../gmm/init/uniform_means2gmm.h"


void InitDescription::parseDescriptionString(std::vector<InitDescription*>& descriptions, std::string descriptionstring, bool first)
{
	std::cout << "Don't know how to parse " << descriptionstring << " yet." << std::endl;
	gmmlab_throw("Not yet implemented.");
}

void initializer::check(idx_type k, idx_type d, idx_type n)
{
	if(k==0)
		gmmlab_throw("Number of means is zero.")
	if(d==0 || n==0)
		gmmlab_throw("Input is empty.");
}

std::vector< InitDescription* > InitDescription::createInitDescriptions(std::string filename)
{
	std::vector<InitDescription*> descriptions;
	std::cout << "Creating init descriptors from file: " << filename << std::endl;
	std::ifstream filestream;
	filestream.open(filename.c_str());

	std::string line;
	
	bool empty = true;
	while(getline(filestream, line))
	{
		bool first = false;
		std::istringstream iss(line);
		std::string entry;
		
		std::string init_id = Description::readString(iss);
		std::string init_description;
		getline(iss, init_description);
		
		// skip comments and empty lines
		if(init_id.empty())
			continue;
		
		empty = false;
		
		// ------------------------------------------------------------------------------------------
		
		/** create description by id string **/
		
		if(init_id == EmptyID::CLASSTAG)
			EmptyID::parseDescriptionString(descriptions, init_description);
// 		
// 		// means
 		if(init_id == UniformMeansID::CLASSTAG)
 			UniformMeansID::parseDescriptionString(descriptions, init_description);
		else if(init_id == AdaptiveMeansID::CLASSTAG)
			AdaptiveMeansID::parseDescriptionString(descriptions, init_description);
// 		else if(init_id == GonzalezMeansID::CLASSTAG)
// 			GonzalezMeansID::parseDescriptionString(descriptions, init_description);
// 		else if(init_id == ResampleAdaptiveMeansID::CLASSTAG)
// 			ResampleAdaptiveMeansID::parseDescriptionString(descriptions, init_description);
// 		
		
		// gmm
		else if(init_id == AdaptiveLloydsMeans2GMMID::CLASSTAG)
			AdaptiveLloydsMeans2GMMID::parseDescriptionString(descriptions, init_description);
		else if(init_id == AdaptiveMeans2GMMID::CLASSTAG)
			AdaptiveMeans2GMMID::parseDescriptionString(descriptions, init_description);
		else if(init_id == AdaptiveSphericalID::CLASSTAG)
			AdaptiveSphericalID::parseDescriptionString(descriptions, init_description);
		else if(init_id == AgglomerativeInitID::CLASSTAG)
			AgglomerativeInitID::parseDescriptionString(descriptions, init_description);
		else if(init_id == AlternatingAdaptiveMeansAndMeans2GMMID::CLASSTAG)
			AlternatingAdaptiveMeansAndMeans2GMMID::parseDescriptionString(descriptions, init_description);
		else if(init_id == AlternatingAdaptiveMeans2GMMThenLloydsID::CLASSTAG)
			AlternatingAdaptiveMeans2GMMThenLloydsID::parseDescriptionString(descriptions, init_description);
		else if(init_id == AlternatingAdaptiveMeans2GMMThenLloydsForGMMID::CLASSTAG)
			AlternatingAdaptiveMeans2GMMThenLloydsForGMMID::parseDescriptionString(descriptions, init_description);
		else if(init_id == Gonzalez2GMMID::CLASSTAG)
			Gonzalez2GMMID::parseDescriptionString(descriptions, init_description);
		else if(init_id == GonzalezForGMMID::CLASSTAG)
			GonzalezForGMMID::parseDescriptionString(descriptions, init_description);
		else if(init_id == GonzalezKwedloID::CLASSTAG)
			GonzalezKwedloID::parseDescriptionString(descriptions, init_description);
		else if(init_id == GonzalezLloydsMeans2GMMID::CLASSTAG)
			GonzalezLloydsMeans2GMMID::parseDescriptionString(descriptions, init_description);
		else if(init_id == GonzalezForGMMThenLloydsForGMMID::CLASSTAG)
			GonzalezForGMMThenLloydsForGMMID::parseDescriptionString(descriptions, init_description);
		else if(init_id == ResampleGMMInitID::CLASSTAG)
			ResampleGMMInitID::parseDescriptionString(descriptions, init_description);
		else if(init_id == RestartInitID::CLASSTAG)
			RestartInitID::parseDescriptionString(descriptions, init_description);
		else if(init_id == UniformLloydsMeans2GMMID::CLASSTAG)
			UniformLloydsMeans2GMMID::parseDescriptionString(descriptions, init_description);
		else if(init_id == UniformSphericalID::CLASSTAG)
			UniformSphericalID::parseDescriptionString(descriptions, init_description);
		else if(init_id == UniformSphericalWithPruningID::CLASSTAG)
			UniformSphericalWithPruningID::parseDescriptionString(descriptions, init_description);
		else if(init_id == UniformMeans2GMMID::CLASSTAG)
			UniformMeans2GMMID::parseDescriptionString(descriptions, init_description);
		
		
		// don't know
		else
		{
			std::cout << "InitDescription doesn't know the id " << init_id << " and cannot create a description using " << init_id << ":" << init_description << std::endl;
			//gmmlab_throw("Unknown init_id.");
		}
		
		
		// ------------------------------------------------------------------------------------------
		
	}
	
	if(empty)
		std::cout << "Warning: init description file is empty." << std::endl;
	
	return(descriptions);
	
}
