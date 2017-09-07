#include "data_description.h"

#include <fstream>

#include "gen_umrc_description.h"
#include "file_description.h"
#include "../settings/settings.h"


std::string DataDescription::tag()
{
	return this->nametag;
}

DataDescription::kIterator DataDescription::first_k() const
{
	return this->kList.begin();
}

DataDescription::kIterator DataDescription::last_k() const
{
	return this->kList.end();
}

std::vector<DataDescription*> DataDescription::createDataDescriptions(std::string const& filename)
{
	std::ifstream filestream;
	std::vector<DataDescription*> descriptors;
	std::cout << "Creating data descriptors from file: " << filename << std::endl;

	filestream.open(filename.c_str());
	std::string line;
	while (getline(filestream, line))
	{
		std::istringstream iss(line);
		
		std::string dataType;
		dataType = readUpperString(iss);
		if(dataType.empty())
			continue; // skip comments and empty lines
		
		if (dataType=="LOAD")
		{
			// read file descriptor from strings formatted as:
			//		load:<filename>:<first coordinate>-<last coordinate>.<block length>:<k range specs>
			
			std::string prefix = readString(iss);

			std::string entry = readString(iss);
			unsigned int from, to, block;
			readRangeSpec(entry, from, to, block);
			unsigned int num = to<from?0:to-from+1;
			
			std::vector<unsigned int> kList = readIntList(iss);
			
			entry = readUpperString(iss);
			bool normalize = (entry=="NORMALIZE");
			
			FileDescription* fdd = new FileDescription(prefix, from, num, block, normalize);
			for (std::size_t i=0; i<kList.size(); ++i)
			{
				fdd->add_k(kList[i]);
			}
			
			descriptors.push_back(fdd);
		}
		else if (dataType=="GEN")
		{
			
			// read data generation descriptor from strings formatted as:
			//		gen:genMode: ...  
				
			
			std::string genMethod = readUpperString(iss);
			
			if(genMethod == "UMRC")
			{

				// uniform_means_and_restricted_covar:
				//  	 gen:umrc:<size>:<dim>:<k>
				//               :<Wexp fp-range-specs>
				//               :<separation fp-range-specs>
				//               :<EWexp fp-range-specs>:<EWmin fp>:<EWmax fp>:<EWPropMin fp>:<EWPropMax fp>
				//               :<ME string>
				//               :<seed range-specs>
				
				unsigned int size = readInt(iss);
				unsigned int d = readInt(iss);
				std::vector<unsigned int> kList = readIntList(iss);
				std::vector<fp_type> wList = readFpTypeList(iss);
				std::vector<fp_type> sepList = readFpTypeList(iss);
				std::vector<fp_type> ewExpList = readFpTypeList(iss);
				fp_type minMinSqrtEW = readFpType(iss);
				fp_type maxMinSqrtEW = readFpType(iss);
				fp_type minSqrtEWProp = readFpType(iss);
				fp_type maxSqrtEWProp = readFpType(iss);
				std::string me = readUpperString(iss);
				std::vector<unsigned int> gList = readIntList(iss);
				
				for(std::size_t g=0; g<gList.size(); ++g)
					for(std::size_t w=0; w<wList.size(); ++w)
						for(idx_type k=0; k<kList.size(); ++k)
							for(std::size_t sep=0; sep<sepList.size(); ++sep)
								for(std::size_t ewExp=0; ewExp<ewExpList.size(); ++ewExp)
									descriptors.push_back(
										new GenUMRCDescription(
											size, d, kList[k], wList[w], 
											sepList[sep],
											ewExpList[ewExp], minMinSqrtEW, maxMinSqrtEW, minSqrtEWProp, maxSqrtEWProp, 
											me, gList[g])
										);
				
			}
			else
			{
					gmmlab_throw("Unknown generation method.");
			}
		}
	}
	filestream.close();
	
	return descriptors;
}
