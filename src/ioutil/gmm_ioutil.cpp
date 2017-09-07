#include "gmm_ioutil.h"

#include "ioutil.h"
#include "../settings/settings.h"

void GMMIOUtil::appendToGMM(std::string const path, Parameters const& desc, std::string const tag)
{
	std::string name = path;
	name.append(GMM_FILE_EXTENSION);

	if (commonSettings().verbose)
		std::cout << "   appending gmm to file " << name << std::endl;
	std::ofstream ofs(name.c_str(), std::ios_base::binary | std::ios_base::out | std::ios_base::app);

	// write tag
	uint32_t taglen = tag.size();
	ofs.write((char*)&taglen, sizeof(uint32_t));
	ofs.write(tag.c_str(), taglen);

	// write dimension
	uint32_t d = desc.means.rows();
	ofs.write((char*)&d, sizeof(uint32_t));

	// write number of components
	uint32_t k = desc.means.cols();
	ofs.write((char*)&k, sizeof(uint32_t));

	double buffer;

	// write weights
	for (idx_type i=0; i<k; ++i)
	{
		buffer = desc.weights[i];
		ofs.write((char*)&buffer, sizeof(double));
	}

	// write means
	for (idx_type i=0; i<k; ++i)
		for (idx_type j=0; j<d; ++j)
		{
			buffer = desc.means(j,i);
			ofs.write((char*)&buffer, sizeof(double));
		}

	// write covariances
	for (idx_type i=0; i<k; ++i)
		for (idx_type j=0; j<d; ++j)
			for (idx_type l=0; l<d; ++l)
			{
				buffer = desc.covariances[i](j,l);
				ofs.write((char*)&buffer, sizeof(double));
			}

	ofs.close();
}

void GMMIOUtil::loadFromGMM(std::string const path, std::vector<Parameters>& gmms, std::vector<std::string>& tags)
{
	if(commonSettings().verbose)
		std::cout << "   reading from file " << path << std::endl;

	std::ifstream ifs(path.c_str(), std::ios_base::binary | std::ios_base::in);
	if(!ifs.is_open())
		gmmlab_throw("Failed to open file.");
	ifs.seekg(0,std::ifstream::end);
	std::streampos end_of_file = ifs.tellg();
	ifs.seekg(0,std::ifstream::beg);

	Parameters desc;
	unsigned int count = 0;
	while (ifs.tellg() < end_of_file)
	{
		// read tag
		uint32_t taglen;
		ifs.read((char*)&taglen, sizeof(uint32_t));
		std::string tag;
		for (uint32_t i=0; i<taglen; ++i)
		{
			char c;
			ifs.read(&c, sizeof(char));
			tag += c;
		}
		tags.push_back(tag);

		// read dimension
		uint32_t d;
		ifs.read((char*)&d, sizeof(uint32_t));


		// read number of components
		uint32_t k;
		ifs.read((char*)&k, sizeof(uint32_t));

		if (commonSettings().verbose)
			std::cout << "      found gmm with tag " << tag << ", dimension " << d << " and k=" << k << std::endl;

		double buffer;

		// read weights
		desc.weights = Vector(k);
		for (idx_type i=0; i<k; ++i)
		{
			ifs.read((char*)&buffer, sizeof(double));
			desc.weights[i] = buffer;
		}

		// read means
		desc.means = Matrix(d,k);
		for (idx_type i=0; i<k; ++i)
			for (idx_type j=0; j<d; ++j)
			{
				ifs.read((char*)&buffer, sizeof(double));
				desc.means(j,i) = buffer;
			}

		// read means
		desc.covariances.clear();
		for (idx_type i=0; i<k; ++i)
		{
			Matrix cov(d,d);
			for (idx_type j=0; j<d; ++j)
				for (idx_type l=0; l<d; ++l)
				{
					ifs.read((char*)&buffer, sizeof(double));
					cov(j,l) = buffer;
				}
			desc.covariances.push_back(cov);
		}

		if (ifs.eof())
			std::cout << "ERROR! - " << path << " corrupted." << std::endl;
		else
		{
			gmms.push_back(desc);
			count++;
		}
	}
	ifs.close();
	if (commonSettings().verbose)
		std::cout << "   " << count << " gmm(s) read." << std::endl;
}

// void GMMIOUtil::loadGMMs(std::string const pattern, std::vector<GMMDesc>& gmms,
// 	std::vector<std::string>& tags)
// {
// 	std::cout << "Loading GMMs from files starting with: " << pattern << std::endl;
// 	idx_type presize = gmms.size();
// 
// 	// get filenames from pattern and sort alphabetically
// 	std::vector<std::string> filelist;
// 	IOUtil::getFilenames(pattern, filelist);
// 	std::sort(filelist.begin(), filelist.end());
// 
// 	if (filelist.empty())
// 		std::cout << "   WARNING - No files found!" << std::endl;
// 	else
// 	{
// 		// check for common file extension .gmm
// 		std::string common = IOUtil::commonFileExtension(filelist);
// 		boost::to_upper(common);
// 		FileFormat ff = UNKNOWN;
// 		if (common == "GMM")
// 			ff = GMM;
// 		else if (common == "LAB")
// 			ff = LAB;
// 		else
// 			std::cout << "   WARNING - Unknown or diverse file extensions!" << std::endl;
// 
// 		// iterate over list and load files
// 		for (std::vector<std::string>::iterator it = filelist.begin(); it!=filelist.end(); ++it)
// 			switch (ff)
// 			{
// 				case GMM:
// 					loadFromGMM(*it, gmms, tags);
// 					break;
// 				case LAB:
// 					break;
// 			}
// 	}
// 	std::cout << "   " << gmms.size()-presize << " gmms read." << std::endl;
// }
