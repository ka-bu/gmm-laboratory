#include "file_description.h"

#include "../ioutil/data_ioutil.h"

FileDescription::FileDescription(std::string const& p, unsigned int fc, unsigned int nc,
	unsigned int bl, bool nrm)
	: prefix(p), firstCoord(fc), numCoords(nc), blockLength(bl), normalize(nrm)
{
	std::stringstream sstream;
	sstream << p.substr(p.find_last_of('/')+1) << "_" << fc << "-" << fc+nc-1 << "." << bl;
	if (nrm)
		sstream << "_normalized";
	this->nametag = sstream.str();
}

void FileDescription::add_k(unsigned int k)
{
	this->kList.push_back(k);
}
	
Parameters FileDescription::retrieve(commonutil::DataSet& input)
{
	DataIOUtil::loadData(input, this->prefix, this->firstCoord, this->numCoords, this->blockLength, this->normalize, true); // noFrameDrop = true
	
	return Parameters();
}
