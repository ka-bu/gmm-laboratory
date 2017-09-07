#include "empty.h"

#include "../../../base.h"
#include "../../../base/parameters.h"

const std::string EmptyID::CLASSTAG = "EmptyInit";

void EmptyID::parseDescriptionString(std::vector< InitDescription* >& descriptions, std::string descriptionstring)
{
	descriptions.push_back(new EmptyID());
}


EmptyID::EmptyID()
{
	std::stringstream sstream;
	sstream << "Empty";
	this->nametag = sstream.str();
}

Parameters EmptyID::compute(commonutil::DataSet const& input, unsigned int k) const
{
	return initializer::emptyGMM();
}


Parameters initializer::emptyGMM()
{
	Parameters gmmdesc;
	gmmdesc.weights.resize(0);
	gmmdesc.means.resize(0,0);
	gmmdesc.covariances.resize(0);	
	return gmmdesc;
}


