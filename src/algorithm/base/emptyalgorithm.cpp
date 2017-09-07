#include "emptyalgorithm.h"




const std::string EmptyAlgorithmAD::CLASSTAG = "EmptyAlgorithm";


EmptyAlgorithmAD::EmptyAlgorithmAD()
{
	this->nametag = CLASSTAG;
}

Algorithm* EmptyAlgorithmAD::create(commonutil::DataSet const& input) const
{
	return new EmptyAlgorithm();
}

void EmptyAlgorithmAD::parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first)
{	
	descriptions.push_back(new EmptyAlgorithmAD());
		
}

// =================================================================================================================================


EmptyAlgorithm::EmptyAlgorithm()
	: Algorithm(commonutil::DataSet(), false)
{
}



void EmptyAlgorithm::run(unsigned int numSteps)
{
}

EmptyAlgorithm* EmptyAlgorithm::toEmptyAlgorithm(Algorithm* a)
{
	return dynamic_cast<EmptyAlgorithm*>(a);
}

