#include "randomalgorithm.h"

RandomAlgorithm::RandomAlgorithm(commonutil::DataSet const& ds, bool v, uint32_t s)
	: Algorithm(ds,v), seed(s), gen(s)
{
	this->seed++;
}

RandomAlgorithm::RandomAlgorithm(commonutil::DataSet const& ds, bool v, std::mt19937& gen)
	: Algorithm(ds,v), seed(0)
{
	this->gen = gen;
}

RandomAlgorithm* RandomAlgorithm::toRandomAlgorithm(Algorithm* a)
{
	return dynamic_cast<RandomAlgorithm*>(a);
}

