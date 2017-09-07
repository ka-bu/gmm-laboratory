#ifndef RANDOMALGORITHM_H
#define RANDOMALGORITHM_H

#include "../../base/parameters.h"
#include "../../base/commonutil.h"

#include <set>
#include <ctime>
#include <vector>
#include "algorithm.h"

/**
* @brief BBKSS algorithm for GMMs
*
* @ingroup learning_algorithms
*/
class RandomAlgorithm : public Algorithm
{
public:
	/**
	* @brief Constructor initializing input and configuring the algorithm.
	* @remarks Remember to init the mixture model.
	*/
	RandomAlgorithm(commonutil::DataSet const& dataset, bool verbose, uint32_t seed);
	
	/**
	* @brief Constructor initializing input and configuring the algorithm.
	* @remarks Remember to init the mixture model. The seed will be set to 0.
	*/
	RandomAlgorithm(commonutil::DataSet const& dataset, bool verbose, std::mt19937& gen);

	virtual ~RandomAlgorithm()
	{
	}

	/**
	* @brief does a dynamic cast of the given GMMAlgorithm to SamplingAlgorithm
	* @return NULL if the GMMAlgorithm is not a SamplingAlgorithm instance
	*/
	static RandomAlgorithm* toRandomAlgorithm(Algorithm* a);

protected:
	uint32_t seed;
	std::mt19937 gen;
};

#endif
