#ifndef SAMPLINGALGORITHM_H
#define SAMPLINGALGORITHM_H

#include "randomalgorithm.h"

#include "../../base.h"
#include "../../base/commonutil.h"

/**
* @brief BBKSS algorithm for GMMs
*
* @ingroup learning_algorithms
*/
class SamplingAlgorithm : public RandomAlgorithm
{
public:
	/**
	* @brief Constructor initializing input and configuring the algorithm.
	* @remarks Remember to init the mixture model.
	*/
	SamplingAlgorithm(commonutil::DataSet const& input, bool verbose, uint32_t seed, fp_type ratio = 0);

	virtual ~SamplingAlgorithm()
	{
	}

	commonutil::DataSet const& getSamples() const;

	/**
	* @brief does a dynamic cast of the given GMMAlgorithm to SamplingAlgorithm
	* @return NULL if the GMMAlgorithm is not a SamplingAlgorithm instance
	*/
	static SamplingAlgorithm* toSamplingAlgorithm(Algorithm* a);

protected:
	commonutil::DataSet samples;
	fp_type ratio;
};

#endif
