#ifndef ALTERNATING_ADAPTIVE_MEANS2GMM_INITALGO_H
#define ALTERNATING_ADAPTIVE_MEANS2GMM_INITALGO_H

#include "../../../base/parameters.h"
#include "../../description/algo_description.h"
#include "../../base/randomalgorithm.h"
#include "../utils/gmmutil.h"


/**
 * @brief Descriptor for Initialization
 */
class AlternatingAdaptiveMeans2GMMInitAlgo : public RandomAlgorithm
{
public:
	AlternatingAdaptiveMeans2GMMInitAlgo(commonutil::DataSet const& input, bool verbose, bool gonzalez_mode, fp_type nonUniformFactor, uint32_t seed = 1);	
	virtual ~AlternatingAdaptiveMeans2GMMInitAlgo()
	{
		
	}
	virtual void run(unsigned int numSteps = 1);

private:
	uint32_t seed;
	bool gonzalez_mode;
	fp_type nonUniformFactor;
};


/**
 * @brief Descriptor for Initialization
 */
class AlternatingAdaptiveMeans2GMMAD : public AlgoDescription
{
public:
	static const std::string CLASSTAG;
	AlternatingAdaptiveMeans2GMMAD(unsigned int steps, bool gonzalez_mode, fp_type nonUniformFactor, uint32_t seed = 1);
	virtual Algorithm* create(commonutil::DataSet const& input) const;
	static void parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first = false);

private:
	uint32_t seed;
	bool gonzalez_mode;
	fp_type nonUniformFactor;
};

#endif
