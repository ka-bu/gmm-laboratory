#ifndef RESTARTLLOYDS_H
#define RESTARTLLOYDS_H

#include "../base/randomalgorithm.h"
#include "../description/algo_description.h"

#include "../../base/parameters.h"
#include "../../base/commonutil.h"

#include <set>
#include <ctime>

/**
 *  
 * @brief  
 *
 * @ingroup learning_algorithms
 */
class RestartLloyds : public RandomAlgorithm
{
public:
	/**
	* @brief Constructor initializing input and configuring the algorithm.
	* @remarks Remember to init the mixture model.
	*/
	RestartLloyds(commonutil::DataSet const& ds, bool verbose, unsigned int subst, uint32_t seed, fp_type lloydsRatio = 0);

	virtual ~RestartLloyds()
	{
	}

	virtual void init(Parameters const& desc);
	virtual void run(unsigned int numSteps = 1);

	/**
	* @brief does a dynamic cast of the given GMMAlgorithm to RestartLloyds
	* @return NULL if the GMMAlgorithm is not a RestartLloyds instance
	*/
	static RestartLloyds* toRestartLloyds(Algorithm* a);

private:
	fp_type lloydsRatio;
	unsigned int subalgoSteps;
	fp_type minCost;
	bool justInitialized;
};

/**
 *@brief: description 
 */
class RestartLloydsAD : public AlgoDescription
{
public:
	static const std::string CLASSTAG;
	RestartLloydsAD(unsigned int steps, uint32_t seed, fp_type lloydsRatio, unsigned int substeps);
	virtual Algorithm* create(const commonutil::DataSet& input) const;
	
	static void parseDescriptionString(std::vector<AlgoDescription*>& descriptions, const std::string descriptionstring, const bool first = false);
private:
	uint32_t seed;
	fp_type lloydsRatio;
	unsigned int subSteps;
	
};

#endif